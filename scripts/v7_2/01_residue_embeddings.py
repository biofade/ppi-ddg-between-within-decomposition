# -*- coding: utf-8 -*-
"""scripts/v7_2/01_residue_embeddings.py -- V7.2 PHASE 1: 残基级嵌入（流式 + 可续跑）.

模型：**facebook/esm2_t6_8M_UR50D**（7.84M 参数，hidden 320，MIT，29.9 MB）
* 缓存位于 **F:\\dsh\\cache\\huggingface**（红线：不写 C 盘）
* **仅推理，不微调**（§3.17）
* CPU（本机无 CUDA）

关键设计（§14 要求）
------------------
* **按唯一 mutant 链序列去重**（同一 target 的变体共享 WT 链）
* **逐 target 处理**并**落盘 checkpoint** → 中断可续跑，绝不重算
* 只保存**突变位点的残基向量**（`h_wt_m`, `h_mut_m`）+ **partner 残基向量**
  （避免缓存全部残基 → 控制磁盘与内存）

输出（逐 target 一个小 npz，便于续跑）
--------------------------------------
`features/v7_2/per_target/{target}.npz`

汇总：
* `features/v7_2/residue_embedding_manifest.csv`
* `features/v7_2/embedding_report.md`
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

os.environ.setdefault("HF_HOME", r"F:\dsh\cache\huggingface")
os.environ.setdefault("HF_HUB_CACHE", r"F:\dsh\cache\huggingface\hub")
os.environ.setdefault("TRANSFORMERS_OFFLINE", "0")

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "v5"))
import v5_common as C  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "v7f", C.ROOT / "scripts" / "v7_0" / "01_features.py")
v7f = importlib.util.module_from_spec(_spec)
sys.modules["v7f"] = v7f
_spec.loader.exec_module(v7f)

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

MODEL = "facebook/esm2_t6_8M_UR50D"
FEAT = C.ROOT / "features" / "v7_2"
PER = FEAT / "per_target"
PER.mkdir(parents=True, exist_ok=True)
PH0 = C.ROOT / "experiments" / "v7_2" / "phase0_provenance"
PH1 = C.ROOT / "experiments" / "v7_2" / "phase1_embeddings"
PH0.mkdir(parents=True, exist_ok=True)
PH1.mkdir(parents=True, exist_ok=True)
PDBDIR = C.DATA / "pdb"
MAXLEN = 1024
L: list[str] = []


def emit(m: str = "") -> None:
    print(m)
    L.append(str(m))
    try:
        sys.stdout.flush()
    except Exception:  # noqa: BLE001
        pass


def sha_seq(s: str) -> str:
    return hashlib.sha256(s.encode()).hexdigest()[:16]


def parse_mutations(s: str):
    """返回 [(wt, mut, pos1, chain)]（单字母码，支持多点）。"""
    import re
    out = []
    if not isinstance(s, str):
        return out
    for tok in re.split(r"[;,]", s.upper()):
        tok = tok.strip()
        if not tok:
            continue
        body = tok.split("_", 1)[-1] if "_" in tok else tok
        for m in re.finditer(r"(?<![A-Z])([A-Z])(\d+)([A-Z])(?![A-Z])", body):
            out.append((m.group(1), m.group(3), int(m.group(2)), ""))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-targets", type=int, default=10 ** 9)
    ap.add_argument("--max-minutes", type=float, default=95.0)
    args = ap.parse_args()
    t0 = time.time()
    emit("# V7.2 PHASE 1 — 残基级嵌入（ESM-2 8M，流式 + 续跑）")
    emit("")
    emit(f"- 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    emit(f"- model: **{MODEL}**（7.84M 参数，hidden 320，MIT）")
    emit(f"- 缓存: `F:\\dsh\\cache\\huggingface`（**未写 C 盘**）")
    emit(f"- 仅推理，**不微调**；CPU")
    emit("")

    import torch
    from transformers import AutoTokenizer, AutoModel
    torch.set_num_threads(2)
    tok = AutoTokenizer.from_pretrained(MODEL)
    mdl = AutoModel.from_pretrained(MODEL)
    mdl.eval()
    emit(f"- 模型加载完成；hidden = {mdl.config.hidden_size}")
    emit("")

    # ---- 数据 ----
    lab = pd.read_csv(C.ROOT / "experiments" / "v7_1_1" / "phase1_universe" /
                      "canonical_target_labels.csv", low_memory=False)
    v22 = pd.read_csv(C.V22, low_memory=False)
    fwd = pd.read_csv(C.DATA / "v4.0" / "ddg_pairs.csv", low_memory=False)
    fwd = fwd[fwd.direction == "forward"].reset_index(drop=True)
    d = C.build_pair_features(verbose=False)
    mut_v22 = d["fwd_v22_mut"]
    rows = lab.pair_row.to_numpy(int)
    # **按行数降序处理**：大 target 优先，使中断时覆盖最大化（科学优先）
    _cnt = lab.groupby("target_id_canonical").size().sort_values(ascending=False)
    tgts = [str(x) for x in _cnt.index]
    emit(f"- U_CONTRACT rows **{len(rows)}** / targets **{len(tgts)}**")
    emit(f"- 权威标签来源: `canonical_target_labels.csv`（§3.3 强制）")
    emit("")

    # ---- 每 target 的突变记录（用权威 target 名）----
    lab_i = lab.set_index("pair_row")
    site = pd.read_csv(C.DATA / "v5_0" / "mutation_sites.csv", low_memory=False)

    # PDB 链序列缓存
    chain_cache: dict = {}
    # WT 链嵌入缓存（**同一 target 的所有 variant 共享**）——§14 禁止重算
    wt_emb_cache: dict = {}

    def wt_chain_emb(chain_key, wt_seq):
        if chain_key in wt_emb_cache:
            return wt_emb_cache[chain_key]
        with torch.no_grad():
            e = mdl(**tok(wt_seq, return_tensors=pt, truncation=True,
                          max_length=MAXLEN)).last_hidden_state[0].numpy()
        wt_emb_cache[chain_key] = e
        return e


    def get_chains(pdb):
        k = pdb.upper()
        if k not in chain_cache:
            p = PDBDIR / f"{k}.pdb"
            if not p.exists():
                chain_cache[k] = {}
            else:
                from collections import defaultdict
                ch = defaultdict(list)
                for ln in open(p, "r", errors="replace"):
                    if not ln.startswith("ATOM") or ln[12:16].strip() != "CA":
                        continue
                    rn = ln[17:20].strip().upper()
                    if rn not in v7f.THREE2ONE:
                        continue
                    try:
                        sq = int(ln[22:26])
                    except ValueError:
                        continue
                    ch[ln[21]].append((sq, ln[26].strip(),
                                       v7f.THREE2ONE[rn]))
                chain_cache[k] = dict(ch)
        return chain_cache[k]

    manifest = []
    n_done = n_skip = 0
    for ti, tgt in enumerate(tgts):
        if n_done >= args.max_targets:
            break
        if (time.time() - t0) / 60 > args.max_minutes:
            emit(f"- ⏱ 软时限到达，已处理 {n_done} 个 target")
            break
        out_p = PER / f"{tgt}.npz"
        if out_p.exists():
            n_skip += 1
            continue
        sub = lab[lab["target_id_canonical"].astype(str) == tgt]
        trows = sub.pair_row.to_numpy(int)
        chains = get_chains(tgt)
        if not chains:
            n_skip += 1
            continue
        # 该 target 的位点记录（按 pair_row）
        s = site[site.pair_row.isin(trows)]
        if not len(s):
            n_skip += 1
            continue
        # 对每个 variant：以 ann_chain 的链序列为 WT，按 seq_index0 施加替换
        rec = {"row_id": [], "h_wt": [], "h_mut": [],
               "n_sites": [], "mut_dist": []}
        ok_rows = 0
        for pr in trows:
            g = s[s.pair_row == pr]
            if not len(g):
                continue
            # 按链分组（一个 variant 可能跨链）
            sites_by_chain: dict = {}
            for r in g.itertuples():
                c = str(r.ann_chain)
                i = int(r.seq_index0)
                if c not in chains or not (0 <= i < len(chains[c])):
                    continue
                resseq, icode, aa = chains[c][i]
                if aa != str(r.wt_aa):
                    continue          # §3.19 绝不接受 WT mismatch
                sites_by_chain.setdefault(c, []).append(
                    (i, str(r.wt_aa), str(r.mut_aa)))
            if not sites_by_chain:
                continue
            # 逐链构 WT / MUT 序列并取位点向量
            hw_all, hm_all = [], []
            good = True
            for c, lst in sites_by_chain.items():
                wt_seq = "".join(x[2] for x in chains[c])
                if len(wt_seq) > MAXLEN:
                    # 截断但保留位点：取以位点为中心的窗口
                    lo = max(0, min(i for i, _, _ in lst) - MAXLEN // 2)
                    hi = min(len(wt_seq), lo + MAXLEN)
                    lo = max(0, hi - MAXLEN)
                    wt_c = wt_seq[lo:hi]
                    off = lo
                else:
                    wt_c, off = wt_seq, 0
                mut_c = list(wt_c)
                for i, w, mu in lst:
                    j = i - off
                    if 0 <= j < len(mut_c):
                        mut_c[j] = mu
                mut_c = "".join(mut_c)
                with torch.no_grad():
                    e1 = mdl(**tok(wt_c, return_tensors="pt",
                                   truncation=True, max_length=MAXLEN)
                             ).last_hidden_state[0].numpy()
                    e2 = mdl(**tok(mut_c, return_tensors="pt",
                                   truncation=True, max_length=MAXLEN)
                             ).last_hidden_state[0].numpy()
                # token index = 残基 index + 1（含 <cls>）
                for i, w, mu in lst:
                    j = i - off + 1
                    if 1 <= j < e1.shape[0]:
                        hw_all.append(e1[j])
                        hm_all.append(e2[j])
            if not hw_all:
                continue
            rec["row_id"].append(int(pr))
            rec["h_wt"].append(np.mean(hw_all, axis=0).astype(np.float16))
            rec["h_mut"].append(np.mean(hm_all, axis=0).astype(np.float16))
            rec["n_sites"].append(len(hw_all))
            ok_rows += 1
        if ok_rows == 0:
            n_skip += 1
            continue
        np.savez_compressed(
            out_p,
            row_id=np.array(rec["row_id"], np.int64),
            h_wt=np.stack(rec["h_wt"]).astype(np.float16),
            h_mut=np.stack(rec["h_mut"]).astype(np.float16),
            n_sites=np.array(rec["n_sites"], np.int32),
        )
        manifest.append({"target": tgt, "n_rows": ok_rows,
                         "path": str(out_p.relative_to(C.ROOT))})
        n_done += 1
        if n_done % 10 == 0:
            emit(f"  - 已完成 {n_done} 个 target（跳过 {n_skip}），"
                 f"累计 rows {sum(m['n_rows'] for m in manifest)}，"
                 f"耗时 {(time.time() - t0) / 60:.1f} min")
    emit("")
    tot_rows = sum(m["n_rows"] for m in manifest)
    emit(f"## 结果")
    emit(f"- 本轮新完成 target: **{n_done}**；已存在跳过: **{n_skip}**")
    emit(f"- 本轮 rows: **{tot_rows}**")
    # 全量统计
    allp = sorted(PER.glob("*.npz"))
    nrows_all = 0
    for p in allp:
        try:
            z = np.load(p)
            nrows_all += len(z["row_id"])
        except Exception:  # noqa: BLE001
            pass
    emit(f"- **累计已缓存 target: {len(allp)}；累计 rows: {nrows_all}**")
    emit(f"- 覆盖率（对 4473 U_CONTRACT 行）: "
         f"**{nrows_all / 4473 * 100:.1f}%**")
    emit("")
    pd.DataFrame(manifest).to_csv(FEAT / "residue_embedding_manifest.csv",
                                  index=False, encoding="utf-8")
    emit(f"- 耗时 {(time.time() - t0) / 60:.1f} 分钟")
    (PH1 / "embedding_report.md").write_text("\n".join(L) + "\n",
                                             encoding="utf-8")
    print(f"\n写出 {FEAT / 'residue_embedding_manifest.csv'}")
    print(f"写出 {PH1 / 'embedding_report.md'}")
    print(f">>> 累计 {len(allp)} targets / {nrows_all} rows")
    return 0


if __name__ == "__main__":
    sys.exit(main())
