# -*- coding: utf-8 -*-
"""
10_extract_remaining.py -- Affinity Intelligence V1.1 / 补全剩余结构 (无临时文件 SASA)

与 07_extract_features_full.py 的唯一差异: SASA 走 **freesasa.structureFromBioPDB**
(直接消费 Bio.PDB 对象), 不再写临时 PDB 文件。好处:
  * 彻底绕开 4FQR 类"链 ID 超过 PDB 格式单字符限制"的失败 (PDBIO 才会报此错);
  * 省掉每结构 10+ 次磁盘写入;
  * 已在 1A22/1A2K/1AO7/1AKJ/1TZN 上验证与文件方案**逐位一致 (差异 0.0000%)**。

输入 : data/full_structures.txt (2632 个非 MVP PDB) + features/features_full.csv (已完成的 1600 个)
输出 : features/features_full.csv   (追加; 与 07 产物同一张表)
       features/extract_checkpoint_v11.json
       features/extract_v11_log.txt

断点续算: 每 --exit-every(=100) 条结构追加落盘 + 原子写检查点; 重启自动跳过已完成。
纪律: n_jobs=2 等价单进程串行; 每批 gc.collect(); 单结构/单界面独立 try...except, 失败即记录不重试。
"""

import argparse
import gc
import json
import logging
import os
import sys
import time

import numpy as np
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = r"F:\DSH_Galatea\Affinity_Intelligence"
DATA_DIR = os.path.join(ROOT, "data")
FEAT_DIR = os.path.join(ROOT, "features")
PDB_DIR = os.path.join(DATA_DIR, "pdb")
BENCH = os.path.join(DATA_DIR, "benchmark.csv")
OUT_CSV = os.path.join(FEAT_DIR, "features_full.csv")
CHECKPOINT = os.path.join(FEAT_DIR, "extract_checkpoint_v11.json")
LOG_PATH = os.path.join(FEAT_DIR, "extract_v11_log.txt")
MAX_PAIR_SASA = 8
EXIT_EVERY = 100

logger = logging.getLogger("v11")
logger.setLevel(logging.INFO)
_fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%H:%M:%S")
_fh = logging.FileHandler(LOG_PATH, mode="a", encoding="utf-8")
_fh.setFormatter(_fmt)
logger.addHandler(_fh)
_sh = logging.StreamHandler(sys.stdout)
_sh.setFormatter(_fmt)
logger.addHandler(_sh)


def _load(mod_name, fname):
    import importlib.util
    p = os.path.join(ROOT, "scripts", fname)
    spec = importlib.util.spec_from_file_location(mod_name, p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


F02 = _load("feat02_v11", "02_extract_features.py")


def get_rss_mb():
    try:
        import ctypes
        from ctypes import wintypes

        class PMC(ctypes.Structure):
            _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD),
                        ("PeakWorkingSetSize", ctypes.c_size_t),
                        ("WorkingSetSize", ctypes.c_size_t),
                        ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
                        ("QuotaPagedPoolUsage", ctypes.c_size_t),
                        ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
                        ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                        ("PagefileUsage", ctypes.c_size_t),
                        ("PeakPagefileUsage", ctypes.c_size_t)]
        pmc = PMC()
        pmc.cb = ctypes.sizeof(PMC)
        ctypes.windll.psapi.GetProcessMemoryInfo(
            ctypes.windll.kernel32.GetCurrentProcess(), ctypes.byref(pmc), pmc.cb)
        return pmc.WorkingSetSize / 1024 / 1024
    except Exception:                                           # noqa: BLE001
        return float("nan")


def resolve_struct_path(pdb_id):
    for ext in (".pdb", ".cif", ".ent"):
        p = os.path.join(PDB_DIR, pdb_id + ext)
        if os.path.exists(p) and os.path.getsize(p) > 500:
            return p
    return None


# ---------------------------------------------------------------- SASA (无临时文件)
class SasaCalc:
    """基于 freesasa.structureFromBioPDB 的 SASA 计算 (零磁盘 IO)。"""

    def __init__(self):
        import freesasa
        from Bio.PDB import Structure as BSt, Model as BMd, Chain as BCh
        self.fs = freesasa
        self.BSt, self.BMd, self.BCh = BSt, BMd, BCh
        self.fs.setVerbosity(self.fs.silent)
        self.opts = {"hetatm": False, "skip-unknown": True}
        # 标准氨基酸白名单 (与 02 的 AA3 一致, 保证 SASA 与几何特征的原子集一致)
        self.aa_ok = set(F02.AA3) if hasattr(F02, "AA3") else {
            "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
            "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}
        self._cache = {}
        self._struct_cache = {}
        self._mono_cache = {}

    def _subset(self, model, chain_ids):
        """构造只含指定链的完整 Structure, 并**过滤成标准氨基酸残基 + 去氢**。

        必须过滤: freesasa 对未知残基 (如 UNK/修饰残基) 会给出半径 -1 而抛错
        ("Radius is <= 0"); 几何特征侧本来就只统计标准氨基酸, 因此过滤后
        两侧原子集一致, SASA 与训练口径仍然可比。
        """
        from Bio.PDB.Residue import Residue as _Res
        st = self.BSt.Structure("sub")
        md = self.BMd.Model(0)
        st.add(md)
        for ch in model:
            if chain_ids is not None and ch.id not in chain_ids:
                continue
            new_ch = self.BCh.Chain(ch.id)
            for res in ch:
                if res.id[0] != " ":                 # 跳过 HETATM/水
                    continue
                if res.get_resname().strip().upper() not in self.aa_ok:
                    continue
                new_res = _Res(res.id, res.get_resname(), res.get_segid())
                for atom in res:
                    el = (atom.element or "").strip().upper()
                    if el in ("H", "D", ""):          # freesasa 默认不含氢
                        continue
                    new_res.add(atom.copy())
                if len(new_res) > 0:
                    new_ch.add(new_res)
            if len(new_ch) > 0:
                md.add(new_ch)
        return st

    def calc(self, model, chain_ids, key, container=None):
        """返回 (总面积, {链: 面积})。

        缓存分两级:
          * _struct_cache: key=链集合 -> freesasa.Structure (构建成本高, 跨行复用);
          * _cache:         key=业务键 -> (面积, 分链面积)。
        container: 含该 model 的 Structure 对象; chain_ids=None 时直接用它 (复杂物)。
        """
        if key in self._cache:
            return self._cache[key]
        ckey = ("ALL",) if chain_ids is None else tuple(sorted(chain_ids))
        st = self._struct_cache.get(ckey)
        if st is None:
            if chain_ids is not None:
                bio = self._subset(model, chain_ids)
            elif container is not None:
                bio = self._subset(model, [c.id for c in model])
            else:
                bio = self._subset(model, [c.id for c in model])
            st = self.fs.structureFromBioPDB(bio, options=self.opts)
            self._struct_cache[ckey] = st
        res = self.fs.calc(st)
        areas = {}
        try:
            aria = res.residueAreas()
            for ch in aria:
                areas[ch] = float(sum(v.total for v in aria[ch].values()))
        except Exception as exc:                                # noqa: BLE001
            logger.debug("residueAreas 失败: %s", exc)
        out = (float(res.totalArea()), areas)
        self._cache[key] = out
        return out

    def monomer(self, model, chain_id):
        """单链 SASA (独立计算, 与旧管线 group_mono(chains={ch}) 等价)。"""
        if chain_id not in self._mono_cache:
            self._mono_cache[chain_id] = self.calc(model, {chain_id}, ("mono", chain_id))[0]
        return self._mono_cache[chain_id]

    def best_key(self, model, chain_ids, key):
        """兼容接口: calc(model, chain_ids, key)[0]"""
        return self.calc(model, chain_ids, key)[0]


def sasa_features_v11(sc, model, structure, lig_chains, rec_chains, pair_contacts=None):
    """V0.3 特征口径的 SASA 计算, 但走零磁盘 IO 路线 (与 02 的数值结果一致)。"""
    pair_contacts = pair_contacts or {}
    out = {}
    total, areas = sc.calc(model, None, ("complex",), container=structure)
    out["sasa_complex"] = total

    model_chains = [c.id for c in model]
    lig_present = [c for c in lig_chains if c in model_chains]
    rec_present = [c for c in rec_chains if c in model_chains]

    def group_area(chs, tag):
        return sc.calc(model, set(chs), ("group", tag, tuple(sorted(chs))))[0]

    try:
        out["sasa_lig_mono"] = group_area(lig_present, "lig")
    except Exception as exc:                                    # noqa: BLE001
        logger.debug("lig 单体 SASA 失败: %s", exc)
        out["sasa_lig_mono"] = float("nan")
    try:
        out["sasa_rec_mono"] = group_area(rec_present, "rec")
    except Exception as exc:                                    # noqa: BLE001
        logger.debug("rec 单体 SASA 失败: %s", exc)
        out["sasa_rec_mono"] = float("nan")

    lig_in = float(sum(areas.get(c, sc.monomer(model, c)) for c in lig_present))
    rec_in = float(sum(areas.get(c, sc.monomer(model, c)) for c in rec_present))
    out["sasa_lig_in_complex"] = lig_in
    out["sasa_rec_in_complex"] = rec_in
    out["dsasa_lig"] = out["sasa_lig_mono"] - lig_in
    out["dsasa_rec"] = out["sasa_rec_mono"] - rec_in
    out["dsasa_total"] = out["dsasa_lig"] + out["dsasa_rec"]
    out["dsasa_frac_lig"] = out["dsasa_lig"] / max(out["sasa_lig_mono"], 1e-6)
    out["dsasa_frac_rec"] = out["dsasa_rec"] / max(out["sasa_rec_mono"], 1e-6)

    # ---- 逐链对 ΔSASA (只对有接触的链对; 上限 MAX_PAIR_SASA) ----
    cand = []
    for (ca, cb), ncnt in pair_contacts.items():
        if ca in lig_present and cb in rec_present:
            cand.append((ca, cb, ncnt))
        elif cb in lig_present and ca in rec_present:
            cand.append((cb, ca, ncnt))
    cand.sort(key=lambda x: (-x[2], x[0], x[1]))
    truncated = len(cand) > MAX_PAIR_SASA
    cand = cand[:MAX_PAIR_SASA]

    dv = []
    for la, rb_, _n in cand:
        try:
            pair_key = ("pair", la, rb_)
            s_pair = sc.calc(model, {la, rb_}, pair_key)[0]
            d = sc.monomer(model, la) + sc.monomer(model, rb_) - s_pair
            if d > 0:
                dv.append(d)
        except Exception as exc:                                # noqa: BLE001
            logger.debug("链对 ΔSASA 失败 %s-%s: %s", la, rb_, exc)

    out["n_chainpairs_with_contact"] = len(cand)
    out["n_chainpairs_sasa_ok"] = len(dv)
    out["pair_sasa_complete"] = int(not truncated)
    if dv:
        out["weakest_chainpair_dsasa"] = float(min(dv))
        out["strongest_chainpair_dsasa"] = float(max(dv))
        out["mean_chainpair_dsasa"] = float(np.mean(dv))
        out["std_chainpair_dsasa"] = float(np.std(dv))
        out["weakest_over_total_dsasa"] = float(min(dv) / max(out["dsasa_total"], 1e-6))
    else:
        for k in ("weakest_chainpair_dsasa", "strongest_chainpair_dsasa",
                  "mean_chainpair_dsasa", "std_chainpair_dsasa",
                  "weakest_over_total_dsasa"):
            out[k] = float("nan")
    return out


def fill_single_pair(feat):
    tot = feat.get("dsasa_total", np.nan)
    feat["weakest_chainpair_dsasa"] = tot
    feat["strongest_chainpair_dsasa"] = tot
    feat["mean_chainpair_dsasa"] = tot
    feat["std_chainpair_dsasa"] = 0.0
    feat["weakest_over_total_dsasa"] = 1.0
    feat["pair_sasa_complete"] = 1
    return feat


# ---------------------------------------------------------------- 检查点 / 落盘
def load_checkpoint(path):
    if not os.path.exists(path):
        return {"done_pdb": [], "failed_pdb": {}, "skipped_pdb": {},
                "started_at": None, "n_rows_written": 0}
    try:
        with open(path, "r", encoding="utf-8") as fh:
            cp = json.load(fh)
        for k, v in (("done_pdb", []), ("failed_pdb", {}), ("skipped_pdb", {})):
            cp.setdefault(k, v)
        return cp
    except Exception as exc:                                    # noqa: BLE001
        logger.error("检查点损坏: %s", exc)
        return {"done_pdb": [], "failed_pdb": {}, "skipped_pdb": {},
                "started_at": None, "n_rows_written": 0, "corrupted": True}


def save_checkpoint(path, cp):
    cp["last_update"] = time.strftime("%Y-%m-%d %H:%M:%S")
    tmp = path + ".tmp"
    try:
        with open(tmp, "w", encoding="utf-8") as fh:
            json.dump(cp, fh, ensure_ascii=False, indent=1)
        os.replace(tmp, path)
    except Exception as exc:                                    # noqa: BLE001
        logger.error("检查点写出失败: %s", exc)


def append_rows(path, rows):
    if not rows:
        return 0
    try:
        df = F02.derive_normalized_features(pd.DataFrame(rows))
    except Exception as exc:                                    # noqa: BLE001
        logger.error("派生特征失败: %s", exc)
        df = pd.DataFrame(rows)
    fresh = not os.path.exists(path) or os.path.getsize(path) < 10
    try:
        if fresh:
            df.to_csv(path, index=False, encoding="utf-8")
            return len(df)
        cols = pd.read_csv(path, nrows=0).columns.tolist()
        missing = [c for c in cols if c not in df.columns]
        extra = [c for c in df.columns if c not in cols]
        if missing or extra:
            logger.warning("列集合变化 (缺 %s / 多 %s); 合并重写", missing, extra)
            old = pd.read_csv(path)
            merged = pd.concat([old, df], ignore_index=True)
            if "interface_key" in merged.columns:
                merged = merged.drop_duplicates(subset=["interface_key"], keep="last")
            merged.to_csv(path, index=False, encoding="utf-8")
            return len(df)
        df[cols].to_csv(path, mode="a", header=False, index=False, encoding="utf-8")
        return len(df)
    except Exception as exc:                                    # noqa: BLE001
        logger.error("落盘失败: %s", exc)
        return 0


# ---------------------------------------------------------------- 主流程
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exit-every", type=int, default=EXIT_EVERY)
    ap.add_argument("--limit", type=int, default=0, help="只跑前 N 个 (验证)")
    args = ap.parse_args()
    exit_every = max(1, args.exit_every)

    t0 = time.time()
    logger.info("=" * 78)
    logger.info("V1.1 补全提取启动 | SASA=structureFromBioPDB(零临时文件) | "
                "输出=%s | 检查点=%s | 落盘间隔=%d", OUT_CSV, CHECKPOINT, exit_every)

    from Bio.PDB import MMCIFParser, PDBParser

    bench = pd.read_csv(BENCH, dtype={"pdb_id": str})
    targets = sorted(bench[bench["is_mvp"] != True]["pdb_id"].unique())   # noqa: E712
    done_csv = set()
    fam_lookup = {}
    if os.path.exists(OUT_CSV):
        try:
            prev = pd.read_csv(OUT_CSV)
            done_csv = set(prev["pdb_id"].astype(str).unique())
            # 家族标签继承: 07/09 已给旧行算好 family_id / split_family, 新行按 PDB 继承,
            # 避免新增行缺列 (否则训练阶段会因 split_family 为空而拒绝)。
            if "family_id" in prev.columns and "split_family" in prev.columns:
                sub = prev.dropna(subset=["family_id"])
                for pid, g in sub.groupby("pdb_id"):
                    fam_lookup[str(pid)] = (g["family_id"].iloc[0],
                                            g["split_family"].iloc[0]
                                            if pd.notna(g["split_family"].iloc[0]) else "")
        except Exception as exc:                                # noqa: BLE001
            logger.warning("读取既有 CSV 失败: %s", exc)
    cp = load_checkpoint(CHECKPOINT)
    done = set(cp["done_pdb"]) | done_csv
    failed = dict(cp["failed_pdb"])
    skipped = dict(cp["skipped_pdb"])
    if cp.get("started_at") is None:
        cp["started_at"] = time.strftime("%Y-%m-%d %H:%M:%S")

    by_pdb = {}
    for rec in bench[bench["is_mvp"] != True].to_dict("records"):         # noqa: E712
        by_pdb.setdefault(rec["pdb_id"], []).append(rec)

    remaining = [p for p in targets if p not in done and p not in skipped]
    if args.limit:
        remaining = remaining[:args.limit]
    logger.info("目标 %d | 已完成 %d (CSV) | 失败记录 %d | 剩余 %d",
                len(targets), len(done_csv), len(failed), len(remaining))
    if not remaining:
        logger.info("无剩余结构; 汇总")
        summarize(targets, cp, t0)
        print("补全完成 (无剩余)")
        return

    parser_pdb, parser_cif = PDBParser(QUIET=True), MMCIFParser(QUIET=True)
    sc = SasaCalc()
    buf, n_ok, n_fail, n_skip = [], 0, 0, 0

    for i, pid in enumerate(remaining, start=1):
        path = resolve_struct_path(pid)
        if path is None:
            skipped[pid] = "no_structure_file"
            n_skip += 1
            if i % exit_every == 0 or i == len(remaining):
                cp["skipped_pdb"] = skipped
                save_checkpoint(CHECKPOINT, cp)
            continue
        recs = by_pdb[pid]
        try:
            parser = parser_cif if path.lower().endswith(".cif") else parser_pdb
            structure = parser.get_structure(pid, path)
            model = structure[0]
            chain_ids = [c.id for c in model]
            ss_map = F02.build_ss_map(path)
        except Exception as exc:                                # noqa: BLE001
            failed[pid] = "parse_error: %s" % str(exc)[:80]
            n_fail += 1
            logger.warning("解析失败(跳过) %s: %s", pid, str(exc)[:70])
            gc.collect()
            continue

        for rec in recs:
            lig = [c.strip() for c in str(rec["ligand_chains"]).split(",") if c.strip()]
            rcp = [c.strip() for c in str(rec["receptor_chains"]).split(",") if c.strip()]
            try:
                lig_p = [c for c in lig if c in chain_ids]
                rcp_p = [c for c in rcp if c in chain_ids]
                if not lig_p or not rcp_p:
                    raise ValueError("链组完全缺失")
                n_missing = len(lig) + len(rcp) - len(lig_p) - len(rcp_p)

                xa, ea, na, ra, rna, ssa = F02.extract_chain_data(model, lig_p, ss_map)
                xb, eb, nb, rb, rnb, ssb = F02.extract_chain_data(model, rcp_p, ss_map)
                if len(xa) == 0 or len(xb) == 0:
                    raise ValueError("链组无有效蛋白原子")

                feat = {"pdb_id": pid, "interface_key": rec["interface_key"],
                        "ligand_chains": rec["ligand_chains"],
                        "receptor_chains": rec["receptor_chains"],
                        "real_dG": float(rec["real_dG"]),
                        "split": rec.get("split", ""),
                        "split_v0.1": rec.get("split_v0.1", rec.get("split", "")),
                        "family_id": fam_lookup.get(pid, ("", ""))[0],
                        "split_family": fam_lookup.get(pid, ("", ""))[1],
                        "subgroup": rec.get("subgroup", ""),
                        "source_dataset": rec.get("source_dataset", ""),
                        "n_missing_chain": int(n_missing),
                        "chain_map_ok": int(n_missing == 0),
                        "self_interface": int(set(lig_p) == set(rcp_p))}
                gf = F02.group_features(xa, ea, na, ra, rna, ssa,
                                        xb, eb, nb, rb, rnb, ssb)
                pc = gf.pop("_pair_contacts", {})
                gf.pop("_pair_hbonds", None)
                feat.update(gf)

                ncp = int(feat.get("n_chainpair_contacts", 1) or 1)
                sc._cache.clear()
                sc._mono_cache.clear()          # 结构与内存缓存跨行复用, 不清 _struct_cache
                if ncp <= 1:
                    sf = sasa_features_v11(sc, model, structure, lig_p, rcp_p,
                                           pair_contacts={})
                    feat.update(sf)
                    fill_single_pair(feat)
                else:
                    feat.update(sasa_features_v11(sc, model, structure, lig_p,
                                                  rcp_p, pair_contacts=pc))

                feat["predicted_dG"] = np.nan
                feat["prodigy_ic"] = np.nan
                feat["feature_version"] = "v1.1-full-freesasa-bioPDB"
                buf.append(feat)
                n_ok += 1
            except Exception as exc:                            # noqa: BLE001
                logger.warning("特征失败(跳过) %s %s|%s: %s", pid,
                               rec["ligand_chains"], rec["receptor_chains"],
                               str(exc)[:80])
                failed["%s|%s|%s" % (pid, rec["ligand_chains"],
                                     rec["receptor_chains"])] = \
                    "feature_error: %s" % str(exc)[:60]
                n_fail += 1

        del structure, model
        done.add(pid)
        # 落盘触发: 结构数达到 exit_every **或** 缓冲行数达到 400 (后者防"100 个结构
        # 恰好都很慢"导致长时间不落盘、崩溃即丢全部进度)
        if i % exit_every == 0 or i == len(remaining) or len(buf) >= 400:
            written = append_rows(OUT_CSV, buf)
            cp["done_pdb"] = sorted(done)
            cp["failed_pdb"] = failed
            cp["skipped_pdb"] = skipped
            cp["n_rows_written"] = int(cp.get("n_rows_written", 0)) + written
            save_checkpoint(CHECKPOINT, cp)
            buf = []
            gc.collect()
            el = time.time() - t0
            logger.info("进度 %d/%d | 成功=%d 失败=%d 跳过=%d | 本批写入=%d | "
                        "内存=%.0fMB | 已用=%.1fmin | ETA=%.1fmin",
                        i, len(remaining), n_ok, n_fail, n_skip, written,
                        get_rss_mb(), el / 60, el / max(i, 1) * (len(remaining) - i) / 60)

    cp["done_pdb"] = sorted(done)
    cp["failed_pdb"] = failed
    cp["skipped_pdb"] = skipped
    save_checkpoint(CHECKPOINT, cp)
    gc.collect()
    summarize(targets, cp, t0)
    print("补全完成，等待下一步指令")


def summarize(targets, cp, t0):
    df = pd.read_csv(OUT_CSV) if os.path.exists(OUT_CSV) else pd.DataFrame()
    print("\n" + "=" * 80)
    print("V1.1 补全提取汇总")
    print("=" * 80)
    print("目标结构      : %d" % len(targets))
    print("已完成结构    : %d" % len(cp.get("done_pdb", [])))
    print("失败结构      : %d" % len(cp.get("failed_pdb", {})))
    print("跳过结构      : %d" % len(cp.get("skipped_pdb", {})))
    print("特征表行数    : %d" % len(df))
    if not df.empty:
        print("覆盖唯一结构  : %d" % df["pdb_id"].nunique())
        print("链路对特征完整: %d / %d"
              % (int((df["pair_sasa_complete"] == 1).sum()), len(df)))
    print("耗时          : %.1f 分钟" % ((time.time() - t0) / 60))
    print("=" * 80)
    logger.info("汇总: 行数=%d 完成=%d 失败=%d 跳过=%d 耗时=%.1fmin",
                len(df), len(cp.get("done_pdb", [])), len(cp.get("failed_pdb", {})),
                len(cp.get("skipped_pdb", {})), (time.time() - t0) / 60)


if __name__ == "__main__":
    main()
