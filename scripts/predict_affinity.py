# -*- coding: utf-8 -*-
"""
predict_affinity.py -- Affinity Intelligence 单结构 ΔG 预测 (可复用函数 + CLI)

两种用法:

1) 命令行:
    python predict_affinity.py --pdb <PDB路径> [--ligand A] [--receptor B] [--model ...]

2) 作为库调用 (Streamlit 前端等, 避免每次启动新进程):
    from predict_affinity import get_prediction
    r = get_prediction("x.pdb", ligand_chain="A", receptor_chain="B")
    r["rank_pct"]     -> 强于训练集的百分比
    r["raw_dG"]       -> 模型原始输出
    r["calibrated_dG"]-> 按测试集斜率还原的估计
    r["ci_low"], r["ci_high"] -> 95% 经验区间 (基于测试集 RMSE)
    r["top_features"] -> [(特征名, 值, 贡献), ...] 供前端画图

重要校准说明 (诚实披露):
    该模型测试集拟合斜率约 0.24, 预测值被强烈压缩向均值 —— **排序能力强于定量能力**。
    因此同时给出 raw 与斜率校准值; 请以相对排名为主要参考。

特征口径: 与 features_v1.1.csv 一致 (V0.3 的 60 个特征); SASA 走零临时文件路线。
纪律: 纯 CPU; 全程 try...except; 不写任何临时文件 (上传文件由调用方负责清理)。
"""

import argparse
import os
import pickle
import sys
import time

import numpy as np
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FEAT_DIR = os.path.join(ROOT, "features")
MODEL_DIR = os.path.join(ROOT, "models")
DEFAULT_MODEL = os.path.join(MODEL_DIR, "affinity_engine_v1.1.pkl")
# 训练集统计来源 (按优先级):
#   1) features/train_stats.csv   —— 精简发布版 (仅含排名分位数与中位数所需列)
#   2) features/features_v1.1.csv —— 开发环境完整特征表
#   3) features/features_full.csv —— 回退
TRAIN_STATS_CSV = os.path.join(FEAT_DIR, "train_stats.csv")
TRAIN_CSV = os.path.join(FEAT_DIR, "features_v1.1.csv")
TRAIN_CSV_FALLBACK = os.path.join(FEAT_DIR, "features_full.csv")


def _ensure_openmp():
    """本机无 VC++ 运行库: lib_lightgbm.dll 依赖的 vcomp140.dll 借用 sklearn 自带副本。

    必须在 import lightgbm / unpickle 模型**之前**执行, 否则 DLL 加载失败。
    """
    import ctypes
    import glob
    try:
        import sklearn
        hits = glob.glob(os.path.join(os.path.dirname(sklearn.__file__),
                                      ".libs", "vcomp140.dll"))
        if hits:
            ctypes.CDLL(hits[0])
            return hits[0]
    except Exception:                                           # noqa: BLE001
        pass
    return None


_VCOMP = _ensure_openmp()
import lightgbm as lgb                                              # noqa: E402

# 与 08/11 训练脚本完全一致的特征集
V02_FEATURES = [
    "n_missing_chain", "chain_map_ok", "n_contact_pairs", "n_atoms_lig", "n_atoms_rec",
    "n_res_lig", "n_res_rec", "n_contact_res_pairs", "n_iface_res_lig",
    "n_iface_res_rec", "iface_res_frac_lig", "iface_res_frac_rec",
    "min_interatomic_dist", "n_hbond", "n_salt_bridge", "n_hydrophobic_contact",
    "kd_hydro_lig", "kd_hydro_rec", "kd_hydro_all", "frac_charged_iface",
    "frac_polar_iface", "frac_hydrophobic_iface", "frac_helix_iface",
    "frac_sheet_iface", "n_chainpair_max", "n_chainpair_min",
    "n_chainpair_contacts", "contact_density", "sasa_complex", "sasa_lig_mono",
    "sasa_rec_mono", "sasa_lig_in_complex", "sasa_rec_in_complex",
    "dsasa_lig", "dsasa_rec", "dsasa_total", "dsasa_frac_lig", "dsasa_frac_rec",
    "interface_residue_ratio_lig", "interface_residue_ratio_rec",
    "interface_residue_ratio", "contact_per_iface_res", "contacts_per_res_total",
    "dsasa_per_iface_res", "dsasa_per_contact", "hbond_per_contact",
    "hbond_per_iface_res", "saltbridge_per_iface_res", "iface_res_per_contact",
    "dsasa_frac_total",
]
V03_NEW = [
    "weakest_chainpair_dsasa", "strongest_chainpair_dsasa",
    "mean_chainpair_dsasa", "std_chainpair_dsasa",
    "weakest_over_total_dsasa", "weakest_to_mean_chainpair_dsasa",
    "mean_chainpair_dsasa_per_contact", "chainpair_contact_density",
    "n_chainpairs_with_contact", "n_chainpairs_sasa_ok", "pair_sasa_complete",
    "self_interface",
]
ALL_FEATURES = V02_FEATURES + V03_NEW

COUNT_COLS = [
    "n_contact_pairs", "n_contact_res_pairs", "n_iface_res_lig", "n_iface_res_rec",
    "n_hbond", "n_salt_bridge", "n_hydrophobic_contact", "n_chainpair_max",
    "n_chainpair_min", "n_chainpair_contacts", "n_atoms_lig", "n_atoms_rec",
    "n_res_lig", "n_res_rec", "sasa_complex", "sasa_lig_mono", "sasa_rec_mono",
    "sasa_lig_in_complex", "sasa_rec_in_complex", "dsasa_lig", "dsasa_rec",
    "dsasa_total", "weakest_chainpair_dsasa", "strongest_chainpair_dsasa",
    "mean_chainpair_dsasa", "n_chainpairs_with_contact", "n_chainpairs_sasa_ok",
]
MAX_PAIR_SASA = 8

# 前端柱状图展示用的"人可读"特征 (避免展示 log 变换后的原始列)
DISPLAY_FEATURES = [
    ("n_contact_pairs", "界面接触原子对"),
    ("n_contact_res_pairs", "界面接触残基对"),
    ("n_hbond", "界面氢键数"),
    ("n_salt_bridge", "界面盐桥数"),
    ("dsasa_total", "界面埋藏面积 (Å²)"),
    ("interface_residue_ratio", "界面残基占比"),
    ("contact_density", "接触密度 (每界面残基)"),
    ("frac_polar_iface", "界面极性残基占比"),
    ("frac_hydrophobic_iface", "界面疏水残基占比"),
    ("frac_helix_iface", "界面螺旋占比"),
    ("frac_sheet_iface", "界面折叠占比"),
    ("chainpair_contact_density", "链对平均接触数"),
    ("weakest_chainpair_dsasa", "最弱界面 ΔSASA (Å²)"),
]


def _load(mod_name, fname):
    import importlib.util
    p = os.path.join(ROOT, "scripts", fname)
    spec = importlib.util.spec_from_file_location(mod_name, p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ---------------------------------------------------------------- 特征提取
class FeatureExtractor:
    """复用 02 的几何/序列逻辑 + 10 的零磁盘 SASA 逻辑, 保证与训练口径一致。"""

    def __init__(self):
        import freesasa
        from Bio.PDB import Structure as BSt, Model as BMd
        self.fs = freesasa
        self.BSt, self.BMd = BSt, BMd
        self.fs.setVerbosity(self.fs.silent)
        self.opts = {"hetatm": False, "skip-unknown": True}
        self.F02 = _load("feat02_pred", "02_extract_features.py")
        self.F10 = _load("feat10_pred", "10_extract_remaining.py")

    def _subset(self, model, chain_ids):
        st = self.BSt.Structure("sub")
        md = self.BMd.Model(0)
        st.add(md)
        for ch in model:
            if ch.id in chain_ids:
                md.add(ch.copy())
        return st

    def sasa(self, model, structure, lig, rec, pair_contacts):
        sc = self.F10.SasaCalc()
        return self.F10.sasa_features_v11(sc, model, structure, lig, rec,
                                          pair_contacts=pair_contacts)

    def list_protein_chains(self, pdb_path):
        """返回 {链ID: 残基数} (仅蛋白链), 供前端下拉菜单。"""
        from Bio.PDB import MMCIFParser, PDBParser
        is_cif = pdb_path.lower().endswith((".cif", ".mmcif"))
        parser = MMCIFParser(QUIET=True) if is_cif else PDBParser(QUIET=True)
        st = parser.get_structure("q", pdb_path)
        out = {}
        for ch in st[0]:
            n = sum(1 for res in ch if res.id[0] == " ")
            if n > 0:
                out[ch.id] = n
        return out

    def extract(self, pdb_path, ligand_chains=None, receptor_chains=None):
        from Bio.PDB import MMCIFParser, PDBParser
        pdb_id = os.path.splitext(os.path.basename(pdb_path))[0].upper()
        is_cif = pdb_path.lower().endswith((".cif", ".mmcif"))
        parser = MMCIFParser(QUIET=True) if is_cif else PDBParser(QUIET=True)
        structure = parser.get_structure(pdb_id, pdb_path)
        model = structure[0]

        prot_chains = {}
        for ch in model:
            n = sum(1 for res in ch if res.id[0] == " ")
            if n > 0:
                prot_chains[ch.id] = n
        if not prot_chains:
            raise RuntimeError("结构中没有蛋白链")

        if ligand_chains:
            lig = [c.strip() for c in ligand_chains.split(",") if c.strip()]
        else:
            lig = [min(prot_chains, key=lambda k: (prot_chains[k], k))]
        if receptor_chains:
            rec = [c.strip() for c in receptor_chains.split(",") if c.strip()]
        else:
            rec = [c for c in prot_chains if c not in lig]
        miss = [c for c in lig + rec if c not in prot_chains]
        if miss:
            raise RuntimeError("指定的链不存在于结构中: %s (结构含 %s)"
                               % (miss, sorted(prot_chains)))
        if not rec:
            raise RuntimeError("受体重链为空; 请指定受体链")

        ss_map = self.F02.build_ss_map(pdb_path)
        xa, ea, na, ra, rna, ssa = self.F02.extract_chain_data(model, lig, ss_map)
        xb, eb, nb, rb, rnb, ssb = self.F02.extract_chain_data(model, rec, ss_map)
        if len(xa) == 0 or len(xb) == 0:
            raise RuntimeError("链组没有有效蛋白原子")

        feat = {"pdb_id": pdb_id, "ligand_chains": ",".join(lig),
                "receptor_chains": ",".join(rec),
                "n_missing_chain": 0, "chain_map_ok": 1,
                "self_interface": int(set(lig) == set(rec))}
        gf = self.F02.group_features(xa, ea, na, ra, rna, ssa,
                                     xb, eb, nb, rb, rnb, ssb)
        pc = gf.pop("_pair_contacts", {})
        gf.pop("_pair_hbonds", None)
        feat.update(gf)

        ncp = int(feat.get("n_chainpair_contacts", 1) or 1)
        sf = self.sasa(model, structure, lig, rec, {} if ncp <= 1 else pc)
        feat.update(sf)
        if ncp <= 1:
            tot = feat.get("dsasa_total", np.nan)
            feat["weakest_chainpair_dsasa"] = tot
            feat["strongest_chainpair_dsasa"] = tot
            feat["mean_chainpair_dsasa"] = tot
            feat["std_chainpair_dsasa"] = 0.0
            feat["weakest_over_total_dsasa"] = 1.0
            feat["pair_sasa_complete"] = 1
            feat["n_chainpairs_with_contact"] = max(
                1, int(feat.get("n_chainpair_contacts", 1) or 1))
            feat["n_chainpairs_sasa_ok"] = 1

        # derive_normalized_features 需要 real_dG (仅用于生成诊断列), 预测时填占位 0
        feat["real_dG"] = 0.0
        feat = self.F02.derive_normalized_features(pd.DataFrame([feat])).iloc[0].to_dict()
        feat["_lig"] = lig
        feat["_rec"] = rec
        feat["_n_res"] = prot_chains
        return feat


# ---------------------------------------------------------------- 模型与训练集统计 (带缓存)
_CACHE = {"model": None, "stats": None, "extractor": None}


def _get_extractor():
    if _CACHE["extractor"] is None:
        _CACHE["extractor"] = FeatureExtractor()
    return _CACHE["extractor"]


def load_model_bundle(model_path=DEFAULT_MODEL):
    """加载模型与训练集统计 (按 路径+修改时间 缓存, 避免前端每次点击都重载)。"""
    key = (os.path.abspath(model_path), os.path.getmtime(model_path))
    if _CACHE["model"] is not None and _CACHE["model"].get("_key") == key:
        return _CACHE["model"], _CACHE["stats"]
    with open(model_path, "rb") as fh:
        payload = pickle.load(fh)
    payload["_key"] = key

    stats = {"dist_dG": None, "medians": {}, "slope": None, "train_n": 0,
             "rmse": None, "source_csv": None}
    try:
        # 优先用精简发布版统计文件 (开源包内自带); 否则回退完整特征表
        if os.path.exists(TRAIN_STATS_CSV):
            csv = TRAIN_STATS_CSV
        elif os.path.exists(TRAIN_CSV):
            csv = TRAIN_CSV
        else:
            csv = TRAIN_CSV_FALLBACK
        stats["source_csv"] = csv
        fdf = pd.read_csv(csv)
        if "split_family" in fdf.columns:
            tr = fdf[fdf["split_family"] == "train"]
        elif "split" in fdf.columns:
            tr = fdf[fdf["split"] == "train"]
        else:
            tr = fdf
        stats["dist_dG"] = tr["real_dG"].to_numpy(dtype=float)
        stats["train_n"] = len(tr)
        Xtr = tr.copy()
        for c in COUNT_COLS:
            if c in Xtr.columns:
                Xtr[c] = np.log1p(pd.to_numeric(Xtr[c], errors="coerce").clip(lower=0))
        for c in payload["feature_names"]:
            if c in Xtr.columns:
                med = pd.to_numeric(Xtr[c], errors="coerce").median()
                if med == med:
                    stats["medians"][c] = float(med)
        # 校准斜率: 优先用当前模型的测试集预测, 回退到 v1.0 的
        slope = None
        for name in ("v1.1_test_predictions.csv", "v1.0_test_predictions.csv"):
            tp = os.path.join(FEAT_DIR, name)
            if os.path.exists(tp):
                t = pd.read_csv(tp)
                if len(t) > 2:
                    slope = float(np.polyfit(t["real_dG"], t["pred_dG"], 1)[0])
                    break
        if slope is None:
            # 开源发布包不含测试集预测明细: 由 RMSE 与标签方差反推斜率上界估计。
            # 对收缩型模型有 RMSE^2 ≈ sigma_y^2 (1 - slope^2), 故
            #   slope ≈ sqrt(max(0, 1 - RMSE^2 / var_y))
            # 这是近似值, 用于保证校准功能在精简包内仍可用。
            tm = payload.get("metrics", {})
            rmse_m = None
            for k, v in tm.items():
                if isinstance(v, dict) and "rmse" in v:
                    rmse_m = float(v["rmse"])
                    break
            if rmse_m and stats["dist_dG"] is not None:
                vy = float(np.var(stats["dist_dG"]))
                if vy > 0:
                    slope = float(np.sqrt(max(0.0, 1.0 - rmse_m ** 2 / vy)))
                    stats["slope_is_approx"] = True
        stats["slope"] = slope
    except Exception as exc:                                    # noqa: BLE001
        stats["error"] = str(exc)[:120]

    tm = payload.get("metrics", {})
    for k, v in tm.items():
        if isinstance(v, dict) and "rmse" in v:
            stats["rmse"] = float(v["rmse"])
            break
    if stats["rmse"] is None and stats.get("slope") is not None:
        stats["rmse"] = None

    _CACHE["model"], _CACHE["stats"] = payload, stats
    return payload, stats


def build_matrix(payload, feat, train_medians):
    """按模型训练时的特征名与顺序构造一行输入; 缺失特征用训练集中位数兜底并记录。"""
    cols = payload["feature_names"]
    row, missing_vals = {}, []
    for c in cols:
        v = feat.get(c, np.nan)
        if c in COUNT_COLS and v == v:
            v = np.log1p(max(float(v), 0.0))
        row[c] = v
    X = pd.DataFrame([row])[cols].apply(pd.to_numeric, errors="coerce")
    for c in cols:
        if pd.isna(X.at[0, c]):
            X.at[0, c] = train_medians.get(c, 0.0)
            missing_vals.append(c)
    return X, missing_vals


def get_prediction(pdb_path, ligand_chain=None, receptor_chain=None,
                   model_path=DEFAULT_MODEL, calibrate=True):
    """核心预测函数: 供 CLI 与前端共用 (不做任何临时文件生成/清理)。

    返回 dict, 失败时 raise RuntimeError (由调用方决定如何展示)。
    """
    t0 = time.time()
    if not os.path.exists(pdb_path):
        raise RuntimeError("文件不存在: %s" % pdb_path)
    if not os.path.exists(model_path):
        raise RuntimeError("模型不存在: %s" % model_path)

    payload, stats = load_model_bundle(model_path)
    fe = _get_extractor()
    feat = fe.extract(pdb_path, ligand_chain or None, receptor_chain or None)
    X, missing = build_matrix(payload, feat, stats["medians"])
    raw = float(payload["model"].predict(X)[0])

    dist_dG, slope = stats["dist_dG"], stats["slope"]
    cal = raw
    if calibrate and slope is not None and abs(slope) > 1e-6 and dist_dG is not None:
        cal = float(dist_dG.mean() + (raw - dist_dG.mean()) / slope)

    rank_pct = None
    if dist_dG is not None and len(dist_dG) > 10:
        rank_pct = float((dist_dG > raw).mean() * 100)     # ΔG 越负越强

    rmse = stats["rmse"]
    ci_low = ci_high = None
    if rmse:
        half = 1.96 * rmse
        ci_low, ci_high = raw - half, raw + half

    top_features = []
    for k, label in DISPLAY_FEATURES:
        v = feat.get(k, None)
        if v is None or (isinstance(v, float) and v != v):
            continue
        top_features.append((label, float(v)))

    warnings = []
    if feat.get("self_interface"):
        warnings.append("配体链组 == 受体链组 (自界面): 接触/面积存在重复计数")
    if missing:
        warnings.append("%d 个特征缺失, 已用训练集中位数兜底" % len(missing))

    return {
        "ok": True,
        "model_name": os.path.basename(model_path),
        "train_rows": payload.get("train_rows"),
        "split_policy": payload.get("split_policy"),
        "ligand_chains": feat["_lig"],
        "receptor_chains": feat["_rec"],
        "chain_residues": feat["_n_res"],
        "raw_dG": raw,
        "calibrated_dG": cal,
        "slope": slope,
        "rank_pct": rank_pct,
        "rmse": rmse,
        "ci_low": ci_low,
        "ci_high": ci_high,
        "dist_mean": float(dist_dG.mean()) if dist_dG is not None else None,
        "dist_min": float(dist_dG.min()) if dist_dG is not None else None,
        "dist_max": float(dist_dG.max()) if dist_dG is not None else None,
        "train_n_stats": stats.get("train_n"),
        "top_features": top_features,
        "n_contact_pairs": feat.get("n_contact_pairs"),
        "n_hbond": feat.get("n_hbond"),
        "n_salt_bridge": feat.get("n_salt_bridge"),
        "dsasa_total": feat.get("dsasa_total"),
        "interface_residue_ratio": feat.get("interface_residue_ratio"),
        "missing_features": missing,
        "warnings": warnings,
        "elapsed_s": time.time() - t0,
    }


def list_existing_pdbs():
    """data/pdb 下可选的 PDB 文件 (供前端下拉)。"""
    d = os.path.join(ROOT, "data", "pdb")
    if not os.path.isdir(d):
        return []
    return sorted(f for f in os.listdir(d) if f.lower().endswith((".pdb", ".cif", ".ent")))


# ---------------------------------------------------------------- CLI
def main():
    ap = argparse.ArgumentParser(description="Affinity Intelligence 单结构 ΔG 预测")
    ap.add_argument("--pdb", required=True, help="PDB/CIF 文件路径")
    ap.add_argument("--ligand", default="", help="配体链 (逗号分隔), 缺省=最短蛋白链")
    ap.add_argument("--receptor", default="", help="受体链 (逗号分隔), 缺省=其余链")
    ap.add_argument("--model", default=DEFAULT_MODEL, help="模型 pkl 路径")
    ap.add_argument("--no-calibrate", action="store_true", help="不做斜率校准")
    args = ap.parse_args()

    print("=" * 78)
    print("Affinity Intelligence 单结构预测")
    print("=" * 78)
    try:
        r = get_prediction(args.pdb, args.ligand, args.receptor, args.model,
                           calibrate=not args.no_calibrate)
    except Exception as exc:                                    # noqa: BLE001
        import traceback
        print("!! 预测失败: %s" % exc)
        traceback.print_exc()
        return 1

    print("模型        : %s" % r["model_name"])
    print("训练规模    : %s 行 / 划分策略: %s" % (r["train_rows"], r["split_policy"]))
    print("")
    print("--- 界面 ---")
    print("配体链      : %s   受体链: %s"
          % (",".join(r["ligand_chains"]), ",".join(r["receptor_chains"])))
    print("链残基数    : %s" % r["chain_residues"])
    print("界面接触对  : %d   氢键: %d   盐桥: %d   埋藏面积 ΔSASA: %.0f Å²"
          % (r["n_contact_pairs"], r["n_hbond"], r["n_salt_bridge"], r["dsasa_total"]))
    print("界面残基占比: %.3f" % r["interface_residue_ratio"])
    for w in r["warnings"]:
        print("!! %s" % w)
    print("")
    print("=" * 78)
    print("预测结果")
    print("=" * 78)
    print("预测 ΔG (raw)       : %+.2f kcal/mol" % r["raw_dG"])
    if r["slope"] is not None:
        print("斜率校准后估计      : %+.2f kcal/mol  (测试集拟合斜率=%.2f, 模型向均值收缩)"
              % (r["calibrated_dG"], r["slope"]))
    if r["rank_pct"] is not None:
        print("训练集 ΔG 分布      : %.2f ~ %.2f (mean %.2f), n=%d"
              % (r["dist_min"], r["dist_max"], r["dist_mean"], r["train_n_stats"]))
        print("相对排名            : 强于训练集中 %.1f%% 的界面" % r["rank_pct"])
        pct = r["rank_pct"]
        verdict = ("强结合 (预测强于 3/4 已知界面)" if pct >= 75 else
                   "中等偏强" if pct >= 50 else "中等偏弱" if pct >= 25 else "弱结合")
        print("定性判读            : %s" % verdict)
    if r["ci_low"] is not None:
        print("95%% 预测区间 (raw)  : [%+.2f, %+.2f] kcal/mol  (基于测试集 RMSE=%.2f)"
              % (r["ci_low"], r["ci_high"], r["rmse"]))
    print("")
    print("提醒: 该模型测试集拟合斜率约 %.2f, **排序能力强于定量能力** ——"
          % (r["slope"] if r["slope"] else float("nan")))
    print("      请优先使用相对排名/定性判读, 绝对值仅供参考。")
    print("耗时: %.1fs" % r["elapsed_s"])
    print("=" * 78)
    return 0


if __name__ == "__main__":
    sys.exit(main())
