# -*- coding: utf-8 -*-
"""scripts/v8_1/145_keyleakage_exact.py -- 精确核验 key-level leakage（直接执行划分代码）。

不重训：只导入 scripts/v6_2/02_split4_ddg.py 的 load_all() / level_units() / make_folds()，
用与正式运行相同的 seed，重建每个 level 的 fold 归属，然后检验：
    同一 obs_key 的所有行是否落在同一 fold？
若四个 level 全部成立，则 key-level leakage 在本 benchmark 的划分中已被构造性排除。
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
R = Path(".").resolve()
sys.path.insert(0, str(R / "scripts" / "v5"))
sys.path.insert(0, str(R / "scripts" / "v6_1"))

spec = importlib.util.spec_from_file_location("s4", R / "scripts/v6_2/02_split4_ddg.py")
s4 = importlib.util.module_from_spec(spec)
sys.modules["s4"] = s4
spec.loader.exec_module(s4)

D = s4.load_all()
obs = D["obs"]
y = D["y"]
print(f"  载入行数 = {len(y):,}   不同 obs_key = {len(np.unique(obs)):,}")

SEEDS = [42, 43, 44]
LEVELS = ["L1_random", "L2_complex", "L3_seq40", "L4_family"]
report = {}
print()
print("=" * 100)
print("同一 obs_key 是否跨 fold（每个 level × 每个 seed）")
print("=" * 100)
print(f"  {'level':<12}{'seed':>6}{'units':>8}{'keys':>8}{'跨fold key':>12}  判定")
allok = True
for lvl in LEVELS:
    units = s4.level_units(D, lvl)
    for sd in SEEDS:
        fl = s4.make_folds(units, sd, 5)
        n = len(y)
        fold_of = np.full(n, -1, dtype=int)
        for k, (tr, te) in enumerate(fl):
            fold_of[te] = k
        assert (fold_of >= 0).all(), "存在未分配行"
        # 每个 obs_key 的 fold 集合
        import pandas as pd
        df = pd.DataFrame({"obs": obs, "fold": fold_of})
        nf = df.groupby("obs").fold.nunique()
        bad = int((nf > 1).sum())
        allok &= (bad == 0)
        report[f"{lvl}_s{sd}"] = dict(units=int(len(np.unique(units))),
                                      keys=int(len(nf)), keys_crossing_folds=bad)
        print(f"  {lvl:<12}{sd:>6}{len(np.unique(units)):>8}{len(nf):>8}{bad:>12}  "
              f"{'OK' if bad == 0 else 'LEAK'}")

print()
print("=" * 100)
if allok:
    print("  结论：四个严格度层级 × 三个 seed，**没有任何 obs_key 跨 fold**。")
    print("  ⇒ key-level leakage 已被划分构造性排除，无需重训。")
    print("  ⇒ Note 6 目前写「未做键级分组」低估了本 benchmark 的划分。")
else:
    print("  结论：存在跨 fold 的 obs_key，key-level leakage 成立，需要键级重训。")
print("=" * 100)

json.dump(dict(seeds=SEEDS, levels=LEVELS, rows=int(len(y)),
               distinct_keys=int(len(np.unique(obs))),
               all_keys_single_fold=bool(allok), detail=report),
          open(R / "v8_1/KEY_LEAKAGE_EXACT.json", "w", encoding="utf-8"),
          indent=2, ensure_ascii=False)
print()
print("wrote v8_1/KEY_LEAKAGE_EXACT.json")
