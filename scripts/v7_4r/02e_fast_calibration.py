# -*- coding: utf-8 -*-
"""scripts/v7_4r/02e_fast_calibration.py -- 在**烧两小时之前**先做 10 分钟裁定.

# 为什么上一版修复不足够（我的错误）

第 1 次 Colab 结果里有两条**互相印证**的异常，我把它们都归因于「起点冲突」，
但其实它们指向同一个更根本的问题：**最小化结果不可复现**。

证据 A：probe 里 `E_pre = 7.69e11` → 最小化 237 s → `E_post = -62874`。
证据 B：identity（**完全相同的输入**）两次结果：
```
+10989, -1298, +2563, +9560, +5995   kJ/mol
```
**符号都变**。在一个 6884 原子体系里，同一输入、同一代码、同一平台，
得到差 1.2e4 kJ/mol 且符号相反的结果 —— 这不是「起点太差」能解释的，
这是**最小化路径不可复现**。

证据 C：identity 与 mutation 的 RMSD 都是 **~1.5 Å**，
与 shell 大小（4~29）无关 ⇒ 这 1.5 Å 是**体系级漂移地板**，不是位点响应。

## 我的上一版修复**没有**解决证据 B

上一版给 `minimize()` 加了 Stage A，但 Stage A **仍然是从同一份未松弛坐标出发、
用同一次 L-BFGS**。如果最终不可复现，Stage A 也不会让它变可复现。
⇒ **上一版很可能仍然过不了收敛门**，让用户再烧两小时是我的失职。

## 因此：先花 ~10 分钟做**直接裁定**

新增 `08d_minimization_determinism.py`（在 probe 之后、任何 20 变体之前），
只处理**同一个 target 的同一个位点**，用增量成本回答 3 个问题：

| 测试 | 做什么 | 回答 |
|---|---|---|
| **T1 最小化器可复现性** | 同一 `System` + 同一坐标，最小化 2 次 | 最小化本身确定吗？ |
| **T2 建系可复现性** | 重跑 `prepare_system` 得新 System，再最小化 | 建系（PDBFixer/参数化）确定吗？ |
| **T3 收敛性** | 记录 energy vs iteration，看是否到平台 | 2000 步够吗？还是根本没收敛？ |

判据（全部为**项目工程阈值**，非领域标准）：
```
T1  |dE| / |E_post| < 1e-4   且  RMSD < 0.01 Å
T2  |dE| / |E_post| < 1e-3   且  RMSD < 0.05 Å
T3  末两次能量变化 < 1e-3 kJ/mol
```
**任一不过 ⇒ STOP**，并把原始诊断打印出来（不再进入 20 变体）。

这个 cell 的成本：~2 次建系 + ~6 次最小化，与用户已经等过的
「identity 5 变体 ≈ 833 s」同量级甚至更低，**但能直接裁定协议可行性**。
"""

from __future__ import annotations

import ast
import json
import sys
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

R = Path(__file__).resolve().parent.parent.parent
NB = R / "colab" / "v7_4r" / "V7_4R_OpenMM_Relaxation_COLAB.ipynb"
nb = json.loads(NB.read_text(encoding="utf-8"))
cells = nb["cells"]


def find(marker):
    for i, c in enumerate(cells):
        if c["cell_type"] == "code" and marker in "".join(c["source"]):
            return i
    raise KeyError(marker)


DET_CELL = r'''
# ============================================================ 08d Minimization determinism gate (FAST, decisive)
# 第 1 次 Colab smoke 暴露：同一输入的 identity 两次最小化差 1.2e4 kJ/mol 且**符号相反**。
# 这必须在跑 20 变体（乃至 3253）之前裁定清楚。
# 本 cell 只用一个 target 的一个位点，成本 ~2 次建系 + ~6 次最小化。
# 判据全部为**项目工程阈值**，非领域标准。
import time

def _minimize_n(system, positions, n_rounds=1, iters=MIN_MAXITER):
    """在**同一 system**上做 n_rounds 次最小化，返回 (positions, energies, n_atoms)。"""
    integ = mm.VerletIntegrator(0.001)
    ctx = mm.Context(system, integ, PLAT)
    ctx.setPositions(positions)
    es = []
    for _ in range(n_rounds):
        mm.LocalEnergyMinimizer.minimize(ctx, tolerance=MIN_TOL, maxIterations=iters)
        st = ctx.getState(getPositions=True, getEnergy=True)
        es.append(st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
    return st.getPositions(), es, system.getNumParticles()

_pdb = fetch_pdb(str(var.iloc[0].target_id))
_sk = {(str(VAL.iloc[0].chain), str(VAL.iloc[0].resseq))}
print("=== 08d 最小化可复现性裁定 ===")
print("  target=%s  位点=%s" % (var.iloc[0].target_id, _sk), flush=True)
print("  platform=%s  MIN_TOL=%s  MIN_MAXITER=%s" % (PLAT_NAME, MIN_TOL, MIN_MAXITER), flush=True)

# ---------- 建系一次 ----------
_t0 = time.time()
_f1, _e1, _s1, _t1 = prepare_system(_pdb)
print("  建系 #1: %s  (%.1fs)  atoms=%s"
      % ("OK" if _e1 is None else _e1, time.time() - _t0,
         _s1.getNumParticles() if _s1 is not None else None), flush=True)
if _e1 is not None:
    raise RuntimeError(f"DETERMINISM_GATE: prepare#1 failed {_e1}")
_n = _s1.getNumParticles()
_pos_start = _f1.positions

# ---------- T3：收敛性（看能量 vs 迭代） ----------
print()
print("--- T3 收敛性（同一 system，逐轮收紧）---", flush=True)
_pos, _es, _ = _minimize_n(_s1, _pos_start, n_rounds=4)
for _i, _e in enumerate(_es):
    print("     round %d: E = %.6e kJ/mol  (%.6e /atom)" % (_i + 1, _e, _e / _n), flush=True)
T3_delta = abs(_es[-1] - _es[-2]) if len(_es) >= 2 else float("nan")
T3_PASS = bool(abs(_es[-1]) / _n < 100.0 and T3_delta < CONVERGENCE_DTOL * 100)
print("     末两次差 = %.6e kJ/mol | per-atom = %.6e | T3_PASS = %s"
      % (T3_delta, abs(_es[-1]) / _n, T3_PASS), flush=True)

# ---------- T1：同一 system 再最小化一次（最小化器可复现性） ----------
print()
print("--- T1 最小化器可复现性（同一 system 同一坐标，独立 Context）---", flush=True)
_t1a = time.time()
_posA, eA, _ = _minimize_n(_s1, _pos_start, n_rounds=1)
_posB, eB, _ = _minimize_n(_s1, _pos_start, n_rounds=1)
_aA = np.array(_posA.value_in_unit(unit.angstrom))
_aB = np.array(_posB.value_in_unit(unit.angstrom))
T1_rmsd = float(np.sqrt(np.mean(np.sum((_aA - _aB) ** 2, axis=1))))
T1_rel = abs(eB[0] - eA[0]) / max(1e-12, abs(eB[0]))
T1_PASS = bool(T1_rel < 1e-4 and T1_rmsd < 0.01)
print("     E_A = %.6e   E_B = %.6e" % (eA[0], eB[0]), flush=True)
print("     |dE| = %.6e  rel = %.3e | RMSD = %.6f A" % (abs(eB[0] - eA[0]), T1_rel, T1_rmsd), flush=True)
print("     T1_PASS = %s   (%.1fs)" % (T1_PASS, time.time() - _t1a), flush=True)

# ---------- T2：重新建系（建系可复现性） ----------
print()
print("--- T2 建系可复现性（重跑 prepare_system，新 System）---", flush=True)
_t2a = time.time()
_f2, _e2, _s2, _t2 = prepare_system(_pdb)
if _e2 is not None:
    print("     建系 #2 失败:", _e2, flush=True)
    T2_PASS = False; T2_rel = float("nan"); T2_rmsd = float("nan")
else:
    _posC, eC, _ = _minimize_n(_s2, _f2.positions, n_rounds=1)
    _aC = np.array(_posC.value_in_unit(unit.angstrom))
    T2_rmsd = float(np.sqrt(np.mean(np.sum((_aA - _aC) ** 2, axis=1))))
    T2_rel = abs(eC[0] - eA[0]) / max(1e-12, abs(eC[0]))
    T2_PASS = bool(T2_rel < 1e-3 and T2_rmsd < 0.05)
    print("     E_A(system1) = %.6e   E_C(system2) = %.6e" % (eA[0], eC[0]), flush=True)
    print("     |dE| = %.6e  rel = %.3e | RMSD = %.6f A" % (abs(eC[0] - eA[0]), T2_rel, T2_rmsd), flush=True)
    print("     T2_PASS = %s   (%.1fs)" % (T2_PASS, time.time() - _t2a), flush=True)

# ---------- 裁定 ----------
DETERMINISM_GATE_PASS = bool(T1_PASS and T2_PASS and T3_PASS)
print()
print("=== MINIMIZATION DETERMINISM GATE ===")
print("   T1 最小化器可复现 : %s" % T1_PASS)
print("   T2 建系可复现     : %s" % T2_PASS)
print("   T3 收敛到合理量级 : %s" % T3_PASS)
print("   DETERMINISM_GATE_PASS =", DETERMINISM_GATE_PASS)
(WORK / "determinism_gate.json").write_text(json.dumps(
    {"T1_minimizer_reproducible": T1_PASS, "T1_rel": T1_rel, "T1_rmsd_A": T1_rmsd,
     "T2_build_reproducible": T2_PASS, "T2_rel": T2_rel, "T2_rmsd_A": T2_rmsd,
     "T3_converged": T3_PASS, "T3_last_delta": T3_delta,
     "E_final_per_atom": abs(_es[-1]) / _n, "n_atoms": _n,
     "platform": PLAT_NAME, "gate_pass": DETERMINISM_GATE_PASS,
     "threshold_note": "prespecified project engineering thresholds, NOT field-standard"},
    indent=2))
if not DETERMINISM_GATE_PASS:
    print()
    print("STOP: 最小化/建系不可复现或未收敛 -> 噪声门必然失败，不要进入 20 变体。")
    print("  上面 T1/T2/T3 的原始数字即为诊断依据。")
    print("  已知可能的成因（待裁定后择一修）：")
    print("    (a) 平台浮点非确定性 -> 改 OpenCL 或 CPU 复现验证")
    print("    (b) 未收敛（迭代不足）-> 提高 iters 或改分阶段阈值")
    print("    (c) 起点严重冲突     -> 需要蛋白质级预松弛")
    raise RuntimeError("MINIMIZATION_DETERMINISM_FAIL —— 依裁定 STOP")
print()
print(">>> 可复现性裁定通过。可以继续进入 20 变体 smoke。")
'''

i09 = find("# ============================================================ 09 Mutation builder")
cells.insert(i09, {"cell_type": "code", "execution_count": None, "metadata": {},
                   "outputs": [], "source": DET_CELL.strip("\n").splitlines(True)})
print("已插入 08d（最小化可复现性裁定）")

NB.write_text(json.dumps(nb, indent=1, ensure_ascii=False), encoding="utf-8")

nb2 = json.loads(NB.read_text(encoding="utf-8"))
code2 = [c for c in nb2["cells"] if c["cell_type"] == "code"]
bad = []
for i, c in enumerate(code2):
    try:
        ast.parse("".join(c["source"]))
    except SyntaxError as e:
        bad.append((i, e.lineno, e.msg))
t = "\n".join("".join(c["source"]) for c in code2)
print(f"cells={len(nb2['cells'])} code={len(code2)} 语法错误={len(bad)}")
for b in bad:
    print("   ", b)
print()
for k, v in (("T1 最小化器可复现", "T1_PASS" in t),
             ("T2 建系可复现", "T2_PASS" in t),
             ("T3 收敛性", "T3_PASS" in t),
             ("总门 DETERMINISM_GATE_PASS", "DETERMINISM_GATE_PASS" in t),
             ("失败即 STOP", "MINIMIZATION_DETERMINISM_FAIL" in t),
             ("写 determinism_gate.json", "determinism_gate.json" in t),
             ("§6 协议数值未变",
              "RADIUS_MOVABLE = 8.0" in t and "RESTRAINT_OUTSIDE = 100.0" in t
              and "PH = 7.0" in t)):
    print(f"  {'✅' if v else '❌'} {k}")
