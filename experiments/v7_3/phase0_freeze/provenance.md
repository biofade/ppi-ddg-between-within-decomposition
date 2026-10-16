# V7.3 PHASE 0 — Provenance

**时间**：2026-10-15 00:15:31

## 命名空间纪律（§3/§58）

```
AUTHORITATIVE TARGET NAMESPACE = experiments/v7_1_1/phase1_universe/canonical_target_labels.csv
JOIN KEY                       = row_id (== pair_row)
FORBIDDEN as target identity   = block_L*.pdb_id, mutation_sites.pdb_id
```

## 未修改的上游

V1–V7.2R-3 全部只读。新命名空间：`scripts/v7_3`、`features/v7_3`、
`experiments/v7_3`、`logs/v7_3`。

## 外部能量引擎

依 §8/§40：**无合法可用的外部能量引擎** → `EXTERNAL_ENERGY_ENGINE = NOT_RUN`；
V7.3 主实验**不依赖**它。

## SASA

`freesasa`（MIT）可用 → 主用；BioPython ShrakeRupley 作交叉校验。

## 术语纪律（§11）

V7.3 主特征是 **fixed-backbone physics-derived perturbation proxies**，
**不是**精确自由能。禁止称其为 `kcal/mol energy` / `binding energy` /
`physical DeltaDeltaG`。
