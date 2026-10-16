# V8 PHASE 0 provenance

* 本目录及全部 V8 产物**只读复用**上游；**不修改** V1–V7 任何输出。
* 权威 target 命名空间：`canonical_target_labels.csv`，以 **`row_id`** 连接。
* **禁止**以 `block_L.pdb_id` / `mutation_sites.pdb_id` 作为 target 身份。
* 冻结 split：`OUTER_FOLDS_V7_2`。
* 平台：**本地 CPU 优先**；GPU（Kaggle）仅在 §51 全部条件满足时使用。
* `V7.4R2_STATUS = PAUSED_RESOURCE_LIMIT`；V8 期间不继续 V7.4R2。
