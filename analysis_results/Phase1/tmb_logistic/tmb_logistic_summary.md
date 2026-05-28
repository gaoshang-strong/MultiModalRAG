# TMB Logistic Regression — Summary

**Script:** `scripts/Phase1/tmb_logistic_all_genes.R`  
**Output:** `tmb_logistic_all_genes.csv`  
**Date:** 2026-05-26

## 方法

对 FDU 队列（n=472，有 TMB 记录）和 TCGA-STAD（n=431）中所有共同基因（n=489）分别拟合单变量 logistic 回归：

```
logit P(gene mutated) = β0 + β1 · log2(TMB)
```

- FDU TMB：TMB_report_numeric（log2 变换，下限 0.1）
- TCGA TMB：每样本 coding 突变数（log2 变换，下限 1）
- CI：Wald interval（`confint.default`）
- FDR：Benjamini-Hochberg，各队列独立校正

## 结果摘要

| | FDU | TCGA |
|---|---|---|
| 建模基因数 | 499 | 489 |
| p < 0.05 | 394 | 421 |
| **两队列均显著且方向一致** | **361** | |
| 两队列 FDR < 0.05 且一致 | 350 | |
| TMB 负相关（OR < 1）基因 | 0 | 0 |

## 关键基因

| 基因 | FDU OR | FDU FDR | TCGA OR | TCGA FDR |
|---|---|---|---|---|
| TRRAP | 2.83 | <0.001 | 2.46 | <0.001 |
| POLD1 | 2.61 | <0.001 | 2.03 | <0.001 |
| SOX9  | 1.97 | <0.001 | 1.46 | 0.008 |
| ATM   | 1.67 | <0.001 | 2.02 | <0.001 |
| TP53  | 1.28 | <0.001 | 1.06 | 0.264 |
| CDH1  | 1.15 | 0.162  | 0.85 | 0.085 |

## 核心结论

### 1. TMB 是全基因组范围的潮水效应

489 个基因中，**没有任何一个基因在两个队列中同时呈 TMB 负相关**。高 TMB 导致所有基因突变率普遍升高，是一种非特异性的背景效应，而非对特定基因有选择性富集。这意味着：

- 在未校正 TMB 的情况下，MSI-H 中几乎所有基因都会显示"富集"（见旧版 MSI Fisher test 结果：301/489 基因显著富集于 MSI-H）
- 这种富集的绝大部分是 **passenger effect**，不代表这些基因与 MSI 有功能关联

### 2. TRRAP、SOX9 与 MSI 的关系重新定性

TRRAP（OR≈2.8）和 SOX9（OR≈2.0）在 TMB logistic 中显著正相关，与其他数百个基因的模式一致，**无法区分于普通的超突变背景效应**。之前基于 Fisher test 发现的 TRRAP 和 SOX9 在 MSI-H 中的强富集（OR=58 和 OR=12），本质上是 MSI-H 极高 TMB 的代理信号，而非这两个基因有 MSI 特异性功能。

### 3. TP53 是真正的离群值

TP53 的 TMB-OR 在 FDU 仅为 1.28（远低于 TRRAP 的 2.83），在 TCGA 接近 1（OR=1.06，p=0.25）。这说明 **TP53 突变的选择压力独立于背景突变率**——TP53 作为经典驱动基因，其突变在低 TMB 的 MSS 背景下同样被强烈正向选择，不依赖超突变环境。

### 4. TP53–TRRAP 互斥的重新解释

综合 TMB logistic 结果与 MSI confounding 分析：

- TP53 对 TMB 依赖弱（OR≈1）；TRRAP 对 TMB 依赖强（OR≈2.8）
- MSI-H 样本 TMB 中位数约为 MSS 的 15 倍，TRRAP 在 MSI-H 中突变率 55%，在 MSS 中仅 2%
- 因此 TP53（主导 MSS）与 TRRAP（主导 MSI-H）的表观互斥，**主要是 TMB/MSI 亚型分层的统计产物**，而非功能性合成致死互斥

MSS 队列内的 ME 趋势（OR=0.267, p=0.053）是目前最接近真实功能性互斥的信号，但样本量不足（仅 10 例 TRRAP+ in MSS），尚需更大队列验证。

## 方法论启示

> **当分析体细胞突变的互斥或共突变关系时，MSI/MSS 的分层混杂本质上是 TMB 混杂。正确的做法是将连续 TMB（log2 变换）作为协变量纳入 logistic 模型，而非使用二元 MSI 标签。前者更基础、更连续、覆盖所有超突变原因（MSI、POLE 突变等）。**
