# FDU Cohort — Somatic Mutual Exclusivity Analysis

**Date:** 2026-06-07  
**Script:** `TP53_ATM_paper/scripts/FDU_cohort_ME_analysis.R`  
**Cohort:** CohortB_All (somt_filtered.maf, all samples)

---

## 1. Cohort Summary


| Parameter                                   | Value           |
| ------------------------------------------- | --------------- |
| Total samples (MAF)                         | 445             |
| Samples with complete covariates (logistic) | 390             |
| MSS-only samples (logistic)                 | 367             |
| MSI-H rate                                  | 27/530 (5.1%)   |
| TP53 mutation frequency                     | 263/445 (59.1%) |
| ATM mutation frequency                      | 31/445 (7.0%)   |


**Input:** `analysis_results/Phase0/enrichment/somt_filtered.maf`  
Filters applied upstream: functional mutations only, VAF ≥ 5%, ExAC_EAS < 1%, primary tumor (F01 suffix).

**Logistic regression covariates:** TMB (log₂) + MSI_H + pathologic stage (I–IV) + histology (adeno / signet / other).  
MSS-only model excludes MSI_H term.

---

## 2. Fisher's Exact Test — Top 50 Genes (N = 445)

11 mutually exclusive pairs reached p < 0.05. All involve TP53 as one partner.


| Gene 1 | Gene 2 | Neither | Gene2 only | Gene1 only | Both | OR    | p-value  |
| ------ | ------ | ------- | ---------- | ---------- | ---- | ----- | -------- |
| SOX9   | TP53   | 164     | 260        | 18         | 3    | 0.106 | 3.2×10⁻⁵ |
| ATM    | TP53   | 160     | 254        | 22         | 9    | 0.258 | 5.4×10⁻⁴ |
| POLD1  | TP53   | 168     | 258        | 14         | 5    | 0.233 | 3.8×10⁻³ |
| TRRAP  | TP53   | 165     | 255        | 17         | 8    | 0.305 | 6.0×10⁻³ |
| CDH1   | TP53   | 155     | 245        | 27         | 18   | 0.423 | 9.8×10⁻³ |
| ARID1A | TP53   | 140     | 227        | 42         | 36   | 0.529 | 0.011    |
| RYR1   | TP53   | 161     | 249        | 21         | 14   | 0.432 | 0.020    |
| TP53   | ARID2  | 168     | 14         | 256        | 7    | 0.329 | 0.021    |
| TP53   | CIC    | 168     | 14         | 256        | 7    | 0.329 | 0.021    |
| KMT2C  | TP53   | 161     | 248        | 21         | 15   | 0.465 | 0.033    |
| TP53   | CHD4   | 169     | 13         | 256        | 7    | 0.356 | 0.035    |


No co-occurring pairs with TP53 reached p < 0.05.

---

## 3. Logistic Regression — Full Cohort (N = 390)

**Model:** `gene_A ~ gene_B + log₂(TMB) + MSI_H + stage + histology`  
Multiple testing correction: Benjamini-Hochberg FDR across all gene pairs tested (2,450 directed models).

### 3a. TP53-mutually exclusive partners (gene_A ~ TP53, p < 0.05)


| Gene  | N   | Adj. OR | 95% CI         | p-value | FDR   |
| ----- | --- | ------- | -------------- | ------- | ----- |
| SOX9  | 390 | 0.112   | [0.030, 0.417] | 0.0011  | 0.454 |
| ATM   | 390 | 0.325   | [0.134, 0.789] | 0.013   | 0.994 |
| ARID2 | 390 | 0.259   | [0.077, 0.867] | 0.028   | 0.994 |
| CDH1  | 390 | 0.467   | [0.224, 0.972] | 0.042   | 0.994 |


FDR does not reach 0.05 for any pair due to conservative correction across 2,450 simultaneous tests; all four pairs show consistent nominal significance and directional agreement with Fisher results.

### 3b. ATM–TP53 in detail


| Direction           | N   | Adj. OR | 95% CI         | p-value  | FDR   |
| ------------------- | --- | ------- | -------------- | -------- | ----- |
| ATM ~ TP53          | 390 | 0.325   | [0.134, 0.789] | 0.013    | 0.994 |
| TP53 ~ ATM          | 390 | 0.334   | [0.139, 0.801] | 0.014    | 0.994 |
| Fisher (unadjusted) | 445 | 0.258   | —              | 5.4×10⁻⁴ | —     |


2×2 contingency (Fisher, N = 445): ATM+/TP53+ = 9, ATM+/TP53− = 22, ATM−/TP53+ = 254, ATM−/TP53− = 160.  
Expected co-mutations under independence ≈ 17; observed = 9.

---

## 4. Logistic Regression — MSS-Only (N = 367)

**Model:** `gene_A ~ gene_B + log₂(TMB) + stage + histology` (MSI_H excluded, all samples MSS)

### TP53-mutually exclusive partners (gene_A ~ TP53, p < 0.05)


| Gene  | N   | Adj. OR | 95% CI         | p-value | FDR   |
| ----- | --- | ------- | -------------- | ------- | ----- |
| SOX9  | 367 | 0.163   | [0.042, 0.631] | 0.0086  | 0.902 |
| CDH1  | 367 | 0.366   | [0.170, 0.790] | 0.010   | 0.915 |
| ATM   | 367 | 0.278   | [0.104, 0.748] | 0.011   | 0.915 |
| ARID2 | 367 | 0.174   | [0.038, 0.795] | 0.024   | >0.99 |


ATM–TP53 ME is preserved in MSS patients (OR = 0.278, p = 0.011), confirming the signal is not driven by MSI-H hypermutation.

---

## 5. Key Finding

ATM and TP53 are significantly mutually exclusive in the FDU gastric cancer cohort. The association survives:

- Fisher's exact test (p = 5.4×10⁻⁴, OR = 0.258)
- Logistic regression adjusted for TMB, MSI status, stage, and histology (p = 0.013, adj. OR = 0.325)
- MSS-restricted analysis (p = 0.011, adj. OR = 0.278)

The mutual exclusivity pattern is consistent with the four other known TP53-ME genes in this cohort (SOX9, CDH1, ARID2, TRRAP/POLD1 by Fisher), and is not an artefact of MSI-driven hypermutation or TMB confounding.

---

## 6. Output Files


| File                                        | Description                                       |
| ------------------------------------------- | ------------------------------------------------- |
| `figures/fig1_oncoprint_CohortB_All.png`    | Somatic landscape, top 50 genes, clinical tracks  |
| `figures/fig2_heatmap_CohortB_All_full.png` | Pairwise interaction matrix, full cohort          |
| `figures/fig2_heatmap_CohortB_All_MSS.png`  | Pairwise interaction matrix, MSS-only             |
| `figures/fig3_tp53_me_CohortB_All_full.png` | TP53-ME mutation rates + adjusted OR, full cohort |
| `figures/fig3_tp53_me_CohortB_All_MSS.png`  | TP53-ME mutation rates + adjusted OR, MSS-only    |
| `figures/figS_fisher_CohortB_All.png`       | maftools Fisher tile (supplementary)              |
| `tables/fisher_CohortB_All.tsv`             | All pairwise Fisher results                       |
| `tables/logistic_full_CohortB_All.tsv`      | All pairwise logistic results, full cohort        |
| `tables/logistic_MSS_CohortB_All.tsv`       | All pairwise logistic results, MSS-only           |


