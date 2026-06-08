# Cross-Cohort Validation — ATM–TP53 Mutual Exclusivity

**Date:** 2026-06-07  
**Script:** `TP53_ATM_paper/scripts/validation_ATM_TP53_ME.R`

---

## 1. Overview

Cross-cohort validation of ATM–TP53 ME in 4 independent gastric cancer cohorts (plus paired MSS sub-strata). Primary pooled estimate uses FDU (full) + TCGA-STAD (MSS) + MSK 2023 + OncoSG 2018.

**Model:** `TP53 ~ ATM + log₂(TMB) + stage + histology` (terms included where available per cohort)  
**MSI handling:** MSI_H is NOT in the model. MSS-only rows are shown as sub-analyses (not used in the primary pooled estimate, except TCGA where MSS-only is the primary stratum — see §2).

---

## 2. MSI Handling Rationale

**TCGA-STAD:** Of 42 ATM+ samples in the full TCGA cohort, **30 (71%) are MSI-H**. MSI-H ATM mutations are hypermutation artifacts from MMR deficiency, not driver events. In MSS-only (N=348), ATM+=12 with 2 co-mutations (expected ≈5.8) — the true ME signal. Full-cohort n_both=14 is dominated by 12 MSI-H artifact cases. → **TCGA MSS-only is the primary TCGA stratum.**

**FDU:** MSI-H rate = 5.1% (low contamination). Full cohort is primary. MSS-only is a sensitivity analysis.

**MSK 2023 / OncoSG 2018:** Full cohort is primary. MSS sub-analysis: OncoSG MSS has ATM+=5 (insufficient for logistic); MSK IMPACT data lacks MSI column.

---

## 3. Cohorts and Results

| Cohort | Stratum | N | ATM+ | TP53+ | Co-mut | Adj. OR | 95% CI | p (logistic) | p (Fisher) | Pool |
|---|---|---|---|---|---|---|---|---|---|---|
| FDU [panel, Chinese] | Full | 402 | 29 | 237 | 9 | 0.303 | [0.127–0.721] | 0.0069 | 0.0027 | ✓ |
| FDU [panel, Chinese] | MSS only | 378 | 22 | 230 | 7 | 0.280 | [0.106–0.744] | 0.011 | 0.0060 | — |
| TCGA-STAD [WES] | **MSS only** | 348 | 12 | 167 | 2 | 0.220 | [0.047–1.037] | 0.056 | 0.037 | ✓ |
| TCGA-STAD [WES] | All (incl. MSI-H) | 431 | 42 | 197 | 14 | 0.467 | [0.218–1.001] | 0.050 | 0.104 | — |
| MSK 2023 [IMPACT341] | Full, STAD only | 278 | 16 | 180 | 4 | 0.136 | [0.039–0.472] | 0.0017 | 0.0018 | ✓ |
| OncoSG 2018 [WES, Asian] | Full | 147 | 9 | 67 | 3 | 0.708 | [0.102–4.906] | 0.727 | 0.510 | ✓ |
| OncoSG 2018 [WES, Asian] | MSS only | 133 | 5 | 66 | 3 | — | — | — (n<5) | 0.680 | — |

Expected co-mutations under independence: FDU ≈17, TCGA(MSS) ≈5.8, MSK ≈10.4, OncoSG ≈4.1.  
All cohorts: OR < 1 (direction consistent). TCGA ALL diluted by MSI-H artifact (OR shifts 0.22→0.47; Fisher p becomes NS).

---

## 4. Meta-Analysis (RE, REML; 4 primary cohorts)

Pooling: FDU (full) + TCGA-STAD (MSS) + MSK 2023 + OncoSG 2018.

| | Value |
|---|---|
| Pooled OR | **0.258** |
| 95% CI | **[0.140–0.477]** |
| p-value | **1.52 × 10⁻⁵** |
| I² | **0.0%** |
| N cohorts | 4 |

---

## 5. Interpretation

**I² = 0.0%** across 4 cohorts spanning Chinese targeted panel (FDU), WES (TCGA, OncoSG), and IMPACT341 (MSK). Zero heterogeneity = pan-population biological constraint, not platform- or cohort-specific.

**TCGA borderline (logistic p=0.056, Fisher p=0.037):** ATM+=12 after MSS filtering limits power. Direction and Fisher test are consistent with ME.

**FDU MSS sensitivity (OR=0.280, p=0.011):** ME persists after removing MSI-H patients.

**TCGA ALL dilution (OR=0.467, p=0.050, Fisher NS):** MSI-H artifact inflates apparent co-mutations; including MSI-H patients shifts OR toward 1 and loses Fisher significance — visual demonstration in forest plot.

---

## 6. Output Files

| File | Description |
|---|---|
| `figures/forest_ATM_TP53.png` | Forest plot: 6 cohort rows + RE pooled; Adj. OR [95% CI] and p columns |
| `figures/summary_tile_ATM.png` | Summary tile: all 7 rows (log₂OR + significance) |
| `figures/fig_tp53rate_by_ATM.png` | TP53 rate in ATM+ vs ATM− per cohort (4 primary) |
| `figures/fig_2x2_tiles.png` | 2×2 contingency tiles with O/E coloring (4 primary) |
| `figures/fig_obs_vs_exp.png` | Observed vs expected co-mutations bubble chart |
| `figures/fig_mutation_freq.png` | ATM% and TP53% mutation frequency bars |
| `figures/fig_combined_validation.png` | Patchwork composite |
| `tables/per_cohort_ATM.tsv` | Full per-cohort logistic + Fisher results |
| `tables/meta_ATM_full.tsv` | RE meta-analysis summary |
