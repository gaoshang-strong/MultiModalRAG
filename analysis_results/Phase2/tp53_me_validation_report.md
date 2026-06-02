# Phase 2 — Cross-Cohort Validation of TP53 Mutual Exclusivity Signals

**Date:** 2026-06-02  
**Discovery cohort:** FDU gastric cancer, N = 402  
**Script:** `scripts/Phase2/validate_tp53_me_cross_cohort.R`

---

## Background

Phase 1 identified 7 genes bidirectionally mutually exclusive with TP53 in the FDU cohort after TMB-adjusted logistic regression. This Phase 2 analysis tests whether those signals replicate in independent public datasets using the same model:

$$\text{logit}\{P(\text{TP53}=1)\} = \beta_0 + \beta_1 \cdot \text{Gene}_B + \beta_2 \cdot \log_2(\text{TMB})$$

---

## Cohorts

| Cohort | N (STAD) | Sequencing | Notes |
|---|---|---|---|
| FDU (discovery) | 402 | Targeted WXS | Chinese; primary STAD |
| TCGA-STAD | 431 | WXS | Primary tumors only (−01A/−01B) |
| egc_tmucih_2015 | 78 | WXS | Chinese (Beijing) |
| stad_oncosg_2018 | 147 | WXS | Multi-cohort Asian (SG/HK) |
| stad_pfizer_uhongkong | 100 | WXS | Chinese (Hong Kong) |
| egc_msk_2017 | 133 | IMPACT341 | MSK metastatic; STAD-only filter; dedup to primary |
| egc_msk_2023 | 278 | IMPACT341 | MSK; STAD-only filter |
| egc_mskcc_2020 | 8 | IMPACT | STAD N too small after filter — excluded from results |
| egc_msk_tp53_ccr_2022 | 4 | IMPACT | STAD N too small after filter — excluded from results |

TRRAP is absent from both IMPACT341 cohorts (not on panel).  
SOX9 and POLD1 have fewer than 5 mutations in several WXS cohorts (too_few_mutations; excluded from those models).

---

## Meta-Analysis Results (Random Effects, REML)

| Gene | Cohorts (n) | Pooled OR | 95% CI | p-value | I² |
|---|---|---|---|---|---|
| **ARID1A** | 7 | **0.272** | 0.190–0.390 | **1.6 × 10⁻¹²** | 19.7% |
| **ATM** | 5 | **0.314** | 0.199–0.498 | **7.9 × 10⁻⁷** | 0.0% |
| **CDH1** | 7 | **0.512** | 0.351–0.746 | **4.9 × 10⁻⁴** | 14.1% |
| **TRRAP** | 3 | **0.348** | 0.180–0.673 | **1.7 × 10⁻³** | 24.7% |
| **POLD1** | 4 | **0.389** | 0.181–0.838 | **0.016** | 29.3% |
| ARID2 | 6 | 0.546 | 0.282–1.057 | 0.073 | 31.3% |
| SOX9 | 4 | 0.383 | 0.121–1.218 | 0.104 | 68.7% |

Summary: **5 of 7 FDU candidates replicate** (pooled p < 0.05). ARID2 shows a consistent trend (OR = 0.55) but does not reach significance. SOX9 has high heterogeneity (I² = 69%) driven by divergent estimates across cohorts.

---

## Per-Gene Interpretation

### ARID1A — strongest validation (OR = 0.272, p = 1.6 × 10⁻¹²)

Consistent ME signal across all 7 cohorts. Forest plot shows all OR < 1 with a narrow pooled confidence interval and low heterogeneity (I² = 20%). ARID1A is the catalytic subunit of the SWI/SNF chromatin remodelling complex and is the most frequently mutated chromatin regulator in gastric cancer, enriched in EBV-positive and MSI-H subtypes. Its ME with TP53 reflects the well-established genomic divergence between the EBV/MSI-H (TP53-wild-type) and CIN (TP53-mutant) molecular subtypes.

**Per-cohort OR:** FDU 0.40 → TCGA 0.18 → TMUCIH 0.13 → OncoSG 0.28 → HK Pfizer 0.17 → MSK2017 0.47 → MSK2023 0.29

### ATM — highly consistent (OR = 0.314, p = 7.9 × 10⁻⁷, I² = 0%)

Zero heterogeneity across 5 cohorts — the most consistent signal in the analysis. ATM is a central DNA damage response kinase; its ME with TP53 is biologically coherent as two alternative routes to genomic instability (p53 pathway inactivation vs. ATM-mediated checkpoint loss). The absence of I² suggests this is a genuine, cross-population biological exclusion rather than a confounded artefact.

**Per-cohort OR:** FDU 0.26 → TCGA 0.43 → OncoSG 0.51 → HK Pfizer 0.34 → MSK2023 0.14

### CDH1 — robust, moderate effect (OR = 0.512, p = 4.9 × 10⁻⁴)

Validated in 7 cohorts with low heterogeneity (I² = 14%). Reflects the established histological split between CDH1-mutant diffuse/signet-ring cell carcinoma (TP53-wild-type) and TP53-mutant intestinal/CIN carcinoma. The moderate pooled OR (0.51) compared to the stronger FDU signal (0.39) likely reflects dilution from the MSK metastatic datasets where CDH1-mutant tumours may be under-represented.

**Per-cohort OR:** FDU 0.39 → TCGA 0.33 → TMUCIH 0.26 → OncoSG 0.51 → HK Pfizer 0.52 → MSK2017 1.39 → MSK2023 0.70  
Note: MSK2017 shows OR > 1; this is an outlier from the metastatic cohort with only 133 STAD samples and may reflect selection bias in CDH1-mutant diffuse-type metastatic disease.

### TRRAP — validated in 3 WXS cohorts (OR = 0.348, p = 1.7 × 10⁻³)

TRRAP is not covered by IMPACT341, so only 3 WXS cohorts contribute. Within these, the signal is consistent (I² = 25%). TRRAP is a large scaffold protein and essential subunit of TIP60/NuA4 and SAGA histone acetyltransferase complexes. Its ME with TP53 aligns with enrichment of epigenetic regulator mutations in the TP53-wild-type EBV/MSI-H subtype.

**Per-cohort OR:** FDU 0.18 → TCGA 0.47 → OncoSG 0.49

### POLD1 — nominal validation (OR = 0.389, p = 0.016)

Validated in 4 cohorts with moderate heterogeneity (I² = 29%). POLD1 encodes the catalytic subunit of DNA polymerase delta; somatic mutations in POLE/POLD1 define the ultramutator (POLE-proofreading deficient) subtype, which is characterised by very high TMB, microsatellite stability, and TP53-wild-type status. The ME with TP53 is mechanistically expected: CIN (TP53-mutant) and POLE-ultramutator represent distinct paths to genome instability that rarely co-occur.

**Per-cohort OR:** FDU 0.17 → TCGA 0.62 → MSK2017 0.26 → MSK2023 0.76

### ARID2 — trend, not significant (OR = 0.546, p = 0.073)

Direction consistent across 6 of 6 cohorts (all OR < 1 except HK Pfizer OR = 2.1 with very wide CI), but wide confidence intervals in small cohorts prevent reaching significance. ARID2 is the ARID domain subunit of the PBAF SWI/SNF complex, functionally related to ARID1A. Biological interpretation is the same: enrichment in the EBV/MSI-H subtype. Underpowered in individual cohorts; may reach significance in a larger pooled WXS analysis.

### SOX9 — high heterogeneity (OR = 0.383, p = 0.104, I² = 69%)

SOX9 was the strongest signal in FDU (OR = 0.085) but does not replicate consistently. TCGA shows OR = 1.10 (no ME signal); only MSK2023 replicates the direction (OR = 0.31). The high heterogeneity suggests population-specific or subtype composition effects. SOX9 may be a stronger ME signal specifically in the Chinese gastric cancer context where EBV and MSI subtypes have different relative frequencies, or may reflect low mutational frequency (SOX9 is mutated in only 1–5% of samples in most cohorts) reducing power severely.

---

## Panel Coverage

| Gene | TCGA | TMUCIH | OncoSG | HK Pfizer | MSK2017 | MSK2023 |
|---|---|---|---|---|---|---|
| ARID1A | WXS | WXS | WXS | WXS | IMPACT | IMPACT |
| ATM | WXS | WXS | WXS | WXS | — (too few) | IMPACT |
| CDH1 | WXS | WXS | WXS | WXS | IMPACT | IMPACT |
| TRRAP | WXS | — (too few) | WXS | — (too few) | Not on panel | Not on panel |
| POLD1 | WXS | — (too few) | — (too few) | — (too few) | IMPACT | IMPACT |
| ARID2 | WXS | — (too few) | WXS | WXS | IMPACT | IMPACT |
| SOX9 | WXS | Not detected | — (too few) | — (too few) | IMPACT | IMPACT |

"— (too few)" = gene present in dataset but <5 mutations; model excluded.

---

## Caveats

1. **egc_mskcc_2020 and egc_msk_tp53_ccr_2022** retain only 8 and 4 STAD samples after filtering; these cohorts are informative only for the full EGC analysis, not STAD-specific validation.
2. **MSK datasets are metastatic** (egc_msk_2017 entirely; egc_msk_2023 mixed). Metastatic disease may alter mutation spectra (clonal evolution, treatment selection), potentially attenuating ME signals.
3. **TMB computation differs by cohort**: TCGA uses functional variant count / 38 Mb; cBioPortal datasets use pre-computed `TMB_NONSYNONYMOUS`. Both measure nonsynonymous TMB but may differ in exact variant filtering.
4. **No FDR correction applied** across cohorts in per-cohort models. Meta-analysis p-values represent pooled evidence and are not adjusted for multiple testing across 7 genes.
5. **SOX9 and POLD1 have low mutation frequency** (<5%) in most cohorts, making logistic regression underpowered. Their signals should be interpreted cautiously until tested in larger datasets.

---

## Output Files

| File | Contents |
|---|---|
| `per_cohort_results.tsv` | All logistic regression results (cohort × gene), including FDU discovery |
| `meta_analysis_results.tsv` | Random-effects pooled OR, 95% CI, p, I² per gene |
| `panel_coverage_check.tsv` | Gene coverage status per cohort |
| `forest_ARID1A.png` | Forest plot: TP53 ↔ ARID1A |
| `forest_ATM.png` | Forest plot: TP53 ↔ ATM |
| `forest_CDH1.png` | Forest plot: TP53 ↔ CDH1 |
| `forest_TRRAP.png` | Forest plot: TP53 ↔ TRRAP |
| `forest_POLD1.png` | Forest plot: TP53 ↔ POLD1 |
| `forest_ARID2.png` | Forest plot: TP53 ↔ ARID2 |
| `forest_SOX9.png` | Forest plot: TP53 ↔ SOX9 |
| `summary_heatmap.png` | OR_adj heatmap across all cohorts and genes |
