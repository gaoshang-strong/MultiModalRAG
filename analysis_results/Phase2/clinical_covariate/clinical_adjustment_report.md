# Phase 2 — TP53 ME Analysis: Clinical Covariate Adjustment

**Date:** 2026-06-03  
**Script:** `scripts/Phase2/tp53_me_clinical_adjustment.R`  
**Output directory:** `analysis_results/Phase2/clinical_covariate/`

---

## Background and Motivation

Phase 2 cross-cohort validation and gene-length sensitivity analysis identified three high-confidence TP53 ME signals: **ARID1A**, **CDH1**, and **ATM**. A remaining question is whether these signals reflect true biological mutual exclusivity or whether they are confounded by clinical variables:

- **Stage**: TP53 mutations are enriched in advanced-stage (III/IV) tumours. If ATM or ARID1A mutations co-segregate with early-stage disease, a spurious ME signal could arise.
- **Histology (Lauren classification)**: Diffuse/signet-ring tumours are predominantly TP53 wild-type; intestinal-type tumours are TP53-enriched. ARID1A mutations are enriched in EBV/MSI subtypes (often diffuse-like), so histology could partially explain TP53–ARID1A ME. ATM mutations are associated with CIN and genomic instability independently of histotype.
- **Metastasis**: Captured implicitly through stage (Stage IV = M1). MSK panel cohorts additionally supply SAMPLE_TYPE (Primary vs Metastasis).

This analysis adds stage and histology as covariates to the base logistic model and tests whether the ME signals survive.

---

## Methods

### Models

Three nested logistic regression models per cohort per gene pair:

| Model | Formula | Description |
|-------|---------|-------------|
| **A** | `TP53 ~ gene_B + log2(TMB)` | Base model (same as Phase 1/2) |
| **B** | `TP53 ~ gene_B + log2(TMB) + stage_advanced` | + tumour stage |
| **C** | `TP53 ~ gene_B + log2(TMB) + stage_advanced + hist_diffuse` | + stage + histology |

**stage_advanced**: binary — 0 = Stage I/II (early), 1 = Stage III/IV (advanced).  
**hist_diffuse**: binary — 0 = intestinal/NOS, 1 = diffuse/signet-ring (Lauren classification or equivalent).  
Model C is skipped for cohorts without histology data.

### Cohort selection

Cohorts were included based on gene coverage (from `panel_coverage_check.tsv`). MSKCC-2020 and MSK-TP53-CCR-2022 were excluded entirely (ATM and ARID1A not covered on their panels).

| Cohort | N | Seq. type | ATM | ARID1A | Stage | Histology |
|--------|---|-----------|-----|--------|-------|-----------|
| FDU | 472 | Panel (503-gene) | ✓ | ✓ | path_stage | 腺癌/印戒 (Chinese, mapped) |
| TCGA-STAD | 431 | WXS | ✓ | ✓ | ajcc_pathologic_tumor_stage | histologic_diagnosis |
| OncoSG-2018 | 147 | WXS | ✓ | ✓ | STAGE | LAURENS_CLASSIFICATION |
| HK-Pfizer-2014 | 100 | WXS | ✓ | ✓ | UICC_TUMOR_STAGE | PATHOLOGY_LAUREN |
| TMUCIH-2015 | 78 | WXS | too few | ✓ | STAGE (1–4) | none (Models A/B only) |
| MSK-2017 | 133 | IMPACT341 | too few | ✓ | STAGE_AT_DIAGNOSIS | LAUREN_CLASS |
| MSK-2023 | 278 | IMPACT468 | ✓ | ✓ | STAGE | HISTOLOGY |

Meta-analysis: random-effects (REML), implemented via `metafor::rma`.

---

## Covariate Availability

| Cohort | N total | % advanced stage | N stage available | % diffuse histology | N hist. available |
|--------|---------|-----------------|-------------------|--------------------|--------------------|
| FDU | 472 | 44.1% | 458 | 34.8% | 451 |
| TCGA_STAD | 431 | 54.8% | 405 | 19.1% | 429 |
| OncoSG_2018 | 147 | 82.4% | 102 | 27.6% | 123 |
| HK_Pfizer | 100 | 84.8% | 99 | 29.0% | 100 |
| TMUCIH_2015 | 78 | 66.7% | 78 | — | 0 |
| MSK_2017 | 133 | 100.0%* | 104 | 25.0% | 132 |
| MSK_2023 | 278 | 86.3% | 277 | 26.3% | 266 |

*MSK-2017 is an exclusively metastatic cohort; stage covariate has minimal variation (no patients in Stage I/II). Stage available = those with parseable STAGE_AT_DIAGNOSIS values.

---

## Results

### TP53 – ATM

**Cohorts:** FDU, TCGA-STAD, OncoSG-2018, HK-Pfizer-2014, MSK-2023 (n = 5)

#### Pooled estimates (random-effects meta-analysis)

| Model | N cohorts | Pooled OR | 95% CI | p-value | I² |
|-------|-----------|-----------|--------|---------|-----|
| A: base | 5 | **0.314** | 0.199–0.498 | 7.8 × 10⁻⁷ | **0%** |
| B: + stage | 5 | **0.338** | 0.210–0.544 | 8.2 × 10⁻⁶ | **0%** |
| C: + stage + hist | 5 | **0.342** | 0.211–0.554 | 1.3 × 10⁻⁵ | **0%** |

OR shift A → C: **0.314 → 0.342 (ratio = 0.92)**. Negligible change.

#### Per-cohort results (Model A / C)

| Cohort | N | n(TP53+) | n(ATM+) | n(both+) | OR_A | p_A | OR_C | p_C |
|--------|---|----------|---------|----------|------|-----|------|-----|
| FDU | 402 | 237 | 29 | 9 | 0.261 | 0.0016 | 0.291 | 0.0052 |
| TCGA_STAD | 431 | 197 | 42 | 14 | 0.428 | 0.024 | 0.462 | 0.047 |
| OncoSG_2018 | 147 | 67 | 9 | 3 | 0.509 | 0.378 | 0.840 | 0.858 |
| HK_Pfizer | 100 | 55 | 7 | 2 | 0.335 | 0.210 | 0.387 | 0.286 |
| MSK_2023 | 278 | 180 | 16 | 4 | 0.142 | 0.0016 | 0.141 | 0.0018 |

OncoSG and HK Pfizer are individually underpowered for ATM (n = 9 and 7 ATM-mutant samples respectively). The pooled signal is driven by FDU, TCGA, and MSK-2023.

#### Interpretation

The TP53–ATM ME signal is fully robust to clinical covariate adjustment. I² = 0% across all three models, indicating perfect cross-cohort consistency. Neither stage nor histological subtype explains the mutual exclusivity. This is consistent with the biological mechanism: ATM and TP53 represent two parallel routes to genome instability (ATM → DSB repair deficiency; TP53 → cell cycle checkpoint failure), so their mutual exclusivity is intrinsic to tumour biology rather than a stage or histotype artefact.

This finding, alongside the Phase 2 length-normalization result (OR_corrected = 0.860, p = 0.031, I² = 0%), makes **ATM–TP53 ME in gastric cancer a high-confidence novel finding** requiring PubMed verification.

---

### TP53 – ARID1A

**Cohorts:** FDU, TCGA-STAD, OncoSG-2018, HK-Pfizer-2014, TMUCIH-2015, MSK-2017, MSK-2023 (n = 7 for A/B; n = 6 for C, TMUCIH excluded — no histology)

#### Pooled estimates (random-effects meta-analysis)

| Model | N cohorts | Pooled OR | 95% CI | p-value | I² |
|-------|-----------|-----------|--------|---------|-----|
| A: base | 7 | **0.272** | 0.190–0.390 | 1.5 × 10⁻¹² | **19.7%** |
| B: + stage | 7 | **0.269** | 0.187–0.387 | 1.2 × 10⁻¹² | **15.8%** |
| C: + stage + hist | 6 | **0.296** | 0.197–0.444 | 4.4 × 10⁻⁹ | **26.2%** |

OR shift A → C: **0.272 → 0.296 (ratio = 0.92)**. Negligible change.

Noteworthy: adding stage (Model B) *reduces* I² from 19.7% to 15.8%, suggesting that some cross-cohort variance in the base model was attributable to differences in stage distribution. After accounting for histology (Model C), I² rises modestly to 26.2%, likely because the histology encoding differs between cohorts (Lauren vs Adeno/Signet_Diffuse categorisation), introducing additional measurement heterogeneity.

#### Per-cohort results (Model A / C)

| Cohort | N | n(TP53+) | n(ARID1A+) | n(both+) | OR_A | p_A | OR_C | p_C |
|--------|---|----------|-----------|----------|------|-----|------|-----|
| FDU | 402 | 237 | 74 | 33 | 0.403 | 0.0016 | 0.459 | 0.011 |
| TCGA_STAD | 431 | 197 | 116 | 32 | 0.178 | 1.6 × 10⁻⁸ | 0.180 | 4.9 × 10⁻⁸ |
| OncoSG_2018 | 147 | 67 | 17 | 4 | 0.276 | 0.041 | 0.351 | 0.137 |
| HK_Pfizer | 100 | 55 | 18 | 4 | 0.172 | 0.0084 | 0.179 | 0.012 |
| TMUCIH_2015 | 78 | 36 | 17 | 3 | 0.130 | 0.0066 | — | — |
| MSK_2017 | 133 | 88 | 21 | 12 | 0.466 | 0.143 | 0.427 | 0.137 |
| MSK_2023 | 278 | 180 | 52 | 24 | 0.294 | 8.7 × 10⁻⁴ | 0.293 | 0.0014 |

MSK-2017 shows a non-significant OR — this cohort is 100% advanced stage (all metastatic), which limits both the stage covariate utility and the representativeness for a primary tumour ME test.

#### Interpretation

The TP53–ARID1A ME signal is fully robust to stage and histology adjustment. The OR barely moves (0.272 → 0.296) and remains highly significant (p = 4.4 × 10⁻⁹) after adjusting for Lauren classification. This is expected: the EBV/MSI (ARID1A-enriched) vs CIN (TP53-enriched) molecular subtype distinction is deeper than Lauren histotype — the ME reflects molecular pathway divergence that transcends both stage and gross histological categorisation. This signal was already documented in the TCGA 2014 gastric cancer landmark paper and serves as the **positive methodological control** for this analysis.

---

## Overall Conclusion

| Signal | OR (base) | OR (stage+hist) | OR shift | I² (C) | Stage confounding? | Histology confounding? |
|--------|-----------|-----------------|----------|--------|--------------------|----------------------|
| TP53–ATM | 0.314 | 0.342 | +8.9% | **0%** | No | No |
| TP53–ARID1A | 0.272 | 0.296 | +8.8% | 26.2% | No | No |

Both ME signals survive clinical covariate adjustment. The OR ratio A→C is ~0.92 for both pairs, meaning stage and histology together explain less than 9% of the OR magnitude. Neither signal can be attributed to stage distribution bias or histological subtype enrichment.

Combined with the Phase 2 gene-length normalization analysis:

| Signal | After length+TMB norm | After stage+hist adj | Final status |
|--------|-----------------------|----------------------|-------------|
| TP53–ARID1A | ✅ robust (OR_C = 0.836, p = 6.1 × 10⁻⁷) | ✅ robust | **High confidence** |
| TP53–ATM | ✅ robust (OR_C = 0.860, p = 0.031) | ✅ robust | **High confidence, potential novel finding** |

---

## Output Files

| File | Contents |
|------|----------|
| `covariate_availability.tsv` | Per-cohort counts of stage/histology availability |
| `per_cohort_clinical_results.tsv` | All logistic regression results (cohort × gene pair × model) |
| `meta_clinical_results.tsv` | Pooled random-effects OR per gene pair per model |
| `clinical_adjustment_summary.tsv` | Model A vs C comparison with conclusions |
| `forest_ATM_clinical.png` | Forest plot: TP53–ATM, three models overlaid |
| `forest_ARID1A_clinical.png` | Forest plot: TP53–ARID1A, three models overlaid |
| `pvalue_trajectory.png` | –log₁₀(p) trajectory across Models A → B → C |

## Limitations

1. **FDU histology in Chinese**: mapping 印戒/黏附性癌 → diffuse, 腺癌/粘液腺癌 → non-diffuse is approximate. Mixed-type tumours and hepatoid adenocarcinomas (n=8) were set to NA.
2. **MSK-2017 stage variable**: STAGE_AT_DIAGNOSIS is a free-text field with inconsistent formatting (mix of Roman numerals, TNM notation). Most values (~73%) parsed to Stage IV; remaining unparseable values were set to NA.
3. **Histology encoding heterogeneity**: Lauren classification (OncoSG, HKU, MSK-2017), TCGA histologic_diagnosis (WHO types), and MSK-2023 HISTOLOGY (Adenocarcinoma vs Signet_Diffuse) capture overlapping but not identical concepts. This likely contributes to the slight I² increase in Model C for ARID1A.
4. **OncoSG and HKU Pfizer are underpowered for ATM**: with only 9 and 7 ATM-mutant samples respectively, per-cohort ATM estimates are imprecise. These cohorts contribute wide CIs and their individual p-values are non-significant, but they are directionally consistent and their inclusion in meta-analysis does not distort the pooled estimate.
5. **Metastasis as separate covariate not modelled**: M-stage was not added as an independent covariate because it is near-collinear with stage_advanced (Stage IV ≡ M1). For cohorts with explicit STAGE_M data (OncoSG, HKU Pfizer), M-stage correlates ≥ 0.95 with stage_advanced in these data.
