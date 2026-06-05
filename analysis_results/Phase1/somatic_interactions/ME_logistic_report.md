# Mutual Exclusivity Analysis — TMB-Adjusted Logistic Regression

**Cohort:** FDU gastric cancer, N = 445 (402 with TMB data)  
**Method:** Binomial logistic regression per gene pair; covariate = log₂(max(TMB, 0.1))  
**Gene set:** Top 50 most-mutated genes (broad run); 2,450 directional models (1,225 pairs × 2 directions)  
**Multiple testing:** Benjamini–Hochberg FDR across all 2,450 models  
**Date:** 2026-06-02

---

## Summary


| Category                              | Count |
| ------------------------------------- | ----- |
| Total directional models tested       | 2,450 |
| ME direction (OR_adj < 1)             | 1,170 |
| Co-occurrence direction (OR_adj > 1)  | 1,295 |
| ME pairs with p_adj < 0.05            | 40    |
| Co-occurrence pairs with p_adj < 0.05 | 92    |
| FDR < 0.05 pairs                      | 0     |


No pair survives FDR < 0.05 after BH correction across 2,450 models, which is expected given the stringency of this multiple-testing burden and the limited sample size. All findings below use nominal p_adj < 0.05 as the threshold.

---

## Bidirectionally Confirmed ME Pairs

A pair is considered robustly mutually exclusive when both directions (A ~ B and B ~ A) return p_adj < 0.05 with concordant OR_adj < 1. This double-check filters out asymmetric model artefacts.


| Gene A    | Gene B     | OR_adj (A~B) | OR_adj (B~A) | p_adj (A~B) | p_adj (B~A) | Also ME by Fisher?                   |
| --------- | ---------- | ------------ | ------------ | ----------- | ----------- | ------------------------------------ |
| **TP53**  | **SOX9**   | 0.085        | 0.117        | 1.5 × 10⁻⁴  | 8.8 × 10⁻⁴  | Yes                                  |
| **TP53**  | **TRRAP**  | 0.184        | 0.351        | 9.1 × 10⁻⁴  | 0.048       | Yes                                  |
| **TP53**  | **ATM**    | 0.261        | 0.309        | 1.6 × 10⁻³  | 5.7 × 10⁻³  | Yes                                  |
| **TP53**  | **POLD1**  | 0.172        | 0.308        | 2.1 × 10⁻³  | 0.044       | Yes                                  |
| **TP53**  | **ARID1A** | 0.403        | 0.526        | 1.6 × 10⁻³  | 0.027       | Yes                                  |
| **TP53**  | **ARID2**  | 0.186        | 0.293        | 2.6 × 10⁻³  | 0.030       | Yes                                  |
| **TP53**  | **CDH1**   | 0.386        | 0.387        | 4.6 × 10⁻³  | 4.6 × 10⁻³  | Yes                                  |
| **PXDNL** | **ATR**    | 0.030        | 0.051        | 4.4 × 10⁻³  | 0.023       | **No — Fisher showed Co-Occurrence** |
| **MUC16** | **CIC**    | 0.110        | 0.133        | 0.020       | 0.044       | **No — Fisher showed Co-Occurrence** |
| **APC**   | **ARID1A** | 0.315        | 0.336        | 0.040       | 0.048       | **No — Fisher showed Co-Occurrence** |


---

## TP53 Mutual Exclusivity — Complete Picture

TP53 is the dominant ME hub. After TMB adjustment, 16 unique genes remain nominally ME with TP53 (p_adj < 0.05, at least one direction). Nine of these are bidirectionally confirmed (see table above plus the additional entries below).

### All TP53-ME partners ranked by p_adj (TP53 ~ X direction)


| Partner | OR_adj | β₁    | p_adj      | FDR_adj | Both directions?                      | Fisher OR |
| ------- | ------ | ----- | ---------- | ------- | ------------------------------------- | --------- |
| SOX9    | 0.085  | −2.47 | 1.5 × 10⁻⁴ | 0.124   | Yes                                   | 0.106     |
| TRRAP   | 0.184  | −1.69 | 9.1 × 10⁻⁴ | 0.249   | Yes                                   | 0.305     |
| RYR1    | 0.247  | −1.40 | 1.6 × 10⁻³ | 0.282   | No (X~TP53 n.s.)                      | 0.432     |
| ARID1A  | 0.403  | −0.91 | 1.6 × 10⁻³ | 0.282   | Yes                                   | 0.529     |
| ATM     | 0.261  | −1.34 | 1.6 × 10⁻³ | 0.282   | Yes                                   | 0.258     |
| POLD1   | 0.172  | −1.76 | 2.1 × 10⁻³ | 0.336   | Yes                                   | 0.233     |
| ARID2   | 0.186  | −1.68 | 2.6 × 10⁻³ | 0.380   | Yes                                   | 0.329     |
| CIC     | 0.199  | −1.62 | 3.2 × 10⁻³ | 0.431   | No (CIC~TP53 p=0.044, included below) | 0.329     |
| CDH1    | 0.386  | −0.95 | 4.6 × 10⁻³ | 0.480   | Yes                                   | 0.423     |
| CHD4    | 0.217  | −1.53 | 5.8 × 10⁻³ | 0.480   | No (CHD4~TP53 n.s.)                   | 0.356     |
| KMT2C   | 0.339  | −1.08 | 8.0 × 10⁻³ | 0.513   | No (KMT2C~TP53 p=0.26)                | 0.465     |
| ABL2    | 0.340  | −1.08 | 0.027      | 0.818   | No (ABL2~TP53 p=0.046)                | 0.426     |
| BRCA2   | 0.319  | −1.14 | 0.028      | 0.818   | No                                    | 0.388     |
| NOTCH3  | 0.315  | −1.16 | 0.033      | 0.840   | No                                    | 0.388     |
| ACVR2A  | 0.328  | −1.12 | 0.034      | 0.849   | No                                    | 0.445     |
| BCORL1  | 0.374  | −0.98 | 0.049      | 0.903   | No                                    | 0.503     |


**Notes on TP53-ME partners:**

- **SOX9** is the strongest single signal (OR_adj = 0.085), robust both before and after TMB adjustment. SOX9 is a transcription factor implicated in lineage specification; its ME with TP53 likely reflects enrichment in the MSI-H/EBV subtype.
- **TRRAP, CHD4, ARID1A, ARID2, KMT2C** are all chromatin remodellers / epigenetic regulators. Their collective ME with TP53 is consistent with enrichment of these mutations in the TP53-wild-type (EBV/MSI-H) gastric cancer molecular subtype.
- **ATM, POLD1, BRCA2** are DNA damage response genes. Their ME with TP53 may reflect mutual exclusivity of two routes to genome instability: p53 pathway inactivation (CIN) vs. DNA repair deficiency (dMMR/hypermutation).
- **CDH1** ME with TP53 (bidirectional, OR_adj ≈ 0.39) reflects the well-established histological divergence: CDH1-mutant diffuse/signet ring cell carcinoma vs. TP53-mutant intestinal/CIN.
- **RYR1, ABL2, NOTCH3, ACVR2A, BCORL1** are nominally significant in one direction only; these should be treated as exploratory signals.

### TP53-ME partners that do NOT survive TMB adjustment

The following gene pairs were ME by Fisher test but become non-significant (or flip to co-occurrence) after adjusting for TMB. These Fisher signals were likely confounded by TMB differences between TP53-mutant and TP53-wild-type tumours.


| Partner | Fisher OR | Fisher p | Logistic OR_adj | Logistic direction   | Interpretation                    |
| ------- | --------- | -------- | --------------- | -------------------- | --------------------------------- |
| NOTCH1  | n/a       | (ME)     | 2.31            | **Co-Occurrence**    | TMB confounder; flipped direction |
| UBR5    | n/a       | (ME)     | 1.82            | **Co-Occurrence**    | TMB confounder                    |
| ANK3    | n/a       | (ME)     | 1.89            | **Co-Occurrence**    | TMB confounder                    |
| GNAS    | n/a       | (ME)     | 1.50            | **Co-Occurrence**    | TMB confounder                    |
| PXDNL   | n/a       | (ME)     | 1.46            | **Co-Occurrence**    | TMB confounder                    |
| PIK3CA  | n/a       | (ME)     | 1.36            | **Co-Occurrence**    | TMB confounder                    |
| ALK     | n/a       | (ME)     | 1.43            | **Co-Occurrence**    | TMB confounder                    |
| CREBBP  | n/a       | (ME)     | 1.08            | Co-Occurrence (n.s.) | TMB confounder                    |


These pairs should be removed from downstream biological interpretation of TP53 mutual exclusivity.

---

## TMB Confounding — Direction-Flip Analysis

A key value of logistic regression over Fisher's test is identifying pairs where apparent ME or co-occurrence is driven by TMB heterogeneity.

### Fisher Co-Occurrence → Logistic ME (p_adj < 0.05)

These pairs appeared to co-occur in the unadjusted test but are actually mutually exclusive after accounting for TMB. The co-occurrence signal was an artifact: both genes happened to accumulate in high-TMB tumours independently, creating a spurious positive association.


| Gene A | Gene B | OR_adj | p_adj      | Fisher OR     | Fisher p |
| ------ | ------ | ------ | ---------- | ------------- | -------- |
| PXDNL  | ATR    | 0.030  | 4.4 × 10⁻³ | 1.01 (Co-occ) | 1.0      |
| ATR    | PXDNL  | 0.051  | 0.023      | —             | —        |
| MUC16  | CREBBP | 0.095  | 0.015      | 1.75 (Co-occ) | 0.35     |
| MUC16  | ACVR2A | 0.110  | 0.016      | 3.07 (Co-occ) | 0.10     |
| MUC16  | CIC    | 0.110  | 0.020      | 1.75 (Co-occ) | 0.35     |
| CIC    | MUC16  | 0.133  | 0.044      | —             | —        |
| CREBBP | MUC16  | 0.136  | 0.046      | —             | —        |
| CHD4   | ALK    | 0.062  | 0.033      | 1.34 (Co-occ) | 0.55     |
| CREBBP | ALK    | 0.066  | 0.036      | 1.27 (Co-occ) | 0.57     |
| ALK    | CREBBP | 0.088  | 0.045      | —             | —        |
| ALK    | KMT2C  | 0.129  | 0.037      | 2.56 (Co-occ) | 0.15     |
| EPHA3  | KMT2C  | 0.115  | 0.037      | 1.36 (Co-occ) | 0.66     |
| APC    | ARID1A | 0.315  | 0.040      | 1.16 (Co-occ) | 0.67     |
| ARID1A | APC    | 0.336  | 0.048      | —             | —        |
| PXDNL  | UBR5   | 0.176  | 0.043      | 3.12 (Co-occ) | 0.064    |


**Notable patterns:**

- **MUC16** emerges as a new ME hub not visible by Fisher: mutually exclusive with CREBBP, ACVR2A, CIC, and PIK3CA after TMB adjustment. MUC16 is a highly mutated mucin gene; its apparent co-occurrence with other genes in Fisher tests is driven by hypermutated tumours.
- **ALK** is ME with CREBBP and KMT2C after TMB adjustment — two chromatin remodellers. ALK mutations in gastric cancer may mark a distinct molecular subtype incompatible with chromatin remodelling gene mutations.
- **APC ↔ ARID1A** mutual exclusivity (bidirectional, OR_adj ≈ 0.32): APC-mutant (Wnt-activated, likely intestinal-type) tumours exclude ARID1A mutations (EBV/MSI-H subtype). This is biologically coherent.
- **PXDNL ↔ ATR** mutual exclusivity (OR_adj = 0.030, the most extreme non-TP53 ME signal): PXDNL is an oxidative enzyme; the mechanism of ME with ATR (replication stress kinase) is unclear and this pair needs biological follow-up.

---

## Comparison: Fisher vs. Logistic Regression


| Result category                                    | Count (directional pairs)                            |
| -------------------------------------------------- | ---------------------------------------------------- |
| ME by both Fisher and logistic (p < 0.05)          | 22                                                   |
| ME by Fisher only (logistic n.s. or co-occurrence) | ~18 (includes TP53 artefact pairs above)             |
| ME by logistic only (Fisher n.s. or co-occurrence) | ~18 (TMB-unmasked, including MUC16, ALK, APC-ARID1A) |
| Agreement: co-occurrence in both                   | 87                                                   |


The overall concordance between Fisher and logistic is high for co-occurrence pairs (most strong co-occurrence signals survive TMB adjustment). For ME pairs, TMB adjustment matters most: it removes ~8 spurious TP53-ME pairs and reveals ~15 new ME pairs masked by TMB confounding.

---

## Prioritised ME Signals for Follow-up

Based on bidirectional confirmation and biological plausibility:


| Priority | Pair          | OR_adj      | Key rationale                                                            |
| -------- | ------------- | ----------- | ------------------------------------------------------------------------ |
| 1        | TP53 ↔ SOX9   | 0.085/0.117 | Strongest overall; bidirectional; SOX9 = EBV/MSI-H lineage marker        |
| 2        | TP53 ↔ ATM    | 0.261/0.309 | DNA repair vs. p53 pathway — two routes to instability                   |
| 3        | TP53 ↔ POLD1  | 0.172/0.308 | Ultramutator (POLE/POLD spectrum) vs. CIN — well-known divergence        |
| 4        | TP53 ↔ CDH1   | 0.386/0.387 | Diffuse vs. intestinal histological split; bidirectional, symmetric OR   |
| 5        | TP53 ↔ ARID1A | 0.403/0.526 | EBV/MSI-H SWI/SNF vs. CIN; bidirectional                                 |
| 6        | APC ↔ ARID1A  | 0.315/0.336 | Wnt-activated intestinal-type vs. epigenetic dysregulation; TMB-unmasked |
| 7        | MUC16 cluster | 0.095–0.198 | TMB-driven artefact removed; MUC16 as subtype marker needs validation    |
| 8        | PXDNL ↔ ATR   | 0.030/0.051 | Most extreme non-TP53 ME; mechanism unclear; exploratory                 |


---

## Caveats

1. **No pair survives FDR < 0.05.** With 2,450 tests and N = 402, the analysis is underpowered for strict FDR control. All results are hypothesis-generating.
2. **MSI status not included in the model.** The current model only adjusts for log₂(TMB). MSI-H tumours have both high TMB and a distinct mutation landscape. Adding MSI as a covariate could further refine the ME signals, particularly for SOX9 and ARID1A.
3. **Histological subtype not included.** Signet ring cell (CDH1-mutant) vs. intestinal type (TP53-mutant) are distinct strata; not adjusting for histology means the CDH1–TP53 ME captures subtype divergence, not an intra-subtype genomic interaction.
4. **Directionality is asymmetric by design.** Logistic regression with a binary outcome is not symmetric: A ~ B and B ~ A can give different β values. A "bidirectionally confirmed" pair is the conservative criterion.
5. **TP53 mutation frequency (60% of cohort) dominates ME signals.** High-frequency genes are more likely to appear ME with rarer genes simply due to distributional constraints; this is partially but not fully corrected by logistic regression.

---

## Output Files Referenced


| File                              | Contents                                                   |
| --------------------------------- | ---------------------------------------------------------- |
| `interactions_broad_logistic.tsv` | Full logistic regression results (2,450 models)            |
| `interactions_broad_logistic.png` | Tile plot: OR_adj colour-coded, ★ = p_adj < 0.05           |
| `interactions_broad_top50.tsv`    | Unadjusted Fisher results for comparison                   |
| `results_summary.md`              | Broader Phase 1 somatic interaction summary (Fisher-based) |


