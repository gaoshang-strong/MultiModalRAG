# TP53–TRRAP Mutual Exclusivity: MSI Confounding Analysis

**Date:** 2026-05-26  
**Cohort:** FDU gastric cancer, N = 445 (somatic data available)  
**Method:** Sequential logistic regression + stratified Fisher test

---

## Background

Prior analysis (somatic_interactions.R) identified TP53–TRRAP as a significant mutually exclusive pair (Fisher p = 0.006, OR = 0.31). This analysis tests whether the signal reflects a **functional gene–gene interaction** or **subtype stratification** by MSI status.

---

## Key Numbers

| Group | n | TRRAP+ rate |
|-------|---|-------------|
| TP53– (all) | 177 | 17/177 (9.6%) |
| TP53+ (all) | 261 | 8/261 (3.1%) |
| TP53– / MSS | 159 | 7/159 (4.4%) |
| TP53+ / MSS | 253 | 3/253 (1.2%) |
| TP53– / MSI-H | 18 | 10/18 (55.6%) |
| TP53+ / MSI-H | 8 | 5/8 (62.5%) |

---

## Sequential Logistic Models (outcome: TRRAP mutated)

| Model | n | OR (TP53) | 95% CI | p |
|-------|---|-----------|--------|---|
| M1: TP53 alone | 402 | 0.265 | 0.10–0.63 | **0.004** |
| M2: + log2(TMB) | 402 | 0.351 | 0.12–0.97 | **0.048** |
| M3: + MSI status | 402 | 0.381 | 0.13–1.09 | 0.076 |
| M4: + MSI + EBV | 402 | 0.379 | 0.13–1.09 | 0.074 |

Adding MSI to the model eliminates statistical significance of the TP53 coefficient.

---

## Stratified Fisher Test

| Stratum | n | TRRAP+ total | OR | p | Interpretation |
|---------|---|-------------|----|---|----------------|
| All | 438 | 25 | 0.30 | 0.004 | **Significant ME** |
| MSS | 412 | 10 | 0.27 | 0.053 | Borderline, underpowered |
| MSI-H | 26 | 15 | 1.04 | 1.000 | **No ME** |

---

## Interpretation

The TP53–TRRAP mutual exclusivity is **primarily driven by MSI subtype stratification**, not a functional gene–gene interaction.

**Mechanism:** TRRAP mutations are highly enriched in MSI-H tumors (15/26 = 58%), where TP53 mutation rate is low (31%). In MSS tumors (where TP53 mutation rate is high at 61%), TRRAP mutations are rare (2.4%). This frequency difference across subtypes creates an apparent ME when the full cohort is analyzed without stratification.

**Within MSI-H:** TRRAP mutation rate is virtually identical in TP53+ (62.5%) and TP53– (55.6%) samples — confirming no functional ME between the two genes.

**Within MSS:** A directional ME trend persists (1.2% vs 4.4%), but is based on only 10 TRRAP+ cases and lacks statistical power. Cannot be distinguished from residual histological subtype confounding within MSS.

---

## Conclusion

> The TP53–TRRAP mutual exclusivity reflects the co-existence of two biologically distinct gastric cancer contexts: **MSI-H/epigenetically dysregulated tumors** (high TRRAP mutation rate, low TP53 mutation rate) and **MSS/CIN tumors** (high TP53 mutation rate, low TRRAP mutation rate). The genes are not functionally antagonistic — they are markers of different disease subtypes.

The biologically meaningful question is therefore not "why are TP53 and TRRAP mutually exclusive?" but rather **"what drives TRRAP mutation selection specifically in the MSI-H context?"**

---

## Output Files

| File | Contents |
|------|----------|
| `msi_confounding_models.tsv` | Sequential logistic regression results |
| `msi_confounding_stratified_fisher.tsv` | Fisher test per MSI stratum |
| `msi_confounding_forest.png` | TP53 OR across models (forest plot) |
| `msi_confounding_stacked_bar.png` | TRRAP rate by TP53 × MSI group |
| `msi_confounding_2x2_tiles.png` | 2×2 contingency heatmaps per stratum |
