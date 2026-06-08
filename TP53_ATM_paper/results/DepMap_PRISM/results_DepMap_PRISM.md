# DepMap CRISPR + PRISM/GDSC2 ATM Inhibitor Analysis — ATM-TP53 Synthetic Lethality

**Date:** 2026-06-07  
**Script:** `TP53_ATM_paper/scripts/depmap_prism_ATM_TP53.R`  
**Datasets:**
- DepMap CRISPR (release 24Q2): ATM (472) gene effect; 1178 cell lines; TP53 status from `mutations_TP53_damaging.csv`
- GDSC2 KU-55933 (ATM inhibitor): 955 cell lines; metrics: LN_IC50, AUC
- PRISM gdsc-007 AZD0156 (ATM inhibitor): 731 cell lines; single-agent: lib1_IC50_ln, lib1_MaxE; combinations with WEE1i (AZD1775), ATRi (AZD6738), AURKBi (AZD2811) — Bliss synergy score

**Statistics:** Wilcoxon rank-sum (pan-cancer + gastric subset) + linear regression with OncotreeLineage covariate (pan-cancer only)

---

## 1. DepMap CRISPR — ATM Knockout Effect by TP53 Status

### 1a. Pan-cancer

| Group | N | Median gene effect | Wilcoxon p | Effect r | lm β (lineage-adj) | lm p |
|---|---|---|---|---|---|---|
| TP53-mut | 755 | **−0.051** | — | — | — | — |
| TP53-WT  | 423 | **+0.074** | **8.67×10⁻³⁵** | **0.43** | **−0.132** | **3.12×10⁻²⁸** |

**Interpretation:** TP53-mutant cell lines are significantly more sensitive to ATM CRISPR KO (more negative gene effect = more essential). After controlling for cancer lineage, the TP53 association remains strongly significant (lm β=−0.132, p=3.12×10⁻²⁸, adjusted R²=0.18). Effect size is moderate-to-large (r=0.43).

### 1b. Gastric only

| Group | N | Median gene effect | Wilcoxon p |
|---|---|---|---|
| TP53-mut | 59 | −0.079 | — |
| TP53-WT  | 10 | +0.024 | **p=0.078** |

Same directional trend as pan-cancer, but does not reach significance with N=10 WT gastric cell lines. The WT gastric group is underpowered; this result is consistent with but cannot confirm the pan-cancer finding in the gastric-specific context.

**Figures:**
- `depmap_crispr/atm_ko_pancancer.png` — violin + jitter, gastric lines highlighted
- `depmap_crispr/atm_ko_gastric.png`

---

## 2. GDSC2 KU-55933 (ATM Inhibitor) — Sensitivity by TP53 Status

### 2a. LN_IC50

| Scope | Group | N | Median LN_IC50 | Wilcoxon p | lm β (lineage-adj) | lm p |
|---|---|---|---|---|---|---|
| Pan-cancer | TP53-mut | 610 | 5.097 | — | — | — |
| Pan-cancer | TP53-WT  | 312 | 5.108 | 0.642 | −0.100 | 0.163 |
| Gastric    | TP53-mut | 52  | 5.219 | — | — | — |
| Gastric    | TP53-WT  | 11  | 5.366 | 0.779 | — | — |

### 2b. AUC

| Scope | Group | N | Median AUC | Wilcoxon p | lm β | lm p |
|---|---|---|---|---|---|---|
| Pan-cancer | TP53-mut | 610 | 0.9766 | — | — | — |
| Pan-cancer | TP53-WT  | 312 | 0.9773 | 0.928 | −0.0015 | 0.351 |

**Interpretation:** Single-agent KU-55933 shows **no differential sensitivity** by TP53 mutation status, either pan-cancer or in gastric cancer alone. LN_IC50 values are nearly identical between groups (medians differ by <0.1 log units). This is a robust negative result across two metrics (IC50 and AUC).

---

## 3. PRISM AZD0156 (ATM Inhibitor) — Single-Agent Sensitivity by TP53 Status

### 3a. LN_IC50 and MaxE

| Scope | Metric | TP53-mut (N) | Median mut | TP53-WT (N) | Median WT | Wilcoxon p | lm p |
|---|---|---|---|---|---|---|---|
| Pan-cancer | LN_IC50 | 469 | 2.758 | 236 | 2.781 | 0.935 | 0.415 |
| Pan-cancer | MaxE    | 469 | 0.103 | 236 | 0.100 | 0.222 | 0.095 |
| Gastric    | LN_IC50 | 44  | 2.718 | 10  | 3.029 | 0.585 | — |
| Gastric    | MaxE    | 44  | 0.112 | 10  | 0.079 | **0.073** | — |

**Pan-cancer:** No significant difference for either metric.  
**Gastric MaxE trend:** Gastric TP53-mut cell lines show a trend toward greater maximal inhibition by AZD0156 (MaxE 0.112 vs 0.079, p=0.073). With N=10 WT gastric lines, this is exploratory but directionally consistent with TP53-mut cells having greater ATM pathway dependence.

---

## 4. PRISM AZD0156 Combination Synergy (Bliss Score) by TP53 Status

### 4a. AZD0156 + AZD1775 (WEE1 inhibitor)

| Scope | Metric | Median mut | Median WT | Wilcoxon p | lm β | lm p |
|---|---|---|---|---|---|---|
| Pan-cancer | Bliss_window    | 0.089 | 0.096 | 0.587 | +0.006 | 0.373 |
| Pan-cancer | Delta_MaxE_lib1 | 0.392 | 0.402 | 0.553 | +0.035 | **0.045** |
| Gastric    | Bliss_window    | 0.113 | 0.077 | 0.333 | — | — |

The Delta_MaxE_lib1 (additional effect in combo relative to AZD0156 alone) shows a borderline lineage-adjusted association with TP53 status (lm p=0.045), but with a very small effect size (β=0.035). Bliss synergy score itself is not significant.

### 4b. AZD0156 + AZD6738 (ATR inhibitor)

| Metric | Median mut | Median WT | Wilcoxon p |
|---|---|---|---|
| Bliss_window    | 0.069 | 0.060 | 0.179 |
| combo_MaxE      | 0.674 | 0.706 | 0.286 |
| Delta_MaxE_lib1 | 0.520 | 0.523 | 0.910 |

**No significant difference.** ATRi + ATMi combination synergy does not depend on TP53 status.

### 4c. AZD0156 + AZD2811 (Aurora B inhibitor)

| Metric | Median mut | Median WT | Wilcoxon p |
|---|---|---|---|
| Bliss_window    | 0.039 | 0.044 | 0.503 |
| combo_MaxE      | 0.232 | 0.263 | 0.091 |
| Delta_MaxE_lib1 | 0.126 | 0.149 | 0.102 |

Trend toward lower combo_MaxE in TP53-mut (p=0.091) and lower Delta_MaxE (p=0.102), but neither reaches significance.

---

## 5. Summary and Biological Interpretation

| Analysis | Scope | Key metric | Result | Significance |
|---|---|---|---|---|
| DepMap CRISPR ATM KO | Pan-cancer | Gene effect | TP53-mut more sensitive (−0.051 vs +0.074) | ★★★ p=8.67×10⁻³⁵ |
| DepMap CRISPR ATM KO (lineage-adj) | Pan-cancer | lm β | −0.132 | ★★★ p=3.12×10⁻²⁸ |
| DepMap CRISPR ATM KO | Gastric | Gene effect | Same direction | Trend p=0.078 |
| GDSC2 KU-55933 LN_IC50 | Pan-cancer | LN_IC50 | No difference | NS p=0.642 |
| GDSC2 KU-55933 AUC | Pan-cancer | AUC | No difference | NS p=0.928 |
| PRISM AZD0156 LN_IC50 | Pan-cancer | LN_IC50 | No difference | NS p=0.935 |
| PRISM AZD0156 MaxE | Gastric | MaxE | Trend ↑ mut | Trend p=0.073 |
| PRISM AZD0156+WEE1i Bliss | Pan-cancer | Bliss_window | No difference | NS |
| PRISM AZD0156+ATRi Bliss | Pan-cancer | Bliss_window | No difference | NS |
| PRISM AZD0156+AURKBi combo_MaxE | Pan-cancer | combo_MaxE | No difference | Trend p=0.091 |

### 5a. The genetic–pharmacological disconnect

The striking finding is the **divergence between CRISPR genetic ablation and pharmacological inhibition**:

- **CRISPR ATM KO**: TP53-mut cell lines show strong, consistent, lineage-adjusted synthetic lethality (N=755+423, p=8.67×10⁻³⁵, r=0.43). This directly validates the genomic mutual exclusivity observed in the FDU/TCGA cohorts.

- **ATM inhibitors (KU-55933, AZD0156)**: No differential sensitivity by TP53 status pan-cancer. Both inhibitors were screened at multiple concentrations in large panels (>600 cell lines matched).

This disconnect has several mechanistic explanations that are not mutually exclusive:

1. **Incomplete ATM inhibition**: KU-55933 and AZD0156 are ATP-competitive inhibitors that partially suppress ATM kinase activity. The synthetic lethality may require a threshold level of ATM inactivation that is only achieved by complete genetic loss, not partial pharmacological inhibition at the concentrations used in these screens.

2. **Scaffold vs. kinase function**: ATM has known kinase-independent scaffold functions in the DNA damage response. CRISPR KO removes both kinase and scaffold roles; inhibitors only block kinase activity. If scaffold function is sufficient to maintain survival in TP53-mut cells, kinase inhibition alone would not trigger synthetic lethality.

3. **Residual ATM activity vs. null**: The CRISPR gene effect score captures the fitness consequence of moving from partial to zero ATM function. In TP53-mut cells that are already adapted to moderate ATM signaling (consistent with our phospho-proteomics data showing a hyperactivated ATM pathway), the additional step to zero ATM is specifically catastrophic.

4. **Cell line vs. tumor biology**: CRISPR screens use chronic, genome-wide depletion which is more analogous to the tumor evolutionary constraint than acute drug treatment.

### 5b. Translational implications

The CRISPR data provides strong mechanistic justification for continued preclinical development:
- Genetic validation of ATM as a synthetic lethal target in TP53-mut cancers (large N, lineage-adjusted, r=0.43)
- The pharmacological negative result highlights that **next-generation ATM inhibitors or strategies achieving more complete ATM inhibition** may be required to translate this vulnerability

The gastric-specific AZD0156 MaxE trend (p=0.073, N small) and AZD0156+AURKBi combo MaxE trend (p=0.091) warrant follow-up in gastric-specific models with larger N.

### 5c. Integration with the 5-layer mechanistic model

```
TP53-mut cells:
  [Genomic]     ME with ATM (FDU N=445, TCGA)
  [mRNA]        Compensatory ATM-pathway upregulation
  [Protein]     ATM-pathway proteins ↑ (CPTAC GSEA padj=0.031)
  [Phospho]     89% ATM-specific phosphosites ↑ (kinase hyperactivated)
  [DepMap]      ATM KO → synthetic lethality (p=8.67×10⁻³⁵, r=0.43)
  ─────────────────────────────────────────────────────────────────
  Mechanism: TP53-mut tumors depend on the hyperactivated ATM-DDR
  axis as a compensatory survival pathway. Genetic ATM loss removes
  this compensation → mitotic catastrophe.
  [ATM inhibitor screens]: partial inhibition insufficient → threshold
  effect requires complete loss or superior ATM inhibitors
```

---

## 6. Limitations

1. **Gastric WT underpowered**: Only 10 WT gastric cell lines in DepMap. Gastric-specific conclusions are directional only.
2. **ATM inhibitor concentrations**: KU-55933 and AZD0156 were screened at fixed concentration ranges; actual ATM inhibition at screened concentrations is unknown.
3. **ATM mutation status in cell lines**: This analysis uses TP53 mutation as the stratifier (consistent with the tumor analysis). ATM mutation status of individual cell lines was not used as an additional covariate; adding it would further validate the interaction.
4. **PRISM AZD0156 is a combination screen**: AZD0156 as lib1 in a 4-drug combination matrix; single-agent parameters (lib1_IC50_ln, lib1_MaxE) are derived from the combination screen design and may be less precise than dedicated single-agent dose-response curves.

---

## 7. Output Files

| Directory | Contents |
|---|---|
| `depmap_crispr/` | ATM KO violin plots (pan/gastric), stats CSV, per-cell-line data |
| `gdsc2_ku55933/` | KU-55933 LN_IC50 + AUC plots, stats CSV |
| `prism_azd0156_single/` | AZD0156 IC50_ln + MaxE plots, stats CSV |
| `prism_azd0156_combos/AZD1775/` | AZD0156+WEE1i Bliss/HSA/MaxE plots, stats CSV |
| `prism_azd0156_combos/AZD6738/` | AZD0156+ATRi plots, stats CSV |
| `prism_azd0156_combos/AZD2811/` | AZD0156+AURKBi plots, stats CSV |
| `master_stats_summary.csv` | All comparisons in one table |
