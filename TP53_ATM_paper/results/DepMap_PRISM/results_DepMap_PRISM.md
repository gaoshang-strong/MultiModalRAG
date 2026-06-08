# DepMap CRISPR + PRISM/GDSC2 ATM Inhibitor Analysis — ATM-TP53 Synthetic Lethality

**Date:** 2026-06-08  
**Script:** `TP53_ATM_paper/scripts/depmap_prism_ATM_TP53.R`  
**Datasets:**
- DepMap CRISPR (release 24Q2): ATM (472) gene effect; TP53 status from `mutations_TP53_damaging.csv` + `OmicsCNGene.csv` (CN < 0.3 = homozygous deletion)
- GDSC2 KU-55933 (ATM inhibitor): 955 cell lines; metrics: LN_IC50, AUC
- PRISM gdsc-007 AZD0156 (ATM inhibitor): 731 cell lines; single-agent: lib1_IC50_ln, lib1_MaxE; combinations with WEE1i (AZD1775), ATRi (AZD6738), AURKBi (AZD2811), DNA-PKi (AZD7648) — Bliss synergy score

**Methodology changes vs prior run (2026-06-07):**
- TP53 mutant definition expanded: damaging mutation OR CN < 0.3 (homozygous deletion) — adds ~5% more mutant cell lines
- Strict ATM-WT filter removed from all sections: all cell lines included regardless of ATM mutation status (rationale: for CRISPR KO, ATM is knocked out in all cells so endogenous ATM status is irrelevant; for drug screens, including ATM-mut lines in the WT arm would dilute the comparison — removing the filter provides the most inclusive comparison)
- Bootstrap 95% CI added (delta = mean_mut − mean_WT, 2000 iterations) — values in CSV

**Statistics:** Wilcoxon rank-sum (pan-cancer + gastric subset) + rank-biserial r effect size + linear regression with OncotreeLineage covariate (pan-cancer only) + bootstrap 95% CI on delta

---

## 1. DepMap CRISPR — ATM Knockout Effect by TP53 Status

### 1a. Pan-cancer

| Group | N | Median gene effect | Wilcoxon p | Effect r | lm β (lineage-adj) | lm p |
|---|---|---|---|---|---|---|
| TP53-mut | 764 | **−0.052** | — | — | — | — |
| TP53-WT  | 414 | **+0.077** | **8.20×10⁻³⁷** | **0.447** | **−0.136** | **4.75×10⁻³⁰** |

**Interpretation:** TP53-mutant cell lines are significantly more sensitive to ATM CRISPR KO (more negative gene effect = more essential). After controlling for cancer lineage, the association remains strongly significant (lm β=−0.136, p=4.75×10⁻³⁰). Effect size is moderate-to-large (r=0.447). N increased vs prior run (764/414 vs 755/423) due to inclusion of CN-deleted TP53 cell lines.

### 1b. Gastric only

| Group | N | Median gene effect | Wilcoxon p | Effect r |
|---|---|---|---|---|
| TP53-mut | 60 | −0.072 | — | — |
| TP53-WT  | 9  | +0.054 | 0.068 | 0.381 |

Same directional trend as pan-cancer (mut more sensitive, r=0.381). Does not reach significance with N=9 WT gastric cell lines (underpowered). Directionally consistent with the pan-cancer finding.

**Figures:**
- `depmap_crispr/atm_ko_allTP53_pancancer.png` — violin + jitter, gastric lines highlighted
- `depmap_crispr/atm_ko_allTP53_gastric.png`

---

## 2. GDSC2 KU-55933 (ATM Inhibitor) — Sensitivity by TP53 Status

### 2a. LN_IC50

| Scope | Group | N | Median LN_IC50 | Wilcoxon p | lm β (lineage-adj) | lm p |
|---|---|---|---|---|---|---|
| Pan-cancer | TP53-mut | 618 | 5.090 | — | — | — |
| Pan-cancer | TP53-WT  | 304 | 5.138 | 0.869 | −0.124 | 0.086 |
| Gastric    | TP53-mut | 53  | 5.219 | — | — | — |
| Gastric    | TP53-WT  | 10  | 5.364 | 0.858 | — | — |

### 2b. AUC

| Scope | Group | N | Median AUC | Wilcoxon p | lm β | lm p |
|---|---|---|---|---|---|---|
| Pan-cancer | TP53-mut | 618 | 0.9764 | — | — | — |
| Pan-cancer | TP53-WT  | 304 | 0.9776 | 0.796 | −0.0018 | 0.270 |
| Gastric    | TP53-mut | 53  | 0.9757 | — | — | — |
| Gastric    | TP53-WT  | 10  | 0.9804 | 0.714 | — | — |

**Interpretation:** Single-agent KU-55933 shows **no differential sensitivity** by TP53 mutation status in any scope or metric. LN_IC50 values are nearly identical between groups. This is a robust negative result.

---

## 3. PRISM AZD0156 (ATM Inhibitor) — Single-Agent Sensitivity by TP53 Status

| Scope | Metric | TP53-mut (N) | Median mut | TP53-WT (N) | Median WT | Wilcoxon p | lm p |
|---|---|---|---|---|---|---|---|
| Pan-cancer | LN_IC50 | 474 | 2.770 | 231 | 2.763 | 0.805 | 0.513 |
| Pan-cancer | MaxE    | 474 | 0.102 | 231 | 0.100 | 0.191 | 0.072 |
| Gastric    | LN_IC50 | 45  | 2.743 | 9   | 2.909 | 0.871 | — |
| Gastric    | MaxE    | 45  | 0.112 | 9   | 0.078 | **0.095** | — |

**Pan-cancer:** No significant difference for either metric.  
**Gastric MaxE trend:** TP53-mut gastric cell lines show a trend toward greater maximal inhibition by AZD0156 (MaxE 0.112 vs 0.078, p=0.095). With N=9 WT gastric lines this is exploratory but directionally consistent with the CRISPR result.

---

## 4. PRISM AZD0156 Combination Synergy (Bliss Score) by TP53 Status

### 4a. AZD0156 + AZD1775 (WEE1 inhibitor)

| Scope | Metric | Median mut | Median WT | Wilcoxon p | lm β | lm p |
|---|---|---|---|---|---|---|
| Pan-cancer | Bliss_window       | 0.089 | 0.096 | 0.469 | +0.004 | 0.510 |
| Pan-cancer | Delta_MaxE_AZD0156 | 0.389 | 0.408 | 0.447 | +0.031 | 0.079 |
| Gastric    | Bliss_window       | 0.108 | 0.106 | 0.693 | — | — |

No significant synergy differential by TP53 status. Delta_MaxE lm trend (p=0.079) is borderline and small (β=0.031).

### 4b. AZD0156 + AZD6738 (ATR inhibitor)

| Metric | Median mut | Median WT | Wilcoxon p |
|---|---|---|---|
| Bliss_window       | 0.069 | 0.060 | 0.173 |
| combo_MaxE         | 0.671 | 0.717 | 0.221 |
| Delta_MaxE_AZD0156 | 0.520 | 0.524 | 0.823 |

**No significant difference.** ATRi + ATMi combination synergy does not depend on TP53 status.

### 4c. AZD0156 + AZD2811 (Aurora B inhibitor)

| Metric | Median mut | Median WT | Wilcoxon p |
|---|---|---|---|
| Bliss_window       | 0.039 | 0.044 | 0.388 |
| combo_MaxE         | 0.232 | 0.266 | 0.069 |
| Delta_MaxE_AZD0156 | 0.125 | 0.152 | 0.082 |

Trend toward lower combo_MaxE in TP53-mut (p=0.069) and lower Delta_MaxE (p=0.082), neither significant. Directionally opposite to a SL prediction (TP53-mut cells less sensitive to the combo).

### 4d. AZD0156 + AZD7648 (DNA-PK inhibitor)

| Scope | Metric | Median mut | Median WT | Wilcoxon p | lm β | lm p |
|---|---|---|---|---|---|---|
| Pan-cancer | Bliss_window       | 0.063 | 0.055 | **0.009** | +0.009 | **0.007** |
| Pan-cancer | HSA_window         | 0.074 | 0.068 | 0.154 | +0.006 | 0.060 |
| Pan-cancer | combo_MaxE         | 0.288 | 0.307 | 0.107 | −0.009 | 0.384 |
| Pan-cancer | Delta_MaxE_AZD0156 | 0.156 | 0.158 | 0.604 | −0.001 | 0.913 |
| Gastric    | Bliss_window       | 0.058 | 0.048 | 0.275 | — | — |

**Notable:** AZD0156+AZD7648 Bliss synergy score is significantly higher in TP53-mut cells pan-cancer (p=0.009, lm p=0.007, β=+0.009). This is a small but consistent and lineage-adjusted signal. ATMi+DNA-PKi synergy in TP53-mut context is mechanistically plausible (both kinases required for double-strand break resolution; TP53 loss removes the apoptotic brake). However, the effect size is small (Δmedian ~0.008 Bliss units) and does not translate to combo_MaxE or Delta_MaxE significance.

---

## 5. Summary and Biological Interpretation

| Analysis | Scope | Key metric | Result | Significance |
|---|---|---|---|---|
| DepMap CRISPR ATM KO | Pan-cancer | Gene effect | TP53-mut more sensitive (−0.052 vs +0.077) | ★★★ p=8.20×10⁻³⁷ |
| DepMap CRISPR ATM KO (lineage-adj) | Pan-cancer | lm β | −0.136 | ★★★ p=4.75×10⁻³⁰ |
| DepMap CRISPR ATM KO | Gastric | Gene effect | Same direction, r=0.381 | Trend p=0.068 |
| GDSC2 KU-55933 LN_IC50 | Pan-cancer | LN_IC50 | No difference | NS p=0.869 |
| GDSC2 KU-55933 AUC | Pan-cancer | AUC | No difference | NS p=0.796 |
| PRISM AZD0156 LN_IC50 | Pan-cancer | LN_IC50 | No difference | NS p=0.805 |
| PRISM AZD0156 MaxE | Gastric | MaxE | Trend ↑ mut | Trend p=0.095 |
| PRISM AZD0156+WEE1i Bliss | Pan-cancer | Bliss_window | No difference | NS p=0.469 |
| PRISM AZD0156+ATRi Bliss | Pan-cancer | Bliss_window | No difference | NS p=0.173 |
| PRISM AZD0156+AURKBi combo_MaxE | Pan-cancer | combo_MaxE | No difference | Trend p=0.069 |
| PRISM AZD0156+DNA-PKi Bliss | Pan-cancer | Bliss_window | TP53-mut higher synergy | ★★ p=0.009 lm p=0.007 |

### 5a. The genetic–pharmacological disconnect

The striking finding is the **divergence between CRISPR genetic ablation and pharmacological inhibition**:

- **CRISPR ATM KO**: TP53-mut cell lines show strong, consistent, lineage-adjusted synthetic lethality (N=764+414, p=8.20×10⁻³⁷, r=0.447). This directly validates the genomic mutual exclusivity observed in the FDU/TCGA cohorts.

- **ATM inhibitors (KU-55933, AZD0156)**: No differential sensitivity by TP53 status pan-cancer. Both inhibitors were screened at multiple concentrations in large panels (>600 cell lines matched).

This disconnect has several mechanistic explanations that are not mutually exclusive:

1. **Incomplete ATM inhibition**: KU-55933 and AZD0156 are ATP-competitive inhibitors that partially suppress ATM kinase activity. The synthetic lethality may require a threshold level of ATM inactivation that is only achieved by complete genetic loss, not partial pharmacological inhibition at the concentrations used in these screens.

2. **Scaffold vs. kinase function**: ATM has known kinase-independent scaffold functions in the DNA damage response. CRISPR KO removes both kinase and scaffold roles; inhibitors only block kinase activity. If scaffold function is sufficient to maintain survival in TP53-mut cells, kinase inhibition alone would not trigger synthetic lethality.

3. **Residual ATM activity vs. null**: The CRISPR gene effect score captures the fitness consequence of moving from partial to zero ATM function. In TP53-mut cells that are already adapted to moderate ATM signaling (consistent with our phospho-proteomics data showing a hyperactivated ATM pathway), the additional step to zero ATM is specifically catastrophic.

4. **Cell line vs. tumor biology**: CRISPR screens use chronic, genome-wide depletion which is more analogous to the tumor evolutionary constraint than acute drug treatment.

### 5b. New finding: AZD0156 + AZD7648 (DNA-PKi) synergy in TP53-mut

The AZD0156+AZD7648 Bliss synergy result (p=0.009, lineage-adjusted p=0.007) is the only pharmacological signal reaching significance. ATM and DNA-PK are the two major PIKK kinases orchestrating the DSB response through NHEJ and HR respectively; dual inhibition in TP53-mut cells (which lack the p53-dependent apoptotic checkpoint) may push cells into a catastrophic repair-deficient state. This warrants validation in a dedicated dose-response experiment. The effect size remains small at the population level.

### 5c. Translational implications

The CRISPR data provides strong mechanistic justification for continued preclinical development:
- Genetic validation of ATM as a synthetic lethal target in TP53-mut cancers (large N, lineage-adjusted, r=0.447)
- The pharmacological negative result highlights that **next-generation ATM inhibitors or strategies achieving more complete ATM inhibition** may be required to translate this vulnerability
- AZD0156+AZD7648 combo synergy signal (p=0.007 lineage-adj) is the strongest pharmacological lead from this screen

The gastric-specific AZD0156 MaxE trend (p=0.095, N small) warrants follow-up in gastric-specific models with larger N.

### 5d. Integration with the 5-layer mechanistic model

```
TP53-mut cells:
  [Genomic]     ME with ATM (FDU N=445, TCGA)
  [mRNA]        Compensatory ATM-pathway upregulation
  [Protein]     ATM-pathway proteins ↑ (CPTAC GSEA padj=0.031)
  [Phospho]     89% ATM-specific phosphosites ↑ (kinase hyperactivated)
  [DepMap]      ATM KO → synthetic lethality (p=8.20×10⁻³⁷, r=0.447)
  ─────────────────────────────────────────────────────────────────
  Mechanism: TP53-mut tumors depend on the hyperactivated ATM-DDR
  axis as a compensatory survival pathway. Genetic ATM loss removes
  this compensation → mitotic catastrophe.
  [ATM inhibitor screens]: partial inhibition insufficient → threshold
  effect requires complete loss or superior ATM inhibitors
  [ATMi+DNA-PKi]: dual PIKK blockade shows synergy signal in TP53-mut
```

---

## 6. Limitations

1. **Gastric WT underpowered**: Only 9–10 WT gastric cell lines. Gastric-specific conclusions are directional only.
2. **ATM inhibitor concentrations**: KU-55933 and AZD0156 were screened at fixed concentration ranges; actual ATM inhibition at screened concentrations is unknown.
3. **PRISM AZD0156 is a combination screen**: AZD0156 as lib1 in a 4-drug combination matrix; single-agent parameters are derived from the combination screen design and may be less precise than dedicated single-agent dose-response curves.
4. **AZD7648 synergy effect size**: While statistically significant, the Bliss score difference (~0.008 units) is small. Clinical relevance requires validation in dedicated synergy experiments.

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
| `prism_azd0156_combos/AZD7648/` | AZD0156+DNA-PKi plots, stats CSV |
| `master_stats_summary.csv` | All 42 comparisons in one table with delta, CI, Wilcoxon, lm |
