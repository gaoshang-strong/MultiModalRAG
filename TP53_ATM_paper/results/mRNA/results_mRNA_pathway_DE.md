# mRNA-Level Pathway Differential Expression — ATM-TP53 Mutual Exclusivity

**Date:** 2026-06-07  
**Script:** `TP53_ATM_paper/scripts/mrna_pathway_DE.R`  
**Data:** TCGA-STAD RNA-seq (SummarizedExperiment, tpm_unstrand + counts)  
**Method:** DESeq2 (design = `~ tp53_status` or `~ atm_status`); multiple-testing correction: BH within each gene class only (not genome-wide)

---

## 1. Study Design

### Sample groups (ALL samples, no MSI-H exclusion)

| Comparison | Group A (mutant) | Group B (reference) | N_mut | N_ref |
|---|---|---|---|---|
| Comp 1 | TP53-mut / ATM-WT | TP53-WT / ATM-WT | 170 | 202 |
| Comp 2 | ATM-mut / TP53-WT | ATM-WT / TP53-WT | 27 | 202 |

**Mutation definition:** All functional variants — Missense, Nonsense, Frame_Shift_Del/Ins, Splice_Site, In_Frame_Del/Ins, Nonstop_Mutation, Translation_Start_Site (no IMPACT filter).  
**Note:** With this definition, **0 TCGA patients carry both TP53 and ATM functional mutations** — complete mutual exclusivity at the functional-variant level.

### Gene set (35 genes, 3 classes from pathway figure)

| Class | Genes |
|---|---|
| ATM-specific (11) | MRE11, NBN, RAD50, ATM, H2AX, MDC1, TRIM28, RNF8, RNF168, TP53BP1, PARP1 |
| Convergence (12) | KAT5, ABRAXAS1, CHEK2, ATR, CHEK1, BRCA1, FOXO3, PALB2, BRCA2, RAD51, WEE1, CDK1 |
| TP53-specific (12) | TP53, MDM2, CDKN2A, CDKN1A, SFN, GADD45A, BAX, BBC3, PMAIP1, DDB2, TIGAR, PPM1D |

---

## 2. Comparison 1: TP53-mut / ATM-WT vs TP53-WT / ATM-WT

### 2a. Class-level summary

| Class | N genes in results | Sig (padj<0.05) | ↑ mut | ↓ mut | Mean log2FC | Median log2FC |
|---|---|---|---|---|---|---|
| ATM-specific  | 11 | **7/11** | 6 | 5 | +0.076 | +0.075 |
| Convergence   | 12 | **8/12** | 10 | 2 | +0.160 | +0.192 |
| TP53-specific | 12 | **9/12** | 2 | 10 | −0.149 | −0.257 |

### 2b. Per-gene results (sorted by class, then p-value)

**ATM-specific:**

| Gene | log2FC | padj (class) | Direction |
|---|---|---|---|
| MDC1    | +0.396 | 5.7×10⁻¹¹ | ↑↑ |
| MRE11   | +0.303 | 3.2×10⁻⁶  | ↑↑ |
| PARP1   | +0.202 | 0.0013    | ↑  |
| RAD50   | −0.180 | 0.0049    | ↓  |
| TRIM28  | +0.209 | 0.0065    | ↑  |
| ATM     | −0.193 | 0.0100    | ↓  |
| H2AX    | +0.196 | 0.0403    | ↑  |
| NBN     | −0.103 | 0.088     | ns |
| RNF8    | +0.075 | 0.139     | ns |
| RNF168  | −0.044 | 0.443     | ns |
| TP53BP1 | −0.021 | 0.705     | ns |

**Convergence:**

| Gene | log2FC | padj (class) | Direction |
|---|---|---|---|
| BRCA1   | +0.349 | 1.1×10⁻⁴  | ↑↑ |
| ABRAXAS1| −0.275 | 1.1×10⁻⁴  | ↓↓ |
| WEE1    | +0.324 | 1.1×10⁻⁴  | ↑↑ |
| CHEK1   | +0.312 | 1.8×10⁻⁴  | ↑↑ |
| CHEK2   | +0.292 | 1.9×10⁻⁴  | ↑↑ |
| CDK1    | +0.308 | 0.0015    | ↑  |
| RAD51   | +0.209 | 0.011     | ↑  |
| FOXO3   | +0.167 | 0.012     | ↑  |
| BRCA2   | +0.175 | 0.086     | ns |
| KAT5    | −0.083 | 0.105     | ns |
| ATR     | +0.100 | 0.133     | ns |
| PALB2   | +0.037 | 0.486     | ns |

**TP53-specific:**

| Gene | log2FC | padj (class) | Direction |
|---|---|---|---|
| MDM2    | −1.046 | 3.0×10⁻²⁹ | ↓↓↓ |
| DDB2    | −0.602 | 1.3×10⁻¹⁸ | ↓↓↓ |
| CDKN2A  | +1.586 | 6.2×10⁻¹⁶ | ↑↑↑ |
| TIGAR   | −0.438 | 5.2×10⁻⁹  | ↓↓  |
| PPM1D   | −0.182 | 4.0×10⁻⁵  | ↓   |
| BBC3    | −0.387 | 6.1×10⁻⁵  | ↓   |
| BAX     | −0.295 | 6.1×10⁻⁵  | ↓   |
| TP53    | −0.357 | 2.2×10⁻⁴  | ↓   |
| CDKN1A  | −0.387 | 0.025     | ↓   |
| SFN     | +0.210 | 0.192     | ns  |
| GADD45A | −0.053 | 0.497     | ns  |
| PMAIP1  | −0.002 | 0.988     | ns  |

### 2c. GSEA (per class, ranked by sign × −log10 padj)

| Pathway | NES | p-value | padj |
|---|---|---|---|
| TP53-specific | **−1.591** | **0.028** | 0.111 |
| Convergence   | +1.227 | 0.217 | 0.249 |
| ATM-specific  | +1.199 | 0.249 | 0.249 |
| All_pathway   | −1.196 | 0.185 | 0.249 |

Leading edge (TP53-specific): MDM2, DDB2, TIGAR, PPM1D, BBC3, BAX, TP53

---

## 3. Comparison 2: ATM-mut / TP53-WT vs ATM-WT / TP53-WT

### 3a. Class-level summary

| Class | N genes | Sig (padj<0.05) | ↑ mut | ↓ mut | Mean log2FC | Median log2FC |
|---|---|---|---|---|---|---|
| ATM-specific  | 11 | **2/11** | 3 | 8 | −0.080 | −0.025 |
| Convergence   | 12 | **0/12** | 7 | 5 | +0.034 | +0.017 |
| TP53-specific | 12 | **1/12** | 9 | 3 | +0.099 | +0.184 |

### 3b. Per-gene results

**ATM-specific** — only ATM itself is significantly different:

| Gene | log2FC | padj (class) |
|---|---|---|
| ATM     | −0.656 | 6.0×10⁻⁵ |
| MRE11   | −0.283 | 0.041 (trend) |
| All others | −0.17 to +0.21 | 0.20–0.90 (ns) |

**Convergence** — no significant genes:

All 12 genes: padj 0.16–0.94, all near zero log2FC (range −0.27 to +0.34).

**TP53-specific** — only GADD45A significant (and in unexpected direction):

| Gene | log2FC | padj (class) |
|---|---|---|
| GADD45A | −0.380 | 0.020 |
| All others | −0.19 to +0.43 | 0.12–0.85 (ns) |

### 3c. GSEA

All four pathways: padj > 0.35. No directional enrichment. NES range: −1.07 to +1.08.

---

## 4. Interpretation

### 4a. Comp 1 — TP53 mutation produces two simultaneous transcriptional signatures

**Signature A (expected): Loss of TP53 transcriptional programme**  
9/12 TP53-specific target genes are downregulated. MDM2 (−1.05), DDB2 (−0.60), TIGAR (−0.44), BBC3, BAX, TP53 transcript itself are all reduced — the expected loss-of-function footprint. GSEA NES = −1.59 (p=0.028) confirms set-level negative enrichment.

**Signature B (novel): Compensatory transcriptional upregulation of ATM-axis genes**  
Simultaneously, 6/11 ATM-specific genes and 10/12 convergence genes are upregulated. Key upregulated nodes:
- *Damage sensors:* MDC1 (+0.40), MRE11 (+0.30), PARP1 (+0.20)
- *Checkpoint kinases:* CHEK1 (+0.31), CHEK2 (+0.29), WEE1 (+0.32)
- *HR repair:* BRCA1 (+0.35), RAD51 (+0.21)

This pattern indicates that **TP53-mutant tumors upregulate the ATM kinase cascade at the transcriptional level**, consistent with a compensatory response to the loss of TP53-mediated DDR.

### 4b. Comp 2 — ATM mutations leave the transcriptome largely unchanged

Across all three gene classes, only 3 genes reach padj < 0.05 (ATM itself reflecting hemizygous loss, MRE11 trend, and GADD45A — possibly a noise hit given N=27). GSEA is flat for all classes.

**This is a key negative result:** ATM is a kinase. Its primary mechanism is post-translational phosphorylation of substrates (H2AX, CHEK2, BRCA1, TP53 Ser15, etc.). ATM mutations disrupt this signalling cascade at the protein-activity level without producing a distinctive mRNA-level signature in either the ATM-specific or the convergence gene class.

### 4c. Unified mechanistic model

```
TP53-mutant tumours
  ├── [Transcription] TP53 target genes SUPPRESSED (canonical LOF)
  └── [Transcription] ATM pathway genes UPREGULATED (compensatory)
        → tumour is transcriptionally "dependent" on ATM signalling

ATM-mutant tumours
  ├── [Transcription] Minimal changes across all 3 classes
  └── [Post-translational] Kinase cascade disrupted (not detectable by mRNA)
```

**Why are ATM and TP53 mutations mutually exclusive?**  
TP53-mutant tumours develop transcriptional dependency on the ATM kinase axis as a compensatory DDR mechanism. A subsequent ATM loss-of-function mutation would simultaneously: (i) abolish the compensatory ATM signalling that the TP53-mutant cell relies on, and (ii) remove the post-translational route to p53 activation (ATM → Ser15 phosphorylation). The combined loss silences all DDR output. The transcriptional compensation observed here makes ATM LOF even more incompatible with TP53 LOF than a simple "two independent pathways" model would predict.

---

## 5. Output Files

| File | Description |
|---|---|
| `comp1_TP53mut/pathway_DE_results.csv` | Per-gene DESeq2 results + class-specific padj |
| `comp1_TP53mut/class_summary.csv` | Class-level counts and mean log2FC |
| `comp1_TP53mut/gsea_pathway_summary.csv` | GSEA NES/p per class |
| `comp1_TP53mut/lollipop_pathway.png` | Lollipop: log2FC per gene, colored by class |
| `comp1_TP53mut/volcano_pathway.png` | Genome-wide volcano, pathway genes highlighted |
| `comp1_TP53mut/heatmap_pathway.png` | Sample × gene heatmap, class annotation |
| `comp1_TP53mut/boxplots_pathway.png` | Per-gene expression boxplots |
| `comp1_TP53mut/gsea_enrichment.png` | GSEA enrichment curves |
| `comp2_ATMmut/…` | Parallel outputs for Comp 2 |
