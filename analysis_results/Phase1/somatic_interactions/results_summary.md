# Phase 1 Somatic Interaction Analysis

**Method:** Pairwise Fisher's exact test via `maftools::somaticInteractions()`  
**Input:** `analysis_results/Phase0/enrichment/somt_filtered.maf`  
**Date:** 2026-05-21  
**Cohort:** N = 445 samples

---

## Overview

Two independent runs were performed to characterise pairwise mutation co-occurrence and mutual exclusivity across the gastric cancer cohort.

| Run | Gene Set | Total Pairs Tested | Sig. Co-occurrences (p<0.05) | Sig. Mutual Exclusivity (p<0.05) |
|-----|----------|--------------------|------------------------------|-----------------------------------|
| A — Broad | Top 50 genes by mutation frequency | 1,225 | 647 | 11 |
| B — Focused | 4 seeds + multi-seed BioGRID neighbors (≥2 seed connections) | 2,278 | 479 | 5 |

Significance thresholds used in maftools plot rendering: p < 0.05 (nominal) and p < 0.01 (highlighted).

---

## Run A — Broad (Top 50 Genes)

### Top Co-occurring Pairs

| Gene 1 | Gene 2 | p-value | Odds Ratio | Co-mut / Total |
|--------|--------|---------|------------|----------------|
| RYR1 | UBR5 | 3.8 × 10⁻¹³ | 29.9 | 16/30 |
| TRRAP | KMT2C | 3.5 × 10⁻¹² | 27.9 | 15/31 |
| UBR5 | TRRAP | 2.0 × 10⁻¹¹ | — | — |
| KMT2C | ACVR2A | 3.5 × 10⁻¹¹ | — | — |
| NRG1 | TRRAP | 8.8 × 10⁻¹¹ | — | — |
| CHD4 | UBR5 | 5.3 × 10⁻¹⁰ | 30.4 | 11/25 |
| CREBBP | ACVR2A | 6.8 × 10⁻¹⁰ | — | — |
| NOTCH1 | ACVR2A | 6.8 × 10⁻¹⁰ | — | — |
| ATR | ACVR2A | 2.0 × 10⁻⁹ | — | — |
| CHD4 | PIK3CA | 2.1 × 10⁻⁹ | 25.6 | 11/28 |
| NOTCH1 | ATR | 3.4 × 10⁻⁹ | — | — |
| KMT2C | ANK3 | 4.4 × 10⁻⁹ | — | — |

The strongest signal is the **RYR1–UBR5** co-occurrence cluster (p = 3.8 × 10⁻¹³), and a broader set of TRRAP/KMT2C/ACVR2A/NRG1 co-mutations emerges as a recurrent chromatin remodelling cluster.

### Mutual Exclusivity (Broad Run)

All significant mutual exclusivity pairs involve **TP53** on one side, reflecting the well-established genomic divergence between TP53-mutant (CIN/intestinal-type) and TP53-wild-type (non-CIN) gastric cancers.

| Gene 1 | Gene 2 | p-value | Odds Ratio | Notes |
|--------|--------|---------|------------|-------|
| SOX9 | TP53 | 3.2 × 10⁻⁵ | 0.11 | Strongest ME signal overall |
| ATM | TP53 | 5.4 × 10⁻⁴ | 0.26 | DNA damage response |
| POLD1 | TP53 | 3.8 × 10⁻³ | 0.23 | DNA polymerase delta; POLE/POLD spectrum |
| TRRAP | TP53 | 6.0 × 10⁻³ | 0.31 | HAT co-activator |
| CDH1 | TP53 | 9.8 × 10⁻³ | 0.42 | Diffuse vs intestinal subtype |
| ARID1A | TP53 | 1.1 × 10⁻² | 0.53 | SWI/SNF chromatin remodelling |
| RYR1 | TP53 | 2.0 × 10⁻² | 0.43 | — |
| ARID2 | TP53 | 2.1 × 10⁻² | 0.33 | SWI/SNF |
| CIC | TP53 | 2.1 × 10⁻² | 0.33 | Transcriptional repressor |
| KMT2C | TP53 | 3.3 × 10⁻² | 0.46 | Histone H3K4 methyltransferase |
| TP53 | CHD4 | 3.5 × 10⁻² | 0.36 | NuRD complex; chromatin remodelling |

---

## Run B — Focused (Seeds + BioGRID Panel)

**Seed genes:** CDH1, TP53, APC, PIK3CA  
**Panel:** Genes connecting to ≥2 seed genes in the BioGRID Tier 3 Neo4j query  
**Gene set present in MAF:** combined from seeds and multi-seed neighbors

### Top Co-occurring Pairs (p < 1 × 10⁻⁷)

| Gene 1 | Gene 2 | p-value | Odds Ratio | Co-mut / Mutated-either |
|--------|--------|---------|------------|------------------------|
| CHD4 | UBR5 | 5.3 × 10⁻¹⁰ | 30.4 | 11/25 |
| CHD4 | PIK3CA | 2.1 × 10⁻⁹ | 25.6 | 11/28 |
| CREBBP | CHD4 | 2.1 × 10⁻⁸ | 27.5 | 9/23 |
| UBR5 | EPHA2 | 2.6 × 10⁻⁸ | 42.3 | 8/23 |
| ANK3 | CHD4 | 3.6 × 10⁻⁸ | 25.3 | 9/24 |
| BUB1B | TOP2A | 6.7 × 10⁻⁸ | 163.2 | 5/8 |
| UBR5 | PIK3CA | 1.2 × 10⁻⁷ | 14.2 | 11/35 |

### Selected Co-occurring Pairs (p < 1 × 10⁻⁵)

| Gene 1 | Gene 2 | p-value | Odds Ratio | Co-mut / Mutated-either |
|--------|--------|---------|------------|------------------------|
| PLK1 | UBR5 | 1.9 × 10⁻⁷ | 114.7 | 6/22 |
| UBR5 | TOP2A | 3.8 × 10⁻⁷ | 35.3 | 7/24 |
| CREBBP | UBR5 | 5.0 × 10⁻⁷ | 16.6 | 9/30 |
| UBR5 | NUP214 | 7.3 × 10⁻⁷ | 57.4 | 6/23 |
| EP300 | UBR5 | 1.1 × 10⁻⁶ | 18.8 | 8/28 |
| PIK3CA | CREBBP | 1.4 × 10⁻⁶ | 14.2 | 9/33 |
| CREBBP | EPHA2 | 4.0 × 10⁻⁶ | 27.1 | 6/21 |
| UBR5 | PTK2 | 1.1 × 10⁻⁵ | Inf | 4/23 |
| ANK3 | UBR5 | 1.1 × 10⁻⁵ | 12.0 | 8/33 |
| HSP90AA1 | TOP2A | 1.4 × 10⁻⁵ | 58.1 | 4/11 |

### Notable Seed-Involving Pairs

| Pair | p-value | Odds Ratio | Event | Notes |
|------|---------|------------|-------|-------|
| UBR5 & APC | 2.7 × 10⁻⁴ | 6.0 | Co-occurrence | UBR5 co-occurs with APC loss |
| APC & CTNNA1 | 5.0 × 10⁻⁴ | 13.7 | Co-occurrence | Both Wnt/adhesion complex members |
| CREBBP & APC | 1.5 × 10⁻³ | 5.7 | Co-occurrence | 7 co-mutations in 48 APC-mutant samples |
| PIK3CA & APC | 1.4 × 10⁻² | 3.4 | Co-occurrence | Wnt + PI3K pathway co-activation |
| CDH1 & AKT1 | 2.5 × 10⁻² | 7.0 | Co-occurrence | 3 co-mutations; CDH1 loss + AKT1 gain |

### Mutual Exclusivity (Focused Run)

All five mutual exclusivity signals involve **TP53**, consistent with the broad run.

| Gene 1 | Gene 2 | p-value | Odds Ratio | Event Ratio | Interpretation |
|--------|--------|---------|------------|-------------|----------------|
| CDH1 | TP53 | 9.8 × 10⁻³ | 0.42 | 18/272 | Diffuse-type (CDH1-mut) vs CIN/intestinal-type (TP53-mut) |
| TP53 | BRCA1 | 1.1 × 10⁻² | 0 | 0/268 | Zero co-mutations; BRCA1 deficiency excludes TP53 mutation |
| TP53 | PTK2 | 2.7 × 10⁻² | 0 | 0/267 | Zero co-mutations; FAK/PTK2 deficiency in TP53-wild-type context |
| TP53 | CHD4 | 3.5 × 10⁻² | 0.36 | 7/269 | NuRD chromatin remodelling mutually exclusive with p53 pathway |
| TP53 | ERBB2 | 4.7 × 10⁻² | 0.36 | 6/268 | ERBB2 amplification/mutation favours TP53-wild-type background |

---

## Key Findings

### 1. UBR5 is the dominant co-occurrence hub

UBR5 (E3 ubiquitin ligase) appears in the most significant co-occurring pairs across both runs. In the focused run it co-occurs significantly with CHD4, PIK3CA, EPHA2, PLK1, TOP2A, NUP214, CREBBP, EP300, PTK2, ANK3, APC, and others (all p < 0.001). The RYR1–UBR5 pair is the single most significant pair in the broad run (p = 3.8 × 10⁻¹³). UBR5 may mark a hypermutated or genomically unstable subgroup prone to accumulating co-occurring driver events.

### 2. CHD4 defines a co-occurrence cluster with chromatin remodellers

CHD4 (NuRD complex ATPase) co-occurs strongly with UBR5, PIK3CA, CREBBP, ANK3, and EP300. The CHD4–UBR5 pair (p = 5.3 × 10⁻¹⁰, OR = 30.4) and CHD4–CREBBP pair (p = 2.1 × 10⁻⁸, OR = 27.5) suggest a chromatin remodelling co-mutation module. CDH4 is paradoxically mutually exclusive with TP53 (p = 0.035), indicating it belongs to the TP53-wild-type disease context.

### 3. CREBBP and EP300 form a co-occurring HAT module

CREBBP and EP300, both histone acetyltransferases and transcriptional co-activators that interact with all three seed genes (Tier 3), co-occur with each other (p = 4.6 × 10⁻⁵) and both co-occur with UBR5, PIK3CA, CHD4, and ANK3. This HAT module may represent a distinct epigenetically dysregulated subtype.

### 4. BUB1B–TOP2A co-occurrence: a mitotic instability pair

BUB1B (spindle assembly checkpoint) and TOP2A (DNA topoisomerase IIα) co-occur with OR = 163.2 (p = 6.7 × 10⁻⁸), the highest odds ratio among high-confidence pairs. Both are markers of chromosomal instability (CIN), and their co-mutation likely reflects a highly aneuploid tumour subgroup. PLK1 co-occurs with both genes independently (p < 0.001 each), reinforcing the mitotic dysregulation theme.

### 5. CDH1 vs TP53 mutual exclusivity — confirmed in both runs

The CDH1–TP53 mutual exclusivity (p = 9.8 × 10⁻³, OR = 0.42) appears robustly in both runs. Of the 272 samples mutated in either gene, only 18 carry mutations in both (6.6%). This reflects the well-established histological divergence: CDH1-mutant tumours are predominantly diffuse-type/signet ring cell, while TP53-mutant tumours are predominantly intestinal-type/CIN. This finding validates the MAF quality and subtype composition of the cohort.

### 6. TP53 mutual exclusivity defines a broad TP53-wild-type module

Across both runs, TP53 is mutually exclusive with SOX9, ATM, POLD1, TRRAP, CDH1, ARID1A, RYR1, ARID2, CIC, KMT2C, CHD4, BRCA1, PTK2, and ERBB2. Many of these genes (ARID1A, ARID2, KMT2C, CHD4, TRRAP) are chromatin remodellers enriched in EBV-positive or MSI-high gastric cancer subtypes, both of which tend to be TP53-wild-type. The BRCA1–TP53 and PTK2–TP53 pairs show zero co-mutations (OR = 0), suggesting near-complete mutual exclusivity rather than just enrichment.

### 7. Seed gene co-occurrence patterns

- **APC** co-occurs with CTNNA1, CREBBP, UBR5, PIK3CA, CTNNB1, and MDM2 — all Wnt/ubiquitin pathway members.
- **PIK3CA** co-occurs with CHD4, CREBBP, UBR5, ANK3, STAT1, DNMT3A, and ERBB2, consistent with a PI3K-pathway co-activation cluster.
- **CDH1** shows a nominally significant co-occurrence with AKT1 (p = 0.025, OR = 7.0), suggesting AKT activation in CDH1-deficient tumours.

---

## Output Files

| File | Contents |
|------|----------|
| `interactions_broad_top50.tsv` | All pairwise Fisher's exact test results, top 50 genes by mutation frequency (1,225 pairs) |
| `interactions_broad_top50.png` | maftools tile plot, broad run |
| `interactions_focused.tsv` | All pairwise results, focused gene set (2,278 pairs) |
| `interactions_focused.png` | maftools tile plot, focused run |
| `interactions_focused_sig.tsv` | Significant pairs only (p < 0.05), focused run (484 pairs) |
