# Phase 1 Neo4j Knowledge Graph Query Results

**Seed genes:** APC, CDH1, TP53  
**Date:** 2026-05-21  
**Source:** Neo4j KB (45k genes, GO, Reactome, BioGRID, UniProt)

---

## Tier 1 — Per-Gene Annotations

### GO Terms (`tier1_per_gene_go_terms.tsv`, 265 entries)

| Gene | GO Terms |
|------|----------|
| TP53 | 129 |
| CDH1 | 70 |
| APC  | 66 |

### Reactome Pathways (`tier1_per_gene_pathways.tsv`, 202 entries)

| Gene | Pathways |
|------|----------|
| TP53 | 124 |
| CDH1 | 47 |
| APC  | 31 |

### UniProt Protein Function (`tier1_per_gene_protein_function.tsv`, 3 entries)

All three genes fully annotated:

- **APC** (P25054) — Tumor suppressor; promotes rapid degradation of CTNNB1; negative regulator of Wnt signaling; associates with microtubules and actin filaments; controls cell migration via HGF-JNK-MMP9 axis.
- **CDH1** (P12830) — Calcium-dependent cell adhesion protein; regulates cell-cell adhesion, mobility, and proliferation; potent invasion suppressor; ligand for integrin αE/β7; recruits DSG2/DSP to desmosomes.
- **TP53** (P04637) — Multifunctional transcription factor inducing cell cycle arrest, DNA repair, or apoptosis; tumor suppressor across many cancer types; regulates circadian clock; activates lincRNA-p21; pro-apoptotic activity modulated by ASPP1/ASPP2/iASPP axis.

---

## Tier 2 — Shared Reactome Pathways

**(`tier2_shared_pathways.tsv`, 185 entries)**

### Shared by all 3 genes

| Pathway ID | Pathway Name |
|------------|--------------|
| R-HSA-109581 | Apoptosis |
| R-HSA-5357801 | Programmed Cell Death |
| R-HSA-162582 | Signal Transduction |
| R-HSA-1643685 | Disease |

### Shared by 2 genes (selected)

| Pathway | Genes |
|---------|-------|
| Apoptotic execution phase | APC, CDH1 |
| Apoptotic cleavage of cellular proteins | APC, CDH1 |
| Deubiquitination | APC, TP53 |
| Post-translational protein modification | APC, TP53 |
| Metabolism of proteins | APC, TP53 |
| Diseases of signal transduction by growth factor receptors | APC, TP53 |
| Immune System | CDH1, TP53 |
| Developmental Biology | CDH1, TP53 |

---

## Tier 3 — BioGRID Panel Neighbors

**(`tier3_panel_biogrid_neighbors.tsv`, 229 entries)**

### Genes connecting to all 3 seeds

| Gene | Full Name | Evidence Rows | Interaction Types |
|------|-----------|---------------|-------------------|
| CTNNB1 | catenin beta 1 | 100 | physical |
| CREBBP | CREB binding lysine acetyltransferase | 48 | physical, genetic+physical |
| HDAC1 | histone deacetylase 1 | 31 | physical |
| PTK2 | protein tyrosine kinase 2 (FAK) | 13 | physical |
| GSK3B | glycogen synthase kinase 3 beta | 12 | physical |
| CDK4 | cyclin dependent kinase 4 | 5 | genetic, physical, genetic+physical |

### Notable 2-seed connectors (top by evidence)

| Gene | Seeds | Evidence Rows | Notes |
|------|-------|---------------|-------|
| MDM2 | CDH1, TP53 | 674 | Master TP53 ubiquitin ligase |
| EP300 | CDH1, TP53 | 111 | HAT transcriptional coactivator |
| USP7 | APC, TP53 | 58 | Deubiquitinase stabilizing p53 |
| BRCA1 | CDH1, TP53 | 24 | DNA repair |
| AXIN1 | APC, TP53 | 22 | Wnt/APC destruction complex |
| HSP90AA1 | CDH1, TP53 | 19 | Molecular chaperone |
| AURKA | APC, TP53 | 18 | Spindle pole kinase |
| PRKDC | CDH1, TP53 | 17 | DNA-PK, double-strand break repair |
| CTNNA1 | APC, CDH1 | 16 | Alpha-catenin, cadherin-actin linker |
| EGFR | CDH1, TP53 | 8 | Receptor tyrosine kinase |
| SMAD3 | CDH1, TP53 | 4 | TGF-β signaling |
| MET | CDH1, TP53 | 3 | HGF receptor, invasion |
| STAT3 | CDH1, TP53 | 3 | JAK-STAT signaling |
| KRAS | CDH1 | 1 | RAS oncogene |
| PIK3CA | APC | 1 | PI3K catalytic subunit |

---

## Tier 4 — Shared BioGRID Hubs (Full Gene Summaries)

**(`tier4_shared_biogrid_hubs.tsv`, 100 entries)**

Extends tier 3 with full NCBI gene summaries. Additional all-3-seed hubs beyond the tier 3 top 6:

| Gene | Full Name | Evidence Rows | Functional Role |
|------|-----------|---------------|-----------------|
| JUP | junction plakoglobin | 18 | Armadillo repeat; component of both desmosomes and adherens junctions |
| SFN | stratifin (14-3-3 sigma) | 10 | Mitotic translation regulator; DNA damage checkpoint |
| CSNK1A1 | casein kinase 1 alpha 1 | 9 | Part of beta-catenin destruction complex |
| RNF43 | ring finger protein 43 | 9 | E3 ubiquitin ligase; negative regulator of Wnt via frizzled ubiquitination |
| TRIM25 | tripartite motif containing 25 | 8 | E3 ubiquitin ligase; antiviral innate immunity |
| CSNK1D | casein kinase 1 delta | 7 | Cell cycle, apoptosis, circadian rhythm, p53 regulation |
| HSPA5 | heat shock protein 5 (GRP78/BiP) | 7 | ER chaperone; UPR regulator; therapeutic target in cancer |
| YWHAZ | 14-3-3 zeta | 7 | Phosphoserine-binding signal transducer; IRS1 interaction |
| CFTR | CF transmembrane conductance regulator | 6 | Chloride channel; ion/water secretion in epithelial tissues |
| YWHAB | 14-3-3 beta | 6 | RAF1/CDC25 interaction; links mitogenic signaling to cell cycle |
| PPP1R13B | ASPP1 | 5 | p53 apoptosis activator; promotes p53 binding to proapoptotic gene promoters |
| UBC | ubiquitin C | 5 | Polyubiquitin precursor; ubiquitination hub |
| EZR | ezrin | 4 | ERM family; plasma membrane–actin cytoskeleton linker; cell adhesion/migration |
| LZTS2 | leucine zipper tumor suppressor 2 | 4 | Represses beta-catenin transcription; negative regulator of Wnt |
| MYO6 | myosin VI | 4 | Reverse-direction actin motor; intracellular vesicle transport |
| PXN | paxillin | 4 | Cytoskeletal scaffold at focal adhesions |
| USP15 | ubiquitin specific peptidase 15 | 4 | Deubiquitinase; SMAD stabilization in TGF-β signaling |
| YWHAE | 14-3-3 epsilon | 4 | CDC25/RAF1 signal transducer; implicated in small cell lung cancer |
| RB1CC1 | RB1 inducible coiled-coil 1 | 3 | Tumor suppressor; autophagy/apoptosis/cell migration coordinator |
| KRT8 | keratin 8 | 3 | Epithelial intermediate filament; structural integrity and signal transduction |
| SPTBN1 | spectrin beta non-erythrocytic 1 | 3 | Plasma membrane–actin cytoskeleton scaffold |

Notable 2-seed hubs with full summaries: **MDM2** (674 rows, p53 E3 ligase), **EP300** (111 rows, HAT coactivator, Rubinstein-Taybi syndrome), **AXIN1** (22 rows, Wnt destruction complex, hepatocellular carcinoma association), **CTNNA1** (16 rows, mechanosensing alpha-catenin).

---

## Tier 5 — Reactome PPI Partners

**(`tier5_reactome_ppi.tsv`, 102 entries; 72 unique partner genes)**

Literature-supported physical interactions from Reactome.

### Coverage by seed gene

| Seed | Reactome PPI Partners |
|------|----------------------|
| TP53 | 62 |
| CDH1 | 36 |
| APC  | 4 (CASP3, BTRC, CTBP1, CTNNB1) |

### Selected TP53 partners (Reactome-curated)

**DNA damage & checkpoint:** ATM, ATR, ATRIP, CHEK1, CHEK2, CDK2, CDK5/CDK5R1, DYRK2, PLK3, HIPK1  
**Ubiquitin/SUMO:** MDM2, MDM4, USP7, USP10, UBE2I, SUMO2, SUMO3  
**Transcription:** EP300, CREBBP, MEN1, RUNX3, SIN3A, MYB, FOS, JUN, E2F4/TFDP1/TFDP2  
**Chromatin:** L3MBTL1, EHMT1, EHMT2, KAT6A, BRPF1, MEAF6, ING5  
**Apoptosis/signaling:** MAPK8, MAPK9, MAPKAPK5, EIF2AK2, BANP, PML  
**Cell cycle:** CDKN2A, GTSE1, PIN1, STK11  

### Selected CDH1 partners (Reactome-curated)

**Catenin complex:** CTNNB1, CTNNA1, CTNND1, JUP  
**Cytoskeleton/adhesion:** VCL, IQGAP1, ACTA1/2, ACTB/C1/G1/G2, PIP5K1C  
**Proteolysis:** MMP3, MMP7, MMP9, ADAM10, ADAM15, FURIN, PCSK7, KLK7, CAPN1/CAPNS1  
**Receptor/signaling:** SRC, FYN, NCSTN, PSEN1, KLRG1, ITGAE, ITGB7  
**Transcription/proliferation:** AURKA, AURKB, CSNK2A1/2/B  

---

## Bonus — Shared GO Biological Process Terms

**(`bonus_shared_go_bp.tsv`, 132 entries)**

### Shared by APC + TP53

| GO Term | Definition |
|---------|------------|
| GO:0006974 — DNA damage response | Change in cell activity in response to DNA damage |
| GO:0008285 — Negative regulation of cell population proliferation | Stops/reduces cell proliferation |

### Shared by APC + CDH1

| GO Term | Definition |
|---------|------------|
| GO:0007155 — Cell adhesion | Attachment of cell to another cell or ECM |
| GO:0016477 — Cell migration | Controlled self-propelled cell movement |

### Shared by APC + CDH1 (apoptosis)

| GO Term | Definition |
|---------|------------|
| GO:0097553 — Apoptotic execution phase (implied via tier2) | Execution of apoptosis program |

---

## Key Findings

### 1. CTNNB1 is the central molecular hub

Beta-catenin physically interacts with all 3 seed genes (100 BioGRID evidence rows), directly bridging the Wnt pathway (APC), epithelial adhesion (CDH1), and transcriptional regulation (TP53). This makes the APC–CTNNB1–CDH1 junction a primary convergence point in gastric carcinogenesis.

### 2. Wnt destruction complex is fully reconstructed

The core complex — APC, AXIN1, GSK3B, CSNK1A1/CSNK1D, BTRC, CTNNB1 — is represented across tiers 2–5, with multiple members connecting to ≥2 seed genes.

### 3. Chromatin regulation as a cross-seed theme

CREBBP, HDAC1, EP300, KDM1A, HDAC2, and EHMT1/2 all bridge multiple seed genes, suggesting that epigenetic dysregulation (acetylation/deacetylation balance) is a shared downstream consequence of APC/CDH1/TP53 mutation.

### 4. MDM2 is the highest-evidence 2-seed hub

674 BioGRID evidence rows linking CDH1 and TP53. The CDH1–MDM2–TP53 axis warrants prioritization for functional analysis, as CDH1 loss may modulate MDM2-mediated p53 degradation.

### 5. TP53 dominates the interaction network

TP53 contributes the majority of annotations across all tiers: 129 GO terms, 124 Reactome pathways, 62 Reactome PPI partners, and the largest BioGRID neighbor set. This is consistent with p53's role as a master regulator and supports its use as an anchor gene in downstream analyses.

### 6. Epithelial integrity genes converge on CDH1

CDH1's Reactome PPI network includes catenin complex members (CTNNB1, CTNNA1, CTNND1, JUP), cytoskeletal linkers (VCL, IQGAP1, actins), and multiple MMPs and ADAMs, reflecting CDH1's dual role in maintaining epithelial structure and suppressing invasion.
