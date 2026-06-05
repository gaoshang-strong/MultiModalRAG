# Analysis Prompt: Why Are These Genes Mutually Exclusively Mutated with TP53 in Gastric Cancer?

## Background and Question

We performed pairwise somatic interaction analysis (Fisher's exact test, maftools::somaticInteractions) on a gastric cancer cohort of **N = 445 patients** sequenced with a targeted DNA panel. The top 50 most frequently mutated genes were tested for all pairwise co-occurrence and mutual exclusivity.

**TP53 is mutated in 59% of the cohort (263/445 samples).**

The following five genes are significantly **mutually exclusively mutated** with TP53:

| Gene | p-value | Odds Ratio | Co-mutations (both mutated) / Total mutated in either |
|------|---------|------------|------------------------------------------------------|
| SOX9 | 3.2 × 10⁻⁵ | 0.11 | — |
| ATM | 5.4 × 10⁻⁴ | 0.26 | — |
| POLD1 | 3.8 × 10⁻³ | 0.23 | — |
| TRRAP | 6.0 × 10⁻³ | 0.31 | — |
| CDH1 | 9.8 × 10⁻³ | 0.42 | 18/272 |

**The core question: What are the biological and/or molecular mechanisms that explain why mutations in POLD1, SOX9, TRRAP, ATM, and CDH1 tend to occur in gastric tumours that do NOT carry TP53 mutations?**

---

## Knowledge Graph Data (Neo4j Queries)

The following data was retrieved by querying these five genes as seeds in a biological knowledge graph (NCBI Gene, UniProtKB Swiss-Prot, Reactome v96, BioGRID 5.0, Gene Ontology).

---

### 1. Per-Gene Protein Functions (UniProt Swiss-Prot)

**ATM** (Q13315 — Serine-protein kinase ATM)
- Serine/threonine kinase activated by DNA double-strand breaks (DSBs); master DNA damage sensor
- Phosphorylates TP53 directly (Ser-15), CHEK2, BRCA1, H2AX, MDM2, MDM4, FANCD2, RBBP8, NBN, RAD50, MRE11
- Also phosphorylates DYRK2 to prevent MDM2-mediated TP53 degradation
- Functions: DNA damage response, apoptosis, cell cycle checkpoint, meiotic recombination, pexophagy
- Keywords: Tumor suppressor, DNA damage, Serine/threonine-protein kinase, Cell cycle

**POLD1** (P28340 — DNA polymerase delta catalytic subunit)
- Catalytic subunit of DNA polymerase delta (Pol-δ); essential for lagging strand synthesis and high-fidelity replication
- Has both DNA polymerase and 3'→5' proofreading exonuclease activity
- Involved in: NER (nucleotide excision repair), BIR (break-induced replication), MMR (mismatch repair), translesion synthesis
- Mutations in POLD1 cause an ultramutator/POLE-like phenotype (proofreading domain mutations → hypermutation)
- Keywords: DNA replication, DNA repair, Exonuclease, Disease variant

**SOX9** (P48436 — Transcription factor SOX-9)
- SRY-related HMG-box transcription factor; binds 5'-ACAAAG-3' motifs in enhancers and super-enhancers
- Master regulator of chondrogenesis and skeletal development
- In epithelial tissues: regulates proliferation/differentiation balance in stem/progenitor cells (lung, kidney branching)
- Suppresses beta-catenin (CTNNB1) signaling and RUNX2 expression in chondrocytes
- Acts downstream of Hedgehog (GLI1/GLI3) and cooperates with SOX5/SOX6
- In cancer: known marker of intestinal stem cell identity; expressed in gastric intestinal metaplasia
- Keywords: Transcription, Activator, Differentiation, Disease variant, Nucleus

**TRRAP** (Q9Y4A5 — Transformation/transcription domain-associated protein)
- Adaptor protein in multiple histone acetyltransferase (HAT) complexes; component of NuA4/TIP60 HAT complex
- Acetylates nucleosomal H4 and H2A; central to epigenetic transcriptional activation
- Required for TP53-, E2F1-, E2F4-mediated transcription activation
- Also required for MYC transcription activation and MYC-driven cell transformation
- Links transcription factors (E1A, MYC, E2F1) to HAT complexes (STAGA)
- May be required for mitotic checkpoint and normal cell cycle progression
- Keywords: Chromatin regulator, Transcription, Activator, Nucleus

**CDH1** (P12830 — Cadherin-1 / E-cadherin)
- Calcium-dependent cell adhesion protein; forms adherens junctions in epithelial cells
- Potent invasion suppressor; loss drives epithelial-to-mesenchymal transition (EMT)
- Anchored to actin cytoskeleton via alpha-, beta-, gamma-catenin complex
- Ligand for integrin αE/β7; recruits DSG2/DSP to desmosomes
- Germline CDH1 loss-of-function mutations cause hereditary diffuse gastric cancer (HDGC)
- Somatic CDH1 mutations define the diffuse-type (signet ring cell) gastric cancer histotype
- Keywords: Cell adhesion, Cell junction, Tumor suppressor, Disease variant

---

### 2. Reactome Pathways (Tier 1 — per gene)

| Gene | Pathway Count | Key Pathway Categories |
|------|---------------|----------------------|
| ATM | 61 | DNA Repair (DSB, HR, NER, MMR, NHEJ), Cell Cycle Checkpoints, TP53 Regulation, Autophagy, Senescence |
| CDH1 | 47 | Cell Adhesion, Adherens Junctions, Apoptosis, EMT, Immune System, Developmental Biology |
| POLD1 | 43 | DNA Replication (S-phase, Lagging/Leading Strand), DNA Repair (BER, NER, MMR, TLS), Telomere Maintenance |
| SOX9 | 17 | Wnt/β-catenin Signaling, Developmental Cell Lineages, Transcriptional Regulation, RUNX2 |
| TRRAP | 11 | Chromatin Organization, HAT (Histone Acetylation), Wnt/TCF Signaling, Deubiquitination |

**Notable ATM-specific pathways:**
- Regulation of TP53 Activity (phosphorylation, degradation, methylation)
- Stabilization of p53
- Transcriptional Regulation by TP53
- TP53 Regulates Transcription of DNA Repair Genes
- DNA Double Strand Break Response
- Homologous Recombination Repair (HRR); defective HRR due to BRCA1/BRCA2/PALB2 loss

---

### 3. Shared Reactome Pathways (Tier 2 — cross-seed)

No pathway is shared by all 5 seeds. Pairwise sharing only:

| Pathway | Seeds | Significance |
|---------|-------|-------------|
| Signal Transduction | CDH1, SOX9, TRRAP | General signaling convergence |
| Cell Cycle | ATM, POLD1 | Replication fidelity and checkpoint |
| DNA Double-Strand Break Repair | ATM, POLD1 | Core DNA repair machinery |
| DNA Repair | ATM, POLD1 | Genome maintenance |
| HDR through Homologous Recombination (HRR) | ATM, POLD1 | HR-directed DSB repair |
| Homology Directed Repair | ATM, POLD1 | DSB repair via HR |
| Developmental Biology | CDH1, SOX9 | Epithelial identity and lineage |
| Developmental Cell Lineages | CDH1, SOX9 | Tissue-specific differentiation |
| Wnt Signaling | SOX9, TRRAP | β-catenin pathway regulation |
| TCF-dependent Wnt Signaling | SOX9, TRRAP | Wnt transcriptional output |
| Gene Expression / Transcription | ATM, SOX9 | Transcriptional regulation |
| RNA Polymerase II Transcription | ATM, SOX9 | General transcription |
| MITF-M regulated melanocyte development | CDH1, SOX9 | Epithelial differentiation |

**Key observation:** ATM and POLD1 share exclusively DNA repair/replication pathways. CDH1 and SOX9 share developmental/epithelial identity pathways. SOX9 and TRRAP share Wnt signaling. There is **no single shared pathway** connecting all 5 genes.

---

### 4. Shared GO Biological Process Terms (Bonus — pairs only, n_seeds ≥ 2)

| GO Term | Seeds | Biological Meaning |
|---------|-------|-------------------|
| Positive regulation of DNA-templated transcription | CDH1, SOX9, TRRAP | Transcriptional activation |
| Regulation of apoptotic process | ATM, SOX9, TRRAP | Apoptosis modulation |
| DNA damage response | ATM, POLD1 | Genome surveillance |
| DNA repair | ATM, POLD1 | Genome maintenance |
| Chromatin remodeling | ATM, SOX9 | Epigenetic regulation |
| Positive regulation of gene expression | ATM, SOX9 | Transcriptional activation |
| Regulation of cell cycle | ATM, TRRAP | Cell cycle control |
| Regulation of cellular response to stress | ATM, TRRAP | Stress response |
| Regulation of gene expression | ATM, CDH1 | Transcriptional control |
| Signal transduction | ATM, SOX9 | Signaling |
| Cellular response to retinoic acid | ATM, SOX9 | Differentiation signaling |
| Positive regulation of transcription by RNA Pol II | ATM, SOX9 | Transcriptional activation |

**Key observation:** Transcriptional regulation and DNA damage response are the two dominant shared functional themes, but no GO term unifies all 5 seeds.

---

### 5. Somatic Panel BioGRID Neighbors — Multi-Seed Connectors (Tier 3)

Panel genes (from the targeted sequencing panel) that are BioGRID 1st-degree interactors of ≥3 seeds:

| Panel Gene | Full Name | Connected Seeds | Evidence Rows | Functional Note |
|-----------|-----------|-----------------|---------------|-----------------|
| CUL3 | cullin 3 | ATM, POLD1, SOX9, TRRAP | 5 | E3 ubiquitin ligase scaffold; ubiquitin-proteasome system |
| BRD4 | bromodomain containing 4 | ATM, POLD1, SOX9, TRRAP | 5 | BET bromodomain; chromatin reader; super-enhancer regulation |
| TP53 | tumor protein p53 | ATM, POLD1, TRRAP | 40 | **The mutually exclusive gene itself** |
| TOP2A | DNA topoisomerase II alpha | ATM, CDH1, POLD1 | 7 | Chromatin topology; replication; CIN marker |
| SMARCA4 | SWI/SNF ATPase subunit | ATM, SOX9, TRRAP | 3 | SWI/SNF chromatin remodeling complex |
| PRKN | parkin E3 ubiquitin ligase | ATM, POLD1, TRRAP | 5 | Mitophagy; proteasomal degradation |
| PRKDC | DNA-PKcs | ATM, CDH1, POLD1 | 7 | NHEJ; DNA double-strand break repair |
| PARP1 | poly(ADP-ribose) polymerase 1 | ATM, POLD1, TRRAP | 4 | BER; DNA damage response; chromatin |
| MRE11 | MRE11 DSB repair nuclease | ATM, CDH1, TRRAP | 7 | MRN complex; DSB sensing and resection |
| KRAS | KRAS proto-oncogene GTPase | CDH1, POLD1, TRRAP | 4 | RAS signaling; oncogenic driver |
| EP300 | EP300 HAT coactivator | CDH1, SOX9, TRRAP | 8 | Histone acetyltransferase; p300; transcriptional coactivation |
| EGFR | epidermal growth factor receptor | CDH1, POLD1, TRRAP | 8 | RTK signaling; proliferation |
| BRCA1 | BRCA1 DNA repair | ATM, CDH1, POLD1 | 14 | HR repair; genome stability; tumour suppressor |
| AURKB | aurora kinase B | CDH1, POLD1, TRRAP | 3 | Mitotic kinase; chromosomal segregation |

**Critical observation: TP53 itself appears as a BioGRID neighbor connecting ATM, POLD1, and TRRAP (40 evidence rows).**

---

### 6. Shared BioGRID Hub Genes — Genes Connecting ≥2 Seeds (Tier 4, top hubs)

| Hub Gene | Connected Seeds | Evidence Rows | Summary |
|----------|-----------------|---------------|---------|
| ZRANB1 | ATM, POLD1, SOX9, TRRAP | 8 | Deubiquitinase; K63-linked polyubiquitin binding; Wnt signaling |
| BRD4 | ATM, POLD1, SOX9, TRRAP | 5 | BET bromodomain; super-enhancer reader; MYC/E2F target genes |
| CUL3 | ATM, POLD1, SOX9, TRRAP | 5 | Cullin-RING E3 ubiquitin ligase; protein degradation scaffold |
| TRIM67 | ATM, POLD1, SOX9, TRRAP | 4 | TRIM E3 ligase; RAS signaling negative regulator |
| TP53 | ATM, POLD1, TRRAP | 40 | **The mutually exclusive gene** |
| PCNA | ATM, POLD1, TRRAP | 16 | DNA sliding clamp; replication processivity; damage bypass |
| BRCA1 | ATM, CDH1, POLD1 | 14 | BASC complex; HR repair; tumour suppressor |
| EP400 | ATM, SOX9, TRRAP | 12 | NuA4/SWR1 HAT complex; H2A.Z deposition; DSB HR repair |
| MDC1 | ATM, POLD1, TRRAP | 12 | γH2AX binding; ATM/MRN recruitment to DSBs |
| MAX | CDH1, POLD1, TRRAP | 11 | MYC heterodimerisation partner; E-box transcription |
| EGFR | CDH1, POLD1, TRRAP | 8 | RTK; MAPK/PI3K signaling |
| EP300 | CDH1, SOX9, TRRAP | 8 | p300 HAT; acetylates H3K27; coactivates TP53, MYC, SOX9 |
| RPA1/RPA2 | ATM, POLD1, TRRAP | 8/7 | ssDNA-binding complex; ATR activation; replication and repair |
| TOP2A | ATM, CDH1, POLD1 | 7 | Topoisomerase; CIN marker |
| PRKDC | ATM, CDH1, POLD1 | 7 | DNA-PKcs; NHEJ |
| MRE11 | ATM, CDH1, TRRAP | 7 | MRN complex; DSB sensing |
| BSG/Basigin | ATM, CDH1, TRRAP | 7 | Plasma membrane; spermatogenesis; tumour progression |
| E2F1 | ATM, CDH1, TRRAP | 5 | Cell cycle transcription factor; S-phase entry; apoptosis |
| KRAS | CDH1, POLD1, TRRAP | 4 | Oncogenic RAS; MAPK/PI3K effector |
| MYC | POLD1, TRRAP | 29 | Proto-oncogene; proliferation; TRRAP-dependent HAT recruitment |
| PARP1 | ATM, POLD1, TRRAP | 4 | DNA damage; BER; chromatin regulation |

---

### 7. Reactome Protein-Protein Interactions (Tier 5 — mechanistic, literature-supported)

**ATM partners (72 total):** BRCA1, BRCA2, CHEK2, MDM2, MDM4, H2AX, MRE11, NBN, RAD50, RPA1, RPA2, ATR, ATRIP, DCLRE1C, EXO1, CDK2, RBBP8, RIF1, RNF168, RNF8, MDC1, KAT5, NSD2, PAXIP1, MAPK8, ABL1, PIDD1, HERC2, H3/H4 histones, MAP1LC3B, PEX5, BLM, DNA2, DYRK2, MLH1, MLH3, SUMO1

**POLD1 partners (22 total):** PCNA, POLD2, POLD3, POLD4, RPA1, RPA2, FEN1, LIG1, LIG3, DNA2, APEX1, POLB, POLA1, PRIM1, PRIM2, XRCC1, RAD18, RBX1, UBE2B, CUL4A, CUL4B, DDB1

**SOX9 partners (3 total):** GATA4, NR5A1 (SF-1), RUNX2

**TRRAP partners (20 total):** CTNNB1 (β-catenin), KAT5 (TIP60/NuA4 HAT), RUVBL1 (AAA+ ATPase), LEF1, TCF7L2, H3/H4 histones

**CDH1 partners (36 total):** CTNNB1, CTNNA1, CTNND1, JUP, VCL, IQGAP1, SRC, FYN, ITGAE, ITGB7, MDM2, MEN1, MMP3, MMP7, MMP9, ADAM10, ADAM15, CASP3, CAPN1, CAPNS1, FURIN, KLK7, KLRG1, NCSTN, PSEN1, PLG, PCSK7, PIP5K1C, ANK3, CANX, actins (ACTA1/2, ACTB/C1/G1/G2)

**Key mechanistic observations:**
- ATM directly phosphorylates TP53 (Ser-15) — ATM is the upstream activator of TP53 in response to DSBs
- TRRAP is required for TP53-mediated transcriptional activation
- CDH1 interacts with MDM2 (the primary E3 ligase that degrades TP53) in Reactome
- POLD1 shares RPA1/RPA2 with ATM in DNA replication and repair contexts
- TRRAP interacts with β-catenin (CTNNB1) and TCF7L2/LEF1, linking it to Wnt signalling

---

## Statistical Context

From the somatic interaction analysis:
- TP53 is mutually exclusive with **all 11 significant mutual exclusivity pairs** identified in the broad run (top 50 genes)
- This pattern is not seen for any other gene — TP53 uniquely defines one pole of the mutual exclusivity landscape
- The mutual exclusivity signal is present at nominal p-value (not all survive BH correction), suggesting the effect is real but moderate
- The extremely high TP53 mutation rate (59%) means the statistical test is sensitive: even small depletion of co-mutations is detectable

---

## Specific Sub-questions to Address

1. **ATM:** ATM directly phosphorylates and activates TP53. If both ATM and TP53 are functional, the ATM→TP53 axis is intact. Does mutation of ATM make TP53 mutation redundant, or does ATM mutation provide an alternative route to p53 pathway inactivation? What does the literature say about ATM vs TP53 mutation patterns in GC molecular subtypes?

2. **POLD1:** POLD1 proofreading-domain mutations cause hypermutation (ultramutator phenotype, similar to POLE). In POLE/POLD1-mutant tumours, why would TP53 mutation be less frequent — is it because hypermutated tumours rely on different transformation mechanisms, or because the high mutational load itself triggers immune pressure that selects against clones accumulating certain drivers?

3. **TRRAP:** TRRAP is required for TP53-dependent transcriptional activation. TRRAP mutation could disrupt the downstream arm of p53 signalling (the transcriptional output) rather than p53 protein itself. Is this a case of functional equivalence — inactivating the effector vs inactivating the regulator?

4. **SOX9:** SOX9 is not a classic DNA damage response or p53 pathway gene. Its mutual exclusivity with TP53 is the most significant (p = 3.2 × 10⁻⁵, OR = 0.11). SOX9 is a marker of intestinal stem cell identity and is expressed in gastric intestinal metaplasia. What gastric cancer molecular subtype is associated with SOX9 mutation, and why would that subtype be TP53-wild-type?

5. **CDH1:** This is the most interpretable case — CDH1 mutations define diffuse-type gastric cancer (Lauren classification), which is characteristically TP53-wild-type, while TP53 mutations define intestinal-type/CIN gastric cancer. How do the TCGA/ACRG molecular subtypes map onto this CDH1 vs TP53 dichotomy?

6. **Are these five genes defining the same TP53-wild-type subgroup or different ones?** The absence of any shared pathway across all 5 genes suggests they may mark different TP53-wild-type tumour subtypes (EBV+, MSI-H, diffuse, POLE-mutant), rather than a single unified molecular mechanism.

---

## Please address the following:

Given the biological functions, pathway memberships, protein-protein interaction networks, and shared molecular themes described above, provide a mechanistic explanation for why each of POLD1, SOX9, TRRAP, ATM, and CDH1 tends to be mutated in gastric cancers that do not carry TP53 mutations. Consider:

- Whether the mutual exclusivity reflects **functional redundancy** (mutation of the ME gene achieves the same oncogenic effect as TP53 mutation)
- Whether it reflects **tumour subtype separation** (different histological or molecular subtypes of gastric cancer that independently tend to be TP53-wild-type)
- Whether it reflects **synthetic lethality or negative epistasis** (co-mutation of both genes is selected against because it is deleterious to the tumour)
- Whether the **ATM→TP53 phosphorylation axis** is particularly relevant for the ATM mutual exclusivity
- How the **TRRAP–TP53 transcriptional co-activation relationship** (TRRAP is required for TP53 transcriptional output) might explain TRRAP mutual exclusivity
- What the **gastric cancer TCGA molecular subtypes** (EBV, MSI, Genomically Stable/GS, Chromosomal Instability/CIN) predict about the distribution of these mutations
