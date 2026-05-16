# FDU Gastric Cancer Cohort — Analysis Plan

Data source: `Target_panel_DNA_seq_FDU/MAF_data/`  
Cohort: ~511 patients, targeted panel DNA-seq, tumor + blood paired, hg38 coordinates.

---

## 1. Somatic Mutation Landscape

**What:** Compute mutation frequency per gene across the cohort. Visualize as an oncoprint (waterfall plot) for the top 20–30 genes.

**Preferred outcome:** A ranked list of the most frequently mutated cancer genes (e.g., TP53, CDH1, ARID1A, RHOA, KRAS) and their co-occurrence patterns across patients.

**Why:** Establishes the baseline mutational landscape of this Chinese gastric cancer cohort and enables comparison to TCGA-STAD. Chinese gastric cancer has known epidemiological and molecular differences from Western cohorts; confirming or contrasting driver gene frequencies is a necessary first step for any downstream analysis.

---

## 2. Significantly Mutated Gene (SMG) Detection

**What:** Apply statistical frameworks (dNdScv or MutSigCV-equivalent) to identify genes mutated more than expected by background mutation rate.

**Preferred outcome:** A set of bona fide driver genes specific to this cohort, with q-values and mutation spectrum context.

**Why:** Not all recurrently mutated genes are drivers. SMG analysis distinguishes selection-driven mutations from passenger noise, producing a credible driver gene list for mechanistic and biomarker follow-up.

---

## 3. Mutational Signature Analysis

**What:** Decompose the somatic mutation catalog into COSMIC v3.3 single-base substitution (SBS) signatures using SigProfilerExtractor or MutationalPatterns.

**Preferred outcome:** Relative contributions of known signatures (e.g., SBS1 aging, SBS3 HRD, SBS17 5-FU exposure, SBS18 oxidative damage, SBS30 germline NTHL1) per patient and per group.

**Why:** Mutational signatures reveal etiology (carcinogen exposure, DNA repair defects, treatment effect). SBS17 is of particular interest given 5-FU use in neoadjuvant regimens. SBS3 (HRD) could identify patients who may benefit from PARP inhibitors. Signature differences between treatment-naive and neoadjuvant groups would directly reflect therapy-induced mutagenesis.

---

## 4. TMB Distribution and Clinical Correlation

**What:** Use `calculated_TMB` from the filtered clinical table. Plot distribution, define high/low cutoffs (e.g., 10 mut/Mb), and correlate TMB with: MSI status, MMR status, HER2, pathological stage, and OS.

**Preferred outcome:** Validated TMB cutoffs with clinical meaning; confirmation or refinement of the TMB ↔ MSI/dMMR relationship in this panel-sequenced cohort.

**Why:** TMB is an approved biomarker for immunotherapy eligibility (FDA, 2020). Panel-derived TMB can differ from WES-derived TMB due to capture bias. Recalibrating and validating TMB in this cohort is essential before using it for clinical decision support or model building.

---

## 5. MSI / dMMR Classification from Somatic Mutations

**What:** Train a classifier (logistic regression, random forest, or gradient boosting) on somatic variant features (indel fraction, mutation count, specific gene mutations) to predict MSI-H / dMMR status, using `y_msi` and `y_dmmr` labels.

**Preferred outcome:** An accurate, panel-calibrated MSI predictor (AUC > 0.90) that does not require IHC or PCR, with feature importance revealing which genomic signals are most predictive.

**Why:** IHC and PCR-based MSI testing can fail on degraded FFPE samples. A sequence-based predictor is a practical fallback and enables retrospective classification of samples with missing IHC data.

---

## 6. Co-mutation and Mutual Exclusivity Analysis

**What:** Test all pairwise gene combinations for co-occurrence or mutual exclusivity using Fisher's exact test with FDR correction. Visualize as a co-mutation heatmap.

**Preferred outcome:** Biologically coherent gene pairs — e.g., mutual exclusivity between KRAS and other RAS/RAF pathway genes, or co-occurrence of CDH1 + RHOA (diffuse-type signature).

**Why:** Co-mutation patterns reveal pathway-level constraints and can stratify patients into molecular subtypes. Mutual exclusivity is a hallmark of synthetic lethality candidates. These patterns inform combination therapy hypotheses.

---

## 7. Molecular Subtyping

**What:** Assign each patient to EBV+, MSI, GS (genomically stable), or CIN subtypes following the TCGA 2014 gastric cancer classification, using available IHC (EBER, MMR), genomic (MSI, copy number proxies), and mutation data.

**Preferred outcome:** Subtype distribution for this Chinese cohort vs. published TCGA proportions; downstream subtype-stratified survival and treatment response analysis.

**Why:** The four TCGA subtypes have distinct prognoses and therapeutic implications. Understanding which subtypes dominate in this Chinese FFPE cohort contextualizes all other findings and aligns the cohort with the international literature.

---

## 8. Treatment Response Prediction (Neoadjuvant Group)

**What:** In the neoadjuvant-treated subset, use pre-treatment genomic features (somatic mutations, TMB, MSI, HER2, PD-L1 CPS, mutational signatures) to predict TRG score (continuous or binarized: good response TRG 0–1 vs. poor response TRG 2–3).

**Preferred outcome:** A genomic predictor of neoadjuvant response with AUC > 0.70, identifying patient subgroups most likely to achieve pathological complete or near-complete response.

**Why:** Neoadjuvant therapy (SOX ± immunotherapy) is standard for locally advanced gastric cancer, but response rates vary widely (~15–30% pCR). A pretreatment genomic predictor would enable patient selection, spare non-responders from toxicity, and guide upfront treatment intensification or de-escalation.

---

## 9. Germline Pathogenic Variant Analysis

**What:** Filter `germ.maf` for pathogenic/likely pathogenic variants (ClinVar CLNSIG) in cancer predisposition genes (BRCA1/2, MLH1, MSH2, MSH6, PMS2, PALB2, ATM, CDH1, TP53). Calculate carrier frequency.

**Preferred outcome:** Prevalence of germline cancer predisposition variants in this gastric cancer cohort; flag patients who may warrant genetic counseling.

**Why:** Hereditary gastric cancer (especially CDH1 germline mutations causing hereditary diffuse gastric cancer) and Lynch syndrome (MLH1/MSH2/MSH6/PMS2) have direct clinical management implications for patients and their families. Population-level carrier frequency in Chinese patients is poorly characterized.

---

## 10. Tumor Clonality and VAF Distribution

**What:** Use allele frequency (`AF`) from somatic VCFs to estimate clonal architecture. Plot VAF distributions per patient; compute clonal vs. subclonal mutation fractions.

**Preferred outcome:** Patient-level clonality scores; association of high subclonal fraction with worse prognosis or poor treatment response.

**Why:** Intra-tumor heterogeneity (high subclonal burden) is associated with treatment resistance and immune evasion. VAF-based clonality estimation is accessible from standard panel data and does not require WGS.

---

## 11. Cross-Cohort Comparison: FDU Panel vs. TCGA-STAD

**What:** Use `combined_TCGA_PANEL_hg38_MAF.csv.gz` to compare mutation frequencies, TMB distributions, driver gene prevalence, and molecular subtype proportions between this Chinese hospital cohort and TCGA-STAD (predominantly Western patients).

**Preferred outcome:** A quantified list of statistically significant differences in mutation landscape (Fisher's exact per gene), TMB, and subtype distribution.

**Why:** Gastric cancer in East Asian populations has distinct molecular features (e.g., higher EBV+, different TP53 spectrum). Identifying these differences validates cohort representativeness and may uncover population-specific therapeutic targets or biomarkers.

---

## 12. Survival Analysis — Genomic Correlates of OS

**What:** Use Cox proportional hazards and Kaplan-Meier analysis to assess association of genomic features (TMB-H/L, MSI status, MMR, HER2, key driver gene mutations, molecular subtype) with overall survival.

**Preferred outcome:** Genomic features independently prognostic for OS after adjusting for stage and treatment; a candidate prognostic signature.

**Why:** OS is the primary clinical endpoint. Identifying genomic features that independently predict survival — beyond stage — adds clinical utility and motivates prospective validation.

---

## 13. HER2 Genomic Context

**What:** Characterize somatic mutation profiles stratified by HER2 IHC status (0/1+/2+/3+). Test whether HER2-high tumors have distinct co-mutations, TMB levels, or mutational signatures.

**Preferred outcome:** Co-mutation landscape specific to HER2+ tumors; any genomic features associated with HER2 amplification that could explain resistance to trastuzumab.

**Why:** HER2+ gastric cancer (IHC 3+ or 2+/FISH+) is treated with trastuzumab, but response rates are ~47% and resistance is common. Genomic co-alterations (e.g., MET, EGFR, KRAS co-amplification) are known resistance mechanisms worth characterizing in this cohort.

---

## 14. Treatment-Naive vs. Neoadjuvant-Treated Mutation Landscape Comparison

**What:** Compare somatic mutation frequency, TMB, mutational signatures, and driver gene prevalence between the two patient groups using the pre-treatment samples.

**Preferred outcome:** Confirmation that the two groups are genomically comparable at baseline (ruling out selection bias), or identification of molecular differences that explain why some patients received neoadjuvant therapy.

**Why:** Selection bias between treatment groups could confound any response or survival analysis. This comparison is a necessary validity check before any cross-group inference.

---

## 15. Germline + Somatic Double-Hit Analysis (Knudson Two-Hit)

**What:** For each patient, cross-reference `germ.maf` and somatic VCFs to find genes carrying both a germline variant (pathogenic or likely pathogenic) and an independent somatic mutation. Also flag somatic variants with high VAF (> 0.7) in known tumor suppressors as a proxy for LOH — i.e., somatic mutation + loss of the remaining wild-type allele.

**Double-hit criteria (any two of the following in the same gene, same patient):**
- Germline pathogenic/likely pathogenic variant (ClinVar CLNSIG)
- Somatic coding mutation
- Allelic imbalance: somatic AF > 0.7 (suggesting LOH of the second allele)
- Homozygous somatic variant (GT = 1/1 in tumor)

**Target genes:** Classical tumor suppressors — TP53, PTEN, APC, CDH1, CDKN2A, RB1, SMAD4, ARID1A — and all MMR genes (MLH1, MSH2, MSH6, PMS2).

**Preferred outcome:** Per-patient double-hit calls; cohort-level frequency of biallelic inactivation per gene; identification of Lynch syndrome cases (germline MMR + somatic second hit).

**Why:** A single heterozygous mutation in a tumor suppressor is rarely sufficient for loss of function. The two-hit model (Knudson 1971) predicts that functional inactivation requires both alleles to be compromised. Treating somatic mutations in isolation — as most MAF-level analyses do — underestimates the true burden of complete tumor suppressor loss. This analysis leverages the unique matched germline+somatic design of this cohort to identify patients with confirmed biallelic inactivation, which has direct implications for prognosis, familial risk (Lynch syndrome), and targeted therapy eligibility (e.g., PARP inhibitors for BRCA/HRD-related double hits).

---

## 16. Pathway-Level Multi-Hit Analysis

**What:** Map all somatic (and germline) mutations onto curated cancer pathways (Wnt/β-catenin, RTK/RAS, PI3K/AKT/mTOR, TGF-β, cell cycle/RB, DNA damage repair, chromatin remodeling, Hippo). For each patient, count the number of distinct genes hit per pathway. Identify patients with 2+ hits in the same pathway ("pathway double/triple hit") even if no single gene is hit twice.

**Triple-hit example:** A patient with ARID1A somatic mutation + SMARCA4 somatic mutation + PBRM1 germline variant = triple hit in the chromatin remodeling pathway.

**Preferred outcome:** Pathway-level multi-hit frequency matrix (patients × pathways); patient subgroups defined by which pathways are maximally disrupted; association of pathway multi-hit burden with MSI, TMB, TRG, and OS.

**Why:** Tumor suppressor genes within the same pathway are functionally redundant — hitting multiple members of the same pathway compounds functional loss even without biallelic inactivation of any single gene. Pathway-level multi-hit burden better captures the degree of pathway disruption than gene-level mutation counts alone. This is especially relevant for pathways like chromatin remodeling (ARID1A, ARID1B, SMARCA4) and Wnt (APC, CTNNB1, RNF43, AXIN1/2), where gastric cancer accumulates multiple hits across pathway members.

---

## 17. Biallelic MMR Inactivation and Lynch Syndrome Stratification

**What:** Specifically for the four MMR genes (MLH1, MSH2, MSH6, PMS2): identify patients with (a) germline pathogenic variant alone (Lynch suspect), (b) germline + somatic second hit (confirmed biallelic loss), or (c) two somatic hits (sporadic biallelic). Cross-validate against clinical MSI/dMMR IHC status.

**Preferred outcome:** Prevalence of germline Lynch syndrome in this Chinese gastric cancer cohort; concordance rate between genomic biallelic MMR inactivation and IHC dMMR calls; identification of patients with dMMR by genomics but pMMR by IHC (discordant cases for clinical follow-up).

**Why:** Lynch syndrome is the most common hereditary gastric cancer predisposition and is actionable (surveillance, cascade testing of relatives, immunotherapy eligibility). IHC can produce false-negative MMR results (especially MLH1 missense mutations that preserve protein expression but abolish function). Genomic biallelic MMR inactivation is orthogonal evidence that resolves IHC discordance and directly informs clinical management. This analysis builds on #9 (germline) and #15 (double-hit) but focuses specifically on the MMR axis given its direct link to MSI, immunotherapy, and hereditary risk.

---

## Summary Table

| # | Analysis | Primary Data | Output |
|---|---|---|---|
| 1 | Mutation landscape / oncoprint | `somt_hg38.maf` | Top mutated genes, frequencies |
| 2 | SMG detection | `somt_hg38.maf` | Driver gene list with q-values |
| 3 | Mutational signatures | `somt_hg38.maf` | Per-patient SBS signature fractions |
| 4 | TMB distribution & correlation | clinical + VCFs | TMB cutoffs, clinical associations |
| 5 | MSI/dMMR classifier | clinical + `somt.maf` | Sequence-based MSI predictor |
| 6 | Co-mutation / mutual exclusivity | `somt_hg38.maf` | Gene pair heatmap |
| 7 | Molecular subtyping | clinical + MAF | EBV/MSI/GS/CIN subtype calls |
| 8 | Treatment response prediction | neoadjuvant clinical + MAF | TRG predictor, AUC |
| 9 | Germline pathogenic variants | `germ.maf` | Carrier frequency, candidate genes |
| 10 | Clonality / VAF distribution | per-patient VCFs | Clonal fraction per patient |
| 11 | Cross-cohort FDU vs. TCGA | combined MAF | Population-level differences |
| 12 | Survival analysis | clinical + MAF | Prognostic genomic features |
| 13 | HER2 genomic context | clinical + MAF | HER2-associated co-mutations |
| 14 | Naive vs. neoadjuvant comparison | clinical + MAF | Baseline genomic equivalence test |
| 15 | Germline + somatic double-hit (Knudson) | `germ.maf` + somatic VCFs | Biallelic inactivation calls per patient/gene |
| 16 | Pathway-level multi-hit analysis | somatic + germline MAFs | Pathway disruption matrix, multi-hit subgroups |
| 17 | Biallelic MMR / Lynch syndrome stratification | `germ.maf` + somatic VCFs + clinical | Lynch prevalence, IHC vs. genomic concordance |
