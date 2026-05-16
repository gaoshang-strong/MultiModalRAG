# Target Panel DNA-seq — FDU Gastric Cancer Cohort

Targeted panel DNA sequencing data from gastric cancer patients collected at Fudan University hospital. The panel covers hundreds of cancer-relevant genes.

---

## Cohort Overview

- **~511 patients** with matched clinical data (530 raw VCF files exist; filtered subset used downstream)
- **Sample ID convention**
  - Tumor tissue: `GZ24XXXXXF01` (`F` = formalin-fixed tissue)
  - Blood / germline: `GZ24XXXXXB01` (`B` = blood)
- Two patient groups based on treatment history:
  - **Treatment-naive**: surgery-first, no neoadjuvant therapy
  - **Neoadjuvant-treated**: received pre-operative chemotherapy ± immunotherapy

---

## Directory Structure

```
Target_panel_DNA_seq_FDU/
└── MAF_data/
    ├── gc_vcfs/
    │   ├── somt_vcf/          # Per-patient somatic variant TSVs (530 files)
    │   ├── germ_vcf/          # Per-patient germline variant TSVs (530 files)
    │   ├── somt.maf           # All-patient merged somatic MAF (hg19)
    │   ├── somt_hg38.maf      # All-patient merged somatic MAF (hg38-lifted)
    │   ├── germ.maf           # All-patient merged germline MAF
    │   └── combined_TCGA_PANEL_hg38_MAF.csv.gz  # FDU panel + TCGA-STAD combined MAF (hg38)
    │
    ├── clinical_with_DNA_seq.csv
    ├── clinical_with_DNA_seq_with_full_TMB.csv
    ├── clinical_with_DNA_seq_with_full_TMB_filtered_samples.csv
    ├── all_genes_detected_with_somatic_variants.csv
    ├── all_genes_detected_with_somatic_variants_with_corrected_gene_names.csv
    ├── 2025.12.29 patients without treatment before surgery.csv
    ├── 2025.12.29 patients with treatment before surgery.csv
    ├── TCGA_STAD.maf
    ├── hg38.chrom.sizes
    └── gencode.v46.annotation.gtf
```

---

## Per-Patient Variant Files (`gc_vcfs/`)

### Somatic variants (`somt_vcf/`)

One file per patient, named `{sample_id}_mutation[.hgvs]_filter.tsv`.

- Produced by paired tumor-vs-blood variant calling (Mutect2 or equivalent)
- Each file has **two genotype columns**: blood normal first, then tumor
- Two naming conventions reflect pipeline versions:
  - Older: `*_mutation_filter.tsv`
  - Newer: `*_mutation.hgvs_filter.tsv` (adds HGVS coordinate column)

### Germline variants (`germ_vcf/`)

One file per patient, named `{sample_id}_mutation[.hgvs]_filter.tsv`.

- Called from blood sample only
- One genotype column (blood)

### Column structure (both somatic and germline TSVs)

Each file is a VCF-derived TSV annotated with ANNOVAR. Key column groups:

| Group | Columns |
|---|---|
| Genomic coordinates | `#CHROM`, `POS`, `REF`, `ALT`, `FILTER` |
| Genotype / allele | `GT`, `AD`, `AF`, `DP` (within FORMAT fields) |
| Functional annotation | `Func.refGene`, `ExonicFunc.refGene`, `AAChange.refGene` (also knownGene, ensGene equivalents) |
| Variant identity | `Gene.refGene`, `cytoBand`, `avsnp150` (dbSNP), `cosmic99` |
| ClinVar | `CLNDN`, `CLNSIG`, `CLNREVSTAT` |
| Population frequency | `ExAC_ALL/AFR/AMR/EAS/FIN/NFE/OTH/SAS`, `1000g2015aug_*`, `esp6500siv2_all` |
| Pathogenicity scores | SIFT, PolyPhen2 (HDIV/HVAR), LRT, MutationTaster, MutationAssessor, FATHMM, CADD, GERP++, phyloP, SiPhy |
| HGVS (newer files) | `transcript`, `gene`, `strand`, `coordinates(gDNA/cDNA/protein)`, `region`, `mut_type`, `aachange` |
| Tumor metrics (somatic) | `TLOD`, `NALOD`, `NLOD` (Mutect2 log-odds scores) |

---

## Aggregated MAF Files

| File | Variants | Notes |
|---|---|---|
| `somt.maf` | ~8,909 | All-patient somatic, hg19 coordinates |
| `somt_hg38.maf` | ~8,909 | hg38-lifted version for cross-cohort analysis |
| `germ.maf` | ~12,996 | All-patient germline |
| `combined_TCGA_PANEL_hg38_MAF.csv.gz` | — | FDU panel + TCGA-STAD merged on hg38; use for comparative analysis |
| `TCGA_STAD.maf` | — | TCGA STAD somatic MAF (reference cohort) |

Standard MAF columns: `Hugo_Symbol`, `Chromosome`, `Start_Position`, `End_Position`, `Reference_Allele`, `Tumor_Seq_Allele2`, `Variant_Classification`, `Variant_Type`, `Tumor_Sample_Barcode`, `tx`, `exon`, `txChange`, `aaChange`.

---

## Clinical Files

### `clinical_with_DNA_seq.csv`

Base clinical table linked to sequencing data. One row per patient.

| Column | Description |
|---|---|
| `sample_id` | Patient ID (e.g. `GZ2400001`) |
| `MMR` | Mismatch repair status: `pMMR` or `dMMR` |
| `HER2` | HER2 IHC score: `0`, `1+`, `2+`, `3+` |
| `TMB` | Tumor mutation burden (reported from clinical report) |
| `MSI` | MSI status: `MSS` or `MSI-H` |
| `y_msi` | Binary label: 1 = MSI-H |
| `y_dmmr` | Binary label: 1 = dMMR |
| `y_any` | Binary label: 1 = MSI-H or dMMR |
| `logTMB` | Natural log of TMB |
| `vcf_path` | Path to per-patient somatic VCF |

### `clinical_with_DNA_seq_with_full_TMB.csv`

Extends the base table with recalculated TMB from the raw variant files:

| Added column | Description |
|---|---|
| `num_somatic_vars` | Count of somatic variants in the VCF |
| `calculated_TMB` | TMB recalculated from panel-adjusted variant count |

### `clinical_with_DNA_seq_with_full_TMB_filtered_samples.csv`

Quality-filtered subset of the above — use this for primary analysis.

### `2025.12.29 patients without treatment before surgery.csv`

Treatment-naive patients (surgery-first). Columns (Chinese headers):

- Demographics: age, sex
- Pathology: pT, pN, pathological stage, lymph node count, primary site, tumor size, differentiation, histological type
- IHC markers: CDX2, MUC6, KI-67 (%), claudin18.2, EBER
- Molecular: MMR, HER2, TMB, MSI (numeric index)
- Outcomes: surgery date, adjuvant treatment regimen, progression date, OS

### `2025.12.29 patients with treatment before surgery.csv`

Neoadjuvant-treated patients. Same base columns plus:

- Pre-op staging: cT, cN, cM, clinical stage
- Neoadjuvant: start date, regimen (e.g. SOX + immunotherapy), number of cycles, response assessment
- **TRG score** (tumor regression grade)
- **CPS** (PD-L1 combined positive score)
- Post-op: adjuvant regimen, progression, OS

### `all_genes_detected_with_somatic_variants[_with_corrected_gene_names].csv`

Simple list of all gene symbols that carry at least one somatic variant across the cohort. The `_corrected` version has standardized/updated HGNC symbols.

---

## Reference Files

| File | Description |
|---|---|
| `hg38.chrom.sizes` | Chromosome sizes for hg38 (used in liftover / bedtools operations) |
| `gencode.v46.annotation.gtf` | GENCODE v46 gene annotation (hg38) |
