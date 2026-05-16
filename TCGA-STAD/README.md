# TCGA-STAD — The Cancer Genome Atlas Stomach Adenocarcinoma

Multi-omic data for the TCGA stomach adenocarcinoma (STAD) cohort, downloaded via `TCGAbiolinks` from the GDC data portal. All data is open-access tier. This serves as the reference/validation cohort alongside the local FDU panel cohort.

---

## Cohort Overview

- **~443 patients** (BCR Biotab clinical supplement)
- **Sample barcode format**: `TCGA-{site}-{patient}-{sample type}` (e.g. `TCGA-BR-4184-01A`)
  - `01A/01B` = primary tumor; `11A` = adjacent normal
- **Sequencing platform**: WXS (whole exome sequencing) for somatic mutations; RNA-seq for expression
- **Genome build**: GRCh38 / hg38
- **Gene annotation**: GENCODE v36 (RNA-seq), GENCODE v46 GTF available in FDU folder

---

## Directory Structure

```
TCGA-STAD/
├── clinical/
│   ├── clinical_patient_stad.csv          # Master patient table (443 patients)
│   ├── clinical_drug_stad.csv             # Drug/chemotherapy records
│   ├── clinical_follow_up_v1.0_stad.csv   # Follow-up visits, vital status updates
│   ├── clinical_radiation_stad.csv        # Radiation therapy records
│   ├── clinical_nte_stad.csv              # New tumor events after initial treatment
│   ├── clinical_omf_v4.0_stad.csv         # Other malignancies
│   └── TCGA-STAD/                         # Raw BCR Biotab .txt files (original source)
│
├── somatic_mutation/
│   ├── somatic_mutation.csv               # All-patient merged MAF (prepared by TCGAbiolinks)
│   └── TCGA-STAD/                         # Per-patient masked somatic MAF .gz files (434 files)
│
├── rnaseq/
│   ├── rnaseq_summarized_experiment.rds   # SummarizedExperiment object (all patients)
│   └── TCGA-STAD/                         # Per-patient STAR augmented gene count TSVs (448 files)
│
├── methylation/
│   └── TCGA-STAD/                         # Per-patient 450k sesame level3 beta TSVs (397 files)
│                                          # Note: RDS not yet prepared (sesameData install needed)
│
├── copy_number/                           # Empty — not downloaded
└── genotype/                              # Empty — not downloaded
```

---

## Clinical Data

### `clinical_patient_stad.csv` — Master patient table

Three header rows (BCR format): row 1 = column names, row 2 = alternate names, row 3 = CDE IDs. Data starts at row 4.

Key columns:

| Column | Description |
|---|---|
| `bcr_patient_barcode` | TCGA patient ID (e.g. `TCGA-BR-4184`) |
| `histological_type` | Histological subtype (adenocarcinoma NOS, diffuse, intestinal, etc.) |
| `neoplasm_histologic_grade` | Tumor grade (G1–G4) |
| `gender` | Patient sex |
| `race` / `ethnicity` | Self-reported demographics |
| `pathologic_T` / `pathologic_N` / `pathologic_M` | TNM staging components |
| `pathologic_stage` | Overall AJCC pathological stage |
| `vital_status` | Alive / Dead |
| `days_to_death` | Days from diagnosis to death (if deceased) |
| `days_to_last_followup` | Days from diagnosis to last known contact |
| `person_neoplasm_cancer_status` | Tumor-free / With tumor at last contact |
| `residual_tumor` | R0/R1/R2 (completeness of resection) |
| `history_of_neoadjuvant_treatment` | Yes/No |
| `primary_lymph_node_presentation_assessment` | Lymph node status |
| `lymph_node_examined_count` / `number_of_lymphnodes_positive_by_he` | Lymph node counts |
| `country_of_procurement` / `city_of_procurement` | Geographic origin of sample |

### `clinical_drug_stad.csv`

One row per drug per patient. Key columns: `drug_name`, `therapy_type`, `days_to_drug_therapy_start`, `days_to_drug_therapy_end`, `measure_of_response`, `number_cycles`.

### `clinical_follow_up_v1.0_stad.csv`

Longitudinal follow-up records. Key columns: `vital_status`, `days_to_last_followup`, `days_to_death`, `person_neoplasm_cancer_status`, `new_tumor_event_after_initial_treatment`, `days_to_new_tumor_event`.

### `clinical_radiation_stad.csv`

Radiation therapy details per patient per treatment course.

### `clinical_nte_stad.csv`

New tumor event records (recurrence, metastasis, second primary) with location and timing.

---

## Somatic Mutation Data

### `somatic_mutation/somatic_mutation.csv`

All-patient merged MAF prepared by `TCGAbiolinks::GDCprepare()`. Standard GDC MAF format with VEP annotations. Key column groups:

| Group | Columns |
|---|---|
| Variant identity | `Hugo_Symbol`, `Chromosome`, `Start_Position`, `End_Position`, `Variant_Classification`, `Variant_Type` |
| Alleles | `Reference_Allele`, `Tumor_Seq_Allele1`, `Tumor_Seq_Allele2`, `HGVSc`, `HGVSp`, `HGVSp_Short` |
| Sample linkage | `Tumor_Sample_Barcode`, `Matched_Norm_Sample_Barcode`, `case_id` |
| Depth / counts | `t_depth`, `t_ref_count`, `t_alt_count`, `n_depth`, `n_ref_count`, `n_alt_count` |
| VEP consequence | `One_Consequence`, `Consequence`, `IMPACT`, `BIOTYPE`, `CANONICAL` |
| Transcript | `Transcript_ID`, `ENSP`, `RefSeq`, `MANE`, `Exon_Number`, `CDS_position`, `Protein_position` |
| Pathogenicity | `SIFT`, `PolyPhen` |
| Population freq | `gnomAD_AF` (+ AFR/AMR/ASJ/EAS/FIN/NFE/OTH/SAS), `gnomAD_non_cancer_AF`, `1000G_AF`, `ESP_AA_AF/ESP_EA_AF` |
| Databases | `COSMIC`, `CLIN_SIG`, `dbSNP_RS`, `hotspot` (boolean) |
| Callers | `callers` (e.g. `mutect2;varscan2`) — ensemble masked MAF |
| RNA support | `RNA_Support`, `RNA_depth`, `RNA_ref_count`, `RNA_alt_count` |
| GDC QC | `GDC_FILTER` |

### Per-patient MAF files (`somatic_mutation/TCGA-STAD/`)

434 gzipped MAF files, one per aliquot. Naming: `{uuid}.wxs.aliquot_ensemble_masked.maf.gz`. These are the GDC-harmonized masked somatic MAFs (germline variants removed). Produced by an ensemble of Mutect2 + VarScan2 callers.

---

## RNA-seq Data

### `rnaseq/rnaseq_summarized_experiment.rds`

`SummarizedExperiment` object containing all patients. Assays available:
- `unstranded` — raw counts (use for DESeq2)
- `stranded_first` / `stranded_second` — strand-specific counts
- `tpm_unstranded` — TPM-normalized expression
- `fpkm_unstranded` — FPKM-normalized expression

Rows = genes (ENSEMBL ID + gene name + gene type), columns = samples.

### Per-patient TSVs (`rnaseq/TCGA-STAD/`)

448 files named `*.rna_seq.augmented_star_gene_counts.tsv`. Each file contains:
- 4 summary rows: `N_unmapped`, `N_multimapping`, `N_noFeature`, `N_ambiguous`
- One row per gene with columns: `gene_id` (Ensembl), `gene_name`, `gene_type`, `unstranded`, `stranded_first`, `stranded_second`, `tpm_unstranded`, `fpkm_unstranded`
- Gene model: GENCODE v36
- ~60,660 gene entries including protein-coding, lncRNA, miRNA, pseudogene, etc.

---

## DNA Methylation Data (450k Array)

### Per-patient beta value files (`methylation/TCGA-STAD/`)

397 files named `*.methylation_array.sesame.level3betas.txt`. Each file contains per-CpG probe beta values (0–1) processed by the sesame pipeline.

**Note**: The merged `SummarizedExperiment` RDS has not yet been generated. To prepare it, `sesameData` must be installed first:
```r
BiocManager::install("sesameData")
query_meth <- GDCquery(project="TCGA-STAD", data.category="DNA Methylation",
                       data.type="Methylation Beta Value",
                       platform="Illumina Human Methylation 450", access="open")
meth_data <- GDCprepare(query_meth, directory="methylation/")
saveRDS(meth_data, "methylation/methylation_summarized_experiment.rds")
```

---

## Data Not Downloaded

| Data type | Status |
|---|---|
| Copy Number Variation (Gene Level) | Not downloaded (`copy_number/` is empty) |
| Genotype / SNP array | Not downloaded (`genotype/` is empty) |
| Gene Fusions (STAR-Fusion) | Not downloaded |

---

## Cross-Cohort Integration

The file `Target_panel_DNA_seq_FDU/MAF_data/gc_vcfs/combined_TCGA_PANEL_hg38_MAF.csv.gz` contains the FDU panel somatic MAF merged with this TCGA-STAD MAF on hg38 coordinates. Use this for direct comparison of mutation frequencies between cohorts.

See `Target_panel_DNA_seq_FDU/README.md` for the FDU cohort data structure.
