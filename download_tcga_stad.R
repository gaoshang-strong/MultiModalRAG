library(TCGAbiolinks)

BASE_DIR <- "/ShangGaoAIProjects/Gastric/genome/TCGA-STAD"
PROJECT  <- "TCGA-STAD"

# ── 1. Clinical ───────────────────────────────────────────────────────────────
message("=== Downloading clinical data ===")
query_clin <- GDCquery(
  project     = PROJECT,
  data.category = "Clinical",
  data.type   = "Clinical Supplement",
  data.format = "BCR Biotab"
)
GDCdownload(query_clin, directory = file.path(BASE_DIR, "clinical"))
clin_list <- GDCprepare(query_clin, directory = file.path(BASE_DIR, "clinical"))
# save each clinical table as its own CSV
for (tbl_name in names(clin_list)) {
  out <- file.path(BASE_DIR, "clinical", paste0(tbl_name, ".csv"))
  write.csv(clin_list[[tbl_name]], out, row.names = FALSE)
}
message("Clinical saved.")

# ── 2. RNA-seq ────────────────────────────────────────────────────────────────
message("=== Querying RNA-seq data ===")
query_rna <- GDCquery(
  project            = PROJECT,
  data.category      = "Transcriptome Profiling",
  data.type          = "Gene Expression Quantification",
  workflow.type      = "STAR - Counts",
  access             = "open"
)

GDCdownload(query_rna, directory = file.path(BASE_DIR, "rnaseq"))
rna_data <- GDCprepare(query_rna, directory = file.path(BASE_DIR, "rnaseq"))
saveRDS(rna_data, file.path(BASE_DIR, "rnaseq", "rnaseq_summarized_experiment.rds"))
message("RNA-seq saved.")

# ── 3. Somatic mutation ───────────────────────────────────────────────────────
message("=== Querying somatic mutation data ===")
query_mut <- GDCquery(
  project       = PROJECT,
  data.category = "Simple Nucleotide Variation",
  data.type     = "Masked Somatic Mutation",
  access        = "open"
)

GDCdownload(query_mut, directory = file.path(BASE_DIR, "somatic_mutation"))
maf <- GDCprepare(query_mut, directory = file.path(BASE_DIR, "somatic_mutation"))
write.csv(maf, file.path(BASE_DIR, "somatic_mutation", "somatic_mutation.csv"), row.names = FALSE)
message("Somatic mutation saved.")


# ── 4. DNA Methylation (450k array) ──────────────────────────────────────────
message("=== Querying DNA methylation data ===")
query_meth <- GDCquery(
  project       = PROJECT,
  data.category = "DNA Methylation",
  data.type     = "Methylation Beta Value",
  platform      = "Illumina Human Methylation 450",
  access        = "open"
)
GDCdownload(query_meth, directory = file.path(BASE_DIR, "methylation"))
meth_data <- GDCprepare(query_meth, directory = file.path(BASE_DIR, "methylation"))
saveRDS(meth_data, file.path(BASE_DIR, "methylation", "methylation_summarized_experiment.rds"))
message("DNA methylation saved.")

# ── 5. Copy Number Alterations ────────────────────────────────────────────────
message("=== Querying copy number alteration data ===")
query_cna <- GDCquery(
  project       = PROJECT,
  data.category = "Copy Number Variation",
  data.type     = "Gene Level Copy Number",
  access        = "open"
)
GDCdownload(query_cna, directory = file.path(BASE_DIR, "cna"))
cna_data <- GDCprepare(query_cna, directory = file.path(BASE_DIR, "cna"))
saveRDS(cna_data, file.path(BASE_DIR, "cna", "cna_summarized_experiment.rds"))
message("CNA saved.")

# ── 6. Gene Fusions ───────────────────────────────────────────────────────────
message("=== Querying gene fusion data ===")
query_fusion <- GDCquery(
  project       = PROJECT,
  data.category = "Transcriptome Profiling",
  data.type     = "Transcript Fusion",
  workflow.type = "STAR-Fusion",
  access        = "open"
)
GDCdownload(query_fusion, directory = file.path(BASE_DIR, "gene_fusion"))
fusion_data <- GDCprepare(query_fusion, directory = file.path(BASE_DIR, "gene_fusion"))
saveRDS(fusion_data, file.path(BASE_DIR, "gene_fusion", "gene_fusion.rds"))
message("Gene fusion saved.")

message("=== All downloads complete ===")
