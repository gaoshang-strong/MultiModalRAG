#!/usr/bin/env Rscript
# Download MSK-IMPACT esophagogastric cancer cohorts from cBioPortal REST API
#
# Studies (MSK targeted panel, EGC — includes gastric + GEJ):
#   egc_msk_2017          n=341 samples / 305 pts   Metastatic EGC (Cancer Discovery 2017)   PMID: 29122777
#   egc_mskcc_2020        n=487 samples             Esophageal/Stomach (MSK, 2020)            PMID: 33795256
#   egc_msk_tp53_ccr_2022 n=237 samples             EGC TP53 (Clin Cancer Res 2022)           PMID: 35377946
#   egc_msk_2023          n=902 patients            Early/average-onset EGC (JNCI 2023)       PMID: 37699004
#
# Usage:
#   /home/sgao30/micromamba/envs/tcga_bioc/bin/Rscript \
#       Other_cBioPartal_stomach/download_msk_cohorts.R

suppressPackageStartupMessages({
  library(httr)
  library(jsonlite)
  library(data.table)
})

BASE      <- "https://www.cbioportal.org/api"
OUT_ROOT  <- "/ShangGaoAIProjects/Gastric/genome/Other_cBioPartal_stomach"
PAGE_SIZE <- 10000L

STUDIES <- list(
  list(id = "egc_msk_2017",          mut_profile = "egc_msk_2017_mutations"),
  list(id = "egc_mskcc_2020",        mut_profile = "egc_mskcc_2020_mutations"),
  list(id = "egc_msk_tp53_ccr_2022", mut_profile = "egc_msk_tp53_ccr_2022_mutations"),
  list(id = "egc_msk_2023",          mut_profile = "egc_msk_2023_mutations")
)

cb_get <- function(endpoint, query = list()) {
  r <- GET(paste0(BASE, endpoint), query = query, timeout(60))
  if (http_error(r)) { cat("  [HTTP", status_code(r), "]\n"); return(NULL) }
  fromJSON(content(r, "text", encoding = "UTF-8"), simplifyDataFrame = TRUE)
}

fetch_mutations <- function(study_id, mut_profile_id, sample_ids) {
  body_base <- list(
    sampleMolecularIdentifiers = lapply(sample_ids, function(sid)
      list(molecularProfileId = mut_profile_id, sampleId = sid))
  )
  r_meta <- POST(paste0(BASE, "/mutations/fetch?projection=META"),
                 body = toJSON(body_base, auto_unbox = TRUE),
                 content_type_json(), timeout(60))
  total <- as.integer(headers(r_meta)$`total-count`)
  cat(sprintf("  Total mutations: %d\n", total))

  pages <- ceiling(total / PAGE_SIZE)
  all_pages <- vector("list", pages)
  for (p in seq_len(pages)) {
    cat(sprintf("  Page %d / %d ...\n", p, pages))
    url <- sprintf("%s/mutations/fetch?projection=SUMMARY&pageSize=%d&pageNumber=%d",
                   BASE, PAGE_SIZE, p - 1L)
    r <- POST(url, body = toJSON(body_base, auto_unbox = TRUE),
              content_type_json(), timeout(120))
    if (http_error(r)) { cat("  [HTTP error on page", p, "]\n"); next }
    all_pages[[p]] <- fromJSON(content(r, "text", encoding = "UTF-8"),
                               simplifyDataFrame = TRUE)
  }
  rbindlist(all_pages, fill = TRUE)
}

fetch_clinical <- function(study_id, level = c("PATIENT", "SAMPLE")) {
  level <- match.arg(level)
  r <- GET(
    sprintf("%s/studies/%s/clinical-data", BASE, study_id),
    query = list(clinicalDataType = level, projection = "SUMMARY"),
    timeout(60)
  )
  if (http_error(r)) return(NULL)
  d <- fromJSON(content(r, "text", encoding = "UTF-8"), simplifyDataFrame = TRUE)
  if (is.null(d) || nrow(d) == 0) return(NULL)
  id_col <- if (level == "PATIENT") "patientId" else "sampleId"
  dcast(as.data.table(d), as.formula(paste(id_col, "~ clinicalAttributeId")),
        value.var = "value", fun.aggregate = function(x) x[1])
}

# ── Main loop ─────────────────────────────────────────────────────────────────
summary_rows <- list()

for (s in STUDIES) {
  sid <- s$id
  cat(sprintf("\n%s\n── %s ──\n%s\n", strrep("=", 60), sid, strrep("=", 60)))

  out_dir <- file.path(OUT_ROOT, sid)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  samp <- cb_get(sprintf("/studies/%s/samples", sid),
                 query = list(projection = "ID", pageSize = 10000))
  if (is.null(samp)) { cat("  [skip] cannot get samples\n"); next }
  sample_ids <- samp$sampleId
  cat(sprintf("  Samples: %d\n", length(sample_ids)))

  mut_file <- file.path(out_dir, "mutations.tsv")
  if (file.exists(mut_file)) {
    cat("  [skip] mutations.tsv already exists\n")
  } else {
    muts <- fetch_mutations(sid, s$mut_profile, sample_ids)
    if (!is.null(muts) && nrow(muts) > 0) {
      fwrite(muts, mut_file, sep = "\t")
      cat(sprintf("  Saved: mutations.tsv (%d rows, %d cols)\n", nrow(muts), ncol(muts)))
    }
  }

  clin_pat_file <- file.path(out_dir, "clinical_patient.tsv")
  if (!file.exists(clin_pat_file)) {
    cp <- fetch_clinical(sid, "PATIENT")
    if (!is.null(cp)) {
      fwrite(cp, clin_pat_file, sep = "\t")
      cat(sprintf("  Saved: clinical_patient.tsv (%d patients, %d attrs)\n", nrow(cp), ncol(cp)))
    }
  } else cat("  [skip] clinical_patient.tsv already exists\n")

  clin_samp_file <- file.path(out_dir, "clinical_sample.tsv")
  if (!file.exists(clin_samp_file)) {
    cs <- fetch_clinical(sid, "SAMPLE")
    if (!is.null(cs)) {
      fwrite(cs, clin_samp_file, sep = "\t")
      cat(sprintf("  Saved: clinical_sample.tsv (%d samples, %d attrs)\n", nrow(cs), ncol(cs)))
    }
  } else cat("  [skip] clinical_sample.tsv already exists\n")

  summary_rows[[sid]] <- data.table(
    study_id  = sid,
    n_samples = length(sample_ids)
  )
}

cat(sprintf("\n%s\nDownload summary\n%s\n", strrep("=", 60), strrep("=", 60)))
print(rbindlist(summary_rows))
cat(sprintf("\nAll outputs in: %s\n", OUT_ROOT))
