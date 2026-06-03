#!/usr/bin/env Rscript
# tp53_me_clinical_adjustment.R
#
# Tests TP53 ME with ATM and ARID1A after controlling for stage and histology.
# Addresses whether the ME signals survive adjustment for tumour biology (Lauren
# classification / histologic type) and disease stage.
#
# Models (per cohort per gene pair):
#   A: TP53 ~ gene_B + log2(TMB)                                  [base]
#   B: TP53 ~ gene_B + log2(TMB) + stage_advanced                 [+stage]
#   C: TP53 ~ gene_B + log2(TMB) + stage_advanced + hist_diffuse  [+stage+histology]
#   (Model C skipped for cohorts without Lauren / histotype data)
#
# Cohort coverage (from panel_coverage_check.tsv):
#   TP53–ATM:   FDU, TCGA-STAD, OncoSG-2018, HK-Pfizer-2014, MSK-2023
#   TP53–ARID1A: + TMUCIH-2015 (Model A/B only), MSK-2017
#   MSKCC-2020 and MSK-TP53-CCR-2022 excluded (genes not covered)
#
# Run:
#   /home/sgao30/micromamba/envs/tcga_bioc/bin/Rscript \
#     scripts/Phase2/tp53_me_clinical_adjustment.R
#
# Output: analysis_results/Phase2/clinical_covariate/

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(tibble)
  library(purrr)
  library(metafor)    # rma for pooling
})

# ── 0. Configuration ───────────────────────────────────────────────────────────

ROOT   <- "/ShangGaoAIProjects/Gastric/genome"
OUTDIR <- file.path(ROOT, "analysis_results/Phase2/clinical_covariate")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

GENE_PAIRS <- list(c("TP53", "ATM"), c("TP53", "ARID1A"))

FUNCTIONAL_CLASSES <- c(
  "Missense_Mutation", "Nonsense_Mutation",
  "Frame_Shift_Del",   "Frame_Shift_Ins",
  "Splice_Site",       "In_Frame_Del",
  "In_Frame_Ins",      "Nonstop_Mutation",
  "Translation_Start_Site"
)
EXOME_MB <- 38.0

ENTREZ_IDS <- c(TP53 = 7157L, ATM = 472L, ARID1A = 8289L)

# ── 1. Helpers ─────────────────────────────────────────────────────────────────

safe_glm <- function(d, formula_str) {
  tryCatch({
    fit <- glm(as.formula(formula_str), data = d, family = binomial)
    cf  <- coef(summary(fit))
    if (nrow(cf) < 2) stop("fewer than 2 coefficients — likely separation")
    b  <- cf[2, "Estimate"]
    se <- cf[2, "Std. Error"]
    p  <- cf[2, "Pr(>|z|)"]
    tibble(OR = exp(b), CI_lo = exp(b - 1.96*se), CI_hi = exp(b + 1.96*se),
           beta = b, se = se, p_value = p, n_complete = nrow(d), status = "ok")
  }, error = \(e) tibble(
    OR = NA_real_, CI_lo = NA_real_, CI_hi = NA_real_,
    beta = NA_real_, se = NA_real_, p_value = NA_real_,
    n_complete = nrow(d), status = str_trunc(conditionMessage(e), 80)
  ))
}

# Run Models A / B / C on a single cohort data frame.
# df must have columns: TP53, <gene_b>, log2_TMB,
#                       stage_advanced (integer or NA),
#                       hist_diffuse   (integer or NA)
run_models <- function(df, gene_b, cohort_id, seq_type) {
  has_stage <- sum(!is.na(df$stage_advanced)) >= 20
  has_hist  <- sum(!is.na(df$hist_diffuse))   >= 20

  map_dfr(c("A", "B", "C"), function(mod) {
    if (mod == "B" && !has_stage) return(NULL)
    if (mod == "C" && (!has_stage || !has_hist)) return(NULL)

    formula_str <- switch(mod,
      A = paste0("TP53 ~ ", gene_b, " + log2_TMB"),
      B = paste0("TP53 ~ ", gene_b, " + log2_TMB + stage_advanced"),
      C = paste0("TP53 ~ ", gene_b, " + log2_TMB + stage_advanced + hist_diffuse")
    )
    label <- switch(mod,
      A = "A: base (TMB only)",
      B = "B: + stage",
      C = "C: + stage + histology"
    )

    d <- switch(mod,
      A = df %>% filter(!is.na(TP53), !is.na(.data[[gene_b]]), !is.na(log2_TMB)),
      B = df %>% filter(!is.na(TP53), !is.na(.data[[gene_b]]), !is.na(log2_TMB),
                        !is.na(stage_advanced)),
      C = df %>% filter(!is.na(TP53), !is.na(.data[[gene_b]]), !is.na(log2_TMB),
                        !is.na(stage_advanced), !is.na(hist_diffuse))
    )

    if (sum(d[[gene_b]]) < 5) {
      return(tibble(cohort = cohort_id, seq_type = seq_type,
                    gene_b = gene_b, model = mod, model_label = label,
                    n = nrow(d), n_tp53 = NA, n_gene = NA, n_both = NA,
                    OR = NA_real_, CI_lo = NA_real_, CI_hi = NA_real_,
                    beta = NA_real_, se = NA_real_, p_value = NA_real_,
                    status = "too_few_mutations"))
    }

    res <- safe_glm(d, formula_str)

    tibble(cohort = cohort_id, seq_type = seq_type,
           gene_b = gene_b, model = mod, model_label = label,
           n      = res$n_complete,
           n_tp53 = sum(d$TP53),
           n_gene = sum(d[[gene_b]]),
           n_both = sum(d$TP53 == 1L & d[[gene_b]] == 1L),
           OR = res$OR, CI_lo = res$CI_lo, CI_hi = res$CI_hi,
           beta = res$beta, se = res$se, p_value = res$p_value,
           status = res$status)
  })
}

# ── 2. Data loaders ────────────────────────────────────────────────────────────

# ---- 2a. FDU (master table already has somatic_ prefix columns) ----
load_fdu <- function() {
  message("Loading FDU...")
  m <- read_tsv(
    file.path(ROOT, "analysis_results/Phase0/filtered_master_sample_table.tsv"),
    show_col_types = FALSE
  )

  m %>% transmute(
    TP53           = as.integer(somatic_TP53   > 0),
    ATM            = as.integer(somatic_ATM    > 0),
    ARID1A         = as.integer(somatic_ARID1A > 0),
    log2_TMB       = log2(pmax(TMB_report_numeric, 0.1)),
    # Stage: IA/IB/IIA/IIB = early(0); IIIA/IIIB/IIIC/IV = advanced(1)
    stage_advanced = case_when(
      path_stage %in% c("IIIA","IIIB","IIIC","IV")        ~ 1L,
      path_stage %in% c("IA","IB","IIA","IIB")             ~ 0L,
      TRUE ~ NA_integer_
    ),
    # Histology: 印戒/黏附性癌 = diffuse(1); 腺癌/粘液腺癌 = non-diffuse(0)
    hist_diffuse   = case_when(
      str_detect(histology, "印戒|黏附")                         ~ 1L,
      str_detect(histology, "腺癌|粘液") & !str_detect(histology, "印戒") ~ 0L,
      TRUE ~ NA_integer_
    )
  ) %>% filter(!is.na(log2_TMB), is.finite(log2_TMB))
}

# ---- 2b. TCGA-STAD ----
load_tcga <- function() {
  message("Loading TCGA-STAD...")
  maf  <- read_csv(file.path(ROOT, "TCGA-STAD/somatic_mutation/somatic_mutation.csv"),
                   show_col_types = FALSE)
  clin <- read_csv(file.path(ROOT, "TCGA-STAD/clinical/clinical_patient_stad.csv"),
                   show_col_types = FALSE)

  # Keep primary tumour samples only
  maf <- maf %>%
    filter(str_detect(Tumor_Sample_Barcode, "-01[AB]-")) %>%
    mutate(sample_id  = Tumor_Sample_Barcode,
           patient_id = substr(Tumor_Sample_Barcode, 1, 12))

  # Build full sample list FIRST (all primary-tumour samples in MAF)
  all_samples <- unique(maf$sample_id)

  tmb_df <- maf %>%
    filter(Variant_Classification %in% FUNCTIONAL_CLASSES) %>%
    count(sample_id, name = "n_mut") %>%
    mutate(TMB = n_mut / EXOME_MB)
  # Floor: samples with 0 coding mutations get TMB = 0.1 (kept in model)
  tmb_df <- tibble(sample_id = all_samples) %>%
    left_join(tmb_df, by = "sample_id") %>%
    mutate(TMB = pmax(replace_na(TMB, 0), 0.1))

  # Gene binary: start from full sample list so WT samples get 0
  gene_long <- maf %>%
    filter(Variant_Classification %in% FUNCTIONAL_CLASSES,
           Hugo_Symbol %in% names(ENTREZ_IDS)) %>%
    distinct(sample_id, patient_id, gene = Hugo_Symbol) %>%
    mutate(val = 1L) %>%
    pivot_wider(names_from = gene, values_from = val, values_fill = 0L)
  for (g in names(ENTREZ_IDS)) if (!g %in% names(gene_long)) gene_long[[g]] <- 0L

  gene_df <- tibble(sample_id  = all_samples,
                    patient_id = substr(all_samples, 1, 12)) %>%
    left_join(gene_long %>% select(sample_id, any_of(names(ENTREZ_IDS))),
              by = "sample_id") %>%
    mutate(across(any_of(names(ENTREZ_IDS)), ~ replace_na(., 0L)))

  # Clinical: filter out metadata rows, derive patient_id
  clin_clean <- clin %>%
    filter(!str_detect(bcr_patient_barcode, "^bcr_|^CDE_")) %>%
    mutate(patient_id = toupper(substr(bcr_patient_barcode, 1, 12)))

  df <- gene_df %>%
    left_join(tmb_df %>% select(sample_id, TMB), by = "sample_id") %>%
    left_join(clin_clean, by = "patient_id") %>%
    transmute(
      TP53, ATM, ARID1A,
      log2_TMB       = log2(TMB),
      stage_advanced = case_when(
        str_detect(ajcc_pathologic_tumor_stage, "Stage III|Stage IV") ~ 1L,
        str_detect(ajcc_pathologic_tumor_stage, "Stage I|Stage II")   ~ 0L,
        TRUE ~ NA_integer_
      ),
      hist_diffuse   = case_when(
        str_detect(histologic_diagnosis, regex("Diffuse|Signet Ring", ignore_case = TRUE)) ~ 1L,
        str_detect(histologic_diagnosis, regex("Intestinal|Tubular|Papillary|Mucinous|NOS",
                                               ignore_case = TRUE)) ~ 0L,
        TRUE ~ NA_integer_
      )
    )

  message("  → ", nrow(df), " samples")
  df
}

# ---- 2c. Generic cBioPortal loader ----
# Returns a tibble with TP53/ATM/ARID1A binary + log2_TMB
# Clinical covariate columns are added by the cohort-specific wrappers below.
load_cbio_base <- function(cohort_id, stad_only = TRUE) {
  base     <- file.path(ROOT, "Other_cBioPartal_stomach", cohort_id)
  mut      <- read_tsv(file.path(base, "mutations.tsv"),        show_col_types = FALSE, comment = "#")
  clin_sp  <- read_tsv(file.path(base, "clinical_sample.tsv"),  show_col_types = FALSE, comment = "#")

  if (stad_only && "CANCER_TYPE_DETAILED" %in% names(clin_sp)) {
    stad_ids <- clin_sp %>%
      filter(str_detect(CANCER_TYPE_DETAILED, regex("Stomach Adenocarcinoma", ignore_case = TRUE))) %>%
      pull(sampleId)
    mut     <- mut     %>% filter(sampleId %in% stad_ids)
    clin_sp <- clin_sp %>% filter(sampleId %in% stad_ids)
    message("  STAD filter: ", nrow(clin_sp), " samples")
  }

  # One sample per patient (prefer Primary)
  if ("patientId" %in% names(mut)) {
    sample_meta <- mut %>% distinct(sampleId, patientId) %>%
      left_join(clin_sp %>% select(sampleId, any_of("SAMPLE_TYPE")), by = "sampleId")
    if ("SAMPLE_TYPE" %in% names(sample_meta)) {
      keep <- sample_meta %>%
        mutate(pref = if_else(SAMPLE_TYPE == "Primary", 1L, 2L)) %>%
        group_by(patientId) %>% slice_min(pref, n = 1, with_ties = FALSE) %>%
        ungroup() %>% pull(sampleId)
    } else {
      keep <- sample_meta %>% group_by(patientId) %>% slice(1) %>%
        ungroup() %>% pull(sampleId)
    }
    mut     <- mut     %>% filter(sampleId %in% keep)
    clin_sp <- clin_sp %>% filter(sampleId %in% keep)
    message("  Dedup: ", length(keep), " samples (1 per patient)")
  }

  # TMB
  tmb_df <- if ("TMB_NONSYNONYMOUS" %in% names(clin_sp)) {
    clin_sp %>% select(sampleId, TMB = TMB_NONSYNONYMOUS) %>% mutate(TMB = as.numeric(TMB))
  } else {
    mut %>% filter(mutationType %in% FUNCTIONAL_CLASSES) %>%
      count(sampleId, name = "n_mut") %>%
      mutate(TMB = n_mut / EXOME_MB) %>%
      select(sampleId, TMB)
  }

  # Gene binary matrix
  gene_map <- enframe(ENTREZ_IDS, name = "gene", value = "entrezGeneId")
  mut_genes <- mut %>%
    filter(mutationType %in% FUNCTIONAL_CLASSES, entrezGeneId %in% ENTREZ_IDS) %>%
    left_join(gene_map, by = "entrezGeneId") %>%
    distinct(sampleId, gene) %>%
    mutate(val = 1L) %>%
    pivot_wider(names_from = gene, values_from = val, values_fill = 0L)

  base_df <- tibble(sampleId = unique(clin_sp$sampleId)) %>%
    left_join(mut_genes, by = "sampleId") %>%
    mutate(across(any_of(names(ENTREZ_IDS)), ~ replace_na(., 0L)))
  for (g in names(ENTREZ_IDS)) if (!g %in% names(base_df)) base_df[[g]] <- 0L

  base_df %>%
    left_join(tmb_df, by = "sampleId") %>%
    filter(!is.na(TMB), TMB > 0) %>%
    mutate(log2_TMB = log2(TMB))
}

# ---- 2d. Per-cohort wrappers: join clinical covariates ----

load_oncosg <- function() {
  message("Loading OncoSG 2018...")
  base <- load_cbio_base("stad_oncosg_2018", stad_only = FALSE)
  clin_pt <- read_tsv(
    file.path(ROOT, "Other_cBioPartal_stomach/stad_oncosg_2018/clinical_patient.tsv"),
    show_col_types = FALSE, comment = "#"
  )
  mut <- read_tsv(
    file.path(ROOT, "Other_cBioPartal_stomach/stad_oncosg_2018/mutations.tsv"),
    show_col_types = FALSE, comment = "#"
  )
  sp_map <- mut %>% distinct(sampleId, patientId) %>% filter(sampleId %in% base$sampleId)

  cov <- clin_pt %>%
    transmute(
      patientId,
      stage_advanced = case_when(
        STAGE %in% c("III","IIIA","IIIB","IIIC","IV") ~ 1L,
        STAGE %in% c("I","IA","IB","II","IIA","IIB")  ~ 0L,
        TRUE ~ NA_integer_
      ),
      hist_diffuse = case_when(
        LAURENS_CLASSIFICATION == "Diffuse"                 ~ 1L,
        LAURENS_CLASSIFICATION %in% c("Intestinal","Mixed") ~ 0L,
        TRUE ~ NA_integer_
      )
    )

  base %>%
    left_join(sp_map, by = "sampleId") %>%
    left_join(cov, by = "patientId") %>%
    select(TP53, ATM, ARID1A, log2_TMB, stage_advanced, hist_diffuse)
}

load_hku_pfizer <- function() {
  message("Loading HK Pfizer 2014...")
  base <- load_cbio_base("stad_pfizer_uhongkong", stad_only = FALSE)
  clin_pt <- read_tsv(
    file.path(ROOT, "Other_cBioPartal_stomach/stad_pfizer_uhongkong/clinical_patient.tsv"),
    show_col_types = FALSE, comment = "#"
  )
  mut <- read_tsv(
    file.path(ROOT, "Other_cBioPartal_stomach/stad_pfizer_uhongkong/mutations.tsv"),
    show_col_types = FALSE, comment = "#"
  )
  sp_map <- mut %>% distinct(sampleId, patientId) %>% filter(sampleId %in% base$sampleId)

  cov <- clin_pt %>%
    transmute(
      patientId,
      stage_advanced = case_when(
        UICC_TUMOR_STAGE %in% c("III","IIIA","IIIB","IV") ~ 1L,
        UICC_TUMOR_STAGE %in% c("I","IA","IB","II")       ~ 0L,
        TRUE ~ NA_integer_
      ),
      hist_diffuse = case_when(
        str_to_lower(PATHOLOGY_LAUREN) == "diffuse"                 ~ 1L,
        str_to_lower(PATHOLOGY_LAUREN) %in% c("intestinal","mixed") ~ 0L,
        TRUE ~ NA_integer_
      )
    )

  base %>%
    left_join(sp_map, by = "sampleId") %>%
    left_join(cov, by = "patientId") %>%
    select(TP53, ATM, ARID1A, log2_TMB, stage_advanced, hist_diffuse)
}

load_tmucih <- function() {
  message("Loading TMUCIH 2015...")
  base <- load_cbio_base("egc_tmucih_2015", stad_only = FALSE)
  clin_pt <- read_tsv(
    file.path(ROOT, "Other_cBioPartal_stomach/egc_tmucih_2015/clinical_patient.tsv"),
    show_col_types = FALSE, comment = "#"
  )
  mut <- read_tsv(
    file.path(ROOT, "Other_cBioPartal_stomach/egc_tmucih_2015/mutations.tsv"),
    show_col_types = FALSE, comment = "#"
  )
  sp_map <- mut %>% distinct(sampleId, patientId) %>% filter(sampleId %in% base$sampleId)

  cov <- clin_pt %>%
    transmute(
      patientId,
      stage_advanced = case_when(
        as.character(STAGE) %in% c("3","4","III","IV") ~ 1L,
        as.character(STAGE) %in% c("1","2","I","II")   ~ 0L,
        TRUE ~ NA_integer_
      ),
      hist_diffuse = NA_integer_  # no Lauren data available
    )

  base %>%
    left_join(sp_map, by = "sampleId") %>%
    left_join(cov, by = "patientId") %>%
    select(TP53, ATM, ARID1A, log2_TMB, stage_advanced, hist_diffuse)
}

load_msk2017 <- function() {
  message("Loading MSK 2017...")
  base    <- load_cbio_base("egc_msk_2017", stad_only = TRUE)
  clin_sp <- read_tsv(
    file.path(ROOT, "Other_cBioPartal_stomach/egc_msk_2017/clinical_sample.tsv"),
    show_col_types = FALSE, comment = "#"
  )
  clin_pt <- read_tsv(
    file.path(ROOT, "Other_cBioPartal_stomach/egc_msk_2017/clinical_patient.tsv"),
    show_col_types = FALSE, comment = "#"
  )
  mut <- read_tsv(
    file.path(ROOT, "Other_cBioPartal_stomach/egc_msk_2017/mutations.tsv"),
    show_col_types = FALSE, comment = "#"
  )
  sp_map <- mut %>% distinct(sampleId, patientId) %>% filter(sampleId %in% base$sampleId)

  # Stage from clinical_sample (STAGE_AT_DIAGNOSIS — messy format)
  stage_sp <- clin_sp %>%
    filter(sampleId %in% base$sampleId) %>%
    transmute(
      sampleId,
      stage_advanced = case_when(
        str_detect(STAGE_AT_DIAGNOSIS, regex("^IV|^Stage IV", ignore_case = TRUE))  ~ 1L,
        str_detect(STAGE_AT_DIAGNOSIS, regex("^III",          ignore_case = TRUE))  ~ 1L,
        str_detect(STAGE_AT_DIAGNOSIS, regex("^I[AB]?$|^II[AB]?$|^Stage I[^I]|^Stage II[^I]",
                                             ignore_case = TRUE))                   ~ 0L,
        TRUE ~ NA_integer_
      )
    )

  # Lauren from clinical_patient
  hist_pt <- clin_pt %>%
    transmute(
      patientId,
      hist_diffuse = case_when(
        str_to_lower(str_trim(LAUREN_CLASS)) %in% c("diffuse","poorly cohesive") ~ 1L,
        str_to_lower(str_trim(LAUREN_CLASS)) %in% c("intestinal","mixed")        ~ 0L,
        TRUE ~ NA_integer_
      )
    )

  base %>%
    left_join(stage_sp %>% select(sampleId, stage_advanced), by = "sampleId") %>%
    left_join(sp_map, by = "sampleId") %>%
    left_join(hist_pt, by = "patientId") %>%
    select(TP53, ATM, ARID1A, log2_TMB, stage_advanced, hist_diffuse)
}

load_msk2023 <- function() {
  message("Loading MSK 2023...")
  base    <- load_cbio_base("egc_msk_2023", stad_only = TRUE)
  clin_pt <- read_tsv(
    file.path(ROOT, "Other_cBioPartal_stomach/egc_msk_2023/clinical_patient.tsv"),
    show_col_types = FALSE, comment = "#"
  )
  mut <- read_tsv(
    file.path(ROOT, "Other_cBioPartal_stomach/egc_msk_2023/mutations.tsv"),
    show_col_types = FALSE, comment = "#"
  )
  sp_map <- mut %>% distinct(sampleId, patientId) %>% filter(sampleId %in% base$sampleId)

  cov <- clin_pt %>%
    transmute(
      patientId,
      stage_advanced = case_when(
        STAGE %in% c("III","IV") ~ 1L,
        STAGE %in% c("I","II")   ~ 0L,
        TRUE ~ NA_integer_
      ),
      hist_diffuse = case_when(
        HISTOLOGY == "Signet_Diffuse" ~ 1L,
        HISTOLOGY == "Adenocarcinoma" ~ 0L,
        TRUE ~ NA_integer_
      )
    )

  base %>%
    left_join(sp_map, by = "sampleId") %>%
    left_join(cov, by = "patientId") %>%
    select(TP53, ATM, ARID1A, log2_TMB, stage_advanced, hist_diffuse)
}

# ── 3. Main analysis ───────────────────────────────────────────────────────────

cohort_list <- list(
  list(id = "FDU",         data = load_fdu(),        seq_type = "panel"),
  list(id = "TCGA_STAD",   data = load_tcga(),        seq_type = "WXS"),
  list(id = "OncoSG_2018", data = load_oncosg(),      seq_type = "WXS"),
  list(id = "HK_Pfizer",   data = load_hku_pfizer(),  seq_type = "WXS"),
  list(id = "TMUCIH_2015", data = load_tmucih(),      seq_type = "WXS"),
  list(id = "MSK_2017",    data = load_msk2017(),     seq_type = "panel"),
  list(id = "MSK_2023",    data = load_msk2023(),     seq_type = "panel")
)

# 3a. Covariate availability summary
cov_summary <- map_dfr(cohort_list, function(co) {
  d <- co$data
  tibble(
    cohort        = co$id,
    seq_type      = co$seq_type,
    n_total       = nrow(d),
    n_tp53_mut    = sum(d$TP53, na.rm = TRUE),
    n_atm_mut     = sum(d$ATM,  na.rm = TRUE),
    n_arid1a_mut  = sum(d$ARID1A, na.rm = TRUE),
    n_stage_avail = sum(!is.na(d$stage_advanced)),
    pct_advanced  = round(100 * mean(d$stage_advanced, na.rm = TRUE), 1),
    n_hist_avail  = sum(!is.na(d$hist_diffuse)),
    pct_diffuse   = round(100 * mean(d$hist_diffuse, na.rm = TRUE), 1)
  )
})
write_tsv(cov_summary, file.path(OUTDIR, "covariate_availability.tsv"))
cat("\n=== Covariate availability ===\n")
print(as.data.frame(cov_summary))

# 3b. Per-cohort logistic regression
per_cohort_rows <- map_dfr(cohort_list, function(co) {
  map_dfr(GENE_PAIRS, function(pair) {
    gene_b <- pair[2]
    if (!gene_b %in% names(co$data)) return(NULL)
    run_models(co$data, gene_b, co$id, co$seq_type)
  })
})
write_tsv(per_cohort_rows, file.path(OUTDIR, "per_cohort_clinical_results.tsv"))

# ── 4. Random-effects meta-analysis ───────────────────────────────────────────

meta_rows <- map_dfr(GENE_PAIRS, function(pair) {
  gb <- pair[2]
  map_dfr(c("A","B","C"), function(mod) {

    sub <- per_cohort_rows %>%
      filter(gene_b == gb, model == mod, status == "ok",
             !is.na(beta), !is.na(se), se > 0, se < 10, n_gene >= 5)

    if (nrow(sub) < 2) {
      return(tibble(gene_b = gb, model = mod,
                    model_label = paste0(mod, ": insufficient cohorts"),
                    n_cohorts = nrow(sub), pooled_OR = NA_real_,
                    CI_lo = NA_real_, CI_hi = NA_real_,
                    p_value = NA_real_, I2 = NA_real_, tau2 = NA_real_,
                    cohorts_included = NA_character_, status = "insufficient"))
    }

    m <- tryCatch(
      rma(yi = sub$beta, sei = sub$se, method = "REML"),
      error = \(e) NULL
    )
    if (is.null(m)) {
      return(tibble(gene_b = gb, model = mod,
                    model_label = sub$model_label[1], n_cohorts = nrow(sub),
                    pooled_OR = NA_real_, CI_lo = NA_real_, CI_hi = NA_real_,
                    p_value = NA_real_, I2 = NA_real_, tau2 = NA_real_,
                    cohorts_included = paste(sub$cohort, collapse=","), status = "meta_error"))
    }

    tibble(
      gene_b            = gb,
      model             = mod,
      model_label       = sub$model_label[1],
      n_cohorts         = m$k,
      pooled_OR         = exp(as.numeric(m$b)),
      CI_lo             = exp(m$ci.lb),
      CI_hi             = exp(m$ci.ub),
      p_value           = as.numeric(m$pval),
      I2                = round(m$I2, 1),
      tau2              = m$tau2,
      cohorts_included  = paste(sub$cohort, collapse = ","),
      status            = "ok"
    )
  })
})
write_tsv(meta_rows, file.path(OUTDIR, "meta_clinical_results.tsv"))

# ── 5. Summary comparison: A vs B vs C ────────────────────────────────────────

summary_df <- meta_rows %>%
  filter(status == "ok") %>%
  select(gene_b, model, pooled_OR, p_value, I2, n_cohorts) %>%
  pivot_wider(
    names_from  = model,
    values_from = c(pooled_OR, p_value, I2, n_cohorts),
    names_glue  = "{.value}_{model}"
  ) %>%
  mutate(
    sig_A = !is.na(p_value_A) & p_value_A < 0.05,
    sig_B = !is.na(p_value_B) & p_value_B < 0.05,
    sig_C = !is.na(p_value_C) & p_value_C < 0.05,
    OR_ratio_A_to_C = round(pooled_OR_A / pooled_OR_C, 3),
    conclusion = case_when(
      is.na(p_value_C)   ~ "C not available",
      sig_A & sig_C      ~ "robust: signal survives clinical adjustment",
      sig_A & !sig_C     ~ "CAUTION: lost after clinical adjustment",
      !sig_A             ~ "not significant in base model",
      TRUE               ~ "other"
    )
  )
write_tsv(summary_df, file.path(OUTDIR, "clinical_adjustment_summary.tsv"))

cat("\n=== Summary (Model A vs C) ===\n")
print(as.data.frame(summary_df %>%
  select(gene_b, pooled_OR_A, p_value_A, pooled_OR_C, p_value_C, OR_ratio_A_to_C, conclusion)))

# ── 6. Forest plots ────────────────────────────────────────────────────────────

MODEL_COLORS <- c(
  "A" = "#1f77b4",  # blue
  "B" = "#ff7f0e",  # orange
  "C" = "#2ca02c"   # green
)

for (gb in c("ATM","ARID1A")) {
  plot_data <- per_cohort_rows %>%
    filter(gene_b == gb, status == "ok", !is.na(OR)) %>%
    mutate(
      label  = paste0(cohort, " [", seq_type, "]"),
      model_f = factor(model, c("A","B","C"),
                       c("A: base","B: +stage","C: +stage+hist"))
    )

  meta_plot <- meta_rows %>%
    filter(gene_b == gb, status == "ok") %>%
    transmute(
      label   = "POOLED (RE meta)",
      OR      = pooled_OR, CI_lo, CI_hi,
      n       = NA_integer_, p_value,
      model   = model,
      model_f = factor(model, c("A","B","C"),
                       c("A: base","B: +stage","C: +stage+hist")),
      is_meta = TRUE
    )

  if (nrow(plot_data) == 0) next

  all_labels <- c(sort(unique(plot_data$label)), "POOLED (RE meta)")

  p <- ggplot(
    plot_data %>% mutate(is_meta = FALSE, label = factor(label, all_labels)),
    aes(x = OR, y = label, color = model_f, shape = model_f)
  ) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = CI_lo, xmax = CI_hi), height = 0.25,
                   position = position_dodge(0.6)) +
    geom_point(size = 2.5, position = position_dodge(0.6)) +
    # Pooled meta diamonds
    geom_errorbarh(
      data = meta_plot %>% mutate(label = factor("POOLED (RE meta)", all_labels)),
      aes(x = OR, xmin = CI_lo, xmax = CI_hi, color = model_f),
      height = 0.35, position = position_dodge(0.6), linewidth = 1
    ) +
    geom_point(
      data = meta_plot %>% mutate(label = factor("POOLED (RE meta)", all_labels)),
      aes(x = OR, color = model_f),
      shape = 18, size = 4, position = position_dodge(0.6)
    ) +
    scale_x_log10(breaks = c(0.1, 0.2, 0.5, 1, 2, 5),
                  labels  = c("0.1","0.2","0.5","1","2","5"),
                  limits  = c(0.05, 6)) +
    scale_color_manual(values = MODEL_COLORS,
                       labels = c("A: base","B: +stage","C: +stage+hist")) +
    scale_shape_manual(values = c(16, 17, 15),
                       labels = c("A: base","B: +stage","C: +stage+hist")) +
    labs(
      title    = paste0("TP53 ME with ", gb, " — clinical covariate adjustment"),
      subtitle = "Model A (base) vs B (+stage) vs C (+stage+histology)",
      x        = "Odds Ratio (log scale, TP53 outcome)",
      y        = NULL,
      color    = "Model", shape = "Model",
      caption  = paste0(
        "OR < 1 = mutual exclusivity with TP53.\n",
        "Diamond = pooled random-effects estimate. Bars = 95% CI.\n",
        "Histology: 0 = intestinal/NOS, 1 = diffuse/signet-ring."
      )
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "right",
          panel.grid.minor = element_blank())

  fname <- file.path(OUTDIR, paste0("forest_", gb, "_clinical.png"))
  ggsave(fname, plot = p, width = 10, height = max(4, nrow(plot_data %>% distinct(label)) * 0.55 + 2),
         dpi = 150)
  message("Saved: ", fname)
}

# ── 7. Side-by-side summary plot (p-value trajectory A→B→C) ──────────────────

traj_data <- meta_rows %>%
  filter(status == "ok") %>%
  mutate(
    model_f = factor(model, c("A","B","C"), c("A: base","B: +stage","C: +stage+hist")),
    neg_log10_p = -log10(p_value),
    label = paste0("TP53 × ", gene_b)
  )

if (nrow(traj_data) > 0) {
  p_traj <- ggplot(traj_data, aes(x = model_f, y = neg_log10_p,
                                   group = gene_b, color = gene_b)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    geom_text(aes(label = paste0("OR=", round(pooled_OR, 2), "\np=", signif(p_value, 2))),
              vjust = -0.6, size = 3.2) +
    scale_color_manual(values = c(ATM = "#d62728", ARID1A = "#9467bd")) +
    labs(
      title    = "TP53 ME signal: effect of clinical covariate adjustment",
      subtitle = "Pooled (random-effects) across cohorts",
      x        = "Model",
      y        = expression(-log[10](p)),
      color    = "Gene B",
      caption  = "Dashed line = p = 0.05 threshold"
    ) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank())

  ggsave(file.path(OUTDIR, "pvalue_trajectory.png"),
         plot = p_traj, width = 7, height = 5, dpi = 150)
}

message("\n=== Done. All outputs in: ", OUTDIR, " ===")
