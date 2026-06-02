#!/usr/bin/env Rscript
# validate_tp53_me_gene_length_sensitivity.R
#
# Sensitivity analysis: does gene CDS length confound the TP53 ME signals?
#
# Three models compared per cohort × gene:
#   Model A (original):  TP53 ~ gene_B_binary   + log2(TMB)
#   Model B (rate):      TP53 ~ mut_rate_B_perMb + log2(TMB)
#   Model C (residual):  TP53 ~ length_residual_B + log2(TMB)
#
#   Model B replaces the binary 0/1 predictor with the number of functional
#   mutations divided by gene CDS length (Mb). A single mutation in TRRAP
#   (11.6 kb CDS) contributes ~1/10 the weight of a mutation in TP53 (1.2 kb).
#
#   Model C uses the Pearson residual of the binary indicator against its
#   Poisson-null expected probability:
#     p_null_i = 1 - exp( -TMB_i * CDS_Mb_G )
#     residual_i = (observed_01 - p_null) / sqrt(p_null * (1 - p_null))
#   This captures "mutation beyond what gene size + TMB alone predict."
#
# CDS lengths are from Ensembl GRCh38 MANE Select canonical transcripts.
# Only WXS cohorts are included (TRRAP not covered by IMPACT341).
#
# Run:
#   /home/sgao30/micromamba/envs/tcga_bioc/bin/Rscript \
#     scripts/Phase2/validate_tp53_me_gene_length_sensitivity.R
#
# Outputs: analysis_results/Phase2/length_norm/
#   comparison_all_models.tsv   — per cohort × gene × model
#   meta_comparison.tsv         — pooled RE estimates per gene per model
#   gene_length_table.tsv       — CDS lengths used
#   sensitivity_scatter.png     — original OR vs. corrected OR
#   sensitivity_barplot.png     — side-by-side model comparison per gene
#   trrap_sensitivity.png       — TRRAP per-cohort model comparison

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(tibble)
  library(purrr)
  library(forcats)
  library(metafor)
})

# ============================================================
# SECTION 0: Configuration
# ============================================================

TARGET_GENES <- c("SOX9", "POLD1", "TRRAP", "ARID2", "ATM", "CDH1", "ARID1A")

ENTREZ_IDS <- c(
  TP53   = 7157L,
  SOX9   = 6662L,
  ATM    =  472L,
  CDH1   =  999L,
  ARID1A = 8289L,
  ARID2  = 196528L,
  POLD1  = 5424L,
  TRRAP  = 8295L
)

FUNCTIONAL_CLASSES <- c(
  "Missense_Mutation", "Nonsense_Mutation",
  "Frame_Shift_Del",   "Frame_Shift_Ins",
  "Splice_Site",       "In_Frame_Del",
  "In_Frame_Ins",      "Nonstop_Mutation",
  "Translation_Start_Site"
)

EXOME_MB <- 38.0

# CDS lengths in Mb — Ensembl GRCh38 MANE Select canonical transcripts.
# These represent the exome-captured coding bases per gene (relevant for WXS TMB).
#
#   TP53   ENST00000269305.9  1,182 bp
#   SOX9   ENST00000245479.4  1,674 bp
#   CDH1   ENST00000261769.10 2,496 bp
#   POLD1  ENST00000320042.8  2,958 bp
#   ARID1A ENST00000374152.8  6,855 bp
#   ARID2  ENST00000420313.7  6,018 bp
#   ATM    ENST00000675843.1  9,171 bp
#   TRRAP  ENST00000371904.8 11,577 bp  ← 9.8× larger than TP53
#
GENE_CDS_MB <- c(
  TP53   = 1182  / 1e6,
  SOX9   = 1674  / 1e6,
  CDH1   = 2496  / 1e6,
  POLD1  = 2958  / 1e6,
  ARID1A = 6855  / 1e6,
  ARID2  = 6018  / 1e6,
  ATM    = 9171  / 1e6,
  TRRAP  = 11577 / 1e6
)

PATHS <- list(
  tcga_maf   = "/ShangGaoAIProjects/Gastric/genome/TCGA-STAD/somatic_mutation/somatic_mutation.csv",
  cbio_base  = "/ShangGaoAIProjects/Gastric/genome/Other_cBioPartal_stomach",
  output_dir = "/ShangGaoAIProjects/Gastric/genome/analysis_results/Phase2/length_norm"
)

# WXS cohorts only (IMPACT341 panels do not cover TRRAP)
WXS_COHORTS <- c(
  "TCGA-STAD",
  "egc_tmucih_2015",
  "stad_oncosg_2018",
  "stad_pfizer_uhongkong"
)

dir.create(PATHS$output_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== Phase 2 sensitivity: gene-length normalization ===\n")
cat("Target genes:", paste(TARGET_GENES, collapse = ", "), "\n")
cat("CDS lengths (kb):\n")
for (g in c("TP53", TARGET_GENES)) {
  cat(sprintf("  %-8s  %.2f kb\n", g, GENE_CDS_MB[g] * 1000))
}
cat("\n")

# ============================================================
# SECTION 1: Data Loading
# ============================================================
# Returns a list with:
#   $data       — tibble: sample_id, TP53, <genes>, TMB (binary 0/1 columns)
#   $count_data — tibble: sample_id, TP53_n, <gene>_n, TMB (raw mutation counts)
#   $cohort     — character
#   $n          — integer

build_matrices <- function(sample_ids, mut_func_df, gene_col, tmb_df) {
  all_genes <- c("TP53", TARGET_GENES)

  # Count matrix: number of functional mutations per sample per gene
  counts <- mut_func_df %>%
    filter(!!sym(gene_col) %in% all_genes) %>%
    count(sample_id, !!sym(gene_col), name = "n_mut") %>%
    pivot_wider(names_from = !!sym(gene_col), values_from = n_mut,
                values_fill = 0L, names_prefix = "n_")

  # Binary matrix
  binary <- mut_func_df %>%
    filter(!!sym(gene_col) %in% all_genes) %>%
    distinct(sample_id, !!sym(gene_col)) %>%
    mutate(val = 1L) %>%
    pivot_wider(names_from = !!sym(gene_col), values_from = val, values_fill = 0L)

  base <- tibble(sample_id = sample_ids)

  mat_bin <- base %>%
    left_join(binary, by = "sample_id") %>%
    mutate(across(any_of(all_genes), ~replace_na(., 0L)))

  mat_cnt <- base %>%
    left_join(counts, by = "sample_id") %>%
    mutate(across(starts_with("n_"), ~replace_na(., 0L)))

  # Attach TMB
  tmb_clean <- tmb_df %>%
    rename(sample_id = 1, TMB = 2) %>%
    mutate(TMB = pmax(as.numeric(TMB), 0.1))

  mat_bin <- mat_bin %>% left_join(tmb_clean, by = "sample_id") %>% filter(!is.na(TMB))
  mat_cnt <- mat_cnt %>% left_join(tmb_clean, by = "sample_id") %>% filter(!is.na(TMB))

  list(binary = mat_bin, counts = mat_cnt)
}

load_tcga <- function(maf_path) {
  cat("Loading TCGA-STAD ...\n")
  maf <- read_csv(maf_path, show_col_types = FALSE) %>%
    mutate(sample_id = Tumor_Sample_Barcode) %>%
    filter(str_detect(sample_id, "-01[AB]-"))

  tmb <- maf %>%
    filter(Variant_Classification %in% FUNCTIONAL_CLASSES) %>%
    count(sample_id, name = "n_mut") %>%
    mutate(TMB = n_mut / EXOME_MB)

  all_samples <- unique(maf$sample_id)
  tmb_full <- tibble(sample_id = all_samples) %>%
    left_join(tmb, by = "sample_id") %>%
    mutate(TMB = replace_na(TMB, 0.1)) %>%
    select(sample_id, TMB)

  mut_func <- maf %>%
    filter(Variant_Classification %in% FUNCTIONAL_CLASSES,
           Hugo_Symbol %in% c("TP53", TARGET_GENES)) %>%
    select(sample_id, gene = Hugo_Symbol)

  mats <- build_matrices(all_samples, mut_func, "gene", tmb_full)
  cat("  →", nrow(mats$binary), "samples\n")
  list(binary = mats$binary, counts = mats$counts, cohort = "TCGA-STAD", n = nrow(mats$binary))
}

load_cbioportal_wxs <- function(cohort_name, base_path) {
  cat("Loading", cohort_name, "...\n")
  mut  <- read_tsv(file.path(base_path, cohort_name, "mutations.tsv"),
                   show_col_types = FALSE, comment = "#")
  clin <- read_tsv(file.path(base_path, cohort_name, "clinical_sample.tsv"),
                   show_col_types = FALSE, comment = "#")

  # Deduplicate: one sample per patient
  if ("patientId" %in% names(mut)) {
    sample_meta <- mut %>%
      distinct(sampleId, patientId) %>%
      left_join(clin %>% select(sampleId, any_of("SAMPLE_TYPE")), by = "sampleId")

    if ("SAMPLE_TYPE" %in% names(sample_meta)) {
      keep <- sample_meta %>%
        mutate(pref = if_else(SAMPLE_TYPE == "Primary", 1L, 2L)) %>%
        group_by(patientId) %>%
        slice_min(pref, n = 1, with_ties = FALSE) %>%
        ungroup() %>% pull(sampleId)
    } else {
      keep <- sample_meta %>%
        group_by(patientId) %>% slice(1) %>% ungroup() %>% pull(sampleId)
    }
    mut  <- mut  %>% filter(sampleId %in% keep)
    clin <- clin %>% filter(sampleId %in% keep)
    cat("  Dedup: kept", length(keep), "samples\n")
  }

  tmb <- if ("TMB_NONSYNONYMOUS" %in% names(clin)) {
    clin %>% select(sample_id = sampleId, TMB = TMB_NONSYNONYMOUS) %>%
      mutate(TMB = pmax(as.numeric(TMB), 0.1))
  } else {
    tmb_c <- mut %>%
      filter(mutationType %in% FUNCTIONAL_CLASSES) %>%
      count(sampleId, name = "n_mut") %>%
      mutate(TMB = n_mut / EXOME_MB)
    tibble(sample_id = unique(clin$sampleId)) %>%
      left_join(tmb_c %>% rename(sample_id = sampleId), by = "sample_id") %>%
      mutate(TMB = replace_na(TMB, 0.1)) %>%
      select(sample_id, TMB)
  }

  gene_map <- tibble(entrezGeneId = ENTREZ_IDS, gene = names(ENTREZ_IDS))
  mut_func <- mut %>%
    filter(mutationType %in% FUNCTIONAL_CLASSES, entrezGeneId %in% ENTREZ_IDS) %>%
    left_join(gene_map, by = "entrezGeneId") %>%
    select(sample_id = sampleId, gene)

  all_samples <- unique(clin$sampleId)
  mats <- build_matrices(all_samples, mut_func, "gene", tmb)
  cat("  →", nrow(mats$binary), "samples\n")
  list(binary = mats$binary, counts = mats$counts, cohort = cohort_name, n = nrow(mats$binary))
}

# ============================================================
# SECTION 2: Gene Length Normalization Functions
# ============================================================

# Model A — original binary predictor
run_model_A <- function(df, gene_b, cohort_name) {
  if (!gene_b %in% names(df)) {
    return(null_result(cohort_name, gene_b, "gene_not_covered", "A"))
  }
  model_df <- df %>%
    select(TP53, gene = all_of(gene_b), TMB) %>%
    filter(!is.na(TP53), !is.na(gene), !is.na(TMB)) %>%
    mutate(log2_TMB = log2(TMB))

  n_gene <- sum(model_df$gene)
  n_tp53 <- sum(model_df$TP53)
  if (n_gene < 5 || n_tp53 < 5) {
    return(null_result(cohort_name, gene_b, "too_few_mutations", "A",
                       n = nrow(model_df), n_tp53 = n_tp53, n_gene = n_gene))
  }
  fit_glm(model_df, "gene", cohort_name, gene_b, "A",
          n_tp53 = n_tp53, n_gene = n_gene, n_both = sum(model_df$TP53 == 1 & model_df$gene == 1))
}

# Model B — mutation rate predictor (n_mutations / CDS_Mb)
# OR interpretation: per 1 mutation/Mb increase in gene B rate, odds of TP53 mutation changes by OR.
# For genes where most samples have 0 or 1 mutation:
#   0 mutations → rate = 0
#   1 mutation  → rate = 1 / CDS_Mb_G  (larger genes give smaller rate)
run_model_B <- function(count_df, binary_df, gene_b, cohort_name) {
  count_col <- paste0("n_", gene_b)
  if (!count_col %in% names(count_df)) {
    return(null_result(cohort_name, gene_b, "gene_not_covered", "B"))
  }
  cds_mb <- GENE_CDS_MB[gene_b]

  model_df <- count_df %>%
    select(sample_id, n_gene = all_of(count_col), TMB) %>%
    left_join(binary_df %>% select(sample_id, TP53), by = "sample_id") %>%
    filter(!is.na(TP53), !is.na(n_gene), !is.na(TMB)) %>%
    mutate(
      log2_TMB  = log2(TMB),
      mut_rate  = n_gene / cds_mb      # mutations per Mb
    )

  n_gene <- sum(model_df$n_gene > 0)
  n_tp53 <- sum(model_df$TP53)
  if (n_gene < 5 || n_tp53 < 5) {
    return(null_result(cohort_name, gene_b, "too_few_mutations", "B",
                       n = nrow(model_df), n_tp53 = n_tp53, n_gene = n_gene))
  }
  # Note: predictor is mut_rate (continuous). OR = per 1 mutation/Mb increase.
  fit_glm(model_df, "mut_rate", cohort_name, gene_b, "B",
          n_tp53 = n_tp53, n_gene = n_gene,
          n_both = sum(model_df$TP53 == 1 & model_df$n_gene > 0))
}

# Model C — Pearson residual predictor
# residual_i = (obs_binary_i - p_null_i) / sqrt(p_null_i * (1 - p_null_i))
# p_null_i   = 1 - exp(-TMB_i * CDS_Mb_G)  [Poisson probability of ≥1 mutation]
# Positive residual = mutated more than expected; negative = less than expected.
# ME signal after length correction: ME-true genes should have negative residuals
# in TP53-mutant samples even after removing the TMB × gene_size expectation.
run_model_C <- function(df, gene_b, cohort_name) {
  if (!gene_b %in% names(df)) {
    return(null_result(cohort_name, gene_b, "gene_not_covered", "C"))
  }
  cds_mb <- GENE_CDS_MB[gene_b]

  model_df <- df %>%
    select(TP53, obs = all_of(gene_b), TMB) %>%
    filter(!is.na(TP53), !is.na(obs), !is.na(TMB)) %>%
    mutate(
      log2_TMB = log2(TMB),
      p_null   = 1 - exp(-TMB * cds_mb),
      p_null   = pmin(pmax(p_null, 1e-6), 1 - 1e-6),   # numerical stability
      residual = (obs - p_null) / sqrt(p_null * (1 - p_null))
    )

  n_gene <- sum(model_df$obs)
  n_tp53 <- sum(model_df$TP53)
  if (n_gene < 5 || n_tp53 < 5) {
    return(null_result(cohort_name, gene_b, "too_few_mutations", "C",
                       n = nrow(model_df), n_tp53 = n_tp53, n_gene = n_gene))
  }
  fit_glm(model_df, "residual", cohort_name, gene_b, "C",
          n_tp53 = n_tp53, n_gene = n_gene,
          n_both = sum(model_df$TP53 == 1 & model_df$obs == 1))
}

# Helper: run glm and return tidy tibble
fit_glm <- function(model_df, pred_col, cohort_name, gene_b, model_label,
                    n_tp53 = NA, n_gene = NA, n_both = NA) {
  tryCatch({
    fit <- glm(
      as.formula(sprintf("TP53 ~ %s + log2_TMB", pred_col)),
      data = model_df, family = binomial()
    )
    cf   <- summary(fit)$coefficients
    beta <- cf[pred_col, "Estimate"]
    se   <- cf[pred_col, "Std. Error"]
    pval <- cf[pred_col, "Pr(>|z|)"]
    tibble(
      cohort   = cohort_name, gene_b = gene_b, model = model_label,
      n        = nrow(model_df), n_tp53 = n_tp53, n_gene = n_gene, n_both = n_both,
      OR_adj   = exp(beta), CI_lower = exp(beta - 1.96 * se), CI_upper = exp(beta + 1.96 * se),
      beta     = beta, se = se, p_value = pval, status = "ok"
    )
  }, error = function(e) {
    null_result(cohort_name, gene_b, paste0("glm_error: ", conditionMessage(e)), model_label,
                n = nrow(model_df), n_tp53 = n_tp53, n_gene = n_gene)
  })
}

null_result <- function(cohort_name, gene_b, status_str, model_label,
                        n = NA, n_tp53 = NA, n_gene = NA) {
  tibble(
    cohort   = cohort_name, gene_b = gene_b, model = model_label,
    n        = n, n_tp53 = n_tp53, n_gene = n_gene, n_both = NA_integer_,
    OR_adj   = NA_real_, CI_lower = NA_real_, CI_upper = NA_real_,
    beta     = NA_real_, se = NA_real_, p_value = NA_real_,
    status   = status_str
  )
}

# ============================================================
# SECTION 3: Load WXS Cohorts
# ============================================================

cohorts <- list()
cohorts[["TCGA-STAD"]]             <- load_tcga(PATHS$tcga_maf)
cohorts[["egc_tmucih_2015"]]       <- load_cbioportal_wxs("egc_tmucih_2015",       PATHS$cbio_base)
cohorts[["stad_oncosg_2018"]]      <- load_cbioportal_wxs("stad_oncosg_2018",      PATHS$cbio_base)
cohorts[["stad_pfizer_uhongkong"]] <- load_cbioportal_wxs("stad_pfizer_uhongkong", PATHS$cbio_base)

cat("\nCohorts loaded:", paste(names(cohorts), collapse = ", "), "\n")
cat("Note: IMPACT341 cohorts (egc_msk_*) excluded — TRRAP not on panel.\n\n")

# ============================================================
# SECTION 4: Run All Three Models
# ============================================================

cat("--- Running models A / B / C ---\n")
results_list <- list()

for (cn in names(cohorts)) {
  co <- cohorts[[cn]]
  for (g in TARGET_GENES) {
    key <- paste(cn, g, sep = "|")
    results_list[[paste0(key, "_A")]] <- run_model_A(co$binary, g, cn)
    results_list[[paste0(key, "_B")]] <- run_model_B(co$counts, co$binary, g, cn)
    results_list[[paste0(key, "_C")]] <- run_model_C(co$binary, g, cn)
  }
}

results <- bind_rows(results_list)

# Model labels for display
MODEL_LABELS <- c(
  A = "A: binary (original)",
  B = "B: mut rate / CDS Mb",
  C = "C: Poisson residual"
)
results <- results %>%
  mutate(model_label = MODEL_LABELS[model])

write_tsv(results, file.path(PATHS$output_dir, "comparison_all_models.tsv"))
cat("Saved: comparison_all_models.tsv\n")

# Quick console summary
cat("\nStatus breakdown:\n")
print(table(results$model, results$status))

# ============================================================
# SECTION 5: Meta-Analysis per Gene per Model
# ============================================================

cat("\n--- Meta-analysis per gene per model ---\n")
meta_list <- list()

for (g in TARGET_GENES) {
  for (m in c("A", "B", "C")) {
    gene_data <- results %>%
      filter(gene_b == g, model == m, status == "ok",
             !is.na(beta), !is.na(se), se > 0, is.finite(beta))

    if (nrow(gene_data) < 2) {
      cat(sprintf("  Skip %s model %s: only %d cohort(s)\n", g, m, nrow(gene_data)))
      next
    }

    tryCatch({
      ma <- rma(yi = beta, sei = se, data = gene_data, method = "REML")
      meta_list[[paste(g, m, sep = "|")]] <- tibble(
        gene_b     = g, model = m, model_label = MODEL_LABELS[m],
        n_cohorts  = nrow(gene_data),
        pooled_OR  = exp(coef(ma)),
        CI_lower   = exp(ma$ci.lb), CI_upper = exp(ma$ci.ub),
        p_value    = ma$pval, I2 = round(ma$I2, 1),
        beta_pool  = coef(ma), se_pool = ma$se
      )
      cat(sprintf("  %-8s model %s  OR=%.3f [%.3f–%.3f]  p=%s  I²=%.0f%%\n",
                  g, m, exp(coef(ma)), exp(ma$ci.lb), exp(ma$ci.ub),
                  format(ma$pval, digits = 2, scientific = TRUE), ma$I2))
    }, error = function(e) {
      cat("  Error:", g, m, conditionMessage(e), "\n")
    })
  }
}

meta_df <- bind_rows(meta_list)
write_tsv(meta_df %>% select(-beta_pool, -se_pool),
          file.path(PATHS$output_dir, "meta_comparison.tsv"))
cat("Saved: meta_comparison.tsv\n")

# ============================================================
# SECTION 6: Gene Length Table
# ============================================================

gene_len_tbl <- tibble(
  gene        = names(GENE_CDS_MB),
  CDS_bp      = as.integer(GENE_CDS_MB * 1e6),
  CDS_kb      = round(GENE_CDS_MB * 1000, 2),
  CDS_Mb      = GENE_CDS_MB,
  ratio_vs_TP53 = round(GENE_CDS_MB / GENE_CDS_MB["TP53"], 1),
  note        = case_when(
    gene == "TRRAP" ~ "largest tested; 9.8x TP53",
    gene == "ATM"   ~ "large; 7.8x TP53",
    gene == "TP53"  ~ "reference",
    TRUE ~ ""
  )
)
write_tsv(gene_len_tbl, file.path(PATHS$output_dir, "gene_length_table.tsv"))
cat("\nGene CDS lengths:\n")
print(gene_len_tbl %>% select(gene, CDS_kb, ratio_vs_TP53, note) %>% as.data.frame())

# ============================================================
# SECTION 7: Plots
# ============================================================

cat("\n--- Generating plots ---\n")

# --- 7a. Scatter: original OR (Model A) vs length-corrected OR (Model C, pooled) ---

scatter_df <- meta_df %>%
  filter(model %in% c("A", "C")) %>%
  select(gene_b, model, pooled_OR, CI_lower, CI_upper, p_value) %>%
  pivot_wider(
    names_from  = model,
    values_from = c(pooled_OR, CI_lower, CI_upper, p_value)
  ) %>%
  filter(!is.na(pooled_OR_A), !is.na(pooled_OR_C)) %>%
  mutate(
    sig_A  = p_value_A  < 0.05,
    sig_C  = p_value_C  < 0.05,
    change = case_when(
      sig_A & sig_C  ~ "significant in both",
      sig_A & !sig_C ~ "lost after correction",
      !sig_A & sig_C ~ "gained after correction",
      TRUE           ~ "not significant"
    ),
    cds_kb = GENE_CDS_MB[gene_b] * 1000
  )

# Diagonal reference lines for ±20% OR change
p_scatter <- ggplot(scatter_df, aes(x = log2(pooled_OR_A), y = log2(pooled_OR_C))) +
  geom_abline(slope = 1, intercept = 0, color = "grey60", linetype = "dashed") +
  geom_abline(slope = 1, intercept =  log2(1.2), color = "grey80", linetype = "dotted") +
  geom_abline(slope = 1, intercept = -log2(1.2), color = "grey80", linetype = "dotted") +
  geom_point(aes(color = change, size = cds_kb), alpha = 0.85) +
  geom_errorbar(
    aes(ymin = log2(CI_lower_C), ymax = log2(CI_upper_C), color = change),
    width = 0.03, alpha = 0.5
  ) +
  geom_errorbarh(
    aes(xmin = log2(CI_lower_A), xmax = log2(CI_upper_A), color = change),
    height = 0.03, alpha = 0.5
  ) +
  geom_text(aes(label = gene_b), vjust = -0.7, size = 3.5, fontface = "bold") +
  scale_color_manual(
    values = c(
      "significant in both"   = "#2166ac",
      "lost after correction" = "#d73027",
      "not significant"       = "#969696"
    ),
    name = "Significance"
  ) +
  scale_size_continuous(name = "Gene CDS (kb)", range = c(2.5, 7)) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.3) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
  labs(
    title    = "Gene length correction: original vs. Poisson-residual model",
    subtitle = "Each point = one gene; size ∝ CDS length; dashed line = identity; dotted = ±20% OR",
    x        = "log₂(OR_adj)  Model A — binary predictor",
    y        = "log₂(OR_adj)  Model C — Poisson residual predictor",
    caption  = "Pooled random-effects OR across WXS cohorts (TCGA, TMUCIH, OncoSG, HK Pfizer)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "right",
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 9, color = "grey40"),
    plot.caption     = element_text(size = 8, color = "grey50")
  )

ggsave(file.path(PATHS$output_dir, "sensitivity_scatter.png"),
       plot = p_scatter, width = 8, height = 6, dpi = 200)
cat("  Saved: sensitivity_scatter.png\n")

# --- 7b. Side-by-side bar: pooled OR per gene per model ---

bar_df <- meta_df %>%
  filter(!is.na(pooled_OR)) %>%
  mutate(
    gene_f  = factor(gene_b, levels = TARGET_GENES),
    log2_OR = log2(pooled_OR),
    sig_lab = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ "ns"
    ),
    model_f = factor(model, levels = c("A", "B", "C"),
                     labels = c("A: binary", "B: rate/Mb", "C: Poisson resid"))
  )

p_bar <- ggplot(bar_df, aes(x = gene_f, y = log2_OR, fill = model_f)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7, alpha = 0.85) +
  geom_errorbar(
    aes(ymin = log2(CI_lower), ymax = log2(CI_upper)),
    position = position_dodge(width = 0.75), width = 0.25, linewidth = 0.5
  ) +
  geom_text(
    aes(label = sig_lab,
        y     = log2(CI_upper) + 0.1),
    position = position_dodge(width = 0.75),
    size = 3.5, vjust = 0
  ) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "black") +
  scale_fill_manual(
    values = c(
      "A: binary"       = "#4575b4",
      "B: rate/Mb"      = "#74add1",
      "C: Poisson resid"= "#313695"
    ),
    name = "Model"
  ) +
  labs(
    title    = "TP53 ME signal: pooled OR across three gene-length models",
    subtitle = "Bars below 0 = mutually exclusive with TP53  |  * p<0.05  ** p<0.01  *** p<0.001\nWXS cohorts only (TCGA + TMUCIH + OncoSG + HK Pfizer)",
    x        = "Gene",
    y        = "log₂(pooled OR_adj)",
    caption  = "Model A = original binary; B = mutation rate per CDS Mb; C = Poisson residual"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title    = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    plot.caption  = element_text(size = 8, color = "grey50"),
    axis.text.x   = element_text(face = "bold", size = 10)
  )

ggsave(file.path(PATHS$output_dir, "sensitivity_barplot.png"),
       plot = p_bar, width = 10, height = 5.5, dpi = 200)
cat("  Saved: sensitivity_barplot.png\n")

# --- 7c. TRRAP per-cohort comparison: 3 models side by side ---

trrap_df <- results %>%
  filter(gene_b == "TRRAP", status == "ok") %>%
  mutate(
    cohort_f = factor(cohort, levels = rev(WXS_COHORTS)),
    model_f  = factor(model, levels = c("A", "B", "C"),
                      labels = c("A: binary", "B: rate/Mb", "C: Poisson resid"))
  )

# Add pooled row
trrap_pooled <- meta_df %>%
  filter(gene_b == "TRRAP") %>%
  transmute(
    cohort   = "Pooled (RE)",
    gene_b   = "TRRAP",
    model    = model,
    OR_adj   = pooled_OR, CI_lower = CI_lower, CI_upper = CI_upper, p_value = p_value,
    n = NA, n_tp53 = NA, n_gene = NA, n_both = NA,
    beta = beta_pool, se = se_pool, status = "pooled",
    model_label = model_label,
    cohort_f = factor("Pooled (RE)", levels = c("Pooled (RE)", rev(WXS_COHORTS))),
    model_f  = factor(model, levels = c("A", "B", "C"),
                      labels = c("A: binary", "B: rate/Mb", "C: Poisson resid"))
  )

trrap_plot_df <- bind_rows(
  trrap_df %>% mutate(cohort_f = factor(as.character(cohort_f),
                                        levels = c("Pooled (RE)", WXS_COHORTS))),
  trrap_pooled
) %>%
  mutate(
    CI_lower_d = pmax(CI_lower, 0.02),
    CI_upper_d = pmin(CI_upper, 50),
    sig_lab    = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ ""
    )
  )

cds_ratio <- round(GENE_CDS_MB["TRRAP"] / GENE_CDS_MB["TP53"], 1)

p_trrap <- ggplot(trrap_plot_df,
                  aes(x = OR_adj, y = model_f, color = model_f, shape = model_f)) +
  facet_wrap(~cohort_f, ncol = 1, strip.position = "left") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_errorbar(
    aes(xmin = CI_lower_d, xmax = CI_upper_d),
    width = 0.25, linewidth = 0.6, orientation = "y"
  ) +
  geom_point(aes(size = ifelse(status == "pooled", 4, 2.5))) +
  geom_text(
    data = trrap_plot_df %>% filter(nchar(sig_lab) > 0),
    aes(x = CI_upper_d * 1.08, label = sig_lab),
    hjust = 0, size = 4, show.legend = FALSE
  ) +
  scale_x_log10(
    breaks = c(0.05, 0.1, 0.2, 0.5, 1, 2, 5),
    labels = c("0.05", "0.1", "0.2", "0.5", "1", "2", "5")
  ) +
  scale_size_identity() +
  scale_color_manual(
    values = c(
      "A: binary"        = "#4575b4",
      "B: rate/Mb"       = "#74add1",
      "C: Poisson resid" = "#313695"
    ), name = "Model"
  ) +
  scale_shape_manual(
    values = c("A: binary" = 16, "B: rate/Mb" = 17, "C: Poisson resid" = 15),
    name = "Model"
  ) +
  labs(
    title    = "TRRAP ↔ TP53: gene-length sensitivity analysis",
    subtitle = sprintf(
      "TRRAP CDS = %.1f kb  (%.1fx larger than TP53 %.1f kb)\nRows = cohorts; each cohort shows 3 models side-by-side",
      GENE_CDS_MB["TRRAP"] * 1000, cds_ratio, GENE_CDS_MB["TP53"] * 1000
    ),
    x       = "Adjusted OR (log scale)   ← ME  |  Co-occur →",
    y       = NULL,
    caption = "WXS cohorts only; * p<0.05  ** p<0.01  *** p<0.001"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position    = "right",
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    strip.text         = element_text(size = 9, face = "bold"),
    strip.placement    = "outside",
    plot.title         = element_text(face = "bold", size = 12),
    plot.subtitle      = element_text(size = 9, color = "grey40"),
    plot.caption       = element_text(size = 8, color = "grey50")
  )

ggsave(file.path(PATHS$output_dir, "trrap_sensitivity.png"),
       plot = p_trrap, width = 9, height = 5.5, dpi = 200)
cat("  Saved: trrap_sensitivity.png\n")

# ============================================================
# SECTION 8: Final Summary
# ============================================================

cat("\n=== Summary: OR shift after gene-length correction ===\n")

summary_tbl <- meta_df %>%
  filter(model %in% c("A", "C")) %>%
  select(gene_b, model, pooled_OR, p_value, I2) %>%
  pivot_wider(
    names_from  = model,
    values_from = c(pooled_OR, p_value, I2)
  ) %>%
  mutate(
    OR_ratio   = round(pooled_OR_C / pooled_OR_A, 3),
    sig_A      = !is.na(p_value_A) & p_value_A < 0.05,
    sig_C      = !is.na(p_value_C) & p_value_C < 0.05,
    conclusion = case_when(
      sig_A & sig_C & OR_ratio > 0.8 & OR_ratio < 1.2 ~ "robust: <20% OR shift",
      sig_A & sig_C                                    ~ "robust: signal survives, OR shifted",
      sig_A & !sig_C                                   ~ "CAUTION: lost after correction",
      TRUE                                             ~ "not significant in either"
    ),
    CDS_kb     = round(GENE_CDS_MB[gene_b] * 1000, 1)
  ) %>%
  arrange(p_value_A) %>%
  select(gene_b, CDS_kb, pooled_OR_A, pooled_OR_C, OR_ratio, p_value_A, p_value_C,
         sig_A, sig_C, conclusion)

cat("\n")
print(as.data.frame(summary_tbl), digits = 3)

write_tsv(summary_tbl, file.path(PATHS$output_dir, "sensitivity_summary.tsv"))

cat("\n=== Complete ===\n")
cat("Output:", PATHS$output_dir, "\n")
cat("  comparison_all_models.tsv\n")
cat("  meta_comparison.tsv\n")
cat("  sensitivity_summary.tsv\n")
cat("  gene_length_table.tsv\n")
cat("  sensitivity_scatter.png\n")
cat("  sensitivity_barplot.png\n")
cat("  trrap_sensitivity.png\n")
cat("\nKey question: if TRRAP's OR remains <1 and significant in Model C,\n")
cat("  the ME signal is not explained by gene length alone.\n")
