#!/usr/bin/env Rscript
# validate_tp53_me_cross_cohort.R
#
# Cross-cohort validation of TP53 mutual exclusivity signals from FDU cohort.
# Runs TMB-adjusted logistic regression in 8 public gastric/EGC datasets,
# then pools estimates via random-effects meta-analysis.
#
# Run:
#   /home/sgao30/micromamba/envs/tcga_bioc/bin/Rscript scripts/Phase2/validate_tp53_me_cross_cohort.R
#
# Inputs:
#   TCGA-STAD/somatic_mutation/somatic_mutation.csv
#   Other_cBioPartal_stomach/{cohort}/mutations.tsv
#   Other_cBioPartal_stomach/{cohort}/clinical_sample.tsv
#
# Outputs (analysis_results/Phase2/):
#   per_cohort_results.tsv       — all models (cohort × gene)
#   meta_analysis_results.tsv    — pooled RE estimates per gene
#   panel_coverage_check.tsv     — gene coverage per cohort
#   forest_{GENE}.png            — forest plot per gene
#   summary_heatmap.png          — OR_adj heatmap (genes × cohorts)

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

# Entrez Gene IDs (used to look up genes in cBioPortal mutations.tsv)
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

# WXS exome size for TCGA TMB computation
EXOME_MB <- 38.0

PATHS <- list(
  tcga_maf      = "/ShangGaoAIProjects/Gastric/genome/TCGA-STAD/somatic_mutation/somatic_mutation.csv",
  cbio_base     = "/ShangGaoAIProjects/Gastric/genome/Other_cBioPartal_stomach",
  output_dir    = "/ShangGaoAIProjects/Gastric/genome/analysis_results/Phase2"
)

dir.create(PATHS$output_dir, recursive = TRUE, showWarnings = FALSE)

# FDU discovery results (bidirectionally confirmed, from interactions_broad_logistic.tsv)
FDU_RESULTS <- tribble(
  ~gene_b,  ~OR_adj, ~beta,    ~se,
  "SOX9",   0.0849,  log(0.0849),  0.6509,
  "POLD1",  0.1722,  log(0.1722),  0.5708,
  "TRRAP",  0.1842,  log(0.1842),  0.5102,
  "ARID2",  0.1855,  log(0.1855),  0.5601,
  "ATM",    0.2611,  log(0.2611),  0.4259,
  "CDH1",   0.3861,  log(0.3861),  0.3358,
  "ARID1A", 0.4029,  log(0.4029),  0.2880
) %>% mutate(
  cohort    = "FDU (discovery)",
  n         = 402L,
  CI_lower  = exp(beta - 1.96 * se),
  CI_upper  = exp(beta + 1.96 * se),
  p_value   = c(1.51e-4, 2.06e-3, 9.14e-4, 2.64e-3, 1.61e-3, 4.60e-3, 1.60e-3),
  n_tp53    = NA_integer_, n_gene = NA_integer_, n_both = NA_integer_,
  status    = "discovery"
)

# Display labels and ordering
COHORT_ORDER <- c(
  "FDU (discovery)",
  "TCGA-STAD",
  "egc_tmucih_2015",
  "stad_oncosg_2018",
  "stad_pfizer_uhongkong",
  "egc_msk_2017",
  "egc_mskcc_2020",
  "egc_msk_tp53_ccr_2022",
  "egc_msk_2023",
  "Pooled (RE)"
)

COHORT_LABELS <- c(
  "FDU (discovery)"        = "FDU  N=402  [WXS, discovery]",
  "TCGA-STAD"              = "TCGA-STAD  N=431  [WXS]",
  "egc_tmucih_2015"        = "TMUCIH 2015  N=78  [WXS]",
  "stad_oncosg_2018"       = "OncoSG 2018  N=147  [WXS]",
  "stad_pfizer_uhongkong"  = "HK Pfizer  N=100  [WXS]",
  "egc_msk_2017"           = "MSK 2017  [IMPACT341, STAD only]",
  "egc_mskcc_2020"         = "MSK 2020  [IMPACT, STAD only]",
  "egc_msk_tp53_ccr_2022"  = "MSK TP53 2022  [IMPACT, STAD only]",
  "egc_msk_2023"           = "MSK 2023  [IMPACT341, STAD only]",
  "Pooled (RE)"            = "Pooled  (random effects)"
)

cat("=== Phase 2: Cross-cohort TP53 ME validation ===\n")
cat("Target genes:", paste(TARGET_GENES, collapse = ", "), "\n\n")

# ============================================================
# SECTION 1: Data Loading Functions
# ============================================================

# Build a list(data=tibble, cohort=str, n=int) with columns:
#   sample_id, TP53, <gene>, ..., TMB
# TMB is the raw value (log2 transform applied at analysis time)

build_binary_matrix <- function(sample_ids, mut_func_df, gene_col,
                                tmb_df, tmb_col = "TMB") {
  all_genes <- c("TP53", TARGET_GENES)

  # Filter on VALUES in gene_col (not column names of the data frame)
  binary <- mut_func_df %>%
    filter(!!sym(gene_col) %in% all_genes) %>%
    distinct(sample_id, !!sym(gene_col)) %>%
    mutate(val = 1L) %>%
    pivot_wider(names_from = !!sym(gene_col), values_from = val, values_fill = 0L)

  # Ensure all expected samples present, fill absent genes with 0
  base <- tibble(sample_id = sample_ids)
  mat <- base %>%
    left_join(binary, by = "sample_id") %>%
    mutate(across(any_of(all_genes), ~replace_na(., 0L)))

  # Attach TMB
  mat <- mat %>%
    left_join(tmb_df %>% rename(sample_id = 1, TMB = all_of(tmb_col)),
              by = "sample_id") %>%
    filter(!is.na(TMB), TMB > 0 | TRUE) %>%   # keep TMB=0 samples; floor applied later
    mutate(TMB = pmax(as.numeric(TMB), 0.1))

  mat
}

# ---- TCGA ----
load_tcga <- function(maf_path) {
  cat("Loading TCGA-STAD ...\n")
  maf <- read_csv(maf_path, show_col_types = FALSE)

  # Keep primary tumor samples only (barcode component 4 = 01A or 01B)
  maf <- maf %>%
    mutate(sample_id = Tumor_Sample_Barcode) %>%
    filter(str_detect(sample_id, "-01[AB]-"))

  # TMB: functional mutations per sample / exome size
  tmb <- maf %>%
    filter(Variant_Classification %in% FUNCTIONAL_CLASSES) %>%
    count(sample_id, name = "n_mut") %>%
    mutate(TMB = n_mut / EXOME_MB)

  # Fill TMB=0 for samples with no functional mutations
  all_samples <- unique(maf$sample_id)
  tmb <- tibble(sample_id = all_samples) %>%
    left_join(tmb, by = "sample_id") %>%
    mutate(TMB = replace_na(TMB, 0.0))

  # Functional variants, target genes only
  mut_func <- maf %>%
    filter(Variant_Classification %in% FUNCTIONAL_CLASSES,
           Hugo_Symbol %in% c("TP53", TARGET_GENES)) %>%
    select(sample_id, gene = Hugo_Symbol)

  mat <- build_binary_matrix(all_samples, mut_func, "gene", tmb, "TMB")

  cat("  →", nrow(mat), "samples after deduplication\n")
  list(data = mat, cohort = "TCGA-STAD", n = nrow(mat))
}

# ---- cBioPortal cohort ----
load_cbioportal <- function(cohort_name, base_path, stad_only = FALSE) {
  cat("Loading", cohort_name, "...\n")

  mut  <- read_tsv(file.path(base_path, cohort_name, "mutations.tsv"),
                   show_col_types = FALSE, comment = "#")
  clin <- read_tsv(file.path(base_path, cohort_name, "clinical_sample.tsv"),
                   show_col_types = FALSE, comment = "#")

  # Filter to STAD only (MSK EGC datasets)
  if (stad_only && "CANCER_TYPE_DETAILED" %in% names(clin)) {
    stad_ids <- clin %>%
      filter(str_detect(CANCER_TYPE_DETAILED,
                        regex("Stomach Adenocarcinoma", ignore_case = TRUE))) %>%
      pull(sampleId)
    clin <- clin %>% filter(sampleId %in% stad_ids)
    mut  <- mut  %>% filter(sampleId %in% stad_ids)
    cat("  STAD filter: kept", nrow(clin), "samples\n")
  }

  # Deduplicate: one sample per patient
  # MSK sample IDs: P-XXXXXXX-TYYY-IMZ → patient = first two components P-XXXXXXX
  if ("patientId" %in% names(mut)) {
    # Build patient → sample table, prefer Primary SAMPLE_TYPE
    sample_meta <- mut %>%
      distinct(sampleId, patientId) %>%
      left_join(clin %>% select(sampleId, any_of("SAMPLE_TYPE")), by = "sampleId")

    if ("SAMPLE_TYPE" %in% names(sample_meta)) {
      # Per patient: prefer Primary; if none, take any
      keep_samples <- sample_meta %>%
        mutate(pref = if_else(SAMPLE_TYPE == "Primary", 1L, 2L)) %>%
        group_by(patientId) %>%
        slice_min(pref, n = 1, with_ties = FALSE) %>%
        ungroup() %>%
        pull(sampleId)
    } else {
      keep_samples <- sample_meta %>%
        group_by(patientId) %>%
        slice(1) %>%
        ungroup() %>%
        pull(sampleId)
    }

    mut  <- mut  %>% filter(sampleId %in% keep_samples)
    clin <- clin %>% filter(sampleId %in% keep_samples)
    cat("  Dedup: kept", length(keep_samples), "samples (1 per patient)\n")
  }

  # TMB: use TMB_NONSYNONYMOUS from clinical_sample if available
  if ("TMB_NONSYNONYMOUS" %in% names(clin)) {
    tmb <- clin %>%
      select(sample_id = sampleId, TMB = TMB_NONSYNONYMOUS) %>%
      mutate(TMB = as.numeric(TMB))
  } else {
    # Fallback: compute from functional mutation count
    tmb_counts <- mut %>%
      filter(mutationType %in% FUNCTIONAL_CLASSES) %>%
      count(sampleId, name = "n_mut") %>%
      mutate(TMB = n_mut / EXOME_MB)
    tmb <- tibble(sample_id = unique(clin$sampleId)) %>%
      left_join(tmb_counts %>% rename(sample_id = sampleId), by = "sample_id") %>%
      mutate(TMB = replace_na(TMB, 0.0))
  }

  # Map entrezGeneId → Hugo Symbol for target genes
  gene_map <- tibble(
    entrezGeneId = ENTREZ_IDS,
    gene = names(ENTREZ_IDS)
  )

  mut_func <- mut %>%
    filter(mutationType %in% FUNCTIONAL_CLASSES,
           entrezGeneId %in% ENTREZ_IDS) %>%
    left_join(gene_map, by = "entrezGeneId") %>%
    select(sample_id = sampleId, gene)

  all_samples <- unique(clin$sampleId)
  mat <- build_binary_matrix(all_samples, mut_func, "gene", tmb, "TMB")

  cat("  →", nrow(mat), "samples\n")
  list(data = mat, cohort = cohort_name, n = nrow(mat))
}

# ============================================================
# SECTION 2: Analysis Function
# ============================================================

run_logistic <- function(cohort_obj, gene_b) {
  df     <- cohort_obj$data
  cohort <- cohort_obj$cohort

  # Check gene coverage (gene present in data AND observed in ≥1 sample)
  if (!gene_b %in% names(df) || sum(df[[gene_b]], na.rm = TRUE) == 0) {
    return(tibble(
      cohort = cohort, gene_b = gene_b,
      n = NA_integer_, n_tp53 = NA_integer_, n_gene = NA_integer_, n_both = NA_integer_,
      OR_adj = NA_real_, CI_lower = NA_real_, CI_upper = NA_real_,
      beta = NA_real_, se = NA_real_, p_value = NA_real_,
      status = "gene_not_covered"
    ))
  }
  if (!"TP53" %in% names(df) || sum(df$TP53, na.rm = TRUE) == 0) {
    return(tibble(
      cohort = cohort, gene_b = gene_b, status = "TP53_absent"
    ))
  }

  model_df <- df %>%
    select(TP53, gene = all_of(gene_b), TMB) %>%
    filter(!is.na(TP53), !is.na(gene), !is.na(TMB)) %>%
    mutate(log2_TMB = log2(TMB))

  n_total <- nrow(model_df)
  n_tp53  <- sum(model_df$TP53)
  n_gene  <- sum(model_df$gene)
  n_both  <- sum(model_df$TP53 == 1L & model_df$gene == 1L)

  # Require at least 5 mutated samples in each margin
  if (n_gene < 5 || n_tp53 < 5) {
    return(tibble(
      cohort = cohort, gene_b = gene_b,
      n = n_total, n_tp53 = n_tp53, n_gene = n_gene, n_both = n_both,
      OR_adj = NA_real_, CI_lower = NA_real_, CI_upper = NA_real_,
      beta = NA_real_, se = NA_real_, p_value = NA_real_,
      status = "too_few_mutations"
    ))
  }

  tryCatch({
    fit <- glm(TP53 ~ gene + log2_TMB, data = model_df, family = binomial())
    cf  <- summary(fit)$coefficients
    beta <- cf["gene", "Estimate"]
    se   <- cf["gene", "Std. Error"]
    pval <- cf["gene", "Pr(>|z|)"]

    tibble(
      cohort  = cohort, gene_b  = gene_b,
      n       = n_total, n_tp53 = n_tp53, n_gene = n_gene, n_both = n_both,
      OR_adj  = exp(beta),
      CI_lower = exp(beta - 1.96 * se),
      CI_upper = exp(beta + 1.96 * se),
      beta    = beta, se = se, p_value = pval,
      status  = "ok"
    )
  }, error = function(e) {
    tibble(
      cohort = cohort, gene_b = gene_b,
      n = n_total, n_tp53 = n_tp53, n_gene = n_gene, n_both = n_both,
      OR_adj = NA_real_, CI_lower = NA_real_, CI_upper = NA_real_,
      beta = NA_real_, se = NA_real_, p_value = NA_real_,
      status = paste0("glm_error: ", conditionMessage(e))
    )
  })
}

# ============================================================
# SECTION 3: Load All Cohorts
# ============================================================

cohorts <- list()

cohorts[["TCGA-STAD"]]             <- load_tcga(PATHS$tcga_maf)
cohorts[["egc_tmucih_2015"]]       <- load_cbioportal("egc_tmucih_2015",       PATHS$cbio_base)
cohorts[["stad_oncosg_2018"]]      <- load_cbioportal("stad_oncosg_2018",      PATHS$cbio_base)
cohorts[["stad_pfizer_uhongkong"]] <- load_cbioportal("stad_pfizer_uhongkong", PATHS$cbio_base)
cohorts[["egc_msk_2017"]]          <- load_cbioportal("egc_msk_2017",          PATHS$cbio_base, stad_only = TRUE)
cohorts[["egc_mskcc_2020"]]        <- load_cbioportal("egc_mskcc_2020",        PATHS$cbio_base, stad_only = TRUE)
cohorts[["egc_msk_tp53_ccr_2022"]] <- load_cbioportal("egc_msk_tp53_ccr_2022", PATHS$cbio_base, stad_only = TRUE)
cohorts[["egc_msk_2023"]]          <- load_cbioportal("egc_msk_2023",          PATHS$cbio_base, stad_only = TRUE)

# Update display labels with actual N
for (cn in names(cohorts)) {
  actual_n <- cohorts[[cn]]$n
  if (cn %in% names(COHORT_LABELS)) {
    COHORT_LABELS[cn] <- str_replace(COHORT_LABELS[cn], "N=\\d+", paste0("N=", actual_n))
  }
}

# ============================================================
# SECTION 4: Run All Models
# ============================================================

cat("\n--- Running logistic regression ---\n")
results_list <- list()
for (cn in names(cohorts)) {
  for (g in TARGET_GENES) {
    results_list[[paste(cn, g, sep = "|")]] <- run_logistic(cohorts[[cn]], g)
  }
}
results <- bind_rows(results_list)

# Coverage summary
coverage <- results %>%
  select(cohort, gene_b, status) %>%
  mutate(covered = status != "gene_not_covered")
write_tsv(coverage, file.path(PATHS$output_dir, "panel_coverage_check.tsv"))
cat("\nPanel coverage check:\n")
print(coverage %>% pivot_wider(names_from = gene_b, values_from = covered))

# Merge with FDU discovery
results_all <- bind_rows(
  FDU_RESULTS %>% select(cohort, gene_b, n, OR_adj, CI_lower, CI_upper,
                          beta, se, p_value, n_tp53, n_gene, n_both, status),
  results
)

write_tsv(results_all, file.path(PATHS$output_dir, "per_cohort_results.tsv"))
cat("\nPer-cohort results saved.\n")

# ============================================================
# SECTION 5: Meta-Analysis
# ============================================================

cat("\n--- Meta-analysis (random effects, REML) ---\n")
meta_list <- list()

for (g in TARGET_GENES) {
  gene_data <- results_all %>%
    filter(gene_b == g, status %in% c("ok", "discovery"),
           !is.na(beta), !is.na(se), se > 0, is.finite(beta))

  if (nrow(gene_data) < 2) {
    cat("  Skipping", g, ": <2 cohorts with valid estimates\n")
    next
  }

  tryCatch({
    ma <- rma(yi = beta, sei = se, data = gene_data, method = "REML")
    meta_list[[g]] <- tibble(
      gene_b       = g,
      n_cohorts    = nrow(gene_data),
      pooled_OR    = exp(coef(ma)),
      CI_lower     = exp(ma$ci.lb),
      CI_upper     = exp(ma$ci.ub),
      p_value      = ma$pval,
      I2           = round(ma$I2, 1),
      tau2         = ma$tau2,
      # store for forest plot
      beta_pooled  = coef(ma),
      se_pooled    = ma$se
    )
    cat(sprintf("  %-8s  pooled OR=%.3f [%.3f–%.3f]  I²=%.1f%%\n",
                g, exp(coef(ma)), exp(ma$ci.lb), exp(ma$ci.ub), ma$I2))
  }, error = function(e) {
    cat("  Meta-analysis error for", g, ":", conditionMessage(e), "\n")
  })
}

meta_df <- bind_rows(meta_list)
write_tsv(meta_df %>% select(-any_of(c("beta_pooled", "se_pooled"))),
          file.path(PATHS$output_dir, "meta_analysis_results.tsv"))
cat("Meta-analysis results saved.\n")

# ============================================================
# SECTION 6: Forest Plots
# ============================================================

cat("\n--- Generating forest plots ---\n")

plot_forest <- function(gene_name, results_df, meta_df) {
  gene_data <- results_df %>%
    filter(gene_b == gene_name,
           status %in% c("ok", "discovery"),
           !is.na(OR_adj), !is.na(CI_lower), !is.na(CI_upper),
           is.finite(OR_adj), is.finite(CI_lower), is.finite(CI_upper))

  pooled_row <- meta_df %>% filter(gene_b == gene_name)
  if (nrow(pooled_row) > 0) {
    gene_data <- bind_rows(
      gene_data,
      tibble(
        cohort   = "Pooled (RE)", gene_b = gene_name,
        n        = sum(gene_data$n, na.rm = TRUE),
        OR_adj   = pooled_row$pooled_OR,
        CI_lower = pooled_row$CI_lower,
        CI_upper = pooled_row$CI_upper,
        p_value  = pooled_row$p_value,
        status   = "pooled",
        n_tp53 = NA_integer_, n_gene = NA_integer_, n_both = NA_integer_
      )
    )
  }

  # Clamp CI for display (avoid extreme axes)
  gene_data <- gene_data %>%
    mutate(
      CI_lower_disp = pmax(CI_lower, 0.01),
      CI_upper_disp = pmin(CI_upper, 100),
      cohort_f = factor(cohort, levels = rev(COHORT_ORDER)),
      direction = case_when(
        status == "pooled" ~ "pooled",
        OR_adj < 1         ~ "ME",
        TRUE               ~ "Co-occur"
      ),
      sig   = !is.na(p_value) & p_value < 0.05,
      label = sapply(as.character(cohort), function(x) {
        if (!is.na(COHORT_LABELS[x])) COHORT_LABELS[x] else x
      })
    )

  I2_str <- if (nrow(pooled_row) > 0) {
    sprintf("I² = %.1f%%  |  Pooled OR = %.3f  [%.3f–%.3f]",
            pooled_row$I2, pooled_row$pooled_OR,
            pooled_row$CI_lower, pooled_row$CI_upper)
  } else { "" }

  p <- ggplot(gene_data, aes(x = OR_adj, y = cohort_f)) +
    # Separator line above pooled
    geom_hline(
      data = gene_data %>% filter(status == "pooled"),
      aes(yintercept = as.numeric(cohort_f) + 0.5),
      color = "grey70", linetype = "solid", linewidth = 0.4
    ) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    # Error bars
    geom_errorbarh(
      aes(xmin = CI_lower_disp, xmax = CI_upper_disp, color = direction),
      height = 0.25, linewidth = 0.65
    ) +
    # Points
    geom_point(
      aes(fill  = direction,
          color = direction,
          size  = ifelse(status == "pooled", 4.5, 2.8),
          shape = ifelse(status == "pooled", 23L, 21L)),
      stroke = 0.8
    ) +
    # Significance star
    geom_text(
      data = gene_data %>% filter(sig, status != "pooled"),
      aes(x = CI_upper_disp * 1.05, label = "*"),
      hjust = 0, size = 4.5, color = "#333333"
    ) +
    scale_shape_identity() +
    scale_size_identity() +
    scale_fill_manual(
      values = c("ME" = "#d73027", "Co-occur" = "#4575b4", "pooled" = "#1a1a1a"),
      guide  = "none"
    ) +
    scale_color_manual(
      values = c("ME" = "#d73027", "Co-occur" = "#4575b4", "pooled" = "#1a1a1a"),
      guide  = "none"
    ) +
    scale_x_log10(
      breaks = c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10),
      labels = c("0.05", "0.1", "0.2", "0.5", "1", "2", "5", "10")
    ) +
    scale_y_discrete(
      labels = function(x) {
        sapply(as.character(x), function(v) {
          lbl <- COHORT_LABELS[v]
          if (!is.na(lbl)) lbl else v
        })
      }
    ) +
    labs(
      title    = sprintf("TP53  ↔  %s   (TMB-adjusted logistic regression)", gene_name),
      subtitle = I2_str,
      x        = "Adjusted OR (log scale)   ← ME  |  Co-occur →",
      y        = NULL,
      caption  = "* p < 0.05;  error bars = 95% CI;  FDU = discovery cohort"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position  = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.title    = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 9, color = "grey40"),
      plot.caption  = element_text(size = 8, color = "grey50"),
      axis.text.y   = element_text(size = 9)
    )

  out_path <- file.path(PATHS$output_dir, sprintf("forest_%s.png", gene_name))
  ggsave(out_path, plot = p, width = 9, height = 5.5, dpi = 200)
  cat("  Saved:", basename(out_path), "\n")
  invisible(p)
}

for (g in TARGET_GENES) {
  plot_forest(g, results_all, meta_df)
}

# ============================================================
# SECTION 7: Summary Heatmap
# ============================================================

cat("\n--- Generating summary heatmap ---\n")

# Add pooled estimates to heatmap
pooled_for_heatmap <- meta_df %>%
  transmute(
    cohort  = "Pooled (RE)",
    gene_b  = gene_b,
    n       = NA_integer_, n_tp53 = NA_integer_, n_gene = NA_integer_, n_both = NA_integer_,
    OR_adj  = pooled_OR,
    CI_lower = CI_lower, CI_upper = CI_upper,
    beta    = beta_pooled, se = se_pooled,
    p_value = p_value,
    status  = "pooled"
  )

heatmap_data <- bind_rows(
  results_all %>% filter(status %in% c("ok", "discovery")),
  pooled_for_heatmap
) %>%
  mutate(
    log2_OR = log2(pmax(pmin(OR_adj, 16), 0.0625)),  # clamp to [-4, +4]
    sig_label = case_when(
      !is.na(p_value) & p_value < 0.001 ~ "***",
      !is.na(p_value) & p_value < 0.01  ~ "**",
      !is.na(p_value) & p_value < 0.05  ~ "*",
      TRUE ~ ""
    ),
    cohort_f = factor(cohort, levels = rev(COHORT_ORDER)),
    gene_f   = factor(gene_b, levels = TARGET_GENES)
  )

cohort_labels_fn <- function(x) {
  sapply(as.character(x), function(v) {
    lbl <- COHORT_LABELS[v]
    if (!is.na(lbl)) lbl else v
  })
}

p_heat <- ggplot(heatmap_data, aes(x = gene_f, y = cohort_f, fill = log2_OR)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sig_label), size = 4.5, vjust = 0.75, color = "white") +
  scale_fill_gradientn(
    colors = c(
      "#313695", "#4575b4", "#74add1", "#abd9e9",
      "#f7f7f7",
      "#fdae61", "#f46d43", "#d73027", "#a50026"
    ),
    limits   = c(-4, 4),
    na.value = "grey88",
    name     = "log₂(OR_adj)",
    guide    = guide_colorbar(barwidth = 8, barheight = 0.8,
                              title.position = "top", title.hjust = 0.5)
  ) +
  scale_y_discrete(labels = cohort_labels_fn) +
  scale_x_discrete(position = "top") +
  labs(
    title    = "TP53 ME candidates — cross-cohort OR_adj summary",
    subtitle = "Red = mutually exclusive  |  Blue = co-occurring  |  Grey = not covered / too few mutations\n* p<0.05  ** p<0.01  *** p<0.001",
    x = NULL, y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x     = element_text(angle = 30, hjust = 0, face = "bold", size = 10),
    axis.text.y     = element_text(size = 9),
    panel.grid      = element_blank(),
    legend.position = "bottom",
    plot.title      = element_text(size = 12, face = "bold"),
    plot.subtitle   = element_text(size = 9, color = "grey40")
  )

ggsave(file.path(PATHS$output_dir, "summary_heatmap.png"),
       plot = p_heat, width = 10, height = 7, dpi = 200)
cat("  Saved: summary_heatmap.png\n")

# ============================================================
# SECTION 8: Print Final Summary
# ============================================================

cat("\n=== Analysis complete ===\n")
cat("Output directory:", PATHS$output_dir, "\n\n")

cat("Meta-analysis summary:\n")
if (nrow(meta_df) > 0) {
  meta_df %>%
    arrange(p_value) %>%
    mutate(
      pooled_OR = round(pooled_OR, 3),
      CI_lower  = round(CI_lower, 3),
      CI_upper  = round(CI_upper, 3),
      p_value   = signif(p_value, 3),
      I2        = paste0(round(I2, 1), "%")
    ) %>%
    select(gene_b, n_cohorts, pooled_OR, CI_lower, CI_upper, p_value, I2) %>%
    { print(as.data.frame(.)); invisible(.) }
}

cat("\nFiles written:\n")
cat("  per_cohort_results.tsv\n")
cat("  meta_analysis_results.tsv\n")
cat("  panel_coverage_check.tsv\n")
for (g in TARGET_GENES) cat(sprintf("  forest_%s.png\n", g))
cat("  summary_heatmap.png\n")
