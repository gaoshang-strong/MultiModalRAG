#!/usr/bin/env Rscript
#
# depmap_prism_ATM_TP53.R
#
# Analyses:
#   1. DepMap CRISPR — ATM KO effect in TP53-mut vs WT cell lines
#   2. GDSC2 KU-55933 — single-agent ATM inhibitor sensitivity vs TP53 status
#   3. PRISM AZD0156 — single-agent ATM inhibitor sensitivity vs TP53 status
#   4. PRISM AZD0156 combos — synergy (Bliss) vs TP53 status
#        Partners: AZD1775 (WEE1i), AZD6738 (ATRi), AZD2811 (AURKBi)
#
# Each analysis run pan-cancer AND gastric-only.
# Stats: Wilcoxon rank-sum + linear regression with lineage covariate.
#
# Usage:
#   micromamba run -n ProjectGeneration \
#     Rscript TP53_ATM_paper/scripts/depmap_prism_ATM_TP53.R

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(broom)
  library(RColorBrewer)
})

# ── Paths ─────────────────────────────────────────────────────────────────────
get_script_dir <- function() {
  args     <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("--file=", "", file_arg)))
  else normalizePath(".")
}
SCRIPT_DIR <- get_script_dir()
REPO_ROOT  <- normalizePath(file.path(SCRIPT_DIR, "../.."))

DEPMAP_DIR <- file.path(REPO_ROOT, "DepMap")
PRISM_DIR  <- file.path(REPO_ROOT, "PRISM")
OUT_BASE   <- file.path(REPO_ROOT, "TP53_ATM_paper/results/DepMap_PRISM")

for (d in file.path(OUT_BASE, c("depmap_crispr","gdsc2_ku55933",
                                 "prism_azd0156_single","prism_azd0156_combos")))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Gastric lineage terms (DepMap OncotreeLineage / OncotreePrimaryDisease)
GASTRIC_LINEAGE <- c("Esophagus/Stomach")
GASTRIC_DISEASE <- c("Gastric Adenocarcinoma","Esophageal Adenocarcinoma",
                     "Esophageal Squamous Cell Carcinoma")
# Gastric TCGA descriptor in GDSC
GASTRIC_TCGA    <- c("STAD","ESCA","STES")

TP53_COLORS <- c(mutant="#E41A1C", wildtype="#377EB8")


# ═══════════════════════════════════════════════════════════════════════════════
# 1. LOAD MASTER DATA
# ═══════════════════════════════════════════════════════════════════════════════

message("Loading DepMap metadata and TP53 status...")
model   <- fread(file.path(DEPMAP_DIR, "Model.csv"))
tp53_df <- fread(file.path(DEPMAP_DIR, "mutations_TP53_damaging.csv"))

# CN homozygous deletion (< 0.3 relative CN) — matches Python skill _data.py
tp53_cn_raw <- fread(file.path(DEPMAP_DIR, "OmicsCNGene.csv"),
                     select = c("V1", "TP53 (7157)"))
setnames(tp53_cn_raw, c("ModelID", "TP53_cn"))
tp53_df <- merge(tp53_df, tp53_cn_raw[, .(ModelID, TP53_cn)],
                 by = "ModelID", all.x = TRUE)
tp53_df[is.na(TP53_cn), TP53_cn := 1.0]  # assume diploid if not profiled
# mutant = damaging mutation OR homozygous CN deletion (CN < 0.3)
tp53_df[, tp53_status := ifelse(TP53 > 0 | TP53_cn < 0.3, "mutant", "wildtype")]

# Master TP53 table with lineage info
tp53_master <- merge(tp53_df[, .(ModelID, tp53_status)],
                     model[, .(ModelID, CellLineName, StrippedCellLineName,
                               OncotreeLineage, OncotreePrimaryDisease,
                               OncotreeSubtype, OncotreeCode)],
                     by = "ModelID", all.x = TRUE)
tp53_master[, is_gastric := OncotreeLineage %in% GASTRIC_LINEAGE |
              OncotreePrimaryDisease %in% GASTRIC_DISEASE]

message(sprintf("  TP53-mut cell lines: %d / %d total",
  sum(tp53_master$tp53_status=="mutant"), nrow(tp53_master)))
message(sprintf("  Gastric cell lines: %d (mut=%d, wt=%d)",
  sum(tp53_master$is_gastric),
  sum(tp53_master$is_gastric & tp53_master$tp53_status=="mutant"),
  sum(tp53_master$is_gastric & tp53_master$tp53_status=="wildtype")))

# Normalised name for matching GDSC/PRISM cell lines
norm_name <- function(x) toupper(gsub("[^A-Z0-9]", "", toupper(x)))
tp53_master[, name_key := norm_name(CellLineName)]


# ═══════════════════════════════════════════════════════════════════════════════
# 2. HELPER FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════════

# Wilcoxon + effect size (rank-biserial r)
wilcox_summary <- function(mut_vals, wt_vals) {
  if (length(mut_vals) < 3 || length(wt_vals) < 3) return(NULL)
  wt_res <- wilcox.test(mut_vals, wt_vals, exact = FALSE)
  n1 <- length(mut_vals); n2 <- length(wt_vals)
  r  <- 1 - 2 * wt_res$statistic / (n1 * n2)   # rank-biserial r
  list(p = wt_res$p.value, r = as.numeric(r),
       n_mut = n1, n_wt = n2,
       median_mut = median(mut_vals, na.rm=TRUE),
       median_wt  = median(wt_vals,  na.rm=TRUE))
}

# Bootstrap 95% CI on delta = mean(mut) - mean(WT) — matches Python skill _stats.py
bootstrap_ci <- function(mut_vals, wt_vals, n_boot = 2000, seed = 42) {
  set.seed(seed)
  deltas <- replicate(n_boot, {
    mean(sample(mut_vals, length(mut_vals), replace = TRUE)) -
    mean(sample(wt_vals,  length(wt_vals),  replace = TRUE))
  })
  list(delta = mean(mut_vals) - mean(wt_vals),
       ci_lo  = unname(quantile(deltas, 0.025)),
       ci_hi  = unname(quantile(deltas, 0.975)))
}

# Linear regression: metric ~ tp53_status + lineage (controls for tissue type)
lm_tp53 <- function(df, metric_col) {
  df <- df[!is.na(df[[metric_col]]) & !is.na(df$tp53_status) &
           !is.na(df$OncotreeLineage), ]
  if (sum(df$tp53_status=="mutant") < 5 || sum(df$tp53_status=="wildtype") < 5)
    return(NULL)
  df$tp53_bin <- as.integer(df$tp53_status == "mutant")
  # Drop lineages with <3 cell lines
  lin_counts <- table(df$OncotreeLineage)
  df <- df[df$OncotreeLineage %in% names(lin_counts[lin_counts >= 3]), ]
  df$OncotreeLineage <- factor(df$OncotreeLineage)
  frm <- as.formula(paste0(metric_col, " ~ tp53_bin + OncotreeLineage"))
  fit <- tryCatch(lm(frm, data=df), error=function(e) NULL)
  if (is.null(fit)) return(NULL)
  s <- summary(fit)
  coef_row <- s$coefficients["tp53_bin", ]
  list(estimate  = coef_row["Estimate"],
       se        = coef_row["Std. Error"],
       t         = coef_row["t value"],
       p         = coef_row["Pr(>|t|)"],
       r2_adj    = s$adj.r.squared,
       n         = nrow(df))
}

# Core violin + jitter plot
plot_violin <- function(df, y_col, y_label, title_str, subtitle_str,
                        outfile, gastric_col = "is_gastric",
                        mut_n = NULL, wt_n = NULL) {
  df <- df[!is.na(df[[y_col]]) & !is.na(df$tp53_status), ]
  df$tp53_status <- factor(df$tp53_status, levels = c("mutant","wildtype"))

  # Subtitle includes sample sizes if not provided
  if (!is.null(mut_n))
    subtitle_str <- paste0(subtitle_str, sprintf("  (mut N=%d, wt N=%d)", mut_n, wt_n))

  p <- ggplot(df, aes_string(x="tp53_status", y=y_col, fill="tp53_status")) +
    geom_violin(trim=TRUE, alpha=0.35, linewidth=0.4) +
    geom_boxplot(width=0.18, outlier.shape=NA, linewidth=0.5, alpha=0.7) +
    geom_jitter(data = df[df[[gastric_col]] == FALSE, ],
                width=0.18, size=0.55, alpha=0.25, color="grey50") +
    geom_jitter(data = df[df[[gastric_col]] == TRUE, ],
                aes(color="Gastric"), width=0.12, size=2.5, alpha=0.85) +
    scale_fill_manual(values=TP53_COLORS, guide="none") +
    scale_color_manual(values=c(Gastric="#FF7F00"), name=NULL) +
    geom_text_repel(data = df[df[[gastric_col]] == TRUE, ],
                    aes_string(label="cell_line_label"),
                    size=2.3, color="#FF7F00", max.overlaps=20,
                    segment.size=0.3, segment.alpha=0.5) +
    labs(title=title_str, subtitle=subtitle_str, x="TP53 status", y=y_label) +
    theme_bw(base_size=12) +
    theme(legend.position="bottom")

  ggsave(outfile, p, width=7, height=6, dpi=150)
  message(sprintf("  Saved: %s", outfile))
}

# Run one full analysis (pan-cancer + gastric) for a given metric
run_analysis <- function(df, metric_col, metric_label, title_base,
                         outdir, prefix, ref_lineage_col="OncotreeLineage") {
  df <- as.data.frame(df)
  df$cell_line_label <- ifelse(df$is_gastric, df$CellLineName, "")

  results <- list()

  for (scope in c("pancancer","gastric")) {
    sub <- if (scope == "gastric") df[df$is_gastric, ] else df
    mut_v <- sub[[metric_col]][sub$tp53_status=="mutant" & !is.na(sub[[metric_col]])]
    wt_v  <- sub[[metric_col]][sub$tp53_status=="wildtype"& !is.na(sub[[metric_col]])]

    wx <- wilcox_summary(mut_v, wt_v)
    lm <- if (scope=="pancancer") lm_tp53(sub, metric_col) else NULL

    if (is.null(wx)) {
      message(sprintf("  [%s | %s] too few samples, skipping", prefix, scope)); next
    }

    boot   <- bootstrap_ci(mut_v, wt_v)
    p_str  <- if (wx$p < 0.001) sprintf("p=%.2e", wx$p) else sprintf("p=%.3f", wx$p)
    r_str  <- sprintf("r=%.2f", wx$r)
    lm_str <- if (!is.null(lm)) sprintf("  |  lm(lineage-adj) β=%.3f p=%.3f",
                                         lm$estimate, lm$p) else ""
    sub_str <- paste0("Wilcoxon ", p_str, " (", r_str, ")",
                      "  n_mut=", wx$n_mut, " n_wt=", wx$n_wt, lm_str)
    scope_label <- if (scope=="gastric") "Gastric only" else "Pan-cancer"
    title_str   <- paste0(title_base, " [", scope_label, "]")

    outfile <- file.path(outdir, paste0(prefix, "_", scope, ".png"))
    plot_violin(sub, metric_col, metric_label, title_str, sub_str,
                outfile, mut_n=wx$n_mut, wt_n=wx$n_wt)

    results[[scope]] <- data.frame(
      scope=scope, metric=metric_col,
      n_mut=wx$n_mut, n_wt=wx$n_wt,
      median_mut=wx$median_mut, median_wt=wx$median_wt,
      delta  = boot$delta,
      ci_lo  = boot$ci_lo,
      ci_hi  = boot$ci_hi,
      wilcox_p=wx$p, effect_r=wx$r,
      lm_beta   = if (!is.null(lm)) lm$estimate else NA,
      lm_p      = if (!is.null(lm)) lm$p        else NA,
      lm_r2adj  = if (!is.null(lm)) lm$r2_adj   else NA
    )
    message(sprintf("  [%s | %s] mut_med=%.3f wt_med=%.3f Wilcoxon %s%s",
      prefix, scope,
      wx$median_mut, wx$median_wt, p_str, lm_str))
  }
  if (length(results) > 0) do.call(rbind, results) else NULL
}


# ═══════════════════════════════════════════════════════════════════════════════
# 3. DEPMAP CRISPR — ATM KO EFFECT
# ═══════════════════════════════════════════════════════════════════════════════
message("\n════════ DepMap CRISPR — ATM KO ════════")

crispr <- fread(file.path(DEPMAP_DIR, "CRISPRGeneEffect.csv"))
setnames(crispr, 1, "ModelID")

# Extract ATM column (format: "ATM (472)")
atm_col <- grep("^ATM \\(", names(crispr), value=TRUE)
atm_df  <- crispr[, .(ModelID, ATM_effect = .SD[[1]]), .SDcols = atm_col]
atm_df  <- merge(atm_df, tp53_master, by="ModelID", all.x=TRUE)
atm_df  <- atm_df[!is.na(tp53_status)]

# All TP53-mut vs TP53-WT (ATM-mut included — ATM KO renders endogenous ATM status irrelevant)
dep_res <- run_analysis(
  atm_df, "ATM_effect",
  "ATM CRISPR Gene Effect (negative = essential)",
  "ATM CRISPR KO — TP53-mut vs TP53-WT",
  file.path(OUT_BASE, "depmap_crispr"), "atm_ko_allTP53")
if (!is.null(dep_res)) dep_res$comparison <- "all_TP53"

write.csv(dep_res, file.path(OUT_BASE, "depmap_crispr", "stats_summary.csv"),
          row.names=FALSE)
write.csv(as.data.frame(atm_df),
          file.path(OUT_BASE, "depmap_crispr", "atm_effect_per_cell_line.csv"),
          row.names=FALSE)


# ═══════════════════════════════════════════════════════════════════════════════
# 4. GDSC2 KU-55933 — SINGLE-AGENT ATM INHIBITOR
# ═══════════════════════════════════════════════════════════════════════════════
message("\n════════ GDSC2 KU-55933 ════════")

# Read GDSC2 (header is on row 5, 0-indexed = skip 4 rows)
gdsc2_raw <- read_excel(file.path(PRISM_DIR, "GDSC2_dose_response.xlsx"), skip=4,
                         col_names=FALSE)
gdsc2_cols <- c("DATASET","NLME_RESULT_ID","NLME_CURVE_ID","COSMIC_ID",
                "CELL_LINE_NAME","SANGER_MODEL_ID","TCGA_DESC","DRUG_ID",
                "DRUG_NAME","TARGET","TARGET_PATHWAY","PUTATIVE_TARGET",
                "GDSC_tissue","MIN_CONC","MAX_CONC","LN_IC50","AUC","RMSE","Z_SCORE")
names(gdsc2_raw) <- gdsc2_cols[seq_len(ncol(gdsc2_raw))]

ku <- as.data.table(gdsc2_raw[gdsc2_raw$DRUG_NAME == "KU-55933", ])
ku[, LN_IC50  := as.numeric(LN_IC50)]
ku[, AUC      := as.numeric(AUC)]
ku[, name_key := norm_name(CELL_LINE_NAME)]
message(sprintf("  KU-55933: %d cell lines", nrow(ku)))

# Link to TP53 status
ku_merged <- merge(ku, tp53_master[, .(name_key, ModelID, tp53_status,
                                        CellLineName, OncotreeLineage,
                                        OncotreePrimaryDisease, is_gastric)],
                   by="name_key", all.x=TRUE)
ku_merged[is.na(is_gastric), is_gastric := TCGA_DESC %in% GASTRIC_TCGA]
ku_merged <- ku_merged[!is.na(tp53_status)]
message(sprintf("  KU-55933: %d cell lines (mut=%d wt=%d)",
  nrow(ku_merged),
  sum(ku_merged$tp53_status=="mutant"),
  sum(ku_merged$tp53_status=="wildtype")))

for (metric in c("LN_IC50","AUC")) {
  mlabel <- if (metric=="LN_IC50") "LN(IC50) [KU-55933]" else "AUC [KU-55933]"
  r <- run_analysis(ku_merged, metric, mlabel,
                    paste0("GDSC2 KU-55933 (ATM inhibitor) — ", metric,
                           " by TP53 Status"),
                    file.path(OUT_BASE,"gdsc2_ku55933"),
                    paste0("ku55933_",tolower(metric)))
  if (!is.null(r)) r$drug <- "KU-55933"
  assign(paste0("ku_res_",metric), r)
}
ku_stats <- rbind(get("ku_res_LN_IC50"), get("ku_res_AUC"))
write.csv(ku_stats, file.path(OUT_BASE,"gdsc2_ku55933","stats_summary.csv"),
          row.names=FALSE)
write.csv(as.data.frame(ku_merged),
          file.path(OUT_BASE,"gdsc2_ku55933","ku55933_per_cell_line.csv"),
          row.names=FALSE)


# ═══════════════════════════════════════════════════════════════════════════════
# 5. PRISM AZD0156 — SINGLE-AGENT
# ═══════════════════════════════════════════════════════════════════════════════
message("\n════════ PRISM AZD0156 — single-agent ════════")

prism <- fread(file.path(PRISM_DIR, "gdsc-007_matrix_results.csv"))

# AZD0156 as lib1: deduplicate per cell line (lib1 metrics identical across
# combo partners for same cell line — take first occurrence)
azd_lib1 <- prism[lib1_name == "AZD0156",
                   .(lib1_IC50_ln = mean(lib1_IC50_ln, na.rm=TRUE),
                     lib1_MaxE    = mean(lib1_MaxE,    na.rm=TRUE)),
                   by = .(CELL_LINE_NAME, TISSUE, CANCER_TYPE)]
azd_lib1[, name_key := norm_name(CELL_LINE_NAME)]
message(sprintf("  AZD0156: %d unique cell lines", nrow(azd_lib1)))

azd_merged <- merge(azd_lib1,
                    tp53_master[, .(name_key, ModelID, tp53_status,
                                    CellLineName, OncotreeLineage,
                                    OncotreePrimaryDisease, is_gastric)],
                    by="name_key", all.x=TRUE)
azd_merged[is.na(is_gastric), is_gastric :=
             grepl("Stomach|Gastric|oesophag|esophag", TISSUE, ignore.case=TRUE)]
azd_merged <- azd_merged[!is.na(tp53_status)]
message(sprintf("  AZD0156 matched: %d cell lines (mut=%d wt=%d)",
  nrow(azd_merged),
  sum(azd_merged$tp53_status=="mutant"),
  sum(azd_merged$tp53_status=="wildtype")))

for (metric in c("lib1_IC50_ln","lib1_MaxE")) {
  mlabel <- if (metric=="lib1_IC50_ln") "LN(IC50) [AZD0156]" else "MaxE (max inhibition) [AZD0156]"
  r <- run_analysis(azd_merged, metric, mlabel,
                    paste0("PRISM AZD0156 (ATM inhibitor) — ", metric,
                           " by TP53 Status"),
                    file.path(OUT_BASE,"prism_azd0156_single"),
                    paste0("azd0156_",sub("lib1_","",metric)))
  if (!is.null(r)) r$drug <- "AZD0156"
  assign(paste0("azd_res_",metric), r)
}
azd_stats <- rbind(get("azd_res_lib1_IC50_ln"), get("azd_res_lib1_MaxE"))
write.csv(azd_stats, file.path(OUT_BASE,"prism_azd0156_single","stats_summary.csv"),
          row.names=FALSE)
write.csv(as.data.frame(azd_merged),
          file.path(OUT_BASE,"prism_azd0156_single","azd0156_per_cell_line.csv"),
          row.names=FALSE)


# ═══════════════════════════════════════════════════════════════════════════════
# 6. PRISM AZD0156 COMBINATIONS — BLISS SYNERGY vs TP53 STATUS
# ═══════════════════════════════════════════════════════════════════════════════
message("\n════════ PRISM AZD0156 combos — Bliss synergy ════════")

COMBOS <- list(
  # AZD0156 as lib1
  AZD1775  = list(partner="AZD1775",  target="WEE1i",   azd_is_lib1=TRUE),
  AZD6738  = list(partner="AZD6738",  target="ATRi",    azd_is_lib1=TRUE),
  AZD2811  = list(partner="AZD2811",  target="AURKBi",  azd_is_lib1=TRUE),
  # AZD7648 (DNA-PKi) as lib1, AZD0156 as lib2 — key combination from Phase4
  AZD7648  = list(partner="AZD7648",  target="DNA-PKi", azd_is_lib1=FALSE)
)

combo_stats_all <- list()

for (nm in names(COMBOS)) {
  partner <- COMBOS[[nm]]$partner
  tgt     <- COMBOS[[nm]]$target
  message(sprintf("  Processing AZD0156 + %s (%s)...", partner, tgt))

  azd_is_lib1 <- COMBOS[[nm]]$azd_is_lib1
  if (azd_is_lib1) {
    sub <- prism[lib1_name=="AZD0156" & lib2_name==partner]
  } else {
    # AZD0156 is lib2; partner (e.g. AZD7648) is lib1
    sub <- prism[lib1_name==partner & lib2_name=="AZD0156"]
  }

  if (nrow(sub) == 0) {
    message(sprintf("  No data for AZD0156 + %s", partner)); next
  }

  # For AZD7648 (lib1), Delta_MaxE_lib1 refers to AZD7648 alone, not AZD0156.
  # Use Delta_MaxE_lib2 instead so "delta" = combo − AZD0156 alone.
  delta_col <- if (azd_is_lib1) "Delta_MaxE_lib1" else "Delta_MaxE_lib2"

  combo_cl <- sub[, .(
    Bliss_window    = mean(Bliss_window,         na.rm=TRUE),
    Bliss_window_SO = mean(Bliss_window_SO,       na.rm=TRUE),
    HSA_window      = mean(HSA_window,            na.rm=TRUE),
    combo_MaxE      = mean(combo_MaxE,            na.rm=TRUE),
    Delta_MaxE_AZD0156 = mean(.SD[[delta_col]],  na.rm=TRUE)
  ), .SDcols = c("Bliss_window","Bliss_window_SO","HSA_window",
                 "combo_MaxE", delta_col),
     by=.(CELL_LINE_NAME, TISSUE, CANCER_TYPE)]
  combo_cl[, name_key := norm_name(CELL_LINE_NAME)]

  combo_merged <- merge(combo_cl,
    tp53_master[, .(name_key, ModelID, tp53_status, CellLineName,
                    OncotreeLineage, OncotreePrimaryDisease, is_gastric)],
    by="name_key", all.x=TRUE)
  combo_merged[is.na(is_gastric), is_gastric :=
    grepl("Stomach|Gastric|oesophag|esophag", TISSUE, ignore.case=TRUE)]
  combo_merged <- combo_merged[!is.na(tp53_status)]
  message(sprintf("  AZD0156 + %s: %d cell lines (mut=%d wt=%d)",
    partner, nrow(combo_merged),
    sum(combo_merged$tp53_status=="mutant"),
    sum(combo_merged$tp53_status=="wildtype")))

  combo_outdir <- file.path(OUT_BASE, "prism_azd0156_combos", nm)
  dir.create(combo_outdir, recursive=TRUE, showWarnings=FALSE)

  metrics_list <- list(
    Bliss_window       = "Bliss Synergy Score (window)",
    HSA_window         = "HSA Synergy Score (window)",
    combo_MaxE         = paste0("Combo MaxE [AZD0156+", partner, "]"),
    Delta_MaxE_AZD0156 = "ΔMaxE (combo − AZD0156 alone)"
  )

  combo_stats_nm <- list()
  for (metric in names(metrics_list)) {
    mlabel <- metrics_list[[metric]]
    r <- run_analysis(combo_merged, metric, mlabel,
                      paste0("AZD0156 + ", partner, " (", tgt, ") — ",
                             metric, " by TP53 Status"),
                      combo_outdir,
                      paste0("azd0156_",nm,"_",metric))
    if (!is.null(r)) { r$drug_combo <- paste0("AZD0156+",partner); combo_stats_nm[[metric]] <- r }
  }
  if (length(combo_stats_nm) > 0) {
    combo_stats_all[[nm]] <- do.call(rbind, combo_stats_nm)
    write.csv(combo_stats_all[[nm]],
              file.path(combo_outdir, "stats_summary.csv"), row.names=FALSE)
  }
  write.csv(as.data.frame(combo_merged),
            file.path(combo_outdir, "combo_per_cell_line.csv"), row.names=FALSE)
}


# ═══════════════════════════════════════════════════════════════════════════════
# 7. SUMMARY TABLE
# ═══════════════════════════════════════════════════════════════════════════════
message("\n════════ Writing master summary ════════")

# Combine all stats (rbindlist with fill=TRUE handles column mismatches)
stats_pieces <- list()
if (!is.null(dep_res))   { dep_res$analysis  <- "DepMap_CRISPR_ATM_KO"; stats_pieces[["dep"]]  <- dep_res }
if (!is.null(ku_stats))  { ku_stats$analysis  <- "GDSC2_KU55933";        stats_pieces[["ku"]]   <- ku_stats }
if (!is.null(azd_stats)) { azd_stats$analysis <- "PRISM_AZD0156";        stats_pieces[["azd"]]  <- azd_stats }
if (length(combo_stats_all) > 0) {
  cs <- do.call(rbind, combo_stats_all)
  cs$analysis <- paste0("PRISM_Combo_", cs$drug_combo)
  stats_pieces[["combo"]] <- cs
}
master_stats <- rbindlist(stats_pieces, fill=TRUE)
write.csv(master_stats,
          file.path(OUT_BASE, "master_stats_summary.csv"), row.names=FALSE)

message("\n=== MASTER SUMMARY ===")
if (!is.null(master_stats)) {
  print(master_stats[, c("analysis","scope","metric","n_mut","n_wt",
                          "median_mut","median_wt","wilcox_p","effect_r",
                          "lm_beta","lm_p")], digits=3)
}

message("\n✓  depmap_prism_ATM_TP53.R complete.")
