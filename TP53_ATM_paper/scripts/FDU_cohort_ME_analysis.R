#!/usr/bin/env Rscript
# =============================================================================
# FDU Cohort — Somatic ME / Co-occurrence Analysis
# TP53-ATM paper, publication-quality figures
# =============================================================================
# Two cohorts:
#   CohortA_TxNaive : somt_filtered.maf ∩ surgery_first
#   CohortB_All     : somt_filtered.maf (all 445 samples)
#
# For each cohort × {full, MSS-only}:
#   fig1_oncoprint   — somatic landscape, top 50 genes + clinical tracks
#   fig2_heatmap     — pairwise interaction matrix (logistic adj.OR)
#   fig3_tp53_me     — TP53-ME genes: mutation rates + adjusted ORs
#   figS_fisher      — maftools somaticInteractions tile (supplementary)
#
# Logistic model (full):  gene_A ~ gene_B + TMB_log + MSI_H + stage + histo
# Logistic model (MSS):   gene_A ~ gene_B + TMB_log + stage + histo
# Multiple testing: BH-FDR; both FDR and uncorrected p retained in tables
#
# Run:
#   /home/sgao30/micromamba/envs/tcga_bioc/bin/Rscript \
#       TP53_ATM_paper/scripts/FDU_cohort_ME_analysis.R
# =============================================================================

suppressPackageStartupMessages({
  library(maftools)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# ── Paths ─────────────────────────────────────────────────────────────────────
ROOT     <- "/ShangGaoAIProjects/Gastric/genome"
FILT_MAF <- file.path(ROOT, "analysis_results/Phase0/enrichment/somt_filtered.maf")
MASTER   <- file.path(ROOT, "analysis_results/Phase0/filtered_master_sample_table.tsv")
OUT      <- file.path(ROOT, "TP53_ATM_paper/results/FDU_cohort")
dir.create(file.path(OUT, "tables"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "figures"), recursive = TRUE, showWarnings = FALSE)

# ── Constants (consistent with Phase 0–1) ─────────────────────────────────────
FUN_CLASSES <- c(
  "Missense_Mutation", "Nonsense_Mutation",
  "Frame_Shift_Del",   "Frame_Shift_Ins",
  "Splice_Site",       "In_Frame_Del", "In_Frame_Ins",
  "Nonstop_Mutation",  "Translation_Start_Site"
)
TOP_N <- 50
MIN_N <- 30   # minimum samples required for logistic regression

# ── Publication theme ─────────────────────────────────────────────────────────
PUB <- theme_classic(base_size = 11) +
  theme(
    axis.text        = element_text(color = "black"),
    axis.title       = element_text(color = "black"),
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 9, color = "grey40"),
    legend.text      = element_text(size = 9),
    legend.title     = element_text(size = 9, face = "bold"),
    strip.text       = element_text(face = "bold", size = 10),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )
ME_COL <- "#C0392B"
CO_COL <- "#2980B9"

# =============================================================================
# 1. LOAD DATA
# =============================================================================
cat("Loading data ...\n")
clin <- fread(MASTER, na.strings = c("", "NA", "X"))

stage_map <- c(IA = "I", IB = "I", IIA = "II", IIB = "II",
               IIIA = "III", IIIB = "III", IIIC = "III", IV = "IV")

clin[, TMB_log  := log2(pmax(suppressWarnings(as.numeric(TMB_report_numeric)), 0.1))]
clin[, MSI_H    := as.integer(MSI %in% c("MSI-H", "MSI"))]
clin[, stage    := factor(stage_map[path_stage], levels = c("I","II","III","IV"))]
clin[, histo    := fcase(
  histology %in% c("腺癌","低分化腺癌","粘液腺癌","肝样腺癌","腺癌和NET"), "adeno",
  histology %in% c("印戒","低黏附性癌","失粘附性癌","黏附性癌","腺癌印戒"), "signet",
  default = "other"
)]
clin[, histo    := factor(histo, levels = c("adeno","signet","other"))]
clin[, MSI_cat  := ifelse(is.na(MSI), NA_character_,
                   ifelse(MSI %in% c("MSI-H","MSI"), "MSI-H", "MSS"))]
clin[, EBV_cat  := ifelse(!is.na(EBER) & EBER == 1, "EBV+", "EBV-")]
clin[, HER2_cat := ifelse(!is.na(HER2) & HER2 == "3+", "HER2+", "HER2-")]
clin[, Stage_cat := ifelse(is.na(stage), "Unknown", as.character(stage))]

naive_ids <- clin[treatment_group == "surgery_first", sample_id]
cat(sprintf("  surgery_first (treatment-naive): %d\n", length(naive_ids)))

maf_raw <- fread(FILT_MAF)
maf_raw[, sample_base := sub("F\\d+$", "", Tumor_Sample_Barcode)]
cat(sprintf("  somt_filtered.maf: %d variants | %d samples\n",
            nrow(maf_raw), uniqueN(maf_raw$sample_base)))

# =============================================================================
# 2. ANALYSIS HELPERS
# =============================================================================

# ── Build maftools object with clinical annotation ─────────────────────────────
build_maf <- function(keep_ids) {
  maf_sub  <- maf_raw[sample_base %in% keep_ids,
                      .SD, .SDcols = setdiff(names(maf_raw), "sample_base")]
  clin_sub <- clin[sample_id %in% keep_ids]
  clin_maf <- data.frame(
    Tumor_Sample_Barcode = paste0(clin_sub$sample_id, "F01"),
    MSI   = clin_sub$MSI_cat,
    EBV   = clin_sub$EBV_cat,
    HER2  = clin_sub$HER2_cat,
    Stage = clin_sub$Stage_cat,
    Sex   = clin_sub$sex,
    stringsAsFactors = FALSE
  )
  read.maf(maf = maf_sub, clinicalData = clin_maf,
           vc_nonSyn = FUN_CLASSES, verbose = FALSE)
}

# ── Build sample × gene binary matrix + clinical covariates ────────────────────
build_model_df <- function(maf_obj, strata = "all") {
  mut_dt   <- unique(as.data.table(maf_obj@data)[, .(Tumor_Sample_Barcode, Hugo_Symbol)])
  mut_dt[, v := 1L]
  mut_wide <- dcast(mut_dt, Tumor_Sample_Barcode ~ Hugo_Symbol, value.var = "v", fill = 0L)
  all_samp <- data.table(Tumor_Sample_Barcode = getSampleSummary(maf_obj)$Tumor_Sample_Barcode)
  mut_wide <- merge(all_samp, mut_wide, by = "Tumor_Sample_Barcode", all.x = TRUE)
  for (col in setdiff(names(mut_wide), "Tumor_Sample_Barcode"))
    set(mut_wide, which(is.na(mut_wide[[col]])), col, 0L)
  mut_wide[, sample_id := sub("F\\d+$", "", Tumor_Sample_Barcode)]
  df <- merge(
    as.data.frame(mut_wide),
    as.data.frame(clin[, .(sample_id, TMB_log, MSI_H, MSI, stage, histo)]),
    by = "sample_id", all = FALSE
  )
  df <- df[complete.cases(df[, c("TMB_log","MSI_H","stage","histo")]), ]
  if (strata == "MSS") df <- df[!is.na(df$MSI_H) & df$MSI_H == 0, ]
  cat(sprintf("    model_df [%s]: %d samples\n", strata, nrow(df)))
  df
}

# ── Single logistic model: gene_a ~ gene_b + covariates ────────────────────────
run_one_glm <- function(df, gene_a, gene_b, fmla_tpl) {
  if (!all(c(gene_a, gene_b) %in% names(df))) return(NULL)
  sub  <- df[!is.na(df[[gene_a]]) & !is.na(df[[gene_b]]), ]
  if (nrow(sub) < MIN_N || var(sub[[gene_a]]) == 0) return(NULL)
  fmla <- as.formula(
    gsub("__B__", paste0("`", gene_b, "`"),
    gsub("__A__", paste0("`", gene_a, "`"), fmla_tpl))
  )
  fit <- tryCatch(suppressWarnings(glm(fmla, sub, family = binomial)), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  ct  <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
  if (is.null(ct)) return(NULL)
  key <- if (paste0("`", gene_b, "`") %in% rownames(ct)) paste0("`", gene_b, "`") else gene_b
  if (!key %in% rownames(ct)) return(NULL)
  b  <- ct[key, ]
  ci <- tryCatch(confint.default(fit)[key, ], error = function(e) c(NA_real_, NA_real_))
  data.frame(
    gene_a = gene_a, gene_b = gene_b, n = nrow(sub),
    beta   = b[1], se = b[2], z = b[3], pvalue = b[4],
    OR     = exp(b[1]), OR_lo95 = exp(ci[1]), OR_hi95 = exp(ci[2]),
    stringsAsFactors = FALSE, row.names = NULL
  )
}

# ── All pairwise (both directions) for a gene set ─────────────────────────────
run_pairwise <- function(df, genes, strata_label) {
  fmla_tpl <- if (strata_label == "MSS")
    "__A__ ~ __B__ + TMB_log + stage + histo" else
    "__A__ ~ __B__ + TMB_log + MSI_H + stage + histo"
  cat(sprintf("    Pairwise logistic [%s, %d genes × %d genes] ...\n",
              strata_label, length(genes), length(genes)))
  pairs <- expand.grid(gene_a = genes, gene_b = genes, stringsAsFactors = FALSE)
  pairs <- pairs[pairs$gene_a != pairs$gene_b, ]
  res   <- do.call(rbind, lapply(seq_len(nrow(pairs)),
    function(i) run_one_glm(df, pairs$gene_a[i], pairs$gene_b[i], fmla_tpl)
  ))
  if (!is.null(res) && nrow(res) > 0) {
    res$FDR   <- p.adjust(res$pvalue, method = "BH")
    res$event <- ifelse(res$beta < 0, "Mutually_Exclusive", "Co_Occurrence")
  }
  res
}

# =============================================================================
# 3. FIGURE FUNCTIONS
# =============================================================================

# ── Fig 1: OncoPrint ──────────────────────────────────────────────────────────
make_fig1 <- function(maf_obj, cohort_label, n_samp) {
  vc_cols <- c(
    Missense_Mutation      = "#3A9D47",
    Nonsense_Mutation      = "#C0392B",
    Frame_Shift_Del        = "#8E44AD",
    Frame_Shift_Ins        = "#E67E22",
    Splice_Site            = "#2980B9",
    In_Frame_Del           = "#F39C12",
    In_Frame_Ins           = "#1ABC9C",
    Nonstop_Mutation       = "#E74C3C",
    Translation_Start_Site = "#95A5A6",
    Multi_Hit              = "#2C3E50"
  )
  ann_cols <- list(
    MSI   = c("MSI-H" = "#1F78B4",  "MSS"  = "#EEEEEE"),
    EBV   = c("EBV+"  = "#33A02C",  "EBV-" = "#EEEEEE"),
    HER2  = c("HER2+" = "#E31A1C",  "HER2-"= "#EEEEEE"),
    Stage = c("I" = "#FEE5D9", "II" = "#FCAE91", "III" = "#FB6A4A",
              "IV" = "#CB181D", "Unknown" = "#DDDDDD"),
    Sex   = c("男" = "#4393C3", "女" = "#D6604D")   # 男/女
  )
  out_png <- file.path(OUT, "figures", sprintf("fig1_oncoprint_%s.png", cohort_label))
  png(out_png, width = 18, height = 11, units = "in", res = 200)
  tryCatch(
    oncoplot(
      maf                    = maf_obj,
      top                    = TOP_N,
      colors                 = vc_cols,
      clinicalFeatures       = c("MSI", "EBV", "HER2", "Stage", "Sex"),
      annotationColor        = ann_cols,
      sortByAnnotation       = TRUE,
      showTumorSampleBarcodes= FALSE,
      drawRowBar             = TRUE,
      drawColBar             = TRUE,
      bgCol                  = "#FFFFFF",
      borderCol              = NA,
      titleText              = sprintf(
        "FDU Cohort — Somatic Landscape  (N = %d, top %d genes, %s)",
        n_samp, TOP_N, cohort_label),
      legendFontSize         = 1.0,
      fontSize               = 0.75
    ),
    error = function(e) cat("  [WARN oncoplot]", e$message, "\n")
  )
  dev.off()
  cat(sprintf("  Fig 1 saved: %s\n", basename(out_png)))
}

# ── Fig 2: Pairwise interaction heatmap ───────────────────────────────────────
make_fig2 <- function(lr_res, top_genes, cohort_label, strata_label, n_samp) {
  if (is.null(lr_res) || nrow(lr_res) == 0) return(invisible(NULL))

  # Collapse to undirected: keep min-p entry per unordered pair
  lr_res$pair_key <- mapply(
    function(a, b) paste(sort(c(a, b)), collapse = "::"),
    lr_res$gene_a, lr_res$gene_b
  )
  lr_u <- lr_res[order(lr_res$pvalue), ]
  lr_u <- lr_u[!duplicated(lr_u$pair_key), ]

  # Build symmetric plot data
  sym <- rbind(
    lr_u[, c("gene_a","gene_b","OR","pvalue","FDR","event")],
    setNames(lr_u[, c("gene_b","gene_a","OR","pvalue","FDR","event")],
             c("gene_a","gene_b","OR","pvalue","FDR","event"))
  )
  sym <- sym[sym$gene_a != sym$gene_b, ]

  # Gene order: TP53 first, then by overall mutation frequency
  freq_ord   <- names(sort(table(c(lr_res$gene_a, lr_res$gene_b)), decreasing = TRUE))
  gene_order <- c("TP53", setdiff(freq_ord, "TP53"))
  gene_order <- gene_order[gene_order %in% top_genes]
  n_g        <- length(gene_order)

  sym$log2OR  <- pmin(pmax(log2(pmax(sym$OR, 1e-4)), -4), 4)
  sym$sig_p   <- sym$pvalue < 0.05
  sym$sig_fdr <- sym$FDR    < 0.05
  sym$gene_a  <- factor(sym$gene_a, levels = gene_order)
  sym$gene_b  <- factor(sym$gene_b, levels = rev(gene_order))

  tp53_x <- which(gene_order == "TP53")
  tp53_y <- which(rev(gene_order) == "TP53")
  sz     <- max(11, n_g * 0.25)

  model_note <- if (strata_label == "MSS")
    "gene_A ~ gene_B + TMB + stage + histo" else
    "gene_A ~ gene_B + TMB + MSI + stage + histo"

  p <- ggplot(sym, aes(gene_a, gene_b)) +
    geom_tile(aes(fill = log2OR), color = "white", linewidth = 0.12) +
    # FDR<0.05: solid star
    geom_point(data = subset(sym, sig_fdr),
               shape = 8, size = 1.8, stroke = 0.8, color = "black") +
    # p<0.05 only: small open circle
    geom_point(data = subset(sym, sig_p & !sig_fdr),
               shape = 1, size = 1.3, stroke = 0.5, color = "grey25") +
    # TP53 column highlight
    {if (length(tp53_x) > 0) list(
      annotate("rect",
               xmin = tp53_x - 0.5, xmax = tp53_x + 0.5,
               ymin = 0.5, ymax = n_g + 0.5,
               fill = NA, color = ME_COL, linewidth = 1.0),
      annotate("rect",
               xmin = 0.5, xmax = n_g + 0.5,
               ymin = tp53_y - 0.5, ymax = tp53_y + 0.5,
               fill = NA, color = ME_COL, linewidth = 1.0)
    ) else NULL} +
    scale_fill_gradient2(
      low = ME_COL, mid = "white", high = CO_COL, midpoint = 0,
      limits = c(-4, 4), oob = squish, na.value = "grey92",
      name  = "log₂(adj.OR)",
      breaks = c(-4, -2, 0, 2, 4),
      labels = c("≤−4\n(ME)", "−2", "0", "+2", "≥4\n(Co-occ)")
    ) +
    PUB +
    theme(
      axis.text.x       = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
      axis.text.y       = element_text(size = 6, color = "black"),
      axis.ticks        = element_blank(),
      axis.line         = element_blank(),
      legend.key.height = unit(1.2, "cm"),
      legend.key.width  = unit(0.4, "cm"),
      plot.margin       = margin(5, 5, 5, 5, "mm")
    ) +
    labs(
      x = NULL, y = NULL,
      title    = sprintf("Pairwise somatic interactions — FDU %s [%s]  N=%d",
                         cohort_label, strata_label, n_samp),
      subtitle = sprintf("Logistic: %s  |  ✶ FDR<0.05  ○ p<0.05", model_note)
    )

  out_png <- file.path(OUT, "figures",
    sprintf("fig2_heatmap_%s_%s.png", cohort_label, strata_label))
  ggsave(out_png, p, width = sz, height = sz * 0.9, dpi = 250, bg = "white")
  cat(sprintf("  Fig 2 saved: %s\n", basename(out_png)))
}

# ── Fig 3: TP53-ME genes — mutation rates + adjusted ORs ──────────────────────
make_fig3 <- function(lr_res, model_df, cohort_label, strata_label) {
  if (is.null(lr_res) || nrow(lr_res) == 0) return(invisible(NULL))
  if (!"TP53" %in% names(model_df)) return(invisible(NULL))

  # partner ~ TP53 direction, ME, p < 0.05
  me <- lr_res[lr_res$gene_b == "TP53" &
               lr_res$event  == "Mutually_Exclusive" &
               lr_res$pvalue < 0.05, ]
  if (nrow(me) == 0) {
    cat(sprintf("  Fig 3: no TP53-ME pairs [%s %s]\n", cohort_label, strata_label))
    return(invisible(NULL))
  }
  me       <- me[order(me$OR), ]               # ascending: most ME first
  partners <- me$gene_a[me$gene_a %in% names(model_df)]
  if (length(partners) == 0) return(invisible(NULL))

  # Mutation rates TP53+ vs TP53−
  rate_list <- lapply(partners, function(g) {
    sub <- model_df[!is.na(model_df$TP53) & !is.na(model_df[[g]]), ]
    data.frame(
      gene  = g,
      group = c("TP53+", "TP53−"),
      n     = c(sum(sub$TP53 == 1), sum(sub$TP53 == 0)),
      n_mut = c(sum(sub$TP53 == 1 & sub[[g]] == 1),
                sum(sub$TP53 == 0 & sub[[g]] == 1)),
      stringsAsFactors = FALSE
    )
  })
  rate_df      <- do.call(rbind, rate_list)
  rate_df$pct  <- rate_df$n_mut / rate_df$n * 100

  # Fisher p-value per partner
  fish_list <- lapply(partners, function(g) {
    sub  <- model_df[!is.na(model_df$TP53) & !is.na(model_df[[g]]), ]
    tbl  <- table(TP53 = sub$TP53, Gene = sub[[g]])
    pval <- tryCatch(fisher.test(tbl)$p.value, error = function(e) NA_real_)
    pmax_val <- max(rate_df$pct[rate_df$gene == g], na.rm = TRUE)
    data.frame(gene = g, fisher_p = pval, pct_max = pmax_val,
               stringsAsFactors = FALSE)
  })
  fish_df      <- do.call(rbind, fish_list)
  fish_df$star <- with(fish_df,
    ifelse(is.na(fisher_p), "",
    ifelse(fisher_p < 0.001, "***",
    ifelse(fisher_p < 0.01,  "**",
    ifelse(fisher_p < 0.05,  "*", "ns")))))

  gene_lvl      <- partners        # already sorted by ascending OR
  rate_df$gene  <- factor(rate_df$gene, levels = gene_lvl)
  fish_df$gene  <- factor(fish_df$gene, levels = gene_lvl)

  me_plot       <- me[me$gene_a %in% gene_lvl, ]
  me_plot$gene  <- factor(me_plot$gene_a, levels = gene_lvl)
  me_plot$sig   <- factor(
    ifelse(me_plot$FDR   < 0.05, "FDR<0.05",
    ifelse(me_plot$pvalue< 0.05, "p<0.05",   "n.s.")),
    levels = c("FDR<0.05","p<0.05","n.s.")
  )

  # Dynamic x-axis limits for OR panel
  or_min <- max(0.01, min(me_plot$OR_lo95, na.rm = TRUE) * 0.6)
  or_max <- min(5,    max(me_plot$OR_hi95, na.rm = TRUE) * 1.5)
  or_max <- max(or_max, 1.5)

  # ── Panel A: mutation rate bar chart ────────────────────────────────────────
  pa <- ggplot(rate_df, aes(x = gene, y = pct, fill = group)) +
    geom_col(position = position_dodge(width = 0.72), width = 0.68, color = NA) +
    geom_text(data = fish_df,
              aes(x = gene, y = pct_max + 1.5, label = star),
              inherit.aes = FALSE, size = 4.5, fontface = "bold", color = "grey20") +
    scale_fill_manual(
      values = c("TP53+" = ME_COL, "TP53−" = CO_COL), name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)),
                       labels = function(x) paste0(x, "%")) +
    PUB +
    theme(
      axis.text.x      = element_text(angle = 35, hjust = 1, face = "bold", size = 10),
      legend.position  = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.background= element_rect(color = "grey70", fill = "white", linewidth = 0.3)
    ) +
    labs(x = NULL, y = "Mutation rate (%)", title = "Mutation rate by TP53 status")

  # ── Panel B: adjusted OR + 95%CI (horizontal dot-CI plot) ───────────────────
  pb <- ggplot(me_plot, aes(x = OR, y = gene, color = sig)) +
    geom_vline(xintercept = 1, linetype = "dashed",
               color = "grey55", linewidth = 0.7) +
    geom_errorbar(aes(xmin = OR_lo95, xmax = OR_hi95),
                  width = 0.28, linewidth = 1.1, orientation = "y") +
    geom_point(size = 4.5) +
    geom_text(
      aes(label = sprintf("OR=%.2f [%.2f–%.2f]\nFDR=%.3f",
                          OR, OR_lo95, OR_hi95, FDR)),
      hjust = -0.1, size = 2.5, color = "grey20", lineheight = 1.15
    ) +
    scale_color_manual(
      values = c("FDR<0.05" = ME_COL, "p<0.05" = "#E67E22", "n.s." = "#7F8C8D"),
      name = NULL, drop = FALSE
    ) +
    scale_x_log10(
      limits = c(or_min, or_max),
      breaks = c(0.05, 0.1, 0.2, 0.5, 1, 2),
      labels = c("0.05","0.1","0.2","0.5","1","2")
    ) +
    PUB +
    theme(
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "bottom"
    ) +
    labs(
      x     = "Adj. OR (log scale) | partner ~ TP53 + TMB + MSI + stage + histo",
      y     = NULL,
      title = "Mutual exclusivity with TP53"
    )

  combined <- pa + pb +
    plot_layout(ncol = 2, widths = c(1.1, 1.6)) +
    plot_annotation(
      title    = sprintf("TP53-mutually exclusive genes — FDU %s [%s]",
                         cohort_label, strata_label),
      subtitle = sprintf("%d partner genes with p < 0.05", length(gene_lvl)),
      theme    = theme(
        plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9, color = "grey40")
      )
    )

  out_png <- file.path(OUT, "figures",
    sprintf("fig3_tp53_me_%s_%s.png", cohort_label, strata_label))
  ggsave(out_png, combined,
         width = 13, height = max(4.5, 2.5 + length(gene_lvl) * 0.6),
         dpi = 250, bg = "white")
  cat(sprintf("  Fig 3 saved: %s\n", basename(out_png)))
}

# =============================================================================
# 4. MAIN ANALYSIS LOOP
# =============================================================================

cohorts <- list(
  CohortA_TxNaive = intersect(unique(maf_raw$sample_base), naive_ids),
  CohortB_All     = unique(maf_raw$sample_base)
)

for (cname in names(cohorts)) {
  keep <- cohorts[[cname]]
  cat(sprintf(
    "\n%s\n%s  (N = %d)\n%s\n",
    strrep("=", 70), cname, length(keep), strrep("=", 70)
  ))

  maf_obj <- build_maf(keep)
  n_samp  <- nrow(getSampleSummary(maf_obj))
  gs      <- getGeneSummary(maf_obj)
  top50   <- head(gs$Hugo_Symbol[order(-gs$MutatedSamples)], TOP_N)
  cat(sprintf("  Top-5: %s\n", paste(head(top50, 5), collapse = ", ")))

  # ── Fisher's exact test (maftools somaticInteractions) ───────────────────────
  cat("  Fisher's exact test (somaticInteractions) ...\n")
  out_fisher_png <- file.path(OUT, "figures",
    sprintf("figS_fisher_%s.png", cname))
  png(out_fisher_png, width = 16, height = 16, units = "in", res = 150)
  fisher_res <- tryCatch(
    somaticInteractions(maf_obj, top = TOP_N,
                        pvalue = c(0.05, 0.01),
                        returnAll = TRUE,
                        showCounts = TRUE,
                        fontSize   = 0.65),
    error = function(e) { cat("  [WARN]", e$message, "\n"); NULL }
  )
  dev.off()
  if (!is.null(fisher_res)) {
    fwrite(as.data.table(fisher_res),
           file.path(OUT, "tables", sprintf("fisher_%s.tsv", cname)),
           sep = "\t")
    n_me_fisher <- sum(!is.na(fisher_res$pValue) &
                         fisher_res$pValue < 0.05 &
                         fisher_res$Event == "Mutually_Exclusive", na.rm = TRUE)
    cat(sprintf("  Fisher: %d ME pairs (p<0.05)\n", n_me_fisher))
  }

  # ── Fig 1: OncoPrint ─────────────────────────────────────────────────────────
  make_fig1(maf_obj, cname, n_samp)

  # ── Full cohort: logistic regression + figs ───────────────────────────────────
  cat("  Full cohort analysis ...\n")
  df_full <- build_model_df(maf_obj, "all")
  lr_full <- run_pairwise(df_full, top50, "all")
  if (!is.null(lr_full) && nrow(lr_full) > 0) {
    fwrite(lr_full,
           file.path(OUT, "tables", sprintf("logistic_full_%s.tsv", cname)),
           sep = "\t")
    cat(sprintf("  Logistic [full]: %d ME pairs (p<0.05), %d (FDR<0.05)\n",
        sum(lr_full$event == "Mutually_Exclusive" & lr_full$pvalue < 0.05),
        sum(lr_full$event == "Mutually_Exclusive" & lr_full$FDR    < 0.05)))
    make_fig2(lr_full, top50, cname, "full", n_samp)
    make_fig3(lr_full, df_full, cname, "full")
  }

  # ── MSS-only: logistic regression + figs ─────────────────────────────────────
  cat("  MSS-only analysis ...\n")
  df_mss <- build_model_df(maf_obj, "MSS")
  lr_mss <- run_pairwise(df_mss, top50, "MSS")
  if (!is.null(lr_mss) && nrow(lr_mss) > 0) {
    fwrite(lr_mss,
           file.path(OUT, "tables", sprintf("logistic_MSS_%s.tsv", cname)),
           sep = "\t")
    n_mss_samp <- sum(!is.na(df_mss$MSI_H) & df_mss$MSI_H == 0)
    make_fig2(lr_mss, top50, cname, "MSS", n_mss_samp)
    make_fig3(lr_mss, df_mss, cname, "MSS")
  }
}

cat(sprintf(
  "\n%s\nDone\nTables  : %s/tables/\nFigures : %s/figures/\n%s\n",
  strrep("=", 70), OUT, OUT, strrep("=", 70)
))
