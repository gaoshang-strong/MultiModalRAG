#!/usr/bin/env Rscript
#
# mrna_pathway_DE.R — DESeq2 differential expression of ATM-TP53 pathway genes
#
# Two comparisons (ALL samples, no MSI-H exclusion):
#   Comp 1: TP53-mut / ATM-WT  vs  TP53-WT / ATM-WT
#   Comp 2: ATM-mut / TP53-WT  vs  ATM-WT  / TP53-WT  (skipped if N < 5)
#
# Multiple-testing correction: BH within each gene class (ATM-specific /
#   Convergence / TP53-specific) — NOT genome-wide FDR.
#
# Usage:
#   /home/sgao30/micromamba/bin/micromamba run -n ProjectGeneration \
#     Rscript TP53_ATM_paper/scripts/mrna_pathway_DE.R
#
# Results saved to:
#   TP53_ATM_paper/results/mRNA/comp1_TP53mut/
#   TP53_ATM_paper/results/mRNA/comp2_ATMmut/

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(fgsea)
  library(dplyr)
  library(RColorBrewer)
  library(data.table)
})

# ── Repo root ─────────────────────────────────────────────────────────────────
get_script_dir <- function() {
  args     <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg))
    dirname(normalizePath(sub("--file=", "", file_arg)))
  else
    normalizePath(".")
}
SCRIPT_DIR <- get_script_dir()
REPO_ROOT  <- normalizePath(file.path(SCRIPT_DIR, "../.."))

RNASEQ_RDS  <- file.path(REPO_ROOT, "TCGA-STAD/rnaseq/rnaseq_summarized_experiment.rds")
SOMATIC_CSV <- file.path(REPO_ROOT, "TCGA-STAD/somatic_mutation/somatic_mutation.csv")
OUTDIR_C1   <- file.path(REPO_ROOT, "TP53_ATM_paper/results/mRNA/comp1_TP53mut")
OUTDIR_C2   <- file.path(REPO_ROOT, "TP53_ATM_paper/results/mRNA/comp2_ATMmut")
dir.create(OUTDIR_C1, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTDIR_C2, recursive = TRUE, showWarnings = FALSE)


# ── Gene sets: three pathway classes from figure ──────────────────────────────
GENES_ATM  <- c("MRE11","NBN","RAD50","ATM","H2AX","MDC1","TRIM28",
                "RNF8","RNF168","TP53BP1","PARP1")
GENES_CONV <- c("KAT5","ABRAXAS1","CHEK2","ATR","CHEK1","BRCA1","FOXO3",
                "PALB2","BRCA2","RAD51","WEE1","CDK1")
GENES_TP53 <- c("TP53","MDM2","CDKN2A","CDKN1A","SFN","GADD45A","BAX",
                "BBC3","PMAIP1","DDB2","TIGAR","PPM1D")
ALL_PATHWAY <- c(GENES_ATM, GENES_CONV, GENES_TP53)

CLASS_ORDER  <- c("ATM-specific", "Convergence", "TP53-specific")
CLASS_COLORS <- c("ATM-specific" = "#1B6BA8",
                  "Convergence"  = "#6A0DAD",
                  "TP53-specific"= "#C0392B")
GENE_CLASS <- c(
  setNames(rep("ATM-specific",  length(GENES_ATM)),  GENES_ATM),
  setNames(rep("Convergence",   length(GENES_CONV)), GENES_CONV),
  setNames(rep("TP53-specific", length(GENES_TP53)), GENES_TP53)
)


# ── Mutation definitions ──────────────────────────────────────────────────────
FUN_CLASSES <- c(
  "Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
  "Frame_Shift_Ins",  "Splice_Site",        "In_Frame_Del",
  "In_Frame_Ins",     "Nonstop_Mutation",   "Translation_Start_Site"
)


# ── Load data ─────────────────────────────────────────────────────────────────
message("Loading RNA-seq SummarizedExperiment...")
se <- readRDS(RNASEQ_RDS)
message(sprintf("  SE: %d genes x %d samples", nrow(se), ncol(se)))

message("Loading somatic mutations...")
maf <- fread(SOMATIC_CSV, data.table = FALSE)

tp53_pts <- maf %>%
  filter(Hugo_Symbol == "TP53",
         Variant_Classification %in% FUN_CLASSES) %>%
  mutate(patient_id = substr(Tumor_Sample_Barcode, 1, 12)) %>%
  distinct(patient_id) %>%
  pull(patient_id)

atm_pts <- maf %>%
  filter(Hugo_Symbol == "ATM",
         Variant_Classification %in% FUN_CLASSES) %>%
  mutate(patient_id = substr(Tumor_Sample_Barcode, 1, 12)) %>%
  distinct(patient_id) %>%
  pull(patient_id)

message(sprintf("  TP53-mut (any functional): N=%d  |  ATM-mut: N=%d",
                length(tp53_pts), length(atm_pts)))


# ── Sample table (primary tumor, one per patient) ─────────────────────────────
barcodes    <- colnames(se)
patient_ids <- substr(barcodes, 1, 12)
sample_type <- substr(barcodes, 14, 15)

sample_tbl <- data.frame(
  barcode     = barcodes,
  patient_id  = patient_ids,
  sample_type = sample_type,
  stringsAsFactors = FALSE
) %>%
  filter(sample_type == "01") %>%
  distinct(patient_id, .keep_all = TRUE) %>%
  mutate(
    tp53_status = ifelse(patient_id %in% tp53_pts, "mutant", "wildtype"),
    atm_status  = ifelse(patient_id %in% atm_pts,  "mutant", "wildtype")
  )

message(sprintf("Primary tumor samples: N=%d  TP53-mut=%d  ATM-mut=%d",
  nrow(sample_tbl),
  sum(sample_tbl$tp53_status == "mutant"),
  sum(sample_tbl$atm_status  == "mutant")
))


# ── Helper: TPM matrix ────────────────────────────────────────────────────────
get_tpm_matrix <- function(se_sub) {
  mat        <- assay(se_sub, "tpm_unstrand")
  gene_names <- rowData(se_sub)$gene_name
  ord        <- order(-rowMeans(mat))
  mat        <- mat[ord, ]
  gene_names <- gene_names[ord]
  keep       <- !duplicated(gene_names)
  mat        <- mat[keep, ]
  rownames(mat) <- gene_names[keep]
  mat
}


# ── Helper: extract pathway genes + within-class BH ──────────────────────────
extract_pathway_results <- function(deseq_res) {
  df <- deseq_res[!is.na(deseq_res$gene_name), ]
  df <- df[order(df$pvalue, na.last = TRUE), ]
  df <- df[!duplicated(df$gene_name), ]

  pw       <- df[df$gene_name %in% ALL_PATHWAY, , drop = FALSE]
  pw$class <- GENE_CLASS[pw$gene_name]

  pw$padj_class <- NA_real_
  for (cl in CLASS_ORDER) {
    idx <- which(pw$class == cl)
    if (length(idx) > 0 && !all(is.na(pw$pvalue[idx])))
      pw$padj_class[idx] <- p.adjust(pw$pvalue[idx], method = "BH")
  }

  pw$class <- factor(pw$class, levels = CLASS_ORDER)
  pw       <- pw[order(pw$class, pw$pvalue, na.last = TRUE), ]

  missing <- setdiff(ALL_PATHWAY, pw$gene_name)
  if (length(missing) > 0)
    message(sprintf("  Not in DESeq2 results: %s", paste(missing, collapse = ", ")))

  pw
}


# ── Helper: GSEA per class ────────────────────────────────────────────────────
run_gsea_3class <- function(deseq_res) {
  df <- deseq_res[!is.na(deseq_res$padj) &
                  !is.na(deseq_res$log2FoldChange) &
                  !is.na(deseq_res$gene_name), ]
  df <- df[order(df$padj), ]
  df <- df[!duplicated(df$gene_name), ]
  r  <- sign(df$log2FoldChange) * -log10(pmax(df$padj, 1e-300))
  names(r) <- df$gene_name
  ranks    <- sort(r, decreasing = TRUE)

  pathways <- list(
    ATM_specific  = GENES_ATM,
    Convergence   = GENES_CONV,
    TP53_specific = GENES_TP53,
    All_pathway   = ALL_PATHWAY
  )
  pathways <- pathways[sapply(pathways, function(g) sum(g %in% names(ranks)) >= 5)]
  if (length(pathways) == 0) { warning("No pathway ≥5 genes in ranks"); return(NULL) }

  set.seed(42)
  fgsea(pathways = pathways, stats = ranks,
        minSize = 5, maxSize = 1000, nPermSimple = 10000)
}


# ── Plots ─────────────────────────────────────────────────────────────────────
plot_lollipop <- function(pw, comp_label, outfile) {
  df <- pw[!is.na(pw$log2FoldChange), ]
  # Order within each class by ascending log2FC (most negative first)
  df <- df[order(df$class, df$log2FoldChange), ]
  df$gene_name <- factor(df$gene_name, levels = df$gene_name)
  df$sig <- with(df, case_when(
    !is.na(padj_class) & padj_class < 0.001 ~ "***",
    !is.na(padj_class) & padj_class < 0.01  ~ "**",
    !is.na(padj_class) & padj_class < 0.05  ~ "*",
    TRUE ~ ""
  ))

  p <- ggplot(df, aes(x = log2FoldChange, y = gene_name, color = class)) +
    geom_vline(xintercept = 0, linewidth = 0.5, color = "grey60", linetype = "dashed") +
    geom_segment(aes(x = 0, xend = log2FoldChange, yend = gene_name),
                 linewidth = 0.8, alpha = 0.7) +
    geom_point(aes(size = -log10(pmax(pvalue, 1e-10))), alpha = 0.9) +
    geom_text(aes(x = log2FoldChange + sign(log2FoldChange + 1e-9) * 0.08,
                  label = sig),
              size = 4, hjust = 0.5, show.legend = FALSE) +
    scale_color_manual(values = CLASS_COLORS, name = "Gene class") +
    scale_size_continuous(name = "-log10(p-value)", range = c(2, 8)) +
    facet_grid(rows = vars(class), scales = "free_y", space = "free_y") +
    labs(
      title    = comp_label,
      subtitle = "Point size = -log10(p);  * = padj < 0.05 within class (BH)",
      x        = "log2 Fold Change  (mutant / wildtype)",
      y        = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      strip.background = element_rect(fill = "grey92"),
      strip.text       = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )

  h <- max(6, nrow(df) * 0.45 + 3)
  ggsave(outfile, p, width = 9, height = h, dpi = 150)
  message(sprintf("  Saved: %s", outfile))
}


plot_volcano_pathway <- function(deseq_res, pw, comp_label, outfile) {
  df <- deseq_res[!is.na(deseq_res$padj) & !is.na(deseq_res$gene_name), ]
  df$neg_log10_padj <- -log10(pmax(df$padj, 1e-50))
  df$class  <- ifelse(df$gene_name %in% ALL_PATHWAY,
                      GENE_CLASS[df$gene_name], "Other")
  df$class  <- factor(df$class, levels = c(CLASS_ORDER, "Other"))
  df$label  <- ifelse(df$gene_name %in% ALL_PATHWAY, df$gene_name, NA_character_)

  p <- ggplot(df[df$class == "Other", ],
              aes(log2FoldChange, neg_log10_padj)) +
    geom_point(color = "#CCCCCC", alpha = 0.25, size = 0.5) +
    geom_point(data = df[df$class != "Other", ],
               aes(color = class), size = 2.8, alpha = 0.85) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = "grey50", linewidth = 0.4) +
    geom_text_repel(data = df[df$class != "Other", ],
                    aes(label = label, color = class),
                    size = 3.1, max.overlaps = 40, show.legend = FALSE,
                    min.segment.length = 0.2) +
    scale_color_manual(values = CLASS_COLORS, name = "Gene class",
                       breaks = CLASS_ORDER, drop = FALSE) +
    labs(
      title    = comp_label,
      subtitle = "y-axis = genome-wide padj (display only); class correction shown in lollipop",
      x        = "log2 Fold Change  (mutant / wildtype)",
      y        = "-log10(padj)  [genome-wide, for context]"
    ) +
    theme_bw(base_size = 12) +
    guides(color = guide_legend(override.aes = list(size = 3.5)))

  ggsave(outfile, p, width = 9, height = 6, dpi = 150)
  message(sprintf("  Saved: %s", outfile))
}


plot_heatmap_3class <- function(tpm_mat, pw, status_vec, var_label, comp_label, outfile) {
  genes <- intersect(as.character(pw$gene_name), rownames(tpm_mat))
  if (length(genes) == 0) return(invisible(NULL))

  # Preserve class-then-log2FC ordering from pw
  gene_order <- intersect(as.character(pw$gene_name), genes)
  mat   <- log2(tpm_mat[gene_order, , drop = FALSE] + 1)
  mat_z <- t(scale(t(mat)))
  mat_z <- pmin(pmax(mat_z, -2.5), 2.5)

  col_order <- order(status_vec == "wildtype", status_vec)
  mat_z     <- mat_z[, col_order, drop = FALSE]

  annot_col <- data.frame(
    Group = factor(status_vec[col_order], levels = c("mutant", "wildtype")),
    row.names = colnames(mat_z)
  )
  colnames(annot_col) <- var_label
  annot_colors_col    <- setNames(
    list(c(mutant = "#E41A1C", wildtype = "#377EB8")),
    var_label
  )

  row_class   <- GENE_CLASS[gene_order]
  annot_row   <- data.frame(Class = factor(row_class, levels = CLASS_ORDER),
                            row.names = gene_order)
  annot_colors_row <- list(Class = CLASS_COLORS)

  h <- max(5, length(genes) * 0.38 + 3)
  png(outfile, width = 12, height = h, units = "in", res = 150)
  pheatmap(
    mat_z,
    annotation_col    = annot_col,
    annotation_row    = annot_row,
    annotation_colors = c(annot_colors_col, annot_colors_row),
    cluster_rows      = FALSE,
    cluster_cols      = FALSE,
    show_colnames     = FALSE,
    color             = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
    breaks            = seq(-2.5, 2.5, length.out = 101),
    main              = comp_label,
    fontsize_row      = 9,
    gaps_row          = cumsum(table(factor(row_class, levels = CLASS_ORDER)))[-3]
  )
  dev.off()
  message(sprintf("  Saved: %s", outfile))
}


plot_boxplots_3class <- function(tpm_mat, pw, status_vec, comp_label, outfile) {
  genes <- intersect(as.character(pw$gene_name), rownames(tpm_mat))
  if (length(genes) == 0) return(invisible(NULL))

  long_df <- lapply(genes, function(g) {
    r   <- pw[as.character(pw$gene_name) == g, ]
    sig <- with(r, case_when(
      !is.na(padj_class) & padj_class < 0.001 ~ "***",
      !is.na(padj_class) & padj_class < 0.01  ~ "**",
      !is.na(padj_class) & padj_class < 0.05  ~ "*",
      TRUE ~ "ns"
    ))
    lfc <- if (!is.na(r$log2FoldChange)) sprintf("%.2f", r$log2FoldChange) else "NA"
    data.frame(
      gene       = g,
      expr       = log2(tpm_mat[g, ] + 1),
      status     = status_vec,
      class      = GENE_CLASS[g],
      gene_label = paste0(g, "\n(FC=", lfc, " ", sig, ")"),
      stringsAsFactors = FALSE
    )
  })
  long_df <- do.call(rbind, long_df)

  # Order facets by class then log2FC
  gene_lbl_order <- unique(long_df$gene_label[
    order(match(long_df$class, CLASS_ORDER), long_df$gene)
  ])
  long_df$gene_label <- factor(long_df$gene_label, levels = gene_lbl_order)

  ncols <- min(6, length(genes))
  nrows <- ceiling(length(genes) / ncols)

  p <- ggplot(long_df, aes(status, expr, fill = status)) +
    geom_boxplot(outlier.size = 0.4, width = 0.6, linewidth = 0.4) +
    scale_fill_manual(values = c(mutant = "#E41A1C", wildtype = "#377EB8")) +
    facet_wrap(~ gene_label, ncol = ncols, scales = "free_y") +
    labs(title = comp_label, x = NULL, y = "log2(TPM + 1)") +
    theme_bw(base_size = 10) +
    theme(legend.position  = "none",
          strip.text       = element_text(face = "bold", size = 8),
          axis.text.x      = element_text(angle = 30, hjust = 1))

  ggsave(outfile, p, width = ncols * 3.2, height = nrows * 3.2, dpi = 150)
  message(sprintf("  Saved: %s", outfile))
}


plot_gsea_enrichment <- function(gsea_res, deseq_res, comp_label, outfile) {
  if (is.null(gsea_res)) return(invisible(NULL))

  df <- deseq_res[!is.na(deseq_res$padj) &
                  !is.na(deseq_res$log2FoldChange) &
                  !is.na(deseq_res$gene_name), ]
  df <- df[order(df$padj), ]; df <- df[!duplicated(df$gene_name), ]
  r  <- sign(df$log2FoldChange) * -log10(pmax(df$padj, 1e-300))
  names(r) <- df$gene_name
  ranks    <- sort(r, decreasing = TRUE)

  pathways_plot <- list(
    ATM_specific  = GENES_ATM,
    Convergence   = GENES_CONV,
    TP53_specific = GENES_TP53
  )

  plots <- lapply(names(pathways_plot), function(nm) {
    g <- pathways_plot[[nm]]
    if (sum(g %in% names(ranks)) < 5) return(NULL)
    row_r <- gsea_res[gsea_res$pathway == nm, ]
    if (nrow(row_r) == 0) return(NULL)
    plotEnrichment(g, ranks) +
      labs(title = nm,
           subtitle = sprintf("NES=%.2f  p=%.3f  padj=%.3f",
             row_r$NES, row_r$pval, row_r$padj)) +
      theme_bw(base_size = 10)
  })
  plots <- Filter(Negate(is.null), plots)
  if (length(plots) == 0) return(invisible(NULL))

  # Arrange side by side
  library(patchwork)
  combined <- Reduce(`+`, plots) + plot_annotation(title = comp_label)
  ggsave(outfile, combined, width = length(plots) * 5, height = 4, dpi = 150)
  message(sprintf("  Saved: %s", outfile))
}


# ── Write results ──────────────────────────────────────────────────────────────
write_results <- function(pw, gsea_res, outdir) {
  write.csv(pw, file.path(outdir, "pathway_DE_results.csv"), row.names = FALSE)

  if (!is.null(gsea_res)) {
    g <- as.data.frame(gsea_res)
    g$leadingEdge <- sapply(g$leadingEdge, paste, collapse = ";")
    write.csv(g, file.path(outdir, "gsea_pathway_summary.csv"), row.names = FALSE)
  }

  summary_tbl <- pw %>%
    group_by(class) %>%
    summarise(
      n_genes        = n(),
      n_sig_padj05   = sum(!is.na(padj_class) & padj_class < 0.05),
      n_up           = sum(!is.na(log2FoldChange) & log2FoldChange > 0),
      n_down         = sum(!is.na(log2FoldChange) & log2FoldChange < 0),
      mean_log2FC    = round(mean(log2FoldChange, na.rm = TRUE), 3),
      median_log2FC  = round(median(log2FoldChange, na.rm = TRUE), 3),
      .groups        = "drop"
    )
  write.csv(summary_tbl, file.path(outdir, "class_summary.csv"), row.names = FALSE)
  message("  class_summary:")
  print(as.data.frame(summary_tbl))
}


# ═══════════════════════════════════════════════════════════════════════════════
# run_comparison() — generic wrapper
# ═══════════════════════════════════════════════════════════════════════════════
run_comparison <- function(design_var, outdir, comp_label) {
  other_var  <- if (design_var == "tp53_status") "atm_status" else "tp53_status"
  var_label  <- if (design_var == "tp53_status") "TP53" else "ATM"

  tbl    <- sample_tbl[sample_tbl[[other_var]] == "wildtype", ]
  n_mut  <- sum(tbl[[design_var]] == "mutant")
  n_wt   <- sum(tbl[[design_var]] == "wildtype")
  message(sprintf("  %s-mut/WT-wt: N=%d  |  WT/WT: N=%d", var_label, n_mut, n_wt))

  if (n_mut < 5) {
    message(sprintf("  Mutant group too small (N=%d) — skipping", n_mut))
    return(invisible(NULL))
  }

  # Subset SE and annotate colData
  se_sub <- se[, tbl$barcode]
  cd     <- as.data.frame(colData(se_sub))
  idx    <- match(colnames(se_sub), tbl$barcode)
  cd[[design_var]] <- tbl[[design_var]][idx]
  cd[[other_var]]  <- tbl[[other_var]][idx]
  colData(se_sub)  <- DataFrame(cd)

  if ("gene_type" %in% colnames(rowData(se_sub)))
    se_sub <- se_sub[rowData(se_sub)$gene_type == "protein_coding", ]

  se_sub[[design_var]] <- factor(se_sub[[design_var]],
                                 levels = c("wildtype", "mutant"))

  # DESeq2
  message("  Running DESeq2...")
  dds  <- DESeqDataSet(se_sub, design = as.formula(paste0("~ ", design_var)))
  keep <- rowSums(counts(dds) >= 10) >= 5
  dds  <- dds[keep, ]
  message(sprintf("  %d genes after count filter", sum(keep)))
  dds <- DESeq(dds, quiet = TRUE)

  res    <- results(dds,
    contrast             = c(design_var, "mutant", "wildtype"),
    independentFiltering = TRUE)
  res_df <- as.data.frame(res)
  res_df$gene_name <- rowData(dds)$gene_name[match(rownames(res_df), rownames(dds))]
  res_df <- res_df[order(res_df$pvalue, na.last = TRUE), ]

  tpm_mat    <- get_tpm_matrix(se_sub)
  status_vec <- as.character(se_sub[[design_var]])

  # Pathway extraction + correction
  pw   <- extract_pathway_results(res_df)
  gsea <- run_gsea_3class(res_df)

  # Plots
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  plot_lollipop(pw, comp_label, file.path(outdir, "lollipop_pathway.png"))
  plot_volcano_pathway(res_df, pw, comp_label,
                       file.path(outdir, "volcano_pathway.png"))
  plot_heatmap_3class(tpm_mat, pw, status_vec, var_label, comp_label,
                      file.path(outdir, "heatmap_pathway.png"))
  plot_boxplots_3class(tpm_mat, pw, status_vec, comp_label,
                       file.path(outdir, "boxplots_pathway.png"))
  plot_gsea_enrichment(gsea, res_df, comp_label,
                       file.path(outdir, "gsea_enrichment.png"))
  write_results(pw, gsea, outdir)

  message(sprintf("  Results → %s\n", outdir))
  invisible(pw)
}


# ═══════════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════════
message("\n════════════════════════════════════════════════════")
message("Comparison 1: TP53-mut / ATM-WT  vs  TP53-WT / ATM-WT")
message("════════════════════════════════════════════════════")
run_comparison(
  design_var = "tp53_status",
  outdir     = OUTDIR_C1,
  comp_label = "Comp 1: TP53-mut / ATM-WT  vs  TP53-WT / ATM-WT  [all samples, N functional variants]"
)

message("\n════════════════════════════════════════════════════")
message("Comparison 2: ATM-mut / TP53-WT  vs  ATM-WT / TP53-WT")
message("════════════════════════════════════════════════════")
run_comparison(
  design_var = "atm_status",
  outdir     = OUTDIR_C2,
  comp_label = "Comp 2: ATM-mut / TP53-WT  vs  ATM-WT / TP53-WT  [all samples, N functional variants]"
)

message("\n✓  mrna_pathway_DE.R complete.")
