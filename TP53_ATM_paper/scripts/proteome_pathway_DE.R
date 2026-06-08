#!/usr/bin/env Rscript
#
# proteome_pathway_DE.R — Limma differential protein abundance, ATM-TP53 pathway
#
# Two comparisons per cohort (ALL samples, MSI as covariate):
#   Comp 1: TP53-mut / ATM-WT  vs  TP53-WT / ATM-WT
#   Comp 2: ATM-mut / TP53-WT  vs  ATM-WT  / TP53-WT  (skipped if N < 5)
#
# Cohorts: TCGA-STAD RPPA (206 proteins + phospho)  |  CPTAC proteome (9613 proteins)
# Correction: BH within each gene class only (not genome-wide)
#
# Usage:
#   micromamba run -n ProjectGeneration \
#     Rscript TP53_ATM_paper/scripts/proteome_pathway_DE.R

suppressPackageStartupMessages({
  library(limma)
  library(fgsea)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(dplyr)
  library(data.table)
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

RPPA_TSV    <- file.path(REPO_ROOT, "TCGA-STAD/RPPA/rppa_matrix.tsv")
TCGA_MAF    <- file.path(REPO_ROOT, "TCGA-STAD/somatic_mutation/somatic_mutation.csv")
CPTAC_PROT  <- Sys.glob(file.path(REPO_ROOT, "CPTAC/proteome/*.tmt18*"))[1]
CPTAC_MAF   <- file.path(REPO_ROOT, "CPTAC/somatic_mutation/cptac_gc_merged.maf.tsv")
CPTAC_IDMAP <- file.path(REPO_ROOT, "CPTAC/proteome/PDC_study_biospecimen_06062026_223409.csv")

OUT_BASE <- file.path(REPO_ROOT, "TP53_ATM_paper/results/Proteome")
for (d in c(file.path(OUT_BASE, "TCGA_RPPA/comp1_TP53mut"),
            file.path(OUT_BASE, "TCGA_RPPA/comp2_ATMmut"),
            file.path(OUT_BASE, "CPTAC/comp1_TP53mut"),
            file.path(OUT_BASE, "CPTAC/comp2_ATMmut")))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)


# ── Gene sets ─────────────────────────────────────────────────────────────────
GENES_ATM  <- c("MRE11","NBN","RAD50","ATM","H2AX","MDC1","TRIM28",
                "RNF8","RNF168","TP53BP1","PARP1")
GENES_CONV <- c("KAT5","ABRAXAS1","CHEK2","ATR","CHEK1","BRCA1","FOXO3",
                "PALB2","BRCA2","RAD51","WEE1","CDK1")
GENES_TP53 <- c("TP53","MDM2","CDKN2A","CDKN1A","SFN","GADD45A","BAX",
                "BBC3","PMAIP1","DDB2","TIGAR","PPM1D")
ALL_PATHWAY <- c(GENES_ATM, GENES_CONV, GENES_TP53)

CLASS_ORDER  <- c("ATM-specific", "Convergence", "TP53-specific")
CLASS_COLORS <- c("ATM-specific"  = "#1B6BA8",
                  "Convergence"   = "#6A0DAD",
                  "TP53-specific" = "#C0392B")
GENE_CLASS <- c(
  setNames(rep("ATM-specific",  length(GENES_ATM)),  GENES_ATM),
  setNames(rep("Convergence",   length(GENES_CONV)), GENES_CONV),
  setNames(rep("TP53-specific", length(GENES_TP53)), GENES_TP53)
)

# Mutation functional classes
FUN_CLASSES <- c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del",
  "Frame_Shift_Ins","Splice_Site","In_Frame_Del","In_Frame_Ins",
  "Nonstop_Mutation","Translation_Start_Site")
CODING_CLASSES <- c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del",
  "Frame_Shift_Ins","Splice_Site","In_Frame_Del","In_Frame_Ins","Nonstop_Mutation")


# ═══════════════════════════════════════════════════════════════════════════════
# DATA LOADERS
# ═══════════════════════════════════════════════════════════════════════════════

load_rppa <- function() {
  message("Loading TCGA RPPA...")
  rppa <- fread(RPPA_TSV)
  # Transpose to proteins × samples
  protein_cols <- setdiff(names(rppa), c("sample_id", "patient_id"))
  mat <- t(as.matrix(rppa[, ..protein_cols]))
  colnames(mat) <- rppa$patient_id
  storage.mode(mat) <- "numeric"
  message(sprintf("  RPPA: %d proteins × %d samples", nrow(mat), ncol(mat)))

  # Mutation status from TCGA MAF
  maf <- fread(TCGA_MAF)
  tp53_pts <- unique(substr(
    maf[Hugo_Symbol=="TP53" & Variant_Classification %in% FUN_CLASSES]$Tumor_Sample_Barcode, 1, 12))
  atm_pts  <- unique(substr(
    maf[Hugo_Symbol=="ATM"  & Variant_Classification %in% FUN_CLASSES]$Tumor_Sample_Barcode, 1, 12))
  # MSI proxy: >500 coding mutations
  msi_pts  <- maf[Variant_Classification %in% CODING_CLASSES,
    .(n = .N), by = .(pid = substr(Tumor_Sample_Barcode, 1, 12))][n > 500]$pid

  sample_tbl <- data.frame(
    sample_id   = colnames(mat),
    tp53_status = ifelse(colnames(mat) %in% tp53_pts, "mutant", "wildtype"),
    atm_status  = ifelse(colnames(mat) %in% atm_pts,  "mutant", "wildtype"),
    msi_h       = as.integer(colnames(mat) %in% msi_pts),
    stringsAsFactors = FALSE
  )
  message(sprintf("  TP53-mut=%d  ATM-mut=%d  MSI-H=%d",
    sum(sample_tbl$tp53_status=="mutant"),
    sum(sample_tbl$atm_status=="mutant"),
    sum(sample_tbl$msi_h)))

  list(mat = mat, sample_tbl = sample_tbl,
       cohort = "TCGA RPPA", unit = "RPPA z-score")
}


load_cptac <- function() {
  message("Loading CPTAC proteome...")
  prot_raw <- fread(CPTAC_PROT)
  # Keep Gene + "CPT* Log Ratio" cols (not Unshared, not NCI7 controls)
  keep_cols <- c("Gene",
    names(prot_raw)[grepl("^CPT.* Log Ratio$", names(prot_raw)) &
                    !grepl("Unshared", names(prot_raw))])
  prot_raw <- prot_raw[, ..keep_cols]
  prot_raw <- prot_raw[Gene != "Mean"]   # drop QC row

  # Duplicate gene names → keep highest mean
  genes     <- prot_raw$Gene
  mat_raw   <- as.matrix(prot_raw[, -"Gene"])
  storage.mode(mat_raw) <- "numeric"
  mean_expr <- rowMeans(mat_raw, na.rm = TRUE)
  ord       <- order(-mean_expr)
  genes     <- genes[ord];  mat_raw <- mat_raw[ord, ]
  keep      <- !duplicated(genes)
  mat_raw   <- mat_raw[keep, ]
  rownames(mat_raw) <- genes[keep]
  colnames(mat_raw) <- sub(" Log Ratio$", "", colnames(mat_raw))

  # Drop proteins with >50% missing; impute rest with row median
  frac_na  <- rowMeans(is.na(mat_raw))
  mat_raw  <- mat_raw[frac_na <= 0.5, ]
  row_med  <- apply(mat_raw, 1, median, na.rm = TRUE)
  for (i in which(rowSums(is.na(mat_raw)) > 0))
    mat_raw[i, is.na(mat_raw[i, ])] <- row_med[i]

  message(sprintf("  CPTAC: %d proteins × %d samples (after NA filter)", nrow(mat_raw), ncol(mat_raw)))

  # ID map: CPT aliquot → Case Submitter ID (C3L-...)
  idmap <- fread(CPTAC_IDMAP)
  aliquot_to_case <- setNames(idmap$`Case Submitter ID`, idmap$`Aliquot Submitter ID`)
  sample_c3l <- aliquot_to_case[colnames(mat_raw)]
  # Drop columns with no C3L mapping (NA)
  keep_mapped <- !is.na(sample_c3l)
  mat_raw    <- mat_raw[, keep_mapped, drop=FALSE]
  sample_c3l <- sample_c3l[keep_mapped]
  colnames(mat_raw) <- sample_c3l   # rename cols to C3L-...

  # Deduplicate: for same C3L case, keep aliquot with highest column mean
  if (any(duplicated(sample_c3l))) {
    col_means <- colMeans(mat_raw, na.rm=TRUE)
    ord_cols  <- order(-col_means)
    mat_raw   <- mat_raw[, ord_cols, drop=FALSE]
    keep_uniq <- !duplicated(colnames(mat_raw))
    mat_raw   <- mat_raw[, keep_uniq, drop=FALSE]
    message(sprintf("  Dedup CPTAC: %d unique C3L cases", ncol(mat_raw)))
  }

  # Mutation status from CPTAC MAF
  maf_c <- fread(CPTAC_MAF)
  setnames(maf_c, make.unique(names(maf_c)))
  maf_c[, cpt_base := sub("[0-9]{4}$", "", Tumor_Sample_Barcode)]
  idmap[, cpt_base := sub("[0-9]{4}$", "", `Aliquot Submitter ID`)]
  base_to_case <- unique(idmap[, .(cpt_base, case_sub = `Case Submitter ID`)])
  maf_m <- merge(maf_c, base_to_case, by = "cpt_base", all.x = TRUE)

  tp53_c3l <- na.omit(unique(maf_m[Hugo_Symbol=="TP53" & Variant_Classification %in% FUN_CLASSES]$case_sub))
  atm_c3l  <- na.omit(unique(maf_m[Hugo_Symbol=="ATM"  & Variant_Classification %in% FUN_CLASSES]$case_sub))
  msi_c3l  <- na.omit(unique(
    maf_m[Variant_Classification %in% CODING_CLASSES,
          .(n=.N), by=case_sub][n > 500]$case_sub))

  sample_tbl <- data.frame(
    sample_id   = colnames(mat_raw),
    tp53_status = ifelse(colnames(mat_raw) %in% tp53_c3l, "mutant", "wildtype"),
    atm_status  = ifelse(colnames(mat_raw) %in% atm_c3l,  "mutant", "wildtype"),
    msi_h       = as.integer(colnames(mat_raw) %in% msi_c3l),
    stringsAsFactors = FALSE
  )
  message(sprintf("  TP53-mut=%d  ATM-mut=%d  MSI-H=%d",
    sum(sample_tbl$tp53_status=="mutant"),
    sum(sample_tbl$atm_status=="mutant"),
    sum(sample_tbl$msi_h)))

  list(mat = mat_raw, sample_tbl = sample_tbl,
       cohort = "CPTAC", unit = "log2 TMT ratio")
}


# ═══════════════════════════════════════════════════════════════════════════════
# LIMMA + PATHWAY HELPERS
# ═══════════════════════════════════════════════════════════════════════════════

run_limma <- function(mat, tbl, design_var) {
  tbl$group <- factor(tbl[[design_var]], levels = c("wildtype", "mutant"))
  tbl$msi_h <- as.integer(tbl$msi_h)

  samples <- intersect(tbl$sample_id, colnames(mat))
  tbl     <- tbl[tbl$sample_id %in% samples, ]
  mat_sub <- mat[, tbl$sample_id, drop = FALSE]

  # Drop rows still all-NA after imputation
  mat_sub <- mat_sub[rowSums(is.na(mat_sub)) == 0, ]

  has_msi_var <- length(unique(tbl$msi_h)) > 1
  design <- if (has_msi_var)
    model.matrix(~ group + msi_h, data = tbl)
  else
    model.matrix(~ group, data = tbl)

  fit  <- lmFit(mat_sub, design)
  fit2 <- eBayes(fit, trend = TRUE)

  res           <- as.data.frame(topTable(fit2, coef = "groupmutant",
                                          number = Inf, sort.by = "P"))
  res$gene_name <- rownames(res)
  message(sprintf("    limma: %d proteins tested  (n_mut=%d  n_wt=%d)",
    nrow(res), sum(tbl$group=="mutant"), sum(tbl$group=="wildtype")))
  res
}


# For RPPA: match pathway by base gene name (strips _PS/PT suffix) → include phospho variants
extract_pathway_results <- function(limma_res, is_rppa = FALSE) {
  df <- limma_res[!is.na(limma_res$gene_name), ]
  df <- df[order(df$P.Value, na.last = TRUE), ]

  if (is_rppa) {
    df$base_gene <- sub("_.*$", "", df$gene_name)
    pw <- df[df$base_gene %in% ALL_PATHWAY, , drop = FALSE]
    pw$class <- GENE_CLASS[pw$base_gene]
  } else {
    df <- df[!duplicated(df$gene_name), ]
    pw <- df[df$gene_name %in% ALL_PATHWAY, , drop = FALSE]
    pw$class <- GENE_CLASS[pw$gene_name]
  }

  pw$padj_class <- NA_real_
  for (cl in CLASS_ORDER) {
    idx <- which(pw$class == cl)
    if (length(idx) > 0 && !all(is.na(pw$P.Value[idx])))
      pw$padj_class[idx] <- p.adjust(pw$P.Value[idx], method = "BH")
  }

  pw$class <- factor(pw$class, levels = CLASS_ORDER)
  pw <- pw[order(pw$class, pw$P.Value, na.last = TRUE), ]

  found   <- if (is_rppa) unique(pw$base_gene) else unique(pw$gene_name)
  missing <- setdiff(ALL_PATHWAY, found)
  if (length(missing) > 0)
    message(sprintf("  Pathway genes not in dataset: %s", paste(missing, collapse=", ")))
  pw
}


run_gsea_3class <- function(limma_res) {
  df <- limma_res[!is.na(limma_res$logFC) & !is.na(limma_res$adj.P.Val) &
                  !is.na(limma_res$gene_name), ]
  df <- df[order(df$adj.P.Val), ]
  # For RPPA: use base gene name
  df$gene_name <- sub("_.*$", "", df$gene_name)
  df <- df[!duplicated(df$gene_name), ]

  r         <- sign(df$logFC) * -log10(pmax(df$adj.P.Val, 1e-300))
  names(r)  <- df$gene_name
  ranks     <- sort(r, decreasing = TRUE)

  pathways  <- list(ATM_specific=GENES_ATM, Convergence=GENES_CONV,
                    TP53_specific=GENES_TP53, All_pathway=ALL_PATHWAY)
  pathways  <- pathways[sapply(pathways, function(g) sum(g %in% names(ranks)) >= 5)]
  if (length(pathways) == 0) return(NULL)

  set.seed(42)
  fgsea(pathways=pathways, stats=ranks, minSize=5, maxSize=1000, nPermSimple=10000)
}


# ═══════════════════════════════════════════════════════════════════════════════
# PLOTS
# ═══════════════════════════════════════════════════════════════════════════════

plot_lollipop <- function(pw, comp_label, unit_label, outfile) {
  df <- pw[!is.na(pw$logFC), ]
  # For RPPA phospho rows, use gene_name; for others use gene_name directly
  df$label <- if ("base_gene" %in% names(df)) df$gene_name else df$gene_name
  df <- df[order(df$class, df$logFC), ]
  df$label <- factor(df$label, levels = df$label)
  df$sig <- with(df, case_when(
    !is.na(padj_class) & padj_class < 0.001 ~ "***",
    !is.na(padj_class) & padj_class < 0.01  ~ "**",
    !is.na(padj_class) & padj_class < 0.05  ~ "*",
    TRUE ~ ""
  ))

  p <- ggplot(df, aes(x = logFC, y = label, color = class)) +
    geom_vline(xintercept = 0, linewidth=0.5, color="grey60", linetype="dashed") +
    geom_segment(aes(x=0, xend=logFC, yend=label), linewidth=0.8, alpha=0.7) +
    geom_point(aes(size = -log10(pmax(P.Value, 1e-10))), alpha=0.9) +
    geom_text(aes(x = logFC + sign(logFC+1e-9)*0.04, label=sig),
              size=4, hjust=0.5, show.legend=FALSE) +
    scale_color_manual(values=CLASS_COLORS, name="Gene class") +
    scale_size_continuous(name="-log10(p)", range=c(2,8)) +
    facet_grid(rows=vars(class), scales="free_y", space="free_y") +
    labs(title=comp_label,
         subtitle="* = padj<0.05 within class (BH)",
         x=paste0("log2 Fold Change  (mutant / wildtype)  [", unit_label, "]"),
         y=NULL) +
    theme_bw(base_size=12) +
    theme(strip.background=element_rect(fill="grey92"),
          strip.text=element_text(face="bold"),
          panel.grid.minor=element_blank())

  h <- max(6, nrow(df)*0.45 + 3)
  ggsave(outfile, p, width=9, height=h, dpi=150)
  message(sprintf("  Saved: %s", outfile))
}


plot_volcano_pathway <- function(limma_res, pw, comp_label, outfile) {
  df <- limma_res[!is.na(limma_res$adj.P.Val) & !is.na(limma_res$gene_name), ]
  df$neg_log10 <- -log10(pmax(df$adj.P.Val, 1e-50))
  base_genes   <- sub("_.*$", "", df$gene_name)
  df$class     <- ifelse(base_genes %in% ALL_PATHWAY,
                         GENE_CLASS[base_genes], "Other")
  df$class     <- factor(df$class, levels=c(CLASS_ORDER, "Other"))
  df$label     <- ifelse(base_genes %in% ALL_PATHWAY, df$gene_name, NA_character_)

  p <- ggplot(df[df$class=="Other", ], aes(logFC, neg_log10)) +
    geom_point(color="#CCCCCC", alpha=0.25, size=0.5) +
    geom_point(data=df[df$class!="Other", ], aes(color=class), size=2.8, alpha=0.85) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey50", linewidth=0.4) +
    geom_text_repel(data=df[df$class!="Other", ],
                    aes(label=label, color=class),
                    size=3.0, max.overlaps=40, show.legend=FALSE) +
    scale_color_manual(values=CLASS_COLORS, name="Gene class",
                       breaks=CLASS_ORDER, drop=FALSE) +
    labs(title=comp_label,
         subtitle="y = genome-wide adj.P.Val (display); class correction in lollipop",
         x="log2 Fold Change  (mutant / wildtype)",
         y="-log10(adj.P.Val)  [genome-wide]") +
    theme_bw(base_size=12) +
    guides(color=guide_legend(override.aes=list(size=3.5)))

  ggsave(outfile, p, width=9, height=6, dpi=150)
  message(sprintf("  Saved: %s", outfile))
}


plot_heatmap_3class <- function(mat, pw, sample_tbl, design_var, comp_label, outfile) {
  gene_col <- if ("base_gene" %in% names(pw)) "gene_name" else "gene_name"
  genes    <- intersect(as.character(pw[[gene_col]]), rownames(mat))
  if (length(genes) == 0) return(invisible(NULL))

  gene_order  <- intersect(as.character(pw[[gene_col]]), genes)
  mat_sub     <- mat[gene_order, sample_tbl$sample_id, drop=FALSE]
  mat_z       <- t(scale(t(mat_sub)))
  mat_z       <- pmin(pmax(mat_z, -2.5), 2.5)

  col_order   <- order(sample_tbl[[design_var]] == "wildtype")
  mat_z       <- mat_z[, col_order, drop=FALSE]
  stbl_ord    <- sample_tbl[col_order, ]

  var_label   <- if (design_var=="tp53_status") "TP53" else "ATM"
  annot_col   <- data.frame(
    Group = factor(stbl_ord[[design_var]], levels=c("mutant","wildtype")),
    MSI_H = factor(stbl_ord$msi_h),
    row.names = stbl_ord$sample_id
  )
  colnames(annot_col)[1] <- var_label
  annot_colors_col <- setNames(
    list(c(mutant="#E41A1C", wildtype="#377EB8"),
         c("0"="white", "1"="#888888")),
    c(var_label, "MSI_H"))

  # Row class annotation
  if ("base_gene" %in% names(pw)) {
    base_for_row <- pw$base_gene[match(gene_order, pw$gene_name)]
  } else {
    base_for_row <- gene_order
  }
  row_class <- GENE_CLASS[base_for_row]
  annot_row <- data.frame(Class=factor(row_class, levels=CLASS_ORDER),
                          row.names=gene_order)

  h <- max(5, length(genes)*0.38 + 3)
  png(outfile, width=12, height=h, units="in", res=150)
  pheatmap(
    mat_z,
    annotation_col    = annot_col,
    annotation_row    = annot_row,
    annotation_colors = c(annot_colors_col, list(Class=CLASS_COLORS)),
    cluster_rows      = FALSE,
    cluster_cols      = FALSE,
    show_colnames     = FALSE,
    color             = colorRampPalette(rev(brewer.pal(9,"RdBu")))(100),
    breaks            = seq(-2.5, 2.5, length.out=101),
    main              = comp_label,
    fontsize_row      = 9,
    gaps_row          = cumsum(table(factor(row_class, levels=CLASS_ORDER)))[-3]
  )
  dev.off()
  message(sprintf("  Saved: %s", outfile))
}


plot_boxplots_3class <- function(mat, pw, sample_tbl, design_var, unit_label, comp_label, outfile) {
  gene_col <- "gene_name"
  genes    <- intersect(as.character(pw[[gene_col]]), rownames(mat))
  if (length(genes) == 0) return(invisible(NULL))

  long_df <- lapply(genes, function(g) {
    r   <- pw[as.character(pw[[gene_col]])==g, ]
    sig <- with(r, case_when(
      !is.na(padj_class) & padj_class<0.001 ~ "***",
      !is.na(padj_class) & padj_class<0.01  ~ "**",
      !is.na(padj_class) & padj_class<0.05  ~ "*",
      TRUE ~ "ns"))
    lfc <- if (!is.na(r$logFC[1])) sprintf("%.2f",r$logFC[1]) else "NA"
    bg  <- if ("base_gene" %in% names(r)) r$base_gene[1] else g
    data.frame(gene=g, expr=mat[g, sample_tbl$sample_id],
               status=sample_tbl[[design_var]], class=GENE_CLASS[bg],
               gene_label=paste0(g,"\n(FC=",lfc," ",sig,")"),
               stringsAsFactors=FALSE)
  })
  long_df <- do.call(rbind, long_df)
  long_df$class      <- factor(long_df$class, levels=CLASS_ORDER)
  lbl_order          <- unique(long_df$gene_label[order(long_df$class, long_df$gene)])
  long_df$gene_label <- factor(long_df$gene_label, levels=lbl_order)

  ncols <- min(6, length(genes))
  nrows <- ceiling(length(genes)/ncols)

  p <- ggplot(long_df, aes(status, expr, fill=status)) +
    geom_boxplot(outlier.size=0.4, width=0.6, linewidth=0.4) +
    scale_fill_manual(values=c(mutant="#E41A1C", wildtype="#377EB8")) +
    facet_wrap(~gene_label, ncol=ncols, scales="free_y") +
    labs(title=comp_label, x=NULL, y=unit_label) +
    theme_bw(base_size=10) +
    theme(legend.position="none",
          strip.text=element_text(face="bold", size=8),
          axis.text.x=element_text(angle=30, hjust=1))

  ggsave(outfile, p, width=ncols*3.2, height=nrows*3.2, dpi=150)
  message(sprintf("  Saved: %s", outfile))
}


write_results <- function(pw, gsea_res, outdir) {
  write.csv(pw, file.path(outdir, "pathway_DE_results.csv"), row.names=FALSE)

  if (!is.null(gsea_res)) {
    g <- as.data.frame(gsea_res)
    g$leadingEdge <- sapply(g$leadingEdge, paste, collapse=";")
    write.csv(g, file.path(outdir, "gsea_pathway_summary.csv"), row.names=FALSE)
  }

  summary_tbl <- pw %>%
    mutate(base = sub("_.*$","",gene_name)) %>%
    group_by(class) %>%
    summarise(
      n_proteins     = n(),
      n_genes        = n_distinct(base),
      n_sig_padj05   = sum(!is.na(padj_class) & padj_class<0.05),
      n_up           = sum(!is.na(logFC) & logFC>0),
      n_down         = sum(!is.na(logFC) & logFC<0),
      mean_logFC     = round(mean(logFC, na.rm=TRUE),3),
      median_logFC   = round(median(logFC, na.rm=TRUE),3),
      .groups        = "drop")
  write.csv(summary_tbl, file.path(outdir, "class_summary.csv"), row.names=FALSE)
  message("  class_summary:")
  print(as.data.frame(summary_tbl))
}


# ═══════════════════════════════════════════════════════════════════════════════
# run_comparison()
# ═══════════════════════════════════════════════════════════════════════════════

run_comparison <- function(cohort_data, design_var, outdir, comp_label, min_n=5) {
  mat       <- cohort_data$mat
  tbl       <- cohort_data$sample_tbl
  unit      <- cohort_data$unit
  is_rppa   <- grepl("RPPA", cohort_data$cohort)

  other_var <- if (design_var=="tp53_status") "atm_status" else "tp53_status"
  var_label <- if (design_var=="tp53_status") "TP53" else "ATM"

  # Filter to comparison-relevant samples
  tbl_sub  <- tbl[tbl[[other_var]]=="wildtype", ]
  n_mut    <- sum(tbl_sub[[design_var]]=="mutant")
  n_wt     <- sum(tbl_sub[[design_var]]=="wildtype")
  message(sprintf("  %s-mut/WT: N=%d  |  WT/WT: N=%d", var_label, n_mut, n_wt))

  if (n_mut < min_n) {
    message(sprintf("  Mutant group too small (N=%d) — skipping", n_mut))
    return(invisible(NULL))
  }
  if (n_mut < 10)
    message(sprintf("  WARNING: N=%d is small — interpret results cautiously", n_mut))

  # Run limma
  limma_res <- run_limma(mat, tbl_sub, design_var)
  pw        <- extract_pathway_results(limma_res, is_rppa=is_rppa)
  gsea      <- run_gsea_3class(limma_res)

  # Subset mat to these samples for plotting
  samp      <- intersect(tbl_sub$sample_id, colnames(mat))
  tbl_plot  <- tbl_sub[tbl_sub$sample_id %in% samp, ]

  dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
  plot_lollipop(pw, comp_label, unit, file.path(outdir,"lollipop_pathway.png"))
  plot_volcano_pathway(limma_res, pw, comp_label, file.path(outdir,"volcano_pathway.png"))
  plot_heatmap_3class(mat, pw, tbl_plot, design_var, comp_label,
                      file.path(outdir,"heatmap_pathway.png"))
  plot_boxplots_3class(mat, pw, tbl_plot, design_var, unit, comp_label,
                       file.path(outdir,"boxplots_pathway.png"))
  write_results(pw, gsea, outdir)
  message(sprintf("  → %s\n", outdir))
}


# ═══════════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════════

rppa_data  <- load_rppa()
cptac_data <- load_cptac()

message("\n════════ TCGA RPPA — Comp 1 ════════")
run_comparison(rppa_data, "tp53_status",
  file.path(OUT_BASE,"TCGA_RPPA/comp1_TP53mut"),
  "TCGA RPPA · Comp 1: TP53-mut/ATM-WT  vs  TP53-WT/ATM-WT")

message("\n════════ TCGA RPPA — Comp 2 ════════")
run_comparison(rppa_data, "atm_status",
  file.path(OUT_BASE,"TCGA_RPPA/comp2_ATMmut"),
  "TCGA RPPA · Comp 2: ATM-mut/TP53-WT  vs  ATM-WT/TP53-WT")

message("\n════════ CPTAC — Comp 1 ════════")
run_comparison(cptac_data, "tp53_status",
  file.path(OUT_BASE,"CPTAC/comp1_TP53mut"),
  "CPTAC · Comp 1: TP53-mut/ATM-WT  vs  TP53-WT/ATM-WT")

message("\n════════ CPTAC — Comp 2 ════════")
run_comparison(cptac_data, "atm_status",
  file.path(OUT_BASE,"CPTAC/comp2_ATMmut"),
  "CPTAC · Comp 2: ATM-mut/TP53-WT  vs  ATM-WT/TP53-WT")

message("\n✓  proteome_pathway_DE.R complete.")
