#!/usr/bin/env Rscript
#
# phospho_pathway_DE.R — Limma differential phosphosite abundance, ATM-TP53 pathway
#
# Two comparisons:
#   Comp 1: TP53-mut / ATM-WT  vs  TP53-WT / ATM-WT
#   Comp 2: ATM-mut / TP53-WT  vs  ATM-WT  / TP53-WT  (N=9, exploratory)
#
# Data: CPTAC Gastric phosphoproteome (45,750 phosphosites, 206 samples)
# Method: limma (eBayes, trend=TRUE); design ~ group + msi_h
#
# Multiple testing:
#   - Genome-wide BH across all phosphosites (reported as adj.P.Val_global)
#   - Within ATM substrate sites only (4 detected sites, reported as padj_atm)
#   - Pathway gene sites: within-class BH (same 3 classes as mRNA/proteome)
#
# IMPORTANT LIMITATION: Key canonical ATM substrate phosphosites are NOT detected
# in this TMT dataset: CHEK2-T68, H2AX-S139 (γH2AX), TRIM28-S824, CHEK1-S345,
# TP53-S15. Results should be interpreted as supportive, not definitive, evidence
# for ATM kinase activity changes.
#
# Usage:
#   micromamba run -n ProjectGeneration \
#     Rscript TP53_ATM_paper/scripts/phospho_pathway_DE.R

suppressPackageStartupMessages({
  library(limma)
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
SCRIPT_DIR  <- get_script_dir()
REPO_ROOT   <- normalizePath(file.path(SCRIPT_DIR, "../.."))

PHOSPHO_TSV <- Sys.glob(file.path(REPO_ROOT, "CPTAC/Phosphoproteome/*phosphosite.tmt18*"))[1]
CPTAC_MAF   <- file.path(REPO_ROOT, "CPTAC/somatic_mutation/cptac_gc_merged.maf.tsv")
CPTAC_IDMAP <- file.path(REPO_ROOT, "CPTAC/proteome/PDC_study_biospecimen_06062026_223409.csv")
OUT_BASE    <- file.path(REPO_ROOT, "TP53_ATM_paper/results/Phospho")

for (d in c(file.path(OUT_BASE, "comp1_TP53mut"),
            file.path(OUT_BASE, "comp2_ATMmut")))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)


# ── Gene / site sets ──────────────────────────────────────────────────────────
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

# ATM direct substrate sites actually detected in this dataset
# Canonical sites MISSING: CHEK2-T68, H2AX-S139, TRIM28-S824, CHEK1-S345, TP53-S15
ATM_SUBSTRATES <- c(
  "ATM_pS1981",    # ATM autophosphorylation (activation marker)
  "NBN_pS343",     # NBS1-S343, canonical ATM substrate
  "BRCA1_pS1524",  # BRCA1-S1524, ATM substrate
  "CHEK1_pS317"    # CHEK1-S317, primarily ATR but within ATM-CHEK1 axis
)

FUN_CLASSES <- c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del",
  "Frame_Shift_Ins","Splice_Site","In_Frame_Del","In_Frame_Ins",
  "Nonstop_Mutation","Translation_Start_Site")
CODING_CLASSES <- c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del",
  "Frame_Shift_Ins","Splice_Site","In_Frame_Del","In_Frame_Ins","Nonstop_Mutation")


# ── Helpers ───────────────────────────────────────────────────────────────────
format_site <- function(s) {
  # "s19" → "pS19";  "t68" → "pT68";  "y123" → "pY123"
  gsub("([sty])(\\d+)", "p\\U\\1\\2", s, perl = TRUE)
}

# Parse "ENSP00000329012.5:s120" → site part "pS120"
parse_site <- function(row_id) {
  site_raw <- sub(".*:", "", row_id)
  format_site(site_raw)
}


# ── Data loader ───────────────────────────────────────────────────────────────
load_phospho <- function() {
  message("Loading CPTAC phosphoproteome...")
  raw <- fread(PHOSPHO_TSV, data.table = TRUE)

  # Identify annotation columns (last 3: Peptide, Gene, Organism)
  lr_cols   <- names(raw)[grepl("^CPT.* Log Ratio$", names(raw)) &
                          !grepl("Unshared", names(raw))]
  gene_col  <- "Gene"
  site_col  <- "Phosphosite"

  # Filter to single-site rows only (skip multisite like s19s22)
  is_single <- grepl("^ENSP[^:]+:[sty]\\d+$", raw[[site_col]])
  raw       <- raw[is_single]
  message(sprintf("  Single-site rows: %d / %d total", nrow(raw), nrow(raw) + sum(!is_single)))

  # Build row labels: GENE_pSXXX
  sites     <- parse_site(raw[[site_col]])
  genes     <- raw[[gene_col]]
  row_ids   <- paste0(genes, "_", sites)

  # Build matrix
  mat_raw <- as.matrix(raw[, ..lr_cols])
  storage.mode(mat_raw) <- "numeric"
  colnames(mat_raw) <- sub(" Log Ratio$", "", lr_cols)  # CPT IDs

  # Deduplicate row labels (same gene+site from different ENSP isoforms)
  row_means <- rowMeans(mat_raw, na.rm = TRUE)
  ord       <- order(-abs(row_means))
  mat_raw   <- mat_raw[ord, ]
  row_ids   <- row_ids[ord]
  keep      <- !duplicated(row_ids)
  mat_raw   <- mat_raw[keep, ]
  rownames(mat_raw) <- row_ids[keep]

  # Drop sites with >50% missing
  frac_na  <- rowMeans(is.na(mat_raw))
  mat_raw  <- mat_raw[frac_na <= 0.5, ]
  message(sprintf("  After NA filter (>50%%): %d sites retained", nrow(mat_raw)))

  # Impute remaining NAs with row median
  row_med <- apply(mat_raw, 1, median, na.rm = TRUE)
  for (i in which(rowSums(is.na(mat_raw)) > 0))
    mat_raw[i, is.na(mat_raw[i, ])] <- row_med[i]

  # ID map: CPT aliquot → C3L case ID
  idmap <- fread(CPTAC_IDMAP)
  aliquot_to_case <- setNames(idmap$`Case Submitter ID`, idmap$`Aliquot Submitter ID`)
  sample_c3l <- aliquot_to_case[colnames(mat_raw)]

  # Drop unmapped columns
  keep_mapped <- !is.na(sample_c3l)
  mat_raw     <- mat_raw[, keep_mapped, drop = FALSE]
  sample_c3l  <- sample_c3l[keep_mapped]
  colnames(mat_raw) <- sample_c3l

  # Deduplicate columns: keep highest-mean aliquot per case
  if (any(duplicated(sample_c3l))) {
    col_means <- colMeans(mat_raw, na.rm = TRUE)
    ord_cols  <- order(-col_means)
    mat_raw   <- mat_raw[, ord_cols, drop = FALSE]
    keep_uniq <- !duplicated(colnames(mat_raw))
    mat_raw   <- mat_raw[, keep_uniq, drop = FALSE]
    message(sprintf("  Dedup: %d unique C3L cases", ncol(mat_raw)))
  }

  # Mutation status
  maf_c <- fread(CPTAC_MAF)
  setnames(maf_c, make.unique(names(maf_c)))
  maf_c[, cpt_base := sub("[0-9]{4}$", "", Tumor_Sample_Barcode)]
  idmap[, cpt_base := sub("[0-9]{4}$", "", `Aliquot Submitter ID`)]
  base_to_case <- unique(idmap[, .(cpt_base, case_sub = `Case Submitter ID`)])
  maf_m <- merge(maf_c, base_to_case, by = "cpt_base", all.x = TRUE)

  tp53_c3l <- na.omit(unique(maf_m[Hugo_Symbol=="TP53" & Variant_Classification %in% FUN_CLASSES]$case_sub))
  atm_c3l  <- na.omit(unique(maf_m[Hugo_Symbol=="ATM"  & Variant_Classification %in% FUN_CLASSES]$case_sub))
  msi_c3l  <- na.omit(unique(
    maf_m[Variant_Classification %in% CODING_CLASSES, .(n=.N), by=case_sub][n > 500]$case_sub))

  sample_tbl <- data.frame(
    sample_id   = colnames(mat_raw),
    tp53_status = ifelse(colnames(mat_raw) %in% tp53_c3l, "mutant", "wildtype"),
    atm_status  = ifelse(colnames(mat_raw) %in% atm_c3l,  "mutant", "wildtype"),
    msi_h       = as.integer(colnames(mat_raw) %in% msi_c3l),
    stringsAsFactors = FALSE
  )
  message(sprintf("  TP53-mut=%d  ATM-mut=%d  MSI-H=%d  total=%d",
    sum(sample_tbl$tp53_status=="mutant"),
    sum(sample_tbl$atm_status=="mutant"),
    sum(sample_tbl$msi_h),
    nrow(sample_tbl)))

  list(mat = mat_raw, sample_tbl = sample_tbl)
}


# ── Limma ─────────────────────────────────────────────────────────────────────
run_limma_phospho <- function(mat, tbl, design_var) {
  tbl$group <- factor(tbl[[design_var]], levels = c("wildtype", "mutant"))
  other_var <- if (design_var == "tp53_status") "atm_status" else "tp53_status"
  tbl       <- tbl[tbl[[other_var]] == "wildtype", ]

  samp    <- intersect(tbl$sample_id, colnames(mat))
  tbl     <- tbl[tbl$sample_id %in% samp, ]
  mat_sub <- mat[, tbl$sample_id, drop = FALSE]
  mat_sub <- mat_sub[rowSums(is.na(mat_sub)) == 0, ]

  has_msi_var <- length(unique(tbl$msi_h)) > 1
  design <- if (has_msi_var)
    model.matrix(~ group + msi_h, data = tbl)
  else
    model.matrix(~ group, data = tbl)

  fit  <- lmFit(mat_sub, design)
  fit2 <- eBayes(fit, trend = TRUE)
  res  <- as.data.frame(topTable(fit2, coef = "groupmutant",
                                 number = Inf, sort.by = "P"))
  res$site_id <- rownames(res)
  res$gene    <- sub("_p[STY].*", "", res$site_id)

  var_label <- if (design_var == "tp53_status") "TP53" else "ATM"
  message(sprintf("  Limma %s-mut: %d sites, n_mut=%d, n_wt=%d",
    var_label, nrow(res), sum(tbl$group=="mutant"), sum(tbl$group=="wildtype")))
  res
}


# ── Extract pathway + ATM substrate sites ─────────────────────────────────────
extract_pathway_sites <- function(limma_res) {
  pw <- limma_res[limma_res$gene %in% ALL_PATHWAY, , drop = FALSE]
  pw$class <- factor(GENE_CLASS[pw$gene], levels = CLASS_ORDER)

  # Within-class BH
  pw$padj_class <- NA_real_
  for (cl in CLASS_ORDER) {
    idx <- which(pw$class == cl)
    if (length(idx) > 0)
      pw$padj_class[idx] <- p.adjust(pw$P.Value[idx], method = "BH")
  }

  # Genome-wide BH (for reference)
  pw$adj.P.Val_global <- limma_res$adj.P.Val[match(rownames(pw), rownames(limma_res))]

  pw <- pw[order(pw$class, pw$P.Value), ]

  missing_genes <- setdiff(ALL_PATHWAY, unique(pw$gene))
  if (length(missing_genes) > 0)
    message(sprintf("  Pathway genes with 0 detected sites: %s",
                    paste(missing_genes, collapse=", ")))
  pw
}


extract_atm_substrate_sites <- function(limma_res) {
  atm_res <- limma_res[limma_res$site_id %in% ATM_SUBSTRATES, , drop = FALSE]
  if (nrow(atm_res) == 0) return(atm_res)

  # BH within these detected ATM substrate sites
  atm_res$padj_atm <- p.adjust(atm_res$P.Value, method = "BH")
  atm_res <- atm_res[order(atm_res$P.Value), ]
  atm_res
}


# ── Plots ──────────────────────────────────────────────────────────────────────
plot_lollipop_pathway <- function(pw, comp_label, outfile) {
  df <- pw[!is.na(pw$logFC), ]
  df <- df[order(df$class, df$logFC), ]
  df$site_label <- factor(df$site_id, levels = df$site_id)
  df$sig <- with(df, case_when(
    !is.na(padj_class) & padj_class < 0.001 ~ "***",
    !is.na(padj_class) & padj_class < 0.01  ~ "**",
    !is.na(padj_class) & padj_class < 0.05  ~ "*",
    TRUE ~ ""
  ))

  p <- ggplot(df, aes(x = logFC, y = site_label, color = class)) +
    geom_vline(xintercept = 0, linewidth = 0.5, color = "grey60", linetype = "dashed") +
    geom_segment(aes(x=0, xend=logFC, yend=site_label), linewidth=0.7, alpha=0.7) +
    geom_point(aes(size = -log10(pmax(P.Value, 1e-10))), alpha=0.9) +
    geom_text(aes(x = logFC + sign(logFC+1e-9)*0.04, label=sig),
              size=3.5, hjust=0.5, show.legend=FALSE) +
    scale_color_manual(values = CLASS_COLORS, name = "Gene class") +
    scale_size_continuous(name = "-log10(p)", range = c(2, 7)) +
    facet_grid(rows = vars(class), scales = "free_y", space = "free_y") +
    labs(title    = comp_label,
         subtitle = "* = padj<0.05 within class (BH) | CPTAC phosphoproteome",
         x = "log2 Fold Change  (mutant / wildtype)  [log2 TMT ratio]",
         y = NULL) +
    theme_bw(base_size = 11) +
    theme(strip.background = element_rect(fill = "grey92"),
          strip.text       = element_text(face = "bold"),
          panel.grid.minor = element_blank(),
          axis.text.y      = element_text(size = 7))

  h <- max(7, nrow(df) * 0.28 + 3)
  ggsave(outfile, p, width = 10, height = h, dpi = 150)
  message(sprintf("  Saved: %s", outfile))
}


plot_atm_substrates <- function(atm_res, mat, sample_tbl, design_var, comp_label, outfile) {
  if (nrow(atm_res) == 0) { message("  No ATM substrate sites to plot"); return(invisible(NULL)) }

  other_var <- if (design_var=="tp53_status") "atm_status" else "tp53_status"
  tbl_sub   <- sample_tbl[sample_tbl[[other_var]] == "wildtype", ]
  samp      <- intersect(tbl_sub$sample_id, colnames(mat))
  tbl_sub   <- tbl_sub[tbl_sub$sample_id %in% samp, ]

  long_df <- lapply(atm_res$site_id, function(s) {
    r   <- atm_res[atm_res$site_id == s, ]
    sig <- with(r, case_when(
      padj_atm < 0.001 ~ "***",
      padj_atm < 0.01  ~ "**",
      padj_atm < 0.05  ~ "*",
      TRUE ~ "ns"))
    lfc <- sprintf("%.2f", r$logFC[1])
    data.frame(
      site   = s,
      expr   = mat[s, tbl_sub$sample_id],
      status = tbl_sub[[design_var]],
      label  = paste0(s, "\n(FC=", lfc, " ", sig, ")"),
      stringsAsFactors = FALSE
    )
  })
  long_df        <- do.call(rbind, long_df)
  lbl_order      <- unique(long_df$label[order(match(long_df$site, atm_res$site_id))])
  long_df$label  <- factor(long_df$label, levels = lbl_order)

  p <- ggplot(long_df, aes(status, expr, fill = status)) +
    geom_boxplot(outlier.size = 0.5, width = 0.55, linewidth = 0.4) +
    geom_jitter(width = 0.12, size = 0.9, alpha = 0.5) +
    scale_fill_manual(values = c(mutant="#E41A1C", wildtype="#377EB8")) +
    facet_wrap(~label, ncol = 4, scales = "free_y") +
    labs(title    = paste0(comp_label, " — ATM substrate sites"),
         subtitle = "BH correction within 4 detected ATM substrates; * padj<0.05",
         x = NULL, y = "log2 TMT ratio") +
    theme_bw(base_size = 11) +
    theme(legend.position = "none",
          strip.text      = element_text(face = "bold", size = 9),
          axis.text.x     = element_text(angle = 30, hjust = 1))

  ggsave(outfile, p, width = 12, height = 5, dpi = 150)
  message(sprintf("  Saved: %s", outfile))
}


plot_heatmap_pathway <- function(mat, pw, sample_tbl, design_var, comp_label, outfile) {
  other_var <- if (design_var=="tp53_status") "atm_status" else "tp53_status"
  tbl_sub   <- sample_tbl[sample_tbl[[other_var]] == "wildtype", ]
  samp      <- intersect(tbl_sub$sample_id, colnames(mat))
  tbl_sub   <- tbl_sub[tbl_sub$sample_id %in% samp, ]

  sites    <- intersect(rownames(pw), rownames(mat))
  if (length(sites) < 2) return(invisible(NULL))

  mat_sub <- mat[sites, tbl_sub$sample_id, drop = FALSE]
  mat_z   <- t(scale(t(mat_sub)))
  mat_z   <- pmin(pmax(mat_z, -2.5), 2.5)

  col_ord <- order(tbl_sub[[design_var]] == "wildtype")
  mat_z   <- mat_z[, col_ord, drop = FALSE]
  tbl_ord <- tbl_sub[col_ord, ]

  var_label <- if (design_var == "tp53_status") "TP53" else "ATM"
  annot_col <- data.frame(
    row.names = tbl_ord$sample_id,
    MSI_H = factor(tbl_ord$msi_h)
  )
  annot_col[[var_label]] <- factor(tbl_ord[[design_var]],
                                   levels = c("mutant","wildtype"))
  annot_col_colors <- setNames(
    list(c(mutant="#E41A1C", wildtype="#377EB8"), c("0"="white","1"="#888888")),
    c(var_label, "MSI_H"))

  row_genes  <- pw$gene[match(sites, rownames(pw))]
  row_class  <- GENE_CLASS[row_genes]
  annot_row  <- data.frame(Class = factor(row_class, levels = CLASS_ORDER),
                            row.names = sites)

  h <- max(6, length(sites) * 0.22 + 3)
  png(outfile, width = 14, height = h, units = "in", res = 150)
  pheatmap(
    mat_z,
    annotation_col    = annot_col,
    annotation_row    = annot_row,
    annotation_colors = c(annot_col_colors, list(Class = CLASS_COLORS)),
    cluster_rows      = TRUE,
    cluster_cols      = FALSE,
    show_colnames     = FALSE,
    color             = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
    breaks            = seq(-2.5, 2.5, length.out = 101),
    main              = comp_label,
    fontsize_row      = 7,
    gaps_col          = cumsum(table(factor(tbl_ord[[design_var]],
                                            levels = c("mutant","wildtype"))))[-2]
  )
  dev.off()
  message(sprintf("  Saved: %s", outfile))
}


plot_volcano_pathway <- function(limma_res, pw, comp_label, outfile) {
  df <- limma_res[!is.na(limma_res$adj.P.Val), ]
  df$neg_log10 <- -log10(pmax(df$adj.P.Val, 1e-50))
  df$class     <- ifelse(df$gene %in% ALL_PATHWAY, GENE_CLASS[df$gene], "Other")
  df$class     <- factor(df$class, levels = c(CLASS_ORDER, "Other"))
  df$label     <- ifelse(df$gene %in% ALL_PATHWAY, df$site_id, NA_character_)

  p <- ggplot(df[df$class == "Other", ], aes(logFC, neg_log10)) +
    geom_point(color = "#CCCCCC", alpha = 0.15, size = 0.4) +
    geom_point(data = df[df$class != "Other", ], aes(color = class),
               size = 2.2, alpha = 0.85) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = "grey50", linewidth = 0.4) +
    geom_text_repel(data = df[df$class != "Other", ],
                    aes(label = label, color = class),
                    size = 2.5, max.overlaps = 50, show.legend = FALSE) +
    scale_color_manual(values = CLASS_COLORS, name = "Gene class",
                       breaks = CLASS_ORDER, drop = FALSE) +
    labs(title    = comp_label,
         subtitle = "Highlighted: 35 ATM-TP53 pathway genes (all detected phosphosites)",
         x = "log2 Fold Change  (mutant / wildtype)",
         y = "-log10(adj.P.Val)  [genome-wide]") +
    theme_bw(base_size = 11)

  ggsave(outfile, p, width = 10, height = 6.5, dpi = 150)
  message(sprintf("  Saved: %s", outfile))
}


write_results <- function(pw, atm_res, outdir, design_var, n_mut, n_wt) {
  write.csv(pw,      file.path(outdir, "pathway_phospho_results.csv"), row.names = FALSE)
  write.csv(atm_res, file.path(outdir, "atm_substrate_results.csv"),   row.names = FALSE)

  var_label <- if (design_var == "tp53_status") "TP53" else "ATM"
  summary_tbl <- pw %>%
    group_by(class) %>%
    summarise(
      n_sites       = n(),
      n_genes       = n_distinct(gene),
      n_sig_class05 = sum(!is.na(padj_class) & padj_class < 0.05),
      n_up          = sum(logFC > 0),
      n_down        = sum(logFC < 0),
      mean_logFC    = round(mean(logFC, na.rm=TRUE), 3),
      .groups       = "drop")
  write.csv(summary_tbl, file.path(outdir, "class_summary.csv"), row.names = FALSE)

  message(sprintf("  %s-mut n=%d vs WT n=%d", var_label, n_mut, n_wt))
  message("  class_summary:")
  print(as.data.frame(summary_tbl))
  if (nrow(atm_res) > 0) {
    message("  ATM substrate sites:")
    print(atm_res[, c("site_id","logFC","P.Value","padj_atm")])
  }
}


# ── run_comparison ─────────────────────────────────────────────────────────────
run_comparison <- function(pdata, design_var, outdir, comp_label) {
  mat <- pdata$mat
  tbl <- pdata$sample_tbl

  other_var <- if (design_var == "tp53_status") "atm_status" else "tp53_status"
  n_mut     <- sum(tbl[tbl[[other_var]]=="wildtype", design_var] == "mutant")
  n_wt      <- sum(tbl[tbl[[other_var]]=="wildtype", design_var] == "wildtype")

  if (n_mut < 5) {
    message(sprintf("  Skipping %s — only %d mutant samples", comp_label, n_mut))
    return(invisible(NULL))
  }
  if (n_mut < 10)
    message(sprintf("  WARNING: n_mut=%d — exploratory only", n_mut))

  limma_res <- run_limma_phospho(mat, tbl, design_var)
  pw        <- extract_pathway_sites(limma_res)
  atm_res   <- extract_atm_substrate_sites(limma_res)

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  plot_lollipop_pathway(pw, comp_label, file.path(outdir, "lollipop_pathway.png"))
  plot_volcano_pathway(limma_res, pw, comp_label, file.path(outdir, "volcano_pathway.png"))
  plot_heatmap_pathway(mat, pw, tbl, design_var, comp_label,
                       file.path(outdir, "heatmap_pathway.png"))
  plot_atm_substrates(atm_res, mat, tbl, design_var, comp_label,
                      file.path(outdir, "atm_substrates_boxplot.png"))
  write_results(pw, atm_res, outdir, design_var, n_mut, n_wt)
  message(sprintf("  → %s\n", outdir))
}


# ── Main ───────────────────────────────────────────────────────────────────────
pdata <- load_phospho()

message("\n════════ Comp 1: TP53-mut / ATM-WT  vs  TP53-WT / ATM-WT ════════")
run_comparison(pdata, "tp53_status",
  file.path(OUT_BASE, "comp1_TP53mut"),
  "CPTAC Phospho · Comp 1: TP53-mut/ATM-WT  vs  TP53-WT/ATM-WT")

message("\n════════ Comp 2: ATM-mut / TP53-WT  vs  ATM-WT / TP53-WT ════════")
run_comparison(pdata, "atm_status",
  file.path(OUT_BASE, "comp2_ATMmut"),
  "CPTAC Phospho · Comp 2: ATM-mut/TP53-WT  vs  ATM-WT/TP53-WT  [N=9, exploratory]")

message("\n✓  phospho_pathway_DE.R complete.")
