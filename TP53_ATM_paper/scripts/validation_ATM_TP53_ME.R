#!/usr/bin/env Rscript
# =============================================================================
# Cross-cohort validation: ATM–TP53 mutual exclusivity
# TP53-ATM paper
# =============================================================================
# Model (consistent across cohorts):
#   TP53 ~ ATM + TMB_log [+ stage_cat] [+ histo_cat]
#   — stage_cat and histo_cat included where data are available
#   — MSI_H is NOT in the model; used only for FDU MSS-only stratum
#
# Cohorts: FDU, TCGA-STAD, OncoSG 2018, HK Pfizer 2014, MSK 2023
# Excluded (same as Phase2): TMUCIH 2015, egc_mskcc_2020, egc_msk_tp53_ccr_2022,
#   egc_msk_2017 (ATM too few after STAD filter)
#
# Min ATM+ to run logistic: 5
#
# Outputs → TP53_ATM_paper/results/validation_cohorts/
#   tables/per_cohort_ATM.tsv     — per-cohort logistic + Fisher
#   tables/meta_ATM_full.tsv      — meta-analysis, full cohort
#   tables/meta_ATM_mss.tsv       — meta-analysis, MSS (FDU only)
#   figures/forest_ATM_TP53.png   — publication-quality forest plot
#   figures/summary_tile_ATM.png  — per-cohort tile summary
#
# Run:
#   /home/sgao30/micromamba/envs/tcga_bioc/bin/Rscript \
#       TP53_ATM_paper/scripts/validation_ATM_TP53_ME.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
  library(metafor)
})

# ── Paths ─────────────────────────────────────────────────────────────────────
ROOT     <- "/ShangGaoAIProjects/Gastric/genome"
OUT      <- file.path(ROOT, "TP53_ATM_paper/results/validation_cohorts")
dir.create(file.path(OUT, "tables"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "figures"), recursive = TRUE, showWarnings = FALSE)

FILT_MAF <- file.path(ROOT, "analysis_results/Phase0/enrichment/somt_filtered.maf")
MASTER   <- file.path(ROOT, "analysis_results/Phase0/filtered_master_sample_table.tsv")
TCGA_MAF <- file.path(ROOT, "TCGA-STAD/somatic_mutation/somatic_mutation.csv")
TCGA_CLIN <- file.path(ROOT, "TCGA-STAD/clinical/clinical_patient_stad.csv")
CBIO     <- file.path(ROOT, "Other_cBioPartal_stomach")
EXOME_MB <- 38.0

FUN_CLASSES <- c(
  "Missense_Mutation", "Nonsense_Mutation",
  "Frame_Shift_Del",   "Frame_Shift_Ins",
  "Splice_Site",       "In_Frame_Del", "In_Frame_Ins",
  "Nonstop_Mutation",  "Translation_Start_Site"
)
ATM_ENTREZ  <- 472L
TP53_ENTREZ <- 7157L
MIN_MUT     <- 5L

PUB <- theme_classic(base_size = 11) +
  theme(
    axis.text    = element_text(color = "black"),
    plot.title   = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

# =============================================================================
# HELPERS
# =============================================================================

# Map stage strings → I / II / III / IV (NA if unrecognised)
harmonize_stage <- function(x) {
  x <- trimws(as.character(x))
  dplyr::case_when(
    grepl("^(Stage\\s*)?I[AB]?$",  x, ignore.case=TRUE) ~ "I",
    grepl("^(Stage\\s*)?II",        x, ignore.case=TRUE) ~ "II",
    grepl("^(Stage\\s*)?III",       x, ignore.case=TRUE) ~ "III",
    grepl("^(Stage\\s*)?IV$",       x, ignore.case=TRUE) ~ "IV",
    TRUE ~ NA_character_
  )
}

# Map Lauren/histology → adeno / signet / other (NA if unknown)
harmonize_histo_lauren <- function(x) {
  x <- trimws(tolower(as.character(x)))
  dplyr::case_when(
    grepl("intestinal",              x) ~ "adeno",
    grepl("diffuse|signet",          x) ~ "signet",
    grepl("mixed",                   x) ~ "other",
    grepl("^adenocarcinoma$",        x) ~ "adeno",
    TRUE ~ NA_character_
  )
}

# Run logistic regression; return one-row data.frame
run_logistic <- function(df, cohort, label, stratum = "full") {
  df  <- as.data.frame(df)
  n_atm  <- sum(df$ATM,  na.rm=TRUE)
  n_tp53 <- sum(df$TP53, na.rm=TRUE)
  n_both <- sum(df$ATM==1 & df$TP53==1, na.rm=TRUE)
  n_tot  <- nrow(df)

  tbl    <- table(ATM=df$ATM, TP53=df$TP53)
  fish   <- tryCatch(fisher.test(tbl), error=function(e) NULL)
  fish_p <- if (!is.null(fish)) fish$p.value  else NA_real_
  fish_or <- if (!is.null(fish)) fish$estimate else NA_real_

  base_row <- data.frame(
    cohort=cohort, label=label, stratum=stratum,
    n=n_tot, n_atm=n_atm, n_tp53=n_tp53, n_both=n_both,
    fish_p=fish_p, fish_or=fish_or,
    OR=NA_real_, OR_lo=NA_real_, OR_hi=NA_real_,
    beta=NA_real_, se=NA_real_, pvalue=NA_real_,
    model=NA_character_, status=NA_character_,
    stringsAsFactors=FALSE
  )

  if (n_atm < MIN_MUT) {
    base_row$status <- "too_few_ATM"; return(base_row)
  }

  # Build formula from available columns
  terms <- "ATM + TMB_log"
  if ("stage_cat" %in% names(df) && sum(!is.na(df$stage_cat)) >= 10)
    terms <- paste0(terms, " + stage_cat")
  if ("histo_cat" %in% names(df) && sum(!is.na(df$histo_cat)) >= 10)
    terms <- paste0(terms, " + histo_cat")
  fmla <- as.formula(paste("TP53 ~", terms))

  # Drop rows with NA in any model variable
  use_cols <- c("TP53","ATM","TMB_log",
                if (grepl("stage_cat", terms)) "stage_cat",
                if (grepl("histo_cat", terms)) "histo_cat")
  df2 <- df[complete.cases(df[, intersect(use_cols, names(df))]), ]

  if (nrow(df2) < 20 || sum(df2$ATM) < MIN_MUT) {
    base_row$status <- "too_few_after_cc"; base_row$model <- deparse(fmla)
    return(base_row)
  }

  # Set factor reference levels
  if ("stage_cat" %in% names(df2))
    df2$stage_cat <- relevel(factor(df2$stage_cat), ref="I")
  if ("histo_cat" %in% names(df2))
    df2$histo_cat <- relevel(factor(df2$histo_cat), ref="adeno")

  fit <- tryCatch(
    suppressWarnings(glm(fmla, df2, family=binomial)),
    error=function(e) NULL
  )
  if (is.null(fit)) {
    base_row$status <- "glm_failed"; base_row$model <- deparse(fmla)
    return(base_row)
  }

  ct <- summary(fit)$coefficients
  b  <- ct["ATM", ]
  ci <- tryCatch(confint.default(fit)["ATM",], error=function(e) c(NA_real_,NA_real_))
  base_row$OR     <- exp(b[1])
  base_row$OR_lo  <- exp(ci[1])
  base_row$OR_hi  <- exp(ci[2])
  base_row$beta   <- b[1]
  base_row$se     <- b[2]
  base_row$pvalue <- b[4]
  base_row$model  <- deparse(fmla)
  base_row$status <- "ok"
  base_row
}

# =============================================================================
# 1. DATA LOADERS
# =============================================================================

# ── FDU (CohortB_All) ─────────────────────────────────────────────────────────
load_fdu <- function() {
  cat("Loading FDU ...\n")
  maf  <- fread(FILT_MAF)
  maf[, sample_base := sub("F\\d+$", "", Tumor_Sample_Barcode)]
  clin <- fread(MASTER, na.strings=c("","NA","X"))

  stage_map <- c(IA="I",IB="I",IIA="II",IIB="II",IIIA="III",IIIB="III",IIIC="III",IV="IV")
  clin[, TMB_log  := log2(pmax(suppressWarnings(as.numeric(TMB_report_numeric)), 0.1))]
  clin[, MSI_H    := as.integer(MSI %in% c("MSI-H","MSI"))]
  clin[, stage_cat := stage_map[as.character(path_stage)]]
  clin[, histo_cat := fcase(
    grepl("腺癌|Adeno",  histology, ignore.case=TRUE), "adeno",
    grepl("印戒|Signet", histology, ignore.case=TRUE), "signet",
    !is.na(histology), "other",
    default=NA_character_
  )]

  mut_dt   <- unique(maf[Variant_Classification %in% FUN_CLASSES &
                           Hugo_Symbol %in% c("TP53","ATM"),
                         .(sample_base, Hugo_Symbol)])
  mut_dt[, v := 1L]
  mut_wide <- dcast(mut_dt, sample_base ~ Hugo_Symbol, value.var="v", fill=0L)
  all_samp <- data.table(sample_base=unique(maf$sample_base))
  mut_wide <- merge(all_samp, mut_wide, by="sample_base", all.x=TRUE)
  for (col in c("ATM","TP53"))
    if (!col %in% names(mut_wide)) mut_wide[, (col):=0L] else
      set(mut_wide, which(is.na(mut_wide[[col]])), col, 0L)

  df <- merge(mut_wide,
              clin[, .(sample_id, TMB_log, MSI_H, stage_cat, histo_cat)],
              by.x="sample_base", by.y="sample_id", all=FALSE)

  df_full <- df[complete.cases(df[, c("TMB_log")])]
  df_mss  <- df_full[MSI_H == 0]
  cat(sprintf("  FDU full: N=%d  ATM+=%d  TP53+=%d\n",
              nrow(df_full), sum(df_full$ATM), sum(df_full$TP53)))
  cat(sprintf("  FDU MSS:  N=%d  ATM+=%d  TP53+=%d\n",
              nrow(df_mss),  sum(df_mss$ATM),  sum(df_mss$TP53)))
  list(full=df_full, mss=df_mss)
}

# ── TCGA-STAD ─────────────────────────────────────────────────────────────────
load_tcga <- function() {
  cat("Loading TCGA-STAD ...\n")
  maf <- fread(TCGA_MAF)
  maf <- maf[grepl("-01[AB]-", Tumor_Sample_Barcode)]
  maf[, sample_id    := Tumor_Sample_Barcode]
  maf[, patient_barcode := sub("-[0-9]{2}[A-Z]-.*", "", sample_id)]

  tmb <- maf[Variant_Classification %in% FUN_CLASSES,
             .(n_mut=.N), by=sample_id]
  all_samp <- data.table(sample_id=unique(maf$sample_id))
  tmb <- merge(all_samp, tmb, by="sample_id", all.x=TRUE)
  tmb[is.na(n_mut), n_mut:=0L]
  tmb[, TMB_log := log2(pmax(n_mut/EXOME_MB, 0.1))]
  tmb[, MSI_H   := as.integer(n_mut > 500)]
  tmb[, patient_barcode := sub("-[0-9]{2}[A-Z]-.*", "", sample_id)]

  # Parse TCGA clinical: row 1 = column names, rows 2-4 = metadata, row 5+ = data
  clin_raw <- tryCatch({
    tmp <- read.csv(TCGA_CLIN, skip=0, header=TRUE, stringsAsFactors=FALSE,
                    na.strings=c("","NA","[Not Available]","[Discrepancy]","[Unknown]"))
    tmp <- tmp[-(1:3), ]      # drop 3 metadata rows
    as.data.table(tmp)
  }, error=function(e) NULL)
  if (!is.null(clin_raw) && "bcr_patient_barcode" %in% names(clin_raw)) {
    # Stage only (TCGA histological_type is unreliable)
    stage_col <- "ajcc_pathologic_tumor_stage"
    if (stage_col %in% names(clin_raw)) {
      clin_sub <- clin_raw[, .(patient_barcode=bcr_patient_barcode,
                               stage_raw=get(stage_col))]
      tmb <- merge(tmb, clin_sub, by="patient_barcode", all.x=TRUE)
      tmb[, stage_cat := harmonize_stage(stage_raw)]
    }
  }

  mut_dt <- unique(maf[Variant_Classification %in% FUN_CLASSES &
                          Hugo_Symbol %in% c("TP53","ATM"),
                        .(sample_id, Hugo_Symbol)])
  mut_dt[, v:=1L]
  mut_wide <- dcast(mut_dt, sample_id ~ Hugo_Symbol, value.var="v", fill=0L)
  mut_wide <- merge(all_samp, mut_wide, by="sample_id", all.x=TRUE)
  for (col in c("ATM","TP53"))
    if (!col %in% names(mut_wide)) mut_wide[,(col):=0L] else
      set(mut_wide, which(is.na(mut_wide[[col]])), col, 0L)

  df <- merge(mut_wide, tmb[, intersect(
                c("sample_id","TMB_log","MSI_H","stage_cat","histo_cat"),
                names(tmb)), with=FALSE],
              by="sample_id", all=FALSE)
  df <- df[complete.cases(df[, c("TMB_log")])]
  df_mss <- df[MSI_H == 0]
  cat(sprintf("  TCGA full: N=%d  ATM+=%d  TP53+=%d  n_both=%d\n",
              nrow(df),     sum(df$ATM),     sum(df$TP53),     sum(df$ATM==1 & df$TP53==1)))
  cat(sprintf("  TCGA MSS:  N=%d  ATM+=%d  TP53+=%d  n_both=%d\n",
              nrow(df_mss), sum(df_mss$ATM), sum(df_mss$TP53), sum(df_mss$ATM==1 & df_mss$TP53==1)))
  list(full=df, mss=df_mss)
}

# ── cBioPortal (generic) ──────────────────────────────────────────────────────
# msi_col / msi_pos_vals: used to define MSI_H (stored but not in model)
# stage_col / histo_col: column name in clinical_patient.tsv (joined via patientId)
load_cbio <- function(cohort_id, label,
                      msi_col=NULL, msi_pos_vals=NULL,
                      stage_col=NULL, histo_col=NULL,
                      histo_type=c("lauren","msk"),
                      stad_only=FALSE, tmb_col="TMB_NONSYNONYMOUS") {
  histo_type <- match.arg(histo_type)
  cat("Loading", cohort_id, "...\n")
  mut  <- fread(file.path(CBIO, cohort_id, "mutations.tsv"))
  clin <- fread(file.path(CBIO, cohort_id, "clinical_sample.tsv"),
                na.strings=c("","NA","N/A","Unknown","[Not Available]"))
  clin_pat <- fread(file.path(CBIO, cohort_id, "clinical_patient.tsv"),
                    na.strings=c("","NA","N/A","Unknown","[Not Available]"))

  # STAD-only filter
  if (stad_only && "CANCER_TYPE_DETAILED" %in% names(clin)) {
    stad_ids <- clin[grepl("(?i)Stomach Adenocarcinoma", CANCER_TYPE_DETAILED), sampleId]
    clin <- clin[sampleId %in% stad_ids]
    mut  <- mut[sampleId  %in% stad_ids]
    cat(sprintf("  STAD filter: %d samples\n", nrow(clin)))
  }

  # Deduplicate: 1 sample per patient (prefer Primary)
  if ("patientId" %in% names(mut)) {
    sm <- unique(mut[, .(sampleId, patientId)])
    if ("SAMPLE_TYPE" %in% names(clin))
      sm <- merge(sm, clin[, .(sampleId, SAMPLE_TYPE)], by="sampleId", all.x=TRUE)
    sm[, pref := if ("SAMPLE_TYPE" %in% names(sm))
      fifelse(SAMPLE_TYPE=="Primary",1L,2L) else 1L]
    setorder(sm, patientId, pref)
    keep <- sm[, .(sampleId=sampleId[1]), by=patientId]$sampleId
    mut  <- mut[sampleId %in% keep]
    clin <- clin[sampleId %in% keep]
    cat(sprintf("  Dedup: %d samples (1/patient)\n", length(keep)))
  }

  # TMB
  if (tmb_col %in% names(clin)) {
    tmb <- clin[, .(sample_id=sampleId,
                    TMB_log=log2(pmax(suppressWarnings(as.numeric(get(tmb_col))),0.1)))]
  } else {
    tc <- mut[mutationType %in% FUN_CLASSES, .(n_mut=.N), by=sampleId]
    tmb <- merge(data.table(sample_id=unique(clin$sampleId)),
                 tc[, .(sample_id=sampleId, TMB_log=log2(pmax(n_mut/EXOME_MB,0.1)))],
                 by="sample_id", all.x=TRUE)
    tmb[is.na(TMB_log), TMB_log:=log2(0.1)]
  }

  # Patient-level clinical (stage, histology, MSI)
  if ("patientId" %in% names(mut)) {
    pid_sid <- unique(mut[sampleId %in% tmb$sample_id, .(sampleId, patientId)])

    pat_cols <- c("patientId",
                  if (!is.null(msi_col)   && msi_col   %in% names(clin_pat)) msi_col,
                  if (!is.null(stage_col) && stage_col %in% names(clin_pat)) stage_col,
                  if (!is.null(histo_col) && histo_col %in% names(clin_pat)) histo_col)
    pat_sub  <- clin_pat[, unique(pat_cols), with=FALSE]
    pid_sid  <- merge(pid_sid, pat_sub, by="patientId", all.x=TRUE)

    # MSI_H
    if (!is.null(msi_col) && msi_col %in% names(pid_sid))
      pid_sid[, MSI_H := as.integer(get(msi_col) %in% msi_pos_vals)]

    # Stage
    if (!is.null(stage_col) && stage_col %in% names(pid_sid))
      pid_sid[, stage_cat := harmonize_stage(get(stage_col))]

    # Histology
    if (!is.null(histo_col) && histo_col %in% names(pid_sid)) {
      if (histo_type == "lauren")
        pid_sid[, histo_cat := harmonize_histo_lauren(get(histo_col))]
      else  # msk: Adenocarcinoma / Signet_Diffuse / Other / Squamous
        pid_sid[, histo_cat := fcase(
          grepl("adenocarcinoma", tolower(get(histo_col))), "adeno",
          grepl("signet|diffuse", tolower(get(histo_col))), "signet",
          !is.na(get(histo_col)), "other",
          default=NA_character_
        )]
    }

    join_cols <- intersect(c("sampleId","MSI_H","stage_cat","histo_cat"), names(pid_sid))
    tmb <- merge(tmb, pid_sid[, join_cols, with=FALSE],
                 by.x="sample_id", by.y="sampleId", all.x=TRUE)
  }

  # Binary ATM / TP53
  all_samp <- data.table(sample_id=unique(clin$sampleId))
  for (gene_sym in c("ATM","TP53")) {
    eid  <- if (gene_sym=="ATM") ATM_ENTREZ else TP53_ENTREZ
    samp <- unique(mut[entrezGeneId==eid & mutationType %in% FUN_CLASSES, sampleId])
    all_samp[, (gene_sym):=as.integer(sample_id %in% samp)]
  }

  df <- merge(all_samp, tmb, by="sample_id", all=FALSE)
  df <- df[complete.cases(df[, c("TMB_log")])]
  cat(sprintf("  %s: N=%d  ATM+=%d  TP53+=%d\n",
              cohort_id, nrow(df), sum(df$ATM), sum(df$TP53)))
  df
}

# =============================================================================
# 2. LOAD ALL COHORTS
# =============================================================================

cat("\n=== Loading cohorts ===\n")
fdu_data <- load_fdu()

tcga_data <- load_tcga()

oncosg_full <- load_cbio("stad_oncosg_2018",
  label      = "OncoSG 2018",
  msi_col    = "MOLECULAR_SUBTYPE", msi_pos_vals = "MSI",
  stage_col  = "STAGE", histo_col = "LAURENS_CLASSIFICATION",
  histo_type = "lauren")
# MSS subset if MSI_H available
oncosg_mss <- if ("MSI_H" %in% names(oncosg_full)) {
  s <- as.data.frame(oncosg_full)[as.data.frame(oncosg_full)$MSI_H==0,]
  cat(sprintf("  OncoSG MSS: N=%d  ATM+=%d\n", nrow(s), sum(s$ATM))); s
} else NULL

msk23_full <- load_cbio("egc_msk_2023",
  label      = "MSK 2023",
  msi_col    = "MSI_TYPE", msi_pos_vals = "Instable",
  stage_col  = "STAGE", histo_col = "HISTOLOGY",
  histo_type = "msk",
  stad_only  = TRUE)
msk23_mss <- if ("MSI_H" %in% names(msk23_full)) {
  s <- as.data.frame(msk23_full)[as.data.frame(msk23_full)$MSI_H==0,]
  cat(sprintf("  MSK2023 MSS: N=%d  ATM+=%d\n", nrow(s), sum(s$ATM))); s
} else NULL

# =============================================================================
# 3. RUN MODELS
# =============================================================================

cat("\n=== Running ATM–TP53 models ===\n")

cohort_meta <- list(
  list(df=fdu_data$full,  cohort="FDU",         label="FDU  [panel, Chinese, discovery]",    stratum="full", pool_meta=TRUE),
  list(df=fdu_data$mss,   cohort="FDU_MSS",     label="FDU (MSS only)  [panel, Chinese]",    stratum="mss",  pool_meta=FALSE),
  list(df=tcga_data$mss,  cohort="TCGA_MSS",    label="TCGA-STAD (MSS)  [WES]",              stratum="mss",  pool_meta=TRUE),
  list(df=tcga_data$full, cohort="TCGA_ALL",    label="TCGA-STAD (all)  [WES, incl. MSI-H]", stratum="full", pool_meta=FALSE),
  list(df=msk23_full,     cohort="MSK_2023",    label="MSK 2023  [IMPACT341, STAD only]",     stratum="full", pool_meta=TRUE),
  list(df=msk23_mss,      cohort="MSK_2023_MSS",label="MSK 2023 (MSS only)  [IMPACT341]",    stratum="mss",  pool_meta=FALSE),
  list(df=oncosg_full,    cohort="OncoSG_2018", label="OncoSG 2018  N=147  [WES, Asian]",    stratum="full", pool_meta=TRUE)
)
# Drop any entries with NULL data (no MSI info → no MSS subset)
cohort_meta <- Filter(function(co) !is.null(co$df), cohort_meta)

res_list <- lapply(cohort_meta, function(co)
  run_logistic(co$df, co$cohort, co$label, co$stratum))
res      <- do.call(rbind, res_list)
rownames(res) <- NULL

fwrite(res, file.path(OUT, "tables", "per_cohort_ATM.tsv"), sep="\t")
cat("\nPer-cohort results:\n")
print(res[, c("cohort","stratum","n","n_atm","n_tp53","n_both",
               "OR","OR_lo","OR_hi","pvalue","fish_p","model","status")])

# =============================================================================
# 4. META-ANALYSIS
# =============================================================================

run_meta <- function(res_sub, label) {
  ma_data <- res_sub[res_sub$status=="ok" & !is.na(res_sub$beta) &
                       !is.na(res_sub$se) & res_sub$se>0, ]
  if (nrow(ma_data) < 2) return(NULL)
  ma <- tryCatch(rma(yi=beta, sei=se, data=ma_data, method="REML"),
                 error=function(e) NULL)
  if (is.null(ma)) return(NULL)
  out <- data.frame(
    stratum=label, n_cohorts=nrow(ma_data),
    pooled_OR=exp(coef(ma)), OR_lo=exp(ma$ci.lb), OR_hi=exp(ma$ci.ub),
    p_value=ma$pval, I2=round(ma$I2,1), tau2=ma$tau2,
    beta_pool=coef(ma), se_pool=ma$se
  )
  cat(sprintf("  [%s] Pooled OR=%.3f [%.3f–%.3f]  p=%.2e  I²=%.1f%%  (%d cohorts)\n",
              label, out$pooled_OR, out$OR_lo, out$OR_hi, out$p_value, out$I2, out$n_cohorts))
  out
}

cat("\n=== Meta-analysis ===\n")
# Add pool_meta flag from cohort_meta
pool_flag <- setNames(
  sapply(cohort_meta, function(co) co$pool_meta),
  sapply(cohort_meta, function(co) co$cohort)
)
res$pool_meta <- pool_flag[res$cohort]

meta_full <- run_meta(res[!is.na(res$pool_meta) & res$pool_meta, ], "primary_4cohorts")

if (!is.null(meta_full))
  fwrite(meta_full, file.path(OUT, "tables", "meta_ATM_full.tsv"), sep="\t")

# =============================================================================
# 5. FOREST PLOT
# =============================================================================

# Display order (top-to-bottom in forest plot):
# FDU full, FDU MSS | TCGA MSS, TCGA ALL | MSK full, MSK MSS | OncoSG full, OncoSG MSS
# Pooled at bottom
cohort_order <- c(
  "FDU", "FDU_MSS",
  "TCGA_MSS", "TCGA_ALL",
  "MSK_2023", "MSK_2023_MSS",
  "OncoSG_2018", "OncoSG_MSS",
  "Pooled (RE)"
)

display_labels <- c(
  "FDU"          = "FDU (all)  [targeted panel, Chinese]  ★ discovery",
  "FDU_MSS"      = "  └ FDU (MSS only)",
  "TCGA_MSS"     = "TCGA-STAD (MSS only)  [WES]",
  "TCGA_ALL"     = "  └ TCGA-STAD (all, incl. MSI-H)",
  "MSK_2023"     = "MSK 2023 (all)  [IMPACT341, STAD only]",
  "MSK_2023_MSS" = "  └ MSK 2023 (MSS only)",
  "OncoSG_2018"  = "OncoSG 2018 (all)  [WES, Asian]",
  "OncoSG_MSS"   = "  └ OncoSG 2018 (MSS only)",
  "Pooled (RE)"  = "Pooled  (random-effects, 4 primary cohorts)"
)

# Include logistic-OK rows; for rows that had too few ATM after complete-cases,
# fall back to Fisher OR (no CI) so the point still appears in the plot.
res_plot <- res[res$status %in% c("ok","too_few_after_cc"), ]
res_plot$OR     <- ifelse(!is.na(res_plot$OR),    res_plot$OR,    res_plot$fish_or)
res_plot$OR_lo  <- ifelse(res_plot$status=="ok",  res_plot$OR_lo, NA_real_)
res_plot$OR_hi  <- ifelse(res_plot$status=="ok",  res_plot$OR_hi, NA_real_)
res_plot$pvalue <- ifelse(!is.na(res_plot$pvalue),res_plot$pvalue,res_plot$fish_p)
res_plot$fisher_only <- res_plot$status != "ok"

plot_df <- res_plot[!is.na(res_plot$OR),
               c("cohort","stratum","n","n_atm","n_both","OR","OR_lo","OR_hi",
                 "pvalue","pool_meta","fisher_only")]
# Add pooled row
if (!is.null(meta_full)) {
  pool_n <- sum(res$n[!is.na(res$pool_meta)&res$pool_meta&res$status=="ok"])
  plot_df <- rbind(plot_df, data.frame(
    cohort="Pooled (RE)", stratum="pooled",
    n=pool_n, n_atm=NA_integer_, n_both=NA_integer_,
    OR=meta_full$pooled_OR, OR_lo=meta_full$OR_lo, OR_hi=meta_full$OR_hi,
    pvalue=meta_full$p_value, pool_meta=NA, fisher_only=FALSE,
    stringsAsFactors=FALSE
  ))
}

plot_df$is_pooled  <- plot_df$cohort == "Pooled (RE)"
plot_df$is_primary <- !is.na(plot_df$pool_meta) & plot_df$pool_meta
plot_df$is_sub     <- !plot_df$is_pooled & !plot_df$is_primary
plot_df$sig        <- !is.na(plot_df$pvalue) & plot_df$pvalue < 0.05
plot_df$OR_lo_d    <- ifelse(!is.na(plot_df$OR_lo), pmax(plot_df$OR_lo, 0.03), plot_df$OR)
plot_df$OR_hi_d    <- ifelse(!is.na(plot_df$OR_hi), pmin(plot_df$OR_hi, 5),   plot_df$OR)
# col_grp: primary / sub / pooled; fisher_only rows get their own shape via fisher_only flag
plot_df$col_grp    <- ifelse(plot_df$is_pooled, "pooled",
                      ifelse(plot_df$is_primary, "primary", "sub"))

# Visual grouping: each paired cohort forms a block separated by a line
# Group breaks happen ABOVE: FDU, TCGA_MSS, MSK_2023, OncoSG_2018, Pooled(RE)
group_tops <- c("FDU","TCGA_MSS","MSK_2023","OncoSG_2018","Pooled (RE)")

valid_order <- cohort_order[cohort_order %in% plot_df$cohort]
plot_df$cohort_f <- factor(plot_df$cohort, levels=rev(valid_order))

# Formatted p-value labels (logistic p for ok rows; Fisher p for fisher_only rows)
fmt_pval <- function(p, fisher_only=FALSE) {
  if (is.na(p)) return("")
  prefix <- if (fisher_only) "Fisher p=" else "p="
  if (p < 0.001)     paste0(prefix, formatC(p, format="e", digits=2))
  else if (p < 0.01) paste0(prefix, sprintf("%.4f", p))
  else               paste0(prefix, sprintf("%.3f",  p))
}
plot_df$pval_lbl <- mapply(fmt_pval, plot_df$pvalue, plot_df$fisher_only)
# Pooled row: add I² to p-value label
pool_i2 <- if (!is.null(meta_full)) meta_full$I2 else 0
plot_df$pval_lbl[plot_df$cohort=="Pooled (RE)"] <-
  sprintf("p=%.2e  I²=%.0f%%", meta_full$p_value, pool_i2)

# OR [95% CI] text column
plot_df$or_ci_lbl <- with(plot_df, ifelse(
  cohort == "Pooled (RE)",
  sprintf("OR=%.2f [%.2f–%.2f]", OR, OR_lo, OR_hi),
  ifelse(!is.na(OR_lo),
    sprintf("%.2f [%.2f–%.2f]", OR, OR_lo, OR_hi),
    ifelse(!is.na(OR), sprintf("%.2f (Fisher)", OR), "")
  )
))

subtitle_str <- if (!is.null(meta_full))
  sprintf("Pooled (RE, 4 cohorts): OR=%.3f [%.3f–%.3f]  p=%.2e  I²=%.1f%%",
          meta_full$pooled_OR, meta_full$OR_lo, meta_full$OR_hi,
          meta_full$p_value, meta_full$I2) else ""

# Compute y positions for separator lines (above group_tops)
sep_df <- plot_df[plot_df$cohort %in% group_tops, c("cohort","cohort_f")]
sep_df$yint <- as.numeric(sep_df$cohort_f) + 0.5

p_forest <- ggplot(plot_df, aes(x=OR, y=cohort_f)) +
  geom_hline(data=sep_df,
             aes(yintercept=yint),
             color="grey60", linewidth=0.45, inherit.aes=FALSE) +
  geom_vline(xintercept=1, linetype="dashed", color="grey50", linewidth=0.6) +
  geom_errorbar(
    aes(xmin=OR_lo_d, xmax=OR_hi_d,
        color=col_grp, linewidth=col_grp),
    width=0.30, orientation="y"
  ) +
  geom_point(
    aes(fill=col_grp, color=col_grp,
        size=col_grp,
        shape=interaction(col_grp, fisher_only, drop=TRUE)),
    stroke=0.8
  ) +
  # Separator line between CI area and text columns
  geom_vline(xintercept=4.8, color="grey75", linewidth=0.35) +
  # OR [95% CI] column
  geom_text(aes(x=7.5, label=or_ci_lbl),
            hjust=0.5, size=3.0, color="grey20", fontface="plain") +
  # p-value column
  geom_text(aes(x=17, label=pval_lbl),
            hjust=0.5, size=3.0, color="grey25") +
  scale_shape_manual(
    values=c("primary.FALSE"=21, "sub.FALSE"=21, "pooled.FALSE"=23,
             "sub.TRUE"=24),   # open triangle for Fisher-only rows
    guide="none") +
  scale_size_manual(
    values=c(primary=3.8, sub=2.8, pooled=5.2), guide="none") +
  scale_linewidth_manual(
    values=c(primary=1.0, sub=0.65, pooled=1.3), guide="none") +
  scale_fill_manual(
    values=c(primary="#C0392B", sub="#F0A0A0", pooled="#1a1a1a"), guide="none") +
  scale_color_manual(
    values=c(primary="#C0392B", sub="#D47070", pooled="#1a1a1a"), guide="none") +
  scale_x_log10(
    limits=c(0.03, 25),
    breaks=c(0.05,0.1,0.2,0.5,1,2),
    labels=c("0.05","0.1","0.2","0.5","1","2")
  ) +
  scale_y_discrete(
    labels=function(x) {
      lbl <- display_labels[as.character(x)]
      ifelse(is.na(lbl), as.character(x), lbl)
    }
  ) +
  coord_cartesian(clip="off") +
  PUB +
  theme(
    panel.grid.major.x = element_line(color="grey90", linewidth=0.3),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y        = element_text(size=9, hjust=1),
    plot.caption       = element_text(size=8, color="grey50"),
    plot.margin        = margin(t=28, r=10, b=5, l=5)
  ) +
  # Column headers (drawn above the top data row using clip="off")
  annotate("text", x=7.5, y=length(valid_order)+0.75,
           label="Adj. OR [95% CI]", fontface="bold", size=3.3, hjust=0.5, color="grey20") +
  annotate("text", x=17,  y=length(valid_order)+0.75,
           label="p (logistic)", fontface="bold", size=3.3, hjust=0.5, color="grey20") +
  annotate("segment", x=4.8, xend=25,
           y=length(valid_order)+0.4, yend=length(valid_order)+0.4,
           color="grey60", linewidth=0.4) +
  labs(
    x        = "Odds Ratio (log scale)    ← Mutually exclusive  |  Co-occurring →",
    y        = NULL,
    title    = "ATM–TP53 mutual exclusivity: cross-cohort validation",
    subtitle = subtitle_str,
    caption  = paste("Model: TP53 ~ ATM + log₂(TMB) + stage + histology (where available)",
                     "  |  Light pink = sub-group rows (not used in pooled)")
  )

ggsave(file.path(OUT, "figures", "forest_ATM_TP53.png"),
       p_forest, width=14, height=7, dpi=250, bg="white")
cat("  Saved: forest_ATM_TP53.png\n")

# ── Tile summary — all cohorts ────────────────────────────────────────────────
tile_df <- res[res$cohort != "Pooled (RE)",
               c("cohort","n","n_atm","n_both","OR","pvalue","status")]
tile_df$log2OR  <- log2(pmax(pmin(tile_df$OR,8),0.125))
tile_df$star    <- with(tile_df, ifelse(
  !is.na(pvalue)&pvalue<0.001,"***",
  ifelse(!is.na(pvalue)&pvalue<0.01,"**",
  ifelse(!is.na(pvalue)&pvalue<0.05,"*",
  ifelse(status=="too_few_ATM","n<5","ns")))))

all_coh_order <- c("FDU","FDU_MSS","TCGA_MSS","TCGA_ALL","MSK_2023","MSK_2023_MSS","OncoSG_2018","OncoSG_MSS")
tile_df$cohort_f <- factor(tile_df$cohort,
                            levels=rev(all_coh_order[all_coh_order %in% tile_df$cohort]))

p_tile <- ggplot(tile_df, aes(x="ATM–TP53", y=cohort_f)) +
  geom_tile(aes(fill=log2OR), color="white", linewidth=1) +
  geom_text(aes(label=paste0(star,"\n",
                  ifelse(is.na(OR),"—",sprintf("OR=%.2f",OR)))),
            size=3.8, fontface="bold", color="white", lineheight=1.2) +
  scale_fill_gradient2(
    low="#C0392B", mid="white", high="#2980B9", midpoint=0,
    limits=c(-3,3), oob=squish, na.value="grey80",
    name="log₂(OR)"
  ) +
  scale_y_discrete(labels=function(x) {
    lbl <- display_labels[as.character(x)]
    ifelse(is.na(lbl),as.character(x),lbl)
  }) +
  PUB +
  theme(axis.ticks=element_blank(), axis.line=element_blank(),
        axis.text.x=element_text(face="bold",size=11)) +
  labs(x=NULL, y=NULL,
       title="ATM–TP53 ME: per-cohort summary",
       subtitle="TP53 ~ ATM + log₂(TMB) + stage + histology (where available)")

ggsave(file.path(OUT,"figures","summary_tile_ATM.png"),
       p_tile, width=6.5, height=6.5, dpi=250, bg="white")
cat("  Saved: summary_tile_ATM.png\n")

# =============================================================================
# 6. ADDITIONAL FIGURES
# =============================================================================

# Helper: Wilson 95% CI for a proportion
wilson_ci <- function(n, x) {
  z <- 1.96; p <- x/n
  d <- 1 + z^2/n
  ctr <- (p + z^2/(2*n)) / d
  mrg <- z * sqrt(p*(1-p)/n + z^2/(4*n^2)) / d
  c(lo=max(0, ctr-mrg), hi=min(1, ctr+mrg))
}

# Short x-axis labels (4 primary cohorts used in meta)
SHORT_LBL <- c(
  "FDU"         = "FDU",
  "TCGA_MSS"    = "TCGA\n(MSS)",
  "MSK_2023"    = "MSK\n2023",
  "OncoSG_2018" = "OncoSG\n2018"
)
COH_ORDER <- c("FDU","TCGA_MSS","MSK_2023","OncoSG_2018")

# Collect per-cohort raw data — primary cohorts used in meta only
main_cohs <- cohort_meta[sapply(cohort_meta, function(x) isTRUE(x$pool_meta))]

# ── Build summary rows ─────────────────────────────────────────────────────────
sum_rows <- lapply(main_cohs, function(co) {
  df  <- as.data.frame(co$df)
  n   <- nrow(df)
  a_p <- df[df$ATM==1,]; a_n <- df[df$ATM==0,]
  n_atm <- nrow(a_p); n_noatm <- nrow(a_n)
  tp53_atm_pos <- sum(a_p$TP53); tp53_atm_neg <- sum(a_n$TP53)
  n_both <- sum(df$ATM==1 & df$TP53==1)
  exp_both <- (sum(df$ATM) / n) * (sum(df$TP53) / n) * n
  ci_pos <- wilson_ci(n_atm,   tp53_atm_pos)
  ci_neg <- wilson_ci(n_noatm, tp53_atm_neg)
  data.frame(
    cohort       = co$cohort,
    short_lbl    = SHORT_LBL[co$cohort],
    n            = n,
    n_atm        = n_atm,
    n_noatm      = n_noatm,
    n_tp53       = sum(df$TP53),
    n_atm_tp53   = n_both,
    n_atm_notp53 = n_atm - n_both,
    n_noatm_tp53 = tp53_atm_neg,
    n_noatm_notp53 = n_noatm - tp53_atm_neg,
    pct_tp53_atm = tp53_atm_pos / n_atm * 100,
    pct_tp53_noatm = tp53_atm_neg / n_noatm * 100,
    pct_tp53_atm_lo = ci_pos["lo"]*100,
    pct_tp53_atm_hi = ci_pos["hi"]*100,
    pct_tp53_noatm_lo = ci_neg["lo"]*100,
    pct_tp53_noatm_hi = ci_neg["hi"]*100,
    pct_atm      = n_atm / n * 100,
    pct_tp53     = sum(df$TP53) / n * 100,
    obs_both     = n_both,
    exp_both     = exp_both,
    fish_p = { fp <- res$fish_p[res$cohort==co$cohort & res$stratum==co$stratum]; if(length(fp)>0) fp[1] else NA_real_ },
    stringsAsFactors=FALSE
  )
})
sum_df <- do.call(rbind, sum_rows)
sum_df$cohort_f <- factor(sum_df$cohort, levels=COH_ORDER[COH_ORDER %in% sum_df$cohort])
sum_df$sig_lbl  <- with(sum_df, ifelse(!is.na(fish_p)&fish_p<0.001,"***",
                          ifelse(!is.na(fish_p)&fish_p<0.01,"**",
                          ifelse(!is.na(fish_p)&fish_p<0.05,"*","ns"))))
sum_df$sig_col  <- with(sum_df, ifelse(!is.na(fish_p)&fish_p<0.05,"#C0392B","grey50"))

# ── Fig A: TP53 mutation rate in ATM+ vs ATM- ─────────────────────────────────
# Wide → long for bar chart
bar_df <- rbind(
  data.frame(sum_df[,c("cohort","cohort_f","short_lbl","sig_lbl","sig_col","fish_p")],
             group="ATM+", pct=sum_df$pct_tp53_atm,
             lo=sum_df$pct_tp53_atm_lo, hi=sum_df$pct_tp53_atm_hi,
             n_grp=sum_df$n_atm, stringsAsFactors=FALSE),
  data.frame(sum_df[,c("cohort","cohort_f","short_lbl","sig_lbl","sig_col","fish_p")],
             group="ATM−", pct=sum_df$pct_tp53_noatm,
             lo=sum_df$pct_tp53_noatm_lo, hi=sum_df$pct_tp53_noatm_hi,
             n_grp=sum_df$n_noatm, stringsAsFactors=FALSE)
)
bar_df$group   <- factor(bar_df$group, levels=c("ATM+","ATM−"))
bar_df$cohort_f <- factor(bar_df$cohort, levels=COH_ORDER[COH_ORDER %in% bar_df$cohort])

# Significance label y-position: just above the ATM+ upper CI
sig_pos <- sum_df[, c("cohort","sig_lbl","sig_col","pct_tp53_atm_hi")]
sig_pos$y <- sig_pos$pct_tp53_atm_hi + 4

p_tp53rate <- ggplot(bar_df, aes(x=cohort_f, y=pct, fill=group)) +
  geom_col(position=position_dodge(0.7), width=0.65, color="white", linewidth=0.4) +
  geom_errorbar(aes(ymin=lo, ymax=hi),
                position=position_dodge(0.7), width=0.22, linewidth=0.7,
                color="grey30") +
  geom_text(data=sig_pos,
            aes(x=factor(cohort, levels=COH_ORDER[COH_ORDER %in% cohort]),
                y=y, label=sig_lbl, color=sig_col),
            inherit.aes=FALSE, size=4.5, fontface="bold") +
  geom_text(data=bar_df[bar_df$group=="ATM+",],
            aes(x=cohort_f, y=-3, label=paste0("n=",n_grp)),
            position=position_nudge(x=-0.17), size=2.8, color="grey40") +
  geom_text(data=bar_df[bar_df$group=="ATM−",],
            aes(x=cohort_f, y=-3, label=paste0("n=",n_grp)),
            position=position_nudge(x=0.17), size=2.8, color="grey40") +
  scale_fill_manual(values=c("ATM+"="#C0392B","ATM−"="#BDC3C7"),
                    name=NULL) +
  scale_color_identity() +
  scale_x_discrete(labels=SHORT_LBL) +
  scale_y_continuous(limits=c(-5, 100), breaks=seq(0,100,20),
                     labels=function(x) ifelse(x<0,"",paste0(x,"%"))) +
  PUB +
  theme(legend.position="top", legend.key.size=unit(0.5,"cm"),
        axis.text.x=element_text(size=10, lineheight=0.9)) +
  labs(x=NULL, y="TP53 mutation rate",
       title="TP53 mutation rate: ATM-mutant vs ATM-wild-type",
       subtitle="Error bars: Wilson 95% CI  |  * Fisher p<0.05  ** p<0.01  *** p<0.001")

ggsave(file.path(OUT,"figures","fig_tp53rate_by_ATM.png"),
       p_tp53rate, width=9, height=5.5, dpi=250, bg="white")
cat("  Saved: fig_tp53rate_by_ATM.png\n")

# ── Fig B: 2×2 Contingency Tiles ──────────────────────────────────────────────
tile2_rows <- lapply(main_cohs, function(co) {
  df <- as.data.frame(co$df)
  n  <- nrow(df)
  p_a <- mean(df$ATM); p_t <- mean(df$TP53)
  cells <- expand.grid(ATM=c(1,0), TP53=c(1,0))
  cells$obs <- c(
    sum(df$ATM==1&df$TP53==1), sum(df$ATM==0&df$TP53==1),
    sum(df$ATM==1&df$TP53==0), sum(df$ATM==0&df$TP53==0)
  )
  cells$exp <- c(p_a*p_t*n, (1-p_a)*p_t*n, p_a*(1-p_t)*n, (1-p_a)*(1-p_t)*n)
  cells$oe  <- cells$obs / cells$exp
  cells$cohort   <- co$cohort
  cells$short_lbl <- SHORT_LBL[co$cohort]
  cells
})
tile2_df <- do.call(rbind, tile2_rows)
tile2_df$cohort_f  <- factor(tile2_df$cohort, levels=COH_ORDER[COH_ORDER %in% tile2_df$cohort])
tile2_df$ATM_lbl   <- ifelse(tile2_df$ATM==1,"ATM+","ATM−")
tile2_df$TP53_lbl  <- ifelse(tile2_df$TP53==1,"TP53+","TP53−")
tile2_df$ATM_lbl   <- factor(tile2_df$ATM_lbl, levels=c("ATM+","ATM−"))
tile2_df$TP53_lbl  <- factor(tile2_df$TP53_lbl, levels=c("TP53+","TP53−"))
tile2_df$cell_lbl  <- with(tile2_df,
  sprintf("n=%d\nO/E=%.2f", obs, oe))
tile2_df$log2oe    <- log2(pmax(pmin(tile2_df$oe, 4), 0.25))
tile2_df$txt_col   <- ifelse(abs(tile2_df$log2oe) > 1, "white", "grey20")

p_2x2 <- ggplot(tile2_df, aes(x=TP53_lbl, y=ATM_lbl, fill=log2oe)) +
  geom_tile(color="white", linewidth=1.2) +
  geom_text(aes(label=cell_lbl, color=txt_col), size=3.1, lineheight=1.2) +
  facet_wrap(~cohort_f, nrow=1,
             labeller=labeller(cohort_f=setNames(SHORT_LBL, names(SHORT_LBL)))) +
  scale_fill_gradient2(
    low="#C0392B", mid="#FDFEFE", high="#2980B9", midpoint=0,
    limits=c(-2,2), oob=squish, na.value="grey85",
    name="log₂(O/E)",
    breaks=c(-2,-1,0,1,2),
    labels=c("≤−2\n(ME)","−1","1:1","+1","≥+2")
  ) +
  scale_color_identity() +
  PUB +
  theme(
    strip.text       = element_text(face="bold", size=9.5),
    strip.background = element_rect(fill="grey93", color=NA),
    axis.text        = element_text(size=9),
    axis.ticks       = element_blank(),
    axis.line        = element_blank(),
    panel.border     = element_rect(color="grey60", fill=NA, linewidth=0.5),
    legend.position  = "right"
  ) +
  labs(x=NULL, y=NULL,
       title="ATM–TP53 contingency: observed counts and O/E ratio",
       subtitle="Red = fewer co-mutations than expected (mutual exclusivity)  |  Blue = co-occurring")

ggsave(file.path(OUT,"figures","fig_2x2_tiles.png"),
       p_2x2, width=13, height=3.8, dpi=250, bg="white")
cat("  Saved: fig_2x2_tiles.png\n")

# ── Fig C: Observed vs Expected co-mutations ───────────────────────────────────
oe_df <- sum_df
oe_df$cohort_f <- factor(oe_df$cohort, levels=COH_ORDER[COH_ORDER %in% oe_df$cohort])

# axis range
max_val <- max(c(oe_df$obs_both, oe_df$exp_both)) * 1.25

p_oe <- ggplot(oe_df, aes(x=exp_both, y=obs_both)) +
  # H0 reference line y = x
  geom_abline(slope=1, intercept=0, linetype="dashed",
              color="grey50", linewidth=0.8) +
  # Shaded ME region (below y=x)
  annotate("text", x=max_val*0.72, y=max_val*0.06,
           label="Fewer than\nexpected\n(ME)", color="#C0392B",
           size=3.5, hjust=0.5, fontface="italic") +
  annotate("text", x=max_val*0.12, y=max_val*0.88,
           label="More than\nexpected\n(co-occurrence)", color="#2980B9",
           size=3.5, hjust=0.5, fontface="italic") +
  # Bubbles
  geom_point(aes(size=n, fill=sig_lbl), shape=21, color="grey20",
             stroke=0.8, alpha=0.9) +
  # Labels
  ggrepel::geom_label_repel(
    aes(label=paste0(SHORT_LBL[cohort], "\nObs=", obs_both,
                     "\nExp=", sprintf("%.1f", exp_both))),
    size=2.9, lineheight=0.9, label.size=0.2, label.padding=unit(0.12,"lines"),
    box.padding=0.5, seed=42, fill="white", color="grey20"
  ) +
  scale_size_continuous(range=c(5,16), name="Total N",
                        breaks=c(100,200,350,450)) +
  scale_fill_manual(
    values=c("***"="#C0392B","**"="#E74C3C","*"="#F39C12","ns"="#AEB6BF"),
    name="Fisher p", labels=c("***"="<0.001","**"="<0.01","*"="<0.05","ns"="≥0.05")
  ) +
  scale_x_continuous(limits=c(0, max_val), expand=c(0.02,0)) +
  scale_y_continuous(limits=c(0, max_val), expand=c(0.02,0)) +
  PUB +
  theme(legend.position="right",
        panel.grid.major=element_line(color="grey93", linewidth=0.4)) +
  labs(x="Expected co-mutations (ATM+ × TP53+ / N)",
       y="Observed co-mutations",
       title="ATM–TP53: observed vs expected co-mutations",
       subtitle="Below dashed line (y = x) = fewer co-mutations than expected by chance")

ggsave(file.path(OUT,"figures","fig_obs_vs_exp.png"),
       p_oe, width=7.5, height=6.5, dpi=250, bg="white")
cat("  Saved: fig_obs_vs_exp.png\n")

# ── Fig D: ATM & TP53 mutation frequency per cohort ───────────────────────────
freq_df <- rbind(
  data.frame(cohort=sum_df$cohort, cohort_f=sum_df$cohort_f,
             gene="TP53", pct=sum_df$pct_tp53, n_mut=sum_df$n_tp53, n=sum_df$n,
             stringsAsFactors=FALSE),
  data.frame(cohort=sum_df$cohort, cohort_f=sum_df$cohort_f,
             gene="ATM",  pct=sum_df$pct_atm,  n_mut=sum_df$n_atm,  n=sum_df$n,
             stringsAsFactors=FALSE)
)
freq_df$cohort_f <- factor(freq_df$cohort, levels=rev(COH_ORDER[COH_ORDER %in% freq_df$cohort]))
freq_df$gene     <- factor(freq_df$gene, levels=c("TP53","ATM"))
freq_df$n_lbl    <- paste0(freq_df$n_mut, "/", freq_df$n)

p_freq <- ggplot(freq_df, aes(x=pct, y=cohort_f, fill=gene)) +
  geom_col(position=position_dodge(0.7), width=0.6,
           color="white", linewidth=0.3) +
  geom_text(aes(x=pct+0.8, label=n_lbl),
            position=position_dodge(0.7), hjust=0, size=3.1, color="grey30") +
  scale_fill_manual(values=c("TP53"="#C0392B","ATM"="#2471A3"), name=NULL) +
  scale_x_continuous(limits=c(0,85), labels=function(x) paste0(x,"%"),
                     expand=c(0,0)) +
  scale_y_discrete(labels=SHORT_LBL) +
  PUB +
  theme(legend.position="top",
        panel.grid.major.x=element_line(color="grey90", linewidth=0.3),
        panel.grid.major.y=element_blank()) +
  labs(x="Mutation frequency", y=NULL,
       title="ATM and TP53 mutation frequency per cohort",
       subtitle="Numbers show mutated samples / total  (TCGA: MSS-only stratum)")

ggsave(file.path(OUT,"figures","fig_mutation_freq.png"),
       p_freq, width=7.5, height=4.5, dpi=250, bg="white")
cat("  Saved: fig_mutation_freq.png\n")

# ── Fig E: Combined multi-panel figure ────────────────────────────────────────
p_combined <- (p_tp53rate | p_freq) /
              (p_2x2) /
              (p_oe | p_forest) +
  plot_annotation(tag_levels="A",
                  theme=theme(plot.title=element_text(size=13,face="bold")))

ggsave(file.path(OUT,"figures","fig_combined_validation.png"),
       p_combined, width=18, height=20, dpi=200, bg="white")
cat("  Saved: fig_combined_validation.png\n")

# =============================================================================
# 7. CONSOLE SUMMARY
# =============================================================================
cat(sprintf("\n%s\nSummary\n%s\n",strrep("=",80),strrep("=",80)))
cat(sprintf("%-20s %-6s %5s %5s %5s %5s  %6s  %8s  %8s\n",
            "Cohort","Strat","N","ATM+","TP53+","both","OR","p","Fisher-p"))
cat(strrep("-",90),"\n")
for (i in seq_len(nrow(res))) {
  r <- res[i,]
  cat(sprintf("%-20s %-6s %5d %5d %5d %5d  %6s  %8s  %8s\n",
              r$cohort, r$stratum, r$n, r$n_atm, r$n_tp53, r$n_both,
              ifelse(is.na(r$OR),"—",sprintf("%.3f",r$OR)),
              ifelse(is.na(r$pvalue),"—",sprintf("%.4f",r$pvalue)),
              ifelse(is.na(r$fish_p),"—",sprintf("%.4f",r$fish_p))))
}
if (!is.null(meta_full))
  cat(sprintf("\nPooled (RE, full cohort): OR=%.3f [%.3f–%.3f]  p=%.2e  I²=%.1f%%\n",
              meta_full$pooled_OR, meta_full$OR_lo, meta_full$OR_hi,
              meta_full$p_value, meta_full$I2))
