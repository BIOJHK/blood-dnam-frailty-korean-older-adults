#############################################################################
# 01. DNA methylation preprocessing & EWAS (Illumina EPIC v1)
# Manuscript section: Genome-wide DNA methylation profiling
#
# Pipeline (as described in the Methods)
#   1) IDAT loading (minfi) -> sex-mismatch QC (sample-level detection QC is reported only)
#   2) Blood cell-type proportions: estimateCellCounts2 + EPIC IDOL (FlowSorted.Blood.EPIC)
#   3) Probe filtering: detection P >= 0.05, SNP-overlapping probes, cross-reactive probes
#   4) Primary data: raw beta -> BMIQ -> M-values
#   5) Primary EWAS  : limma, M ~ FI (continuous) + age + sex + cell proportions, BH-FDR
#   6) Secondary EWAS: frailty status (3-level categorical) + same covariates
#   7) Data for prediction / exploratory domain analyses: ssNoob -> BMIQ
#############################################################################

suppressPackageStartupMessages({
  library(minfi)
  library(IlluminaHumanMethylationEPICmanifest)
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  library(FlowSorted.Blood.EPIC)
  library(wateRmelon)   # BMIQ
  library(limma)
  library(readxl)
})


## ===========================================================================
## 0. Settings — edit paths and phenotype column names here only
## ===========================================================================
cfg <- list(
  idat_dir    = "idat_directory",
  sheet_file  = "clinical_sheet_file",
  snp_file    = "filtering_option",    # SNP probe lists (sheets: CpG_rs / SBE_rs / Probe_rs / stem_cg)
  cross_files = "cross_active_probe",  # cross-reactive probe list(s) (column name 'probe')
  out_dir     = "Results_directory",

  ## QC criteria — as stated in the Methods
  detp_cut         = 0.05,  # Methods: probes with detection P >= 0.05 are excluded
  probe_fail_rule  = "any", # a probe is excluded if detection P >= 0.05 in any sample
  snp_maf          = 0,     # Methods: probes overlapping known SNPs excluded (dropLociWithSnps, MAF > 0)

  ## Steps not described in the Methods — OFF by default (results stay as reported)
  sample_detp_drop = FALSE, # TRUE: drop samples with failed-probe fraction > sample_fail_frac (OFF: report only)
  sample_fail_frac = 0.20,
  drop_XY          = FALSE, # TRUE: remove all probes on chrX / chrY

  ## Phenotype column names (sample sheet / pData)
  col_FI      = "FI",
  col_age     = "Age",
  col_sex     = "Sex",               # 1 = M, 0 = F
  col_frail   = "Frailty_status",
  frail_levels = c("Non-frail", "Pre-frail", "frail"),   # first level = reference

  ## Cell-type proportions
  cell_types = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu"),
  cell_ref   = "Neu",   # proportions sum to ~1, so one is left out of the model (avoids collinearity)

  beta_eps = 1e-6       # beta bounded to [eps, 1 - eps] for BMIQ / logit transformation
)
dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)
out <- function(f) file.path(cfg$out_dir, f)

message("minfi version: ", packageVersion("minfi"))  # check against the version reported in the Methods


## ===========================================================================
## 1. Sample sheet & IDAT loading
## ===========================================================================
chip <- as.data.frame(read_xlsx(cfg$sheet_file, skip = 6))
colnames(chip) <- as.character(chip[1, ])   # sheet layout: the row after row 7 holds the real header
chip <- chip[-1, , drop = FALSE]
chip$Basename <- file.path(cfg$idat_dir,
                           paste0(chip$Sentrix_ID, "_", chip$Sentrix_Position))

RGSet <- read.metharray.exp(targets = chip, force = TRUE)
controlStripPlot(RGSet, controls = "BISULFITE CONVERSION II")
message("Loaded samples: ", ncol(RGSet))


## ===========================================================================
## 2. Sample QC
## ===========================================================================
## 2-1. Sex mismatch (predicted vs recorded sex) -----------------------------
pred_sex     <- getSex(mapToGenome(preprocessRaw(RGSet)))
recorded_sex <- factor(pData(RGSet)[[cfg$col_sex]], levels = c("1", "0"), labels = c("M", "F"))

sex_tab <- data.frame(Sample_Name  = pData(RGSet)$Sample_Name,
                      recorded     = recorded_sex,
                      predicted    = pred_sex$predictedSex,
                      xMed = pred_sex$xMed, yMed = pred_sex$yMed,
                      row.names    = colnames(RGSet))
sex_tab$match <- !is.na(sex_tab$recorded) &
                 as.character(sex_tab$recorded) == sex_tab$predicted   # missing recorded sex -> treated as mismatch
write.table(sex_tab, out("QC_sex_check.txt"), sep = "\t", quote = FALSE, col.names = NA)

RGSet_qc <- RGSet[, sex_tab$match]
message("Sex QC: ", ncol(RGSet), " -> ", ncol(RGSet_qc), " samples")

## 2-2. Sample-level detection P ---------------------------------------------
detP   <- detectionP(RGSet_qc, type = "m+u")
failed <- detP >= cfg$detp_cut
samp_fail <- colMeans(failed)
bad_samples <- names(samp_fail)[samp_fail > cfg$sample_fail_frac]
write.table(data.frame(Sample_Name = pData(RGSet_qc)$Sample_Name, frac_failed = samp_fail),
            out("QC_sample_detP.txt"), sep = "\t", quote = FALSE, col.names = NA)
if (length(bad_samples) > 0 && !cfg$sample_detp_drop) {
  warning(length(bad_samples), " sample(s) exceed ", cfg$sample_fail_frac,
          " failed-probe fraction but are KEPT (as in the Methods): ", paste(bad_samples, collapse = ", "))
}
if (length(bad_samples) > 0 && cfg$sample_detp_drop) {
  message("Samples removed (detP): ", paste(bad_samples, collapse = ", "))
  RGSet_qc <- RGSet_qc[, !colnames(RGSet_qc) %in% bad_samples]
  detP     <- detP[, colnames(RGSet_qc)]
  failed   <- failed[, colnames(RGSet_qc)]
}
message("Final analytic samples: ", ncol(RGSet_qc))


## ===========================================================================
## 3. Blood cell-type proportions (estimateCellCounts2, EPIC IDOL reference)
## ===========================================================================
data(IDOLOptimizedCpGs)

cc <- estimateCellCounts2(
  RGSet_qc,                                  # un-normalized RGChannelSet
  compositeCellType = "Blood",
  processMethod     = "preprocessNoob",
  probeSelect       = "IDOL",
  cellTypes         = cfg$cell_types,
  referencePlatform = "IlluminaHumanMethylationEPIC",
  IDOLOptimizedCpGs = IDOLOptimizedCpGs,     # IDOL CpG set for EPIC
  returnAll         = FALSE,
  meanPlot          = FALSE
)
cell_prop <- as.data.frame(cc$prop)[colnames(RGSet_qc), cfg$cell_types]

write.table(data.frame(Sample_Name = pData(RGSet_qc)$Sample_Name, cell_prop),
            out("cell_proportions_estimateCellCounts2_IDOL.txt"),
            sep = "\t", quote = FALSE, col.names = NA)


## ===========================================================================
## 4. Probe filtering — one keep_probes list, applied to both datasets
## ===========================================================================
MSet_raw   <- preprocessRaw(RGSet_qc)
anno       <- getAnnotation(MSet_raw)
probe_type <- setNames(anno$Type, rownames(anno))   # "I" / "II"
all_probes <- rownames(MSet_raw)

read_ids <- function(file, sheet = NULL, col = NULL, col_names = TRUE) {
  x <- as.data.frame(read_xlsx(file, sheet = sheet, col_names = col_names))
  unique(na.omit(as.character(if (is.null(col)) x[[1]] else x[[col]])))
}

## 4-1. Detection P ----------------------------------------------------------
detP <- detP[all_probes, colnames(MSet_raw)]
fail_detP <- if (identical(cfg$probe_fail_rule, "any")) {
  rowSums(detP >= cfg$detp_cut) > 0
} else {
  rowMeans(detP >= cfg$detp_cut) > cfg$probe_fail_rule
}

## 4-2. SNP-overlapping probes: published lists + minfi dropLociWithSnps -----
snp_list <- unique(c(
  read_ids(cfg$snp_file, sheet = "CpG_rs",   col = "cgID"),
  read_ids(cfg$snp_file, sheet = "SBE_rs",   col = "cgID"),
  read_ids(cfg$snp_file, sheet = "Probe_rs", col = "cgID"),
  read_ids(cfg$snp_file, sheet = "stem_cg",  col_names = FALSE)
))
gmset     <- mapToGenome(MSet_raw)
snp_minfi <- setdiff(rownames(gmset),
                     rownames(dropLociWithSnps(gmset, snps = c("CpG", "SBE"), maf = cfg$snp_maf)))

## 4-3. Cross-reactive probes ------------------------------------------------
cross_reactive <- unique(unlist(lapply(cfg$cross_files, read_ids, col = "probe")))

## 4-4. (Optional) chrX / chrY probes ----------------------------------------
xy_probes <- if (cfg$drop_XY) rownames(anno)[anno$chr %in% c("chrX", "chrY")] else character(0)

## Apply sequentially and log the number removed at each step (for Methods / supplementary table)
filters <- list(
  detection_P    = all_probes[fail_detP],
  SNP_list       = snp_list,
  SNP_minfi      = snp_minfi,
  cross_reactive = cross_reactive,
  sex_chromosome = xy_probes
)
keep_probes <- all_probes
filter_log  <- data.frame(step = "start", removed = 0, remaining = length(keep_probes))
for (nm in names(filters)) {
  before      <- length(keep_probes)
  keep_probes <- setdiff(keep_probes, filters[[nm]])
  filter_log  <- rbind(filter_log, data.frame(step = nm,
                                              removed   = before - length(keep_probes),
                                              remaining = length(keep_probes)))
}
print(filter_log)
write.table(filter_log, out("QC_probe_filtering.txt"), sep = "\t", quote = FALSE, row.names = FALSE)


## ===========================================================================
## 5. BMIQ normalization (Type I / II probe correction)
## ===========================================================================
run_bmiq <- function(beta, probe_type) {
  design <- ifelse(probe_type[rownames(beta)] == "I", 1L, 2L)
  stopifnot(!anyNA(design))
  if (anyNA(beta)) {
    na_rows <- rowSums(is.na(beta)) > 0
    message("Dropping ", sum(na_rows), " probes with NA beta before BMIQ")
    beta <- beta[!na_rows, , drop = FALSE]; design <- design[!na_rows]
  }
  norm <- vapply(seq_len(ncol(beta)), function(i) {
    set.seed(2026)
    wateRmelon::BMIQ(beta[, i], design, plots = FALSE, pri = FALSE, sampleID = i)$nbeta
  }, numeric(nrow(beta)))
  dimnames(norm) <- dimnames(beta)
  pmin(pmax(norm, cfg$beta_eps), 1 - cfg$beta_eps)
}
beta2m <- function(b) log2(b / (1 - b))


## ===========================================================================
## 6. Primary EWAS data: raw beta (M / (M + U)) -> BMIQ -> M-values
## ===========================================================================
beta_raw     <- getBeta(MSet_raw[keep_probes, ], offset = 0, betaThreshold = cfg$beta_eps)
beta_primary <- run_bmiq(beta_raw, probe_type)
M_primary    <- beta2m(beta_primary)

saveRDS(beta_primary, out("beta_raw_BMIQ_filtered.rds"))
saveRDS(M_primary,    out("M_raw_BMIQ_filtered.rds"))


## ===========================================================================
## 7. Phenotypes + covariates
## ===========================================================================
pheno <- as.data.frame(pData(RGSet_qc))

## Check that frailty labels match frail_levels exactly (a case mismatch would silently become NA)
frail_raw <- unique(na.omit(as.character(pheno[[cfg$col_frail]])))
if (length(setdiff(frail_raw, cfg$frail_levels)) > 0) {
  stop("Frailty values do not match frail_levels: ",
       paste(setdiff(frail_raw, cfg$frail_levels), collapse = ", "),
       "  (labels in data: ", paste(frail_raw, collapse = ", "), ")")
}

pheno <- data.frame(
  Sample_Name = pheno$Sample_Name,
  FI      = as.numeric(pheno[[cfg$col_FI]]),
  Age     = as.numeric(pheno[[cfg$col_age]]),
  Sex     = factor(pheno[[cfg$col_sex]], levels = c("0", "1"), labels = c("F", "M")),
  Frailty = factor(pheno[[cfg$col_frail]], levels = cfg$frail_levels),
  cell_prop[rownames(pheno), ],
  row.names = rownames(pheno)
)
stopifnot(identical(rownames(pheno), colnames(M_primary)))

cell_cov <- setdiff(cfg$cell_types, cfg$cell_ref)   # Neu = reference (left out of the model)


## ===========================================================================
## 8. EWAS functions (limma moderated statistics + BH-FDR)
## ===========================================================================
ann_cols <- c("chr", "pos", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_Island")

genomic_lambda <- function(p) median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1)

run_ewas <- function(M, pheno, formula) {
  vars <- all.vars(formula)
  ok   <- complete.cases(pheno[, vars])
  message(sum(ok), " / ", nrow(pheno), " samples with complete data for: ", deparse(formula))
  design <- model.matrix(formula, data = pheno[ok, ])
  fit    <- eBayes(lmFit(M[, ok], design))
  list(fit = fit, design = design, n = sum(ok))
}

extract_tab <- function(fit, coef) {
  tt <- topTable(fit, coef = coef, number = Inf, adjust.method = "BH", sort.by = "P")
  cbind(CpG = rownames(tt), tt, as.data.frame(anno[rownames(tt), ann_cols]))
}


## ===========================================================================
## 9. Primary EWAS: M ~ FI + Age + Sex + cell proportions
## ===========================================================================
f_primary <- as.formula(paste("~ FI + Age + Sex +", paste(cell_cov, collapse = " + ")))
ewas1     <- run_ewas(M_primary, pheno, f_primary)
res_FI    <- extract_tab(ewas1$fit, "FI")

message(sprintf("Primary EWAS (n=%d): FDR<0.05 = %d, FDR<0.1 = %d CpGs, lambda = %.3f",
                ewas1$n, sum(res_FI$adj.P.Val < 0.05), sum(res_FI$adj.P.Val < 0.1),
                genomic_lambda(res_FI$P.Value)))
write.table(res_FI, out("EWAS_primary_FI_continuous.txt"), sep = "\t", quote = FALSE, row.names = FALSE)


## ===========================================================================
## 10. Secondary EWAS: frailty status (3-level categorical) + same covariates
## ===========================================================================
f_secondary <- as.formula(paste("~ Frailty + Age + Sex +", paste(cell_cov, collapse = " + ")))
ewas2       <- run_ewas(M_primary, pheno, f_secondary)
frail_coefs <- grep("^Frailty", colnames(ewas2$design), value = TRUE)   # vs reference (Non-frail)

for (cf in frail_coefs) {
  tab <- extract_tab(ewas2$fit, cf)
  message(sprintf("Secondary EWAS %s: FDR<0.05 = %d, FDR<0.1 = %d CpGs",
                  cf, sum(tab$adj.P.Val < 0.05), sum(tab$adj.P.Val < 0.1)))
  write.table(tab, out(paste0("EWAS_secondary_", cf, "_vs_", cfg$frail_levels[1], ".txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
}
## Overall difference across the 3 groups (moderated F-test)
res_F <- extract_tab(ewas2$fit, frail_coefs)
write.table(res_F, out("EWAS_secondary_Frailty_overall_Ftest.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)


## ===========================================================================
## 11. Data for prediction / exploratory domain analyses: ssNoob -> BMIQ
##     (same samples, same keep_probes)
## ===========================================================================
MSet_ssnoob <- preprocessNoob(RGSet_qc, dyeMethod = "single")   # single-sample Noob
beta_ssnoob <- getBeta(MSet_ssnoob, betaThreshold = cfg$beta_eps)[keep_probes, colnames(RGSet_qc)]
beta_pred   <- run_bmiq(beta_ssnoob, probe_type)

saveRDS(beta_pred,         out("beta_ssNoob_BMIQ_filtered.rds"))
saveRDS(beta2m(beta_pred), out("M_ssNoob_BMIQ_filtered.rds"))
saveRDS(pheno,             out("pheno_with_cellprop.rds"))
## -> Downstream nested 5-fold CV marker discovery / MSE, RMSE, MAE use beta_pred (or M-values) and pheno.
##    (Note: CpG selection inside CV must be done within the training fold only to avoid information leakage.)
