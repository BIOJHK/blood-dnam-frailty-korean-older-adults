#############################################################################
# 02. FI EWAS sensitivity analysis (additional covariates: smoking / diabetes / polypharmacy)
# Manuscript section: Sensitivity analyses
#
# Inputs: outputs of 01_preprocessing_EWAS.R (out_dir)
#   - M_raw_BMIQ_filtered.rds                        : primary EWAS M-values (raw -> BMIQ)
#   - cell_proportions_estimateCellCounts2_IDOL.txt  : estimateCellCounts2 + EPIC IDOL cell proportions
#   - EWAS_primary_FI_continuous.txt                 : primary EWAS results -> reference CpG set
#
# Model: M ~ FI + Age + Sex + CD8T + CD4T + NK + Bcell + Mono (+ additional covariate)
#        (Neu = reference; same covariate set as the primary EWAS in 01)
#############################################################################

suppressPackageStartupMessages({
  library(openxlsx)
  library(limma)
})


## ===========================================================================
## 0. Settings — edit paths here only
## ===========================================================================
v2_dir <- "Results_directory"   # out_dir of 01_preprocessing_EWAS.R

cfg <- list(
  mval_file     = file.path(v2_dir, "M_raw_BMIQ_filtered.rds"),
  cell_file     = file.path(v2_dir, "cell_proportions_estimateCellCounts2_IDOL.txt"),
  primary_file  = file.path(v2_dir, "EWAS_primary_FI_continuous.txt"),
  clinical_file = "clinical_information_file",   # clinical information (xlsx)
  paper_rds     = "paper_sensitivity_rds",       # submitted sensitivity results — read-only; comparison skipped if absent
  out_dir       = file.path(v2_dir, "sensitivity"),

  id_col   = "Sentrix_id",                         # sample ID column in the clinical sheet (= M-value column names)
  fdr_cut  = 0.1,                                  # significance & reference CpG threshold
  base_cov = c("Age", "Sex"),
  cell_cov = c("CD8T", "CD4T", "NK", "Bcell", "Mono"),

  ## Additional covariates: na_codes = codes treated as missing, as_factor = categorical or not
  sens = list(
    smoking      = list(var = "Smoking",      na_codes = c("", "9"),   as_factor = TRUE),
    diabetes     = list(var = "Diabetes",     na_codes = character(0), as_factor = FALSE),
    polypharmacy = list(var = "Polypharmacy", na_codes = character(0), as_factor = FALSE)
  )
)
dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)
out <- function(f) file.path(cfg$out_dir, f)


## ===========================================================================
## 1. Read & align data
## ===========================================================================
Mval  <- readRDS(cfg$mval_file)
cell  <- read.delim(cfg$cell_file, row.names = 1, check.names = FALSE)
pheno <- read.xlsx(cfg$clinical_file)
pheno[[cfg$id_col]] <- trimws(as.character(pheno[[cfg$id_col]]))

idx <- match(colnames(Mval), pheno[[cfg$id_col]])
if (anyNA(idx)) stop("Samples missing from clinical information: ", paste(colnames(Mval)[is.na(idx)], collapse = ", "))
ph <- pheno[idx, ]
rownames(ph) <- ph[[cfg$id_col]]

missing_cell <- setdiff(colnames(Mval), rownames(cell))
if (length(missing_cell)) stop("Samples missing from cell proportion file: ", paste(missing_cell, collapse = ", "))

## If the clinical sheet already has cell-proportion columns, keep them as *_clinical and replace with IDOL values
old_cols <- intersect(cfg$cell_cov, colnames(ph))
if (length(old_cols)) {
  names(ph)[match(old_cols, names(ph))] <- paste0(old_cols, "_clinical")
  ## Check whether the existing values come from the same source as IDOL (r ~ 1 means identical)
  print(round(sapply(old_cols, function(k)
    cor(as.numeric(ph[[paste0(k, "_clinical")]]), cell[rownames(ph), k], use = "complete.obs")), 3))
}
ph[cfg$cell_cov] <- cell[rownames(ph), cfg$cell_cov]

base_vars <- c("FI", cfg$base_cov, cfg$cell_cov)
stopifnot(identical(rownames(ph), colnames(Mval)),
          !anyNA(ph[base_vars]))


## ===========================================================================
## 2. Helper functions
## ===========================================================================
fit_fi <- function(M, ph, extra = NULL) {
  d <- model.matrix(reformulate(c(base_vars, extra)), data = ph)
  stopifnot(nrow(d) == ncol(M))
  topTable(eBayes(lmFit(M, d)), coef = "FI", number = Inf,
           adjust.method = "BH", sort.by = "none")
}

clean_var <- function(x, na_codes, as_factor) {
  x <- trimws(as.character(x))
  x[x %in% na_codes] <- NA
  if (as_factor) factor(x) else type.convert(x, as.is = TRUE)
}

## r: results to compare, b: reference results, ids: CpGs to compare
conc_n <- function(r, b, ids) {
  ids <- intersect(ids, intersect(rownames(r), rownames(b)))
  n   <- length(ids)
  fmt <- function(k) sprintf("%d/%d (%.1f%%)", k, n, 100 * k / n)
  data.frame(
    n_CpG          = n,
    same_direction = fmt(sum(sign(r[ids, "logFC"]) == sign(b[ids, "logFC"]), na.rm = TRUE)),
    nomP05         = fmt(sum(r[ids, "P.Value"]   < 0.05,        na.rm = TRUE)),
    nomP001        = fmt(sum(r[ids, "P.Value"]   < 0.001,       na.rm = TRUE)),
    FDR            = fmt(sum(r[ids, "adj.P.Val"] < cfg$fdr_cut, na.rm = TRUE)),
    logFC_ratio    = round(median(r[ids, "logFC"] / b[ids, "logFC"], na.rm = TRUE), 3),
    cor_t          = round(cor(r[ids, "t"], b[ids, "t"], use = "complete.obs"), 3)
  )
}


## ===========================================================================
## 3. Model fitting: base (all samples) + each additional-covariate model & base on the same samples
## ===========================================================================
res    <- list(base = fit_fi(Mval, ph))
n_used <- c(base = ncol(Mval))

for (m in names(cfg$sens)) {
  s <- cfg$sens[[m]]
  ph[[s$var]] <- clean_var(ph[[s$var]], s$na_codes, s$as_factor)
  keep <- !is.na(ph[[s$var]])
  message(m, ": ", sum(keep), " / ", length(keep), " samples")
  print(table(ph[[s$var]], useNA = "ifany"))

  res[[paste0("base_", m)]] <- fit_fi(Mval[, keep], ph[keep, ])           # base model, same samples
  res[[m]]                  <- fit_fi(Mval[, keep], ph[keep, ], s$var)    # additional-covariate model
  n_used[m] <- sum(keep)
}


## ===========================================================================
## 4. Number of significant CpGs per model (FDR < fdr_cut)
## ===========================================================================
n_sig <- data.frame(model   = names(res),
                    n_sig   = sapply(res, function(t) sum(t$adj.P.Val < cfg$fdr_cut, na.rm = TRUE)),
                    row.names = NULL)
print(n_sig)


## ===========================================================================
## 5. Reference CpG set: FDR < fdr_cut in the primary EWAS (01)
## ===========================================================================
primary <- read.delim(cfg$primary_file, check.names = FALSE)
rownames(primary) <- primary$CpG
ids <- primary$CpG[primary$adj.P.Val < cfg$fdr_cut]
message("Reference CpGs (primary EWAS FDR < ", cfg$fdr_cut, "): ", length(ids),
        "  <- check against the number reported in the manuscript (930)")


## ===========================================================================
## 6. Concordance
## ===========================================================================
sens_models <- names(cfg$sens)

## 6-1. Each model vs the base model on the same samples
summ_vs_base <- do.call(rbind, c(
  list(cbind(model = "base", conc_n(res$base, res$base, ids))),
  lapply(sens_models, function(m)
    cbind(model = m, conc_n(res[[m]], res[[paste0("base_", m)]], ids)))
))

## 6-2. Each model vs the primary EWAS (the base row should have cor_t ~ 1: reproducibility check)
summ_vs_primary <- do.call(rbind, lapply(c("base", sens_models), function(m)
  cbind(model = m, conc_n(res[[m]], primary, ids))))

print(summ_vs_base)
print(summ_vs_primary)

## 6-3. (Optional) comparison with the submitted sensitivity results — the file is only read
if (file.exists(cfg$paper_rds)) {
  paper  <- readRDS(cfg$paper_rds)
  common <- intersect(sens_models, names(paper))

  summ_paper_vs_primary <- do.call(rbind, lapply(common, function(m)
    cbind(model = m, conc_n(paper[[m]], primary, ids))))
  summ_rerun_vs_paper <- do.call(rbind, lapply(common, function(m)
    cbind(model = m, conc_n(res[[m]], paper[[m]], ids))))
  n_sig_paper <- sapply(paper[common], function(t) sum(t$adj.P.Val < cfg$fdr_cut, na.rm = TRUE))

  print(summ_paper_vs_primary)
  print(summ_rerun_vs_paper)
  print(n_sig_paper)

  write.table(summ_paper_vs_primary, out("concordance_paper_vs_primary.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(summ_rerun_vs_paper, out("concordance_rerun_vs_paper.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
} else {
  message("paper_rds not found -> comparison with submitted results skipped")
}


## ===========================================================================
## 7. Save (never overwrite the submitted rds)
## ===========================================================================
res_file <- out("limma_FI_sensitivity_IDOL.rds")
if (file.exists(cfg$paper_rds) &&
    normalizePath(res_file, mustWork = FALSE) == normalizePath(cfg$paper_rds)) {
  stop("Output file path is identical to the submitted rds. Change out_dir.")
}
saveRDS(res, res_file)

write.table(data.frame(n_sig, n_samples = n_used[sub("^base_", "", n_sig$model)]),
            out("n_significant.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summ_vs_base,    out("concordance_vs_samesample_base.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summ_vs_primary, out("concordance_vs_primary.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
