#############################################################################
# 07. FI prediction: nested cross-validation (limma screening + elastic net)
# Manuscript section: FI prediction (nested cross-validation)
#
# Uses the generic implementation in nested_cv/nested_cv.R without modification
#   (Zenodo record 24211256, v1.0.0, CC BY 4.0; see nested_cv/README.md for the method).
#   - Within every training fold: missingness / variance filtering, median imputation,
#     limma screening of the top_n CpGs (M ~ FI + covariates), elastic net on beta values
#   - Outer test samples are never used for screening or tuning
#   - Metrics: Pearson R, RMSE, MAE, predictive R2, CCC (out-of-fold predictions)
#
# Inputs: outputs of 01_preprocessing_EWAS.R (out_dir)
#   - beta_ssNoob_BMIQ_filtered.rds : ssNoob -> BMIQ beta, all QC-passed CpGs (CpG x sample)
#   - pheno_with_cellprop.rds       : FI, Age, Sex (F/M), cell proportions
#
# Run from the repository root (the path to nested_cv.R is relative to it).
#############################################################################

## ===========================================================================
## 0. Settings — edit paths and CV settings here only
## ===========================================================================
v2_dir <- "Results_directory"   # out_dir of 01_preprocessing_EWAS.R

cfg <- list(
  beta_file      = file.path(v2_dir, "beta_ssNoob_BMIQ_filtered.rds"),
  pheno_file     = file.path(v2_dir, "pheno_with_cellprop.rds"),
  nested_cv_file = file.path("nested_cv", "nested_cv.R"),
  out_dir        = file.path(v2_dir, "nested_cv_FI"),   # must not exist yet (no overwriting)

  outcome     = "FI",
  ## Covariates adjust limma screening only; they are not glmnet predictors.
  ## Neu is the reference cell type (proportions sum to ~1), as in the EWAS.
  covariates  = c("Sex", "Age", "CD8T", "CD4T", "NK", "Bcell", "Mono"),
  categorical = "Sex",
  balance     = c("Sex", "Age"),   # fold allocation balances FI, sex and age

  ## Match these to the Methods
  top_n       = 800L,
  outer_folds = 10L,
  inner_folds = 5L,
  repeats     = 1L,
  seed        = 1L
)

source(cfg$nested_cv_file)


## ===========================================================================
## 1. Read data & build metadata
## ===========================================================================
beta  <- readRDS(cfg$beta_file)        # CpG x sample, values in [0, 1]
pheno <- readRDS(cfg$pheno_file)
stopifnot(identical(rownames(pheno), colnames(beta)))

meta <- data.frame(sample_id = rownames(pheno),
                   pheno[, c(cfg$outcome, cfg$covariates)],
                   check.names = FALSE, stringsAsFactors = FALSE)

## nested_cv.R rejects missing outcomes / covariates instead of dropping samples silently
incomplete <- meta$sample_id[!complete.cases(meta)]
if (length(incomplete))
  stop(length(incomplete), " sample(s) with missing FI / covariates: ", paste(incomplete, collapse = ", "),
       "  -> decide on exclusion explicitly and report it in the Methods")
message("Samples: ", nrow(meta), " | CpGs: ", nrow(beta))


## ===========================================================================
## 2. Nested CV
## ===========================================================================
result <- run_nested_cv(
  beta, meta,
  outcome     = cfg$outcome,
  covariates  = cfg$covariates,
  categorical = cfg$categorical,
  balance     = cfg$balance,
  top_n       = cfg$top_n,
  outer_folds = cfg$outer_folds,
  inner_folds = cfg$inner_folds,
  repeats     = cfg$repeats,
  seed        = cfg$seed,
  output_dir  = cfg$out_dir
)

print(result$metrics)
print(result$model_summary[, c("outer_fold", "alpha", "lambda_fraction", "n_nonzero",
                               "fraction_at_grid_boundary")])
message("CpGs selected in >= 50% of outer models (exploratory): ",
        nrow(result$exploratory_markers_ge_50pct))
