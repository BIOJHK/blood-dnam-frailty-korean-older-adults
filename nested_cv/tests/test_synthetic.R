# Artificial in-memory data only; no actual study data is read.
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- dirname(dirname(sub("^--file=", "", script_arg[1L])))
source(file.path(root, "nested_cv.R"))
check_packages()
set.seed(871)
n <- 100L; p <- 900L
meta <- data.frame(sample_id = sprintf("synthetic_%03d", seq_len(n)),
  outcome = rnorm(n), sex = rep(c("A", "B"), n / 2), age = runif(n, 30, 80),
  cell_A = runif(n, .1, .3), cell_B = runif(n, .2, .4))
latent <- matrix(rnorm(p * n, sd = .8), p, n)
latent[1:20, ] <- latent[1:20, ] + rep(meta$outcome * .8, each = 20L)
beta <- plogis(latent)
dimnames(beta) <- list(sprintf("probe_%04d", seq_len(p)), meta$sample_id)
beta[30, 1] <- NA_real_
beta[31, ] <- .4

# A held-out perturbation must not affect preprocessing or screening.
checked <- validate_inputs(beta, meta, "outcome", c("sex", "age", "cell_A", "cell_B"), "sex")
tr <- seq_len(80L)
a <- prepare_training(beta, tr, .05, 1e-6)
rank_a <- select_markers(a, checked, tr, "outcome", c("sex", "age", "cell_A", "cell_B"), 800L)
changed_beta <- beta
changed_beta[, 81:100] <- matrix(runif(p * 20L), p)
changed_meta <- checked
changed_meta$outcome[81:100] <- 1000
b <- prepare_training(changed_beta, tr, .05, 1e-6)
rank_b <- select_markers(b, changed_meta, tr, "outcome", c("sex", "age", "cell_A", "cell_B"), 800L)
stopifnot(identical(a, b), identical(rank_a, rank_b), !"probe_0031" %in% rank_a$CpG)

out <- tempfile("generic_nested_cv_test_")
res <- run_nested_cv(beta, meta[n:1, ],
  covariates = c("sex", "age", "cell_A", "cell_B"), categorical = "sex",
  top_n = 800L, outer_folds = 10L, inner_folds = 5L,
  balance_candidates = 10L, output_dir = out)
stopifnot(nrow(res$oof_predictions) == n, !anyDuplicated(res$oof_predictions$sample_id),
  all(is.finite(res$oof_predictions$predicted)), nrow(res$model_summary) == 10L,
  all(res$model_summary$n_screened == 800L), all(res$marker_frequency$total_models == 10L),
  all(res$marker_frequency$frequency >= 0 & res$marker_frequency$frequency <= 1),
  all(res$fold_assignments$inner_fold[res$fold_assignments$role == "test"] |> is.na()),
  all(table(res$fold_assignments$sample_id[res$fold_assignments$role == "test"]) == 1L))

# Verify reported metrics independently, and reject ID mismatch / bad beta range.
z <- res$oof_averaged
stopifnot(abs(res$metrics$RMSE - sqrt(mean((z$observed - z$predicted)^2))) < 1e-12)
for (o in seq_len(10L)) {
  grid <- read.csv(file.path(out, sprintf("repeat_01_outer_%02d", o), "inner_tuning.csv"))
  stopifnot(abs(min(grid$inner_MSE) - res$model_summary$inner_MSE[o]) < 1e-12)
}
# Repeated-CV averaging, no-covariate mode and seed reproducibility.
small_run <- function(destination) run_nested_cv(beta, meta,
  covariates = character(), categorical = character(), balance = character(),
  top_n = 30L, outer_folds = 3L, inner_folds = 2L, repeats = 2L,
  alphas = 1, nlambda = 5L, balance_candidates = 2L, output_dir = destination)
small_a <- small_run(tempfile("generic_repeat_a_"))
small_b <- small_run(tempfile("generic_repeat_b_"))
stopifnot(identical(small_a, small_b), nrow(small_a$oof_predictions) == 2L * n,
  nrow(small_a$oof_averaged) == n, nrow(small_a$metrics_by_repeat) == 2L)
expected <- tapply(small_a$oof_predictions$predicted, small_a$oof_predictions$sample_id, mean)
stopifnot(max(abs(small_a$oof_averaged$predicted - expected[small_a$oof_averaged$sample_id])) < 1e-12)
bad <- meta; bad$sample_id[1] <- "unknown"
stopifnot(inherits(try(validate_inputs(beta, bad, "outcome", c("sex", "age"), "sex"), silent = TRUE), "try-error"))
bad_beta <- beta; bad_beta[1, 1] <- 1.2
stopifnot(inherits(try(validate_inputs(bad_beta, meta, "outcome", c("sex", "age"), "sex"), silent = TRUE), "try-error"))
cat("PASS: synthetic nested CV, OOF coverage, training-only screening, metrics and input guards.\n")
cat("Temporary synthetic outputs:", out, "\n")
print(res$metrics)
print(utils::sessionInfo())
