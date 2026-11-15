# Generic nested cross-validation for continuous outcomes from normalized beta values.
# No cohort-specific data, sample identifiers, or configuration is included.

check_packages <- function() {
  for (p in c("limma", "glmnet"))
    if (!requireNamespace(p, quietly = TRUE)) stop("Install required package: ", p)
}

validate_inputs <- function(beta, metadata, outcome, covariates, categorical) {
  if (!is.matrix(beta) || !is.numeric(beta)) stop("beta must be a numeric CpG-by-sample matrix.")
  for (ids in list(rownames(beta), colnames(beta)))
    if (is.null(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids))
      stop("beta row and column names must be nonempty and unique.")
  if (!is.data.frame(metadata) || !"sample_id" %in% names(metadata))
    stop("metadata must contain sample_id.")
  metadata$sample_id <- as.character(metadata$sample_id)
  if (anyNA(metadata$sample_id) || any(!nzchar(metadata$sample_id)) || anyDuplicated(metadata$sample_id))
    stop("metadata sample_id must be nonempty and unique.")
  if (!setequal(colnames(beta), metadata$sample_id)) stop("Sample IDs must match exactly between inputs.")
  metadata <- metadata[match(colnames(beta), metadata$sample_id), , drop = FALSE]
  if (outcome %in% covariates || anyDuplicated(covariates)) stop("Invalid covariate list.")
  required <- c(outcome, covariates)
  if (!all(required %in% names(metadata))) stop("Missing metadata columns.")
  if (!all(categorical %in% covariates)) stop("categorical must be a subset of covariates.")
  if (anyNA(metadata[, required, drop = FALSE])) stop("Missing outcomes/covariates are not supported.")
  if (!is.numeric(metadata[[outcome]]) || any(!is.finite(metadata[[outcome]])) ||
      stats::sd(metadata[[outcome]]) == 0) stop("Outcome must be finite numeric and variable.")
  for (nm in covariates) {
    if (nm %in% categorical) {
      metadata[[nm]] <- factor(metadata[[nm]])
      if (nlevels(metadata[[nm]]) < 2L) stop("Constant categorical covariate: ", nm)
    } else if (!is.numeric(metadata[[nm]]) || any(!is.finite(metadata[[nm]]))) {
      stop("Continuous covariates must be finite numeric; list factors in categorical: ", nm)
    }
  }
  if (any(is.infinite(beta)) || any(beta < 0 | beta > 1, na.rm = TRUE))
    stop("beta values must be within [0,1]; NA values are allowed.")
  metadata
}

# Candidate folds balance sample counts and standardized means/proportions of
# outcome and specified covariates. Only allocation balance, never model error,
# is used to choose folds. Inner folds are built from outer training samples only.
balanced_folds <- function(metadata, outcome, balance, k, candidates = 100L) {
  n <- nrow(metadata)
  if (k < 2L || n < 2L * k) stop("Need at least two observations per fold.")
  cols <- unique(c(outcome, balance))
  blocks <- lapply(cols, function(nm) {
    v <- metadata[[nm]]
    if (is.factor(v)) {
      sapply(levels(v), function(level) as.numeric(v == level))
    } else matrix(v, ncol = 1L)
  })
  z <- do.call(cbind, blocks)
  keep <- apply(z, 2L, stats::sd) > 0
  z <- scale(z[, keep, drop = FALSE])
  best_score <- Inf
  best <- NULL
  for (i in seq_len(candidates)) {
    f <- sample(rep(seq_len(k), length.out = n))
    score <- sum(vapply(seq_len(k), function(j) sum(colMeans(z[f == j, , drop = FALSE])^2), numeric(1)))
    if (score < best_score) { best <- f; best_score <- score }
  }
  best
}

# All missingness filtering, imputation and variance filtering use training data.
prepare_training <- function(beta, train, max_missing, eps) {
  b <- beta[, train, drop = FALSE]
  keep <- rowMeans(is.na(b)) <= max_missing
  b <- b[keep, , drop = FALSE]
  if (!nrow(b)) stop("No probes remain after training-only missingness filtering.")
  med <- apply(b, 1L, stats::median, na.rm = TRUE)
  missing <- which(is.na(b), arr.ind = TRUE)
  if (nrow(missing)) b[missing] <- med[missing[, 1L]]
  keep <- apply(b, 1L, stats::sd) > 1e-12 & is.finite(med)
  b <- b[keep, , drop = FALSE]
  med <- med[keep]
  if (nrow(b) < 2L) stop("Fewer than two variable probes in a training fold.")
  clipped <- pmin(pmax(b, eps), 1 - eps)
  list(beta = b, m = log2(clipped / (1 - clipped)), medians = med)
}

select_markers <- function(prepared, metadata, train, outcome, covariates, top_n) {
  d <- data.frame(.outcome = metadata[[outcome]][train])
  # Synthetic column names allow arbitrary user column names without formula parsing.
  for (j in seq_along(covariates)) d[[paste0(".cov", j)]] <- metadata[[covariates[j]]][train]
  design <- stats::model.matrix(~ ., data = d)
  if (qr(design)$rank < ncol(design))
    stop("Rank-deficient training design: check constant/rare covariates and compositional collinearity.")
  if (nrow(design) - ncol(design) < 3L) stop("Insufficient residual degrees of freedom for limma.")
  fit <- limma::eBayes(limma::lmFit(prepared$m, design), trend = FALSE, robust = FALSE)
  p <- fit$p.value[, ".outcome"]
  ids <- names(p)
  if (is.null(ids)) ids <- rownames(prepared$m)
  ordering <- order(p, ids, na.last = NA)
  if (length(ordering) < top_n)
    stop("Fewer than top_n eligible CpGs in a training fold. Reduce top_n explicitly.")
  i <- head(ordering, top_n)
  data.frame(CpG = ids[i], P_value = unname(p[i]), stringsAsFactors = FALSE)
}

prediction_matrix <- function(beta, samples, markers, medians) {
  b <- beta[markers, samples, drop = FALSE]
  missing <- which(is.na(b), arr.ind = TRUE)
  if (nrow(missing)) b[missing] <- medians[markers][missing[, 1L]]
  t(b)
}

# Each training set estimates its own lambda_max. A common lambda/lambda_max
# grid aligns paths despite fold-specific selected predictors and scaling.
fit_path <- function(x, y, alpha, fractions) {
  pilot <- glmnet::glmnet(x, y, family = "gaussian", alpha = alpha,
                          standardize = TRUE, intercept = TRUE, nlambda = 5L)
  lambda_max <- pilot$lambda[1L]
  if (!is.finite(lambda_max) || lambda_max <= 0) stop("Invalid lambda_max.")
  lambda <- lambda_max * fractions
  fit <- glmnet::glmnet(x, y, family = "gaussian", alpha = alpha,
                        standardize = TRUE, intercept = TRUE, lambda = lambda)
  if (fit$jerr != 0L || length(fit$lambda) != length(lambda)) stop("glmnet path did not fully converge.")
  list(fit = fit, lambda = lambda)
}

regression_metrics <- function(observed, predicted) {
  mse <- mean((observed - predicted)^2)
  vx <- mean((observed - mean(observed))^2)
  vy <- mean((predicted - mean(predicted))^2)
  cv <- mean((observed - mean(observed)) * (predicted - mean(predicted)))
  data.frame(n = length(observed),
    Pearson_R = if (vy > 0) stats::cor(observed, predicted) else NA_real_,
    RMSE = sqrt(mse), MAE = mean(abs(observed - predicted)),
    predictive_R2 = 1 - mse / vx,
    CCC = 2 * cv / (vx + vy + (mean(observed) - mean(predicted))^2))
}

run_nested_cv <- function(beta, metadata, outcome = "outcome",
                          covariates = c("sex", "age"), categorical = "sex",
                          balance = c("sex", "age"), top_n = 800L,
                          outer_folds = 10L, inner_folds = 5L, repeats = 1L,
                          alphas = c(0.05, 0.1, 0.25, 0.5, 0.75, 1),
                          nlambda = 100L, min_fraction = 0.001,
                          seed = 1L, balance_candidates = 100L,
                          max_missing = 0.05, beta_epsilon = 1e-6,
                          coefficient_tolerance = 1e-10,
                          output_dir = "nested_cv_output") {
  check_packages()
  metadata <- validate_inputs(beta, metadata, outcome, covariates, categorical)
  integer_args <- c(top_n, outer_folds, inner_folds, repeats, nlambda, balance_candidates)
  if (any(!is.finite(integer_args)) || any(integer_args != floor(integer_args)) || any(integer_args < 1))
    stop("Count parameters must be positive integers.")
  if (top_n < 2 || outer_folds < 2 || inner_folds < 2 || nlambda < 2) stop("Invalid counts.")
  if (!all(balance %in% covariates)) stop("balance must be a subset of covariates.")
  if (!length(alphas) || any(!is.finite(alphas)) || any(alphas <= 0 | alphas > 1) || anyDuplicated(alphas))
    stop("alphas must be unique values in (0,1]; ridge alpha=0 is not supported.")
  if (!is.finite(min_fraction) || min_fraction <= 0 || min_fraction >= 1 ||
      !is.finite(max_missing) || max_missing < 0 || max_missing >= 1 ||
      !is.finite(beta_epsilon) || beta_epsilon <= 0 || beta_epsilon >= 0.5 ||
      !is.finite(coefficient_tolerance) || coefficient_tolerance < 0) stop("Invalid numeric settings.")
  if (dir.exists(output_dir) || file.exists(output_dir)) stop("Output path already exists; use a new directory.")
  dir.create(output_dir, recursive = TRUE)
  # Disable early path truncation, restoring the user's glmnet controls on exit.
  old_control <- glmnet::glmnet.control()
  on.exit(do.call(glmnet::glmnet.control, old_control), add = TRUE)
  glmnet::glmnet.control(fdev = 0, devmax = 1)
  set.seed(seed)
  fractions <- exp(seq(log(1), log(min_fraction), length.out = nlambda))
  settings <- as.list(environment())
  settings <- settings[c("outcome", "covariates", "categorical", "balance", "top_n", "outer_folds",
                         "inner_folds", "repeats", "alphas", "nlambda", "min_fraction", "seed",
                         "balance_candidates", "max_missing", "beta_epsilon", "coefficient_tolerance")]
  dput(settings, file = file.path(output_dir, "run_settings.R"))
  writeLines(capture.output(utils::sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
  y <- metadata[[outcome]]
  n <- ncol(beta)
  all_predictions <- list()
  all_folds <- list()
  all_selections <- list()
  summaries <- list()
  model_index <- 0L
  for (r in seq_len(repeats)) {
    outer <- balanced_folds(metadata, outcome, balance, outer_folds, balance_candidates)
    for (o in seq_len(outer_folds)) {
      model_index <- model_index + 1L
      message("Repeat ", r, "/", repeats, ", outer fold ", o, "/", outer_folds)
      test <- which(outer == o)
      train <- which(outer != o)
      inner <- balanced_folds(metadata[train, , drop = FALSE], outcome, balance, inner_folds, balance_candidates)
      all_folds[[model_index]] <- data.frame(repetition = r, outer_fold = o,
        sample_id = metadata$sample_id,
        role = ifelse(outer == o, "test", "train"), inner_fold = NA_integer_)
      all_folds[[model_index]]$inner_fold[train] <- inner
      sse <- matrix(0, nrow = length(alphas), ncol = nlambda)
      inner_rankings <- vector("list", inner_folds)
      for (i in seq_len(inner_folds)) {
        tr <- train[inner != i]; va <- train[inner == i]
        prep <- prepare_training(beta, tr, max_missing, beta_epsilon)
        ranking <- select_markers(prep, metadata, tr, outcome, covariates, top_n)
        inner_rankings[[i]] <- ranking
        xtr <- t(prep$beta[ranking$CpG, , drop = FALSE])
        xva <- prediction_matrix(beta, va, ranking$CpG, prep$medians)
        for (a in seq_along(alphas)) {
          path <- fit_path(xtr, y[tr], alphas[a], fractions)
          pred <- predict(path$fit, newx = xva, s = path$lambda)
          sse[a, ] <- sse[a, ] + colSums((pred - y[va])^2)
        }
      }
      mse <- sse / length(train)
      # Ties: prefer stronger regularization fraction, then larger alpha.
      grid <- expand.grid(alpha_index = seq_along(alphas), fraction_index = seq_along(fractions))
      grid$alpha <- alphas[grid$alpha_index]
      grid$lambda_fraction <- fractions[grid$fraction_index]
      grid$inner_MSE <- as.vector(mse)
      best <- order(grid$inner_MSE, -grid$lambda_fraction, -grid$alpha)[1L]
      a <- grid$alpha_index[best]; j <- grid$fraction_index[best]
      prep <- prepare_training(beta, train, max_missing, beta_epsilon)
      ranking <- select_markers(prep, metadata, train, outcome, covariates, top_n)
      xtr <- t(prep$beta[ranking$CpG, , drop = FALSE])
      xte <- prediction_matrix(beta, test, ranking$CpG, prep$medians)
      path <- fit_path(xtr, y[train], alphas[a], fractions)
      prediction <- as.numeric(predict(path$fit, newx = xte, s = path$lambda[j]))
      coefficients <- as.matrix(stats::coef(path$fit, s = path$lambda[j]))
      selected <- rownames(coefficients) != "(Intercept)" & abs(coefficients[, 1L]) > coefficient_tolerance
      all_selections[[model_index]] <- data.frame(repetition = rep(r, sum(selected)),
        outer_fold = rep(o, sum(selected)), CpG = rownames(coefficients)[selected],
        coefficient = coefficients[selected, 1L], row.names = NULL)
      all_predictions[[model_index]] <- data.frame(repetition = r, outer_fold = o,
        sample_id = metadata$sample_id[test], observed = y[test], predicted = prediction)
      summaries[[model_index]] <- data.frame(repetition = r, outer_fold = o,
        alpha = alphas[a], lambda_fraction = fractions[j], lambda_min = path$lambda[j],
        inner_MSE = grid$inner_MSE[best], n_train = length(train), n_test = length(test),
        n_screened = nrow(ranking), n_nonzero = sum(selected),
        fraction_at_grid_boundary = j %in% c(1L, nlambda))
      model_dir <- file.path(output_dir, sprintf("repeat_%02d_outer_%02d", r, o))
      dir.create(model_dir)
      utils::write.csv(grid, file.path(model_dir, "inner_tuning.csv"), row.names = FALSE)
      utils::write.csv(ranking, file.path(model_dir, "outer_training_top_markers.csv"), row.names = FALSE)
      for (i in seq_len(inner_folds)) utils::write.csv(inner_rankings[[i]],
        file.path(model_dir, sprintf("inner_%02d_training_top_markers.csv", i)), row.names = FALSE)
      saveRDS(list(model = path$fit, lambda_min = path$lambda[j], alpha = alphas[a],
                   markers = ranking$CpG, imputation_medians = prep$medians[ranking$CpG]),
              file.path(model_dir, "outer_model.rds"))
    }
  }
  oof <- do.call(rbind, all_predictions)
  if (nrow(oof) != n * repeats || any(table(oof$sample_id) != repeats) ||
      anyDuplicated(oof[, c("repetition", "sample_id")]) || any(!is.finite(oof$predicted)))
    stop("OOF coverage verification failed.")
  averaged <- stats::aggregate(predicted ~ sample_id, oof, mean)
  averaged <- averaged[match(metadata$sample_id, averaged$sample_id), , drop = FALSE]
  averaged$observed <- y
  metrics <- regression_metrics(y, averaged$predicted)
  by_repeat <- do.call(rbind, lapply(seq_len(repeats), function(r) {
    d <- oof[oof$repetition == r, ]
    cbind(repetition = r, regression_metrics(d$observed, d$predicted))
  }))
  selected <- do.call(rbind, all_selections)
  counts <- table(factor(selected$CpG, levels = rownames(beta)))
  frequency <- data.frame(CpG = names(counts), selected_models = as.integer(counts),
                          total_models = model_index, frequency = as.integer(counts) / model_index)
  frequency <- frequency[order(-frequency$frequency, frequency$CpG), ]
  outputs <- list(oof_predictions = oof, oof_averaged = averaged, metrics = metrics,
    metrics_by_repeat = by_repeat, fold_assignments = do.call(rbind, all_folds),
    model_summary = do.call(rbind, summaries), selected_coefficients = selected,
    marker_frequency = frequency, exploratory_markers_ge_50pct = frequency[frequency$frequency >= .5, ])
  for (nm in names(outputs)) utils::write.csv(outputs[[nm]],
    file.path(output_dir, paste0(nm, ".csv")), row.names = FALSE)
  message("Completed. OOF results: ", normalizePath(output_dir))
  invisible(outputs)
}
