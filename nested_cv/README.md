# Nested cross-validation from normalized methylation beta values

Generic R code for continuous-outcome prediction using training-fold-specific
limma screening and elastic net regression. This release contains **no real
participant data and no study-specific configuration**. It is a reusable method
implementation, not a reproduction package for a particular paper or figure.

## What the code does

Defaults: top **800** CpGs, **10-fold outer CV**, **5-fold inner CV**, one outer
CV repetition, and minimum-inner-MSE (`lambda.min` rule) tuning. The number of
repetitions and all major analysis settings can be changed through arguments.

For every outer split:

1. Hold out the outer test set completely from fitting and tuning.
2. Within each inner split, use only the inner training set to filter missing or
   constant probes, estimate imputation medians, convert beta values to M-values,
   and run limma: `M ~ outcome + covariates`.
3. Select the top 800 eligible CpGs by the **unadjusted, two-sided moderated
   t-test P-value for the outcome**. P-value ties are ordered by CpG ID. This is
   a ranking rule, not a significance threshold.
4. Fit Gaussian `glmnet` models on the selected **beta values**, using training-only
   standardization. Evaluate the alpha/regularization grid on inner validation
   samples. Covariates adjust limma selection; they are **not predictors in glmnet**.
5. Choose the alpha and relative lambda with the lowest pooled inner-validation
   mean squared error, without using outer test outcomes.
6. Repeat preprocessing and marker selection on the entire outer training set,
   refit its selected model, and predict the held-out outer test samples once.

### Meaning of `lambda.min` here

A single call to `cv.glmnet` after selecting CpGs on the full outer training set
would not repeat feature selection inside each inner fold. This implementation
therefore uses a **custom inner CV loop**. It implements the minimum-error rule;
it does not obtain `cv.glmnet(...)$lambda.min` from a fixed feature matrix.

Because feature sets differ across folds, paths are aligned by
`lambda_fraction = lambda / lambda_max`. Each training fold estimates its own
`lambda_max` using `glmnet`. Defaults are 100 log-spaced fractions from 1 to 0.001
and alpha values 0.05, 0.1, 0.25, 0.5, 0.75, 1. The inner score pools squared
errors over all inner validation observations. Exact score ties prefer a larger
fraction and then larger alpha. The selected fraction is multiplied by the
outer-training `lambda_max` to give the saved numeric `lambda_min`.

This distinction must be described accurately if reporting the method in a
publication. The code does not use the one-standard-error rule. A selected grid
boundary is flagged in `model_summary.csv`; outer performance must not be used
to choose a preferred grid after inspecting results.

## Requirements and installation

Use R with the CRAN package `glmnet` and the Bioconductor package `limma`:

```r
install.packages("glmnet")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("limma", ask = FALSE)
```

The code is single-process and portable across Windows, macOS and Linux. The
large CpG matrix and its training subsets must fit in RAM. Package versions and
R/platform details are automatically recorded in every run's `sessionInfo.txt`.
See `VALIDATION.md` for the tested environment and synthetic checks.

## Inputs

### 1. `beta.rds`

An R numeric matrix saved with `saveRDS(beta, "beta.rds")`:

- Rows: unique CpG identifiers in `rownames(beta)`.
- Columns: unique sample identifiers in `colnames(beta)`.
- Values: normalized beta values between 0 and 1, including endpoints; `NA` is
  permitted. Infinite values are rejected.
- Include **all eligible QC-passed CpGs before outcome-based screening**. Do not
  supply a full-cohort EWAS-selected set of 800 probes.

### 2. `metadata.csv`

One row per sample, with a unique `sample_id`, a finite numeric continuous
`outcome`, and the requested covariates. The generic command-line wrapper uses
`sex` (categorical) and `age` (numeric). Labels of sex categories are unrestricted.
Input sample sets must match exactly; metadata is reordered to match beta columns.
Missing outcomes or covariates are rejected rather than silently dropping samples.

Cell proportions, if used, must be additional numeric columns. Include them in
`covariates`. If the proportions sum to one, **omit one cell type** because the
limma model includes an intercept. The code rejects rank-deficient training
designs, including absent factor levels, constant covariates and compositional
collinearity; it does not silently remove requested adjustment terms.

This code starts **after normalization**. It does not perform IDAT reading,
Noob/BMIQ normalization, probe annotation/QC, batch correction, or cell-type
deconvolution. Cell proportions may be estimated upstream with a documented
reference-based method such as FlowSorted.Blood.EPIC. These upstream procedures
and versions remain the data provider's responsibility. Outcome-dependent or
cohort-fitted upstream transformations may themselves require training-fold
handling; this code cannot undo leakage already present in its inputs.

## Run

From the folder containing these scripts:

```sh
Rscript run_analysis.R beta.rds metadata.csv results_new
```

The output directory must not already exist, preventing accidental overwriting.
For custom columns and cell-adjusted limma, call the function directly:

```r
source("nested_cv.R")
beta <- readRDS("beta.rds")
meta <- read.csv("metadata.csv", check.names = FALSE)

result <- run_nested_cv(
  beta, meta,
  outcome = "outcome",
  covariates = c("sex", "age", "cell_A", "cell_B"),
  categorical = "sex",
  balance = c("sex", "age"),
  top_n = 800,
  outer_folds = 10,
  inner_folds = 5,
  repeats = 1,
  seed = 1,
  output_dir = "results_new"
)
```

`cell_A` and `cell_B` are generic placeholders, not required biological cell
types. Replace them with the appropriate independent cell columns in your input.
For no covariates, pass `covariates = character()`, `categorical = character()`
and `balance = character()`. A different categorical covariate must be explicitly
listed in `categorical`; numerical codes are otherwise treated as continuous.

### Fold allocation and missing values

Each outer and inner allocation considers 100 randomly generated, nearly
equal-sized fold assignments. It chooses the assignment with the smallest sum
of squared fold means of standardized outcome and balance covariates (dummy
indicators for categorical covariates). Thus it approximately balances outcome,
age and sex by default; it is not exact matching or guaranteed stratification.
No model-performance criterion is used in fold construction. Inner allocation
uses only outer training samples. Set `balance_candidates` to change this number.

Within each training fold, probes with more than `max_missing = 0.05` missing
values are excluded, remaining missing beta values are replaced by training
medians, and probes with training SD <= 1e-12 are excluded. Those same medians
are applied to validation/test samples. For limma only, beta values are clipped
to `[1e-6, 1 - 1e-6]` before `log2(beta / (1 - beta))`. The elastic net uses the
unclipped, imputed beta values. Fewer than `top_n` eligible probes causes an error.

## Outputs

| File | Meaning |
|---|---|
| `oof_predictions.csv` | One held-out prediction per sample per repetition |
| `oof_averaged.csv` | Predictions averaged across repetitions for each sample |
| `metrics.csv` | Pearson R, RMSE, MAE, predictive R-squared and CCC from averaged OOF predictions |
| `metrics_by_repeat.csv` | The same metrics separately for each repetition |
| `fold_assignments.csv` | Outer training/test membership and inner-fold membership |
| `model_summary.csv` | Selected alpha, relative and absolute lambda, inner MSE, model sizes |
| `selected_coefficients.csv` | Nonzero CpG coefficients in each outer model |
| `marker_frequency.csv` | Selection counts and proportions for all input CpGs |
| `exploratory_markers_ge_50pct.csv` | CpGs with nonzero coefficients in at least half of outer models |
| `run_settings.R` | Settings for that user's run, automatically generated |
| `sessionInfo.txt` | R and package versions |
| `repeat_XX_outer_XX/` | Inner tuning scores, each training set's top-CpG ranking, fitted outer model |

Predictive R-squared is `1 - sum((observed-predicted)^2) /
sum((observed-mean(observed))^2)` and can be negative; it is not squared Pearson R.
CCC uses empirical moments with denominator n for variance/covariance. Pearson R
is undefined (NA) for constant predictions. There is no AUC calculation because
this generic workflow has a continuous outcome and no universal group cutoffs.

A CpG is counted as selected if its absolute coefficient exceeds 1e-10 in an
outer model. The frequency denominator is `outer_folds * repeats`, not the
number of inner models. The >=50% list is **exploratory** and may contain any
number of CpGs, including zero. It is not a separately validated fixed panel.
Repeated CV is optional; repeated predictions are averaged per person before
the primary metric calculation, never treated as independent participants.

The program checks that every sample receives exactly one OOF prediction per
repetition. Run outputs contain input identifiers, outcomes and trained models;
they are generated locally and are **not uploaded automatically**. The released
source archive does not include any run outputs or study settings.

## Synthetic self-test (no dataset download)

```sh
Rscript tests/test_synthetic.R
```

The test generates entirely artificial data in memory, executes the default
800-CpG / outer-10 / inner-5 structure and full tuning grid, checks complete OOF
coverage and minimum-error selection, and checks training-only preprocessing and
marker selection. Smaller additional runs test repeated-CV averaging and seed
reproducibility. Test outputs are temporary. Synthetic results are not clinical evidence.

## GitHub and Zenodo

These files can be placed in a GitHub repository. Choose a software license and
set the real authors, title, description and release version before publishing;
this template does not assign authorship or a license on the owner's behalf.
For a DOI, archive the release in Zenodo (manually upload the source ZIP, or
connect a GitHub repository and enable release archiving). Cite the DOI of the
specific archived version used. A GitHub URL alone is not a DOI.

Describe this deposit as a generic implementation. Do not claim that it alone
reproduces a publication's numerical results without the corresponding data and
settings. Data availability and code availability are separate statements.

## Method references

- limma: https://bioconductor.org/packages/limma/
- glmnet: https://glmnet.stanford.edu/
- glmnet reference: https://stat.ethz.ch/CRAN/web/packages/glmnet/refman/glmnet.html
- Zenodo deposit instructions: https://help.zenodo.org/docs/deposit/create-new-upload/
