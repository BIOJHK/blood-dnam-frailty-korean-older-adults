# Usage: Rscript run_analysis.R beta.rds metadata.csv new_output_directory
# Edit the generic column mappings below if your input uses different names.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: Rscript run_analysis.R beta.rds metadata.csv new_output_directory")
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- dirname(sub("^--file=", "", script_arg[1L]))
source(file.path(script_dir, "nested_cv.R"))
beta <- readRDS(args[1L])
metadata <- read.csv(args[2L], check.names = FALSE, stringsAsFactors = FALSE)
# For cell-adjusted limma selection, add your cell columns to covariates below.
# If cell proportions sum to one, omit one cell type to avoid collinearity.
run_nested_cv(beta, metadata,
  outcome = "outcome", covariates = c("sex", "age"), categorical = "sex",
  balance = c("sex", "age"), top_n = 800L, outer_folds = 10L,
  inner_folds = 5L, repeats = 1L, output_dir = args[3L])
