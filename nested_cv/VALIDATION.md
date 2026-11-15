# Validation

Validated on 2026-11-14 with entirely synthetic data; no participant data was
read or bundled. Test script: `tests/test_synthetic.R`.

Environment: Windows 11, R 4.5.3, limma 3.66.0, glmnet 5.2, Matrix 1.7-5.

Checks:

- Full outer 10-fold / inner 5-fold workflow, top 800 out of 900 artificial probes,
  six alpha values and 100 relative lambda values.
- Every artificial sample receives exactly one outer held-out prediction.
- Training-only screening and imputation are unchanged when held-out beta values
  and outcomes are perturbed.
- Constant training probes are excluded from screening.
- Selected tuning parameters attain the minimum saved inner MSE.
- Reported OOF RMSE agrees with a direct independent calculation.
- Repeated CV averages predictions per person correctly.
- Repeating a run with the same seed yields identical returned results.
- No-covariate mode works; mismatched sample IDs and invalid beta ranges fail.

These are software checks, not evidence of predictive validity in any cohort.
No synthetic performance numbers are presented as research findings. New data
may expose design problems such as insufficient samples, rare factor levels,
rank-deficient covariates or fewer than 800 eligible probes; the code raises
errors rather than silently changing the analysis in these cases.
