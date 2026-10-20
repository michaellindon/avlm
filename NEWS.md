# avlm (development version)

## New Features

- `av.lm()` gains a `cluster` argument for CR1 cluster-robust standard errors (the
  sandwich scaled by `G / (G - 1)`, without Stata/fixest's extra `(n - 1) / (n - k)`). When `cluster` is supplied, the
  anytime-valid p-values, F-test, and confidence intervals use `G - 1` (clusters minus
  one) as the denominator degrees of freedom, rather than `n - p` (observations minus
  coefficients). `cluster` cannot be combined with `vcov_estimator`. The method is described
  in Lindon (2026), "Anytime-Valid Inference in Linear Models with Clustered Outcomes: Exact
  Multivariate Tests and Cluster-Robust Confidence Sequences", <https://arxiv.org/abs/2610.25613>.
- `cluster` accepts a vector of cluster identifiers, a one-sided formula (`cluster = ~ cyl`),
  or a column name (`cluster = "cyl"`). Formulas and column names are looked up in the
  model's `data` and aligned with the observations used in the fit.
- With `cluster`, the anytime-valid sample size is the number of clusters `G`, so `g` is on
  the cluster scale, and the joint F-test uses the cluster-robust e-value with exponent `G/2`.
  The F-test p-value is 1 when the cluster-robust covariance of the tested coefficients is
  singular (e.g. fewer clusters than tested coefficients).

## Bug Fixes

- The heteroskedasticity-robust (`vcov_estimator`) F-test now uses the Wald statistic of the
  tested coefficients, inverting their covariance rather than taking a block of the full
  precision matrix, and stores it as `Q / d`. Previously the evidence was overstated when
  testing more than one coefficient.

# avlm 0.1.0

## New Features

- First release of the **avlm** package, which implements anytime-valid inference for linear models and ANOVA.
- Provides a drop-in replacement for `summary.lm()` and `summary.aov()` with sequential p-values and confidence intervals.
- Includes support for:
  - `av()`: applies the anytime-valid transformation to `lm` and `aov` objects.
  - `confint()`: computes anytime-valid confidence intervals for regression coefficients.

## Infrastructure

- Full documentation with `roxygen2`.
- Unit tests using `testthat`.

## Notes

- Based on the methodology described in  
  *"Anytime-Valid Linear Models and Regression Adjusted Causal Inference in Randomized Experiments"*  
  available at: <https://arxiv.org/abs/2210.08589>


