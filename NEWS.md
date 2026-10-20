# avlm 0.2.0

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
- `vcov()` on an `avlm` object returns the covariance used for inference (classical, HC or
  CR1).
- `optimal_g()` gains a `nu` argument, so it can tune `g` for cluster-robust inference
  (`optimal_g(G, nu = G - 1, alpha = 0.05)`).
- The summary prints the value of `g`.

## Bug Fixes

- The heteroskedasticity-robust (`vcov_estimator`) F-test now uses the Wald statistic of the
  tested coefficients, inverting their covariance rather than taking a block of the full
  precision matrix, and stores it as `Q / d`. Previously the evidence was overstated when
  testing more than one coefficient.
- `confint()` no longer misplaces intervals when the fit has aliased (`NA`) coefficients;
  previously every coefficient after an aliased one got `NA` or the wrong standard error.
  Heteroskedasticity- and cluster-robust summaries of such fits no longer fail with a
  singular-matrix error.
- Heteroskedasticity- and cluster-robust standard errors now account for `weights`
  (rows and residuals are scaled by `sqrt(w)`; zero weights are dropped). Previously the
  weights were ignored, giving wrong standard errors for weighted fits.
- Intercept-only models now work with `vcov_estimator` and `cluster`.
- `cluster` now accepts one observation per cluster, which gives the HC0 sandwich scaled by
  `n / (n - 1)`.
- `optimal_g()` now searches below `g = 1`, where the optimum lies for small samples.
- Anytime-valid ANOVA p-values use the number of observations as the sample size; models
  without an intercept were off by one.
- `av()` checks that `g` is a single positive number and gives informative errors for
  `glm`, multivariate (`mlm`) and multi-stratum (`Error()`) fits instead of failing later or
  returning meaningless output.

## Internal

- The classical, heteroskedasticity-robust and cluster-robust code paths share one
  covariance routine, and unused helpers were removed.

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


