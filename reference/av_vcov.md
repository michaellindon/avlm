# Coefficient Covariance and Degrees of Freedom for Anytime-Valid Inference

Computes, for an `avlm` object, the coefficient estimates, their
covariance matrix and the two quantities the anytime-valid e-values
need: the denominator degrees of freedom `nu` and the sample size `n` on
which the precision `g` acts.

## Usage

``` r
av_vcov(object, summ = NULL)
```

## Arguments

- object:

  An `avlm` object.

- summ:

  Optionally, `summary.lm(object)`, to avoid computing it twice.

## Value

A list with elements `coefficients` (non-aliased estimates), `vcov`,
`nu`, `n`, `type` (`"classical"`, an `HC` type, or `"CR1"`) and `G`
(number of clusters, or `NULL`).

## Details

- Classical: \\s^2 (X'X)^{-1}\\, \\\nu = n - k\\ (Lindon et al., 2026,
  Section 4).

- Heteroskedasticity-robust (`HC0`–`HC3`): the sandwich \\(X'X)^{-1}
  X'\Omega X (X'X)^{-1}\\, \\\nu = n - k\\ (Lindon et al., 2026, Section
  5).

- Cluster-robust (`CR1`): the sandwich of the summed cluster scores
  scaled by \\G/(G-1)\\, \\\nu = G - 1\\ and sample size \\G\\ (Lindon,
  2026, eq. CR1).

Aliased (`NA`) coefficients are dropped. For weighted fits, the rows of
\\X\\ and the residuals are scaled by \\\sqrt{w}\\ and observations with
zero weight are dropped.
