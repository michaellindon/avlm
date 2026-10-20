# Resolve a Cluster Argument into a Factor Aligned with the Fitted Model

Resolve a Cluster Argument into a Factor Aligned with the Fitted Model

## Usage

``` r
resolve_cluster_factor(model, cluster)
```

## Arguments

- model:

  A fitted `lm` object.

- cluster:

  `NULL`, a vector of cluster identifiers with one entry per observation
  used in the fit, a one-sided formula such as `~ cyl`, or a single
  column name such as `"cyl"`. Formulas and column names are looked up
  in the `data` the model was fit with and subset to the observations
  used in the fit.

## Value

A factor of cluster identifiers, or `NULL` if `cluster` is `NULL`.
