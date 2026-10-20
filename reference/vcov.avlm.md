# Covariance Matrix for Anytime-Valid lm (avlm) Objects

Returns the covariance matrix used for anytime-valid inference: the
classical one, or the heteroskedasticity-robust or cluster-robust one
requested in
[`av.lm`](https://michaellindon.github.io/avlm/reference/av.lm.md).

## Usage

``` r
# S3 method for class 'avlm'
vcov(object, complete = TRUE, ...)
```

## Arguments

- object:

  An `avlm` object.

- complete:

  Logical; if `TRUE` (the default), aliased coefficients are included
  with `NA` rows and columns, as in
  [`vcov`](https://rdrr.io/r/stats/vcov.html).

- ...:

  Unused.

## Value

A covariance matrix.

## Examples

``` r
fit <- av(lm(mpg ~ wt + hp, data = mtcars), cluster = ~ carb)
vcov(fit)
#>              (Intercept)           wt            hp
#> (Intercept)  5.450355831 -1.782104584  4.628435e-03
#> wt          -1.782104584  0.647418361 -3.024426e-03
#> hp           0.004628435 -0.003024426  3.960511e-05
```
