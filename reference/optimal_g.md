# Computes the value of g such that width of the \\1-\alpha\\ confidence interval at sample size n is minimized

Minimizes the radius \\R(g, n, 1, \alpha)\\ of the asymptotic t
confidence sequence (Lindon et al., 2026, Section 5.1) over \\g\\.

## Usage

``` r
optimal_g(n, number_of_coefficients, alpha, nu = n - number_of_coefficients)
```

## Arguments

- n:

  A positive sample size integer. For cluster-robust inference, the
  number of clusters.

- number_of_coefficients:

  A positive integer of coefficients in the full model

- alpha:

  A positive numeric scalar in (0,1) for nominal Type I error.

- nu:

  Denominator degrees of freedom. Defaults to
  `n - number_of_coefficients`; for cluster-robust inference
  (`av(..., cluster = )`) use `nu = n - 1` with `n` the number of
  clusters.

## Value

A positive numeric scalar representing the optimal \\g\\ that minimizes
the CI width.

## Examples

``` r
n <- 10000
alpha <- 0.05
g_star <- optimal_g(n, 5, alpha)
cat("The optimal g is:", g_star, "\n")
#> The optimal g is: 1217.196 

# cluster-robust inference with 200 clusters
optimal_g(200, alpha = 0.05, nu = 199)
#> [1] 23.79657
```
