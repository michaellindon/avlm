# Convert a Linear Model (lm) Object to Anytime-Valid lm (avlm)

Converts an `lm` object into an anytime-valid version by storing a
precision parameter `g` as an attribute and updating the object's class.
The resulting object supports anytime-valid inference.

## Usage

``` r
# S3 method for class 'lm'
av(model, g = 1, vcov_estimator = NULL, cluster = NULL, ...)
```

## Arguments

- model:

  An `lm` object from a linear model fit.

- g:

  A positive precision parameter for the anytime-valid mixture (see
  Lindon et al., 2026, Section 4.3.1). It acts like a prior sample size:
  larger values favour smaller effects and give narrower intervals at
  large \\n\\, at the cost of wider intervals at small \\n\\. Use
  [`optimal_g`](https://michaellindon.github.io/avlm/reference/optimal_g.md)
  to choose it for a target sample size. Default is 1.

- vcov_estimator:

  Optional character string specifying the type of
  heteroskedasticity-robust standard errors to use. Must be one of
  `"HC0"`, `"HC1"`, `"HC2"`, or `"HC3"`. If `NULL` (default), no robust
  variance estimation is applied. Cannot be combined with `cluster`.

- cluster:

  Optional cluster specification for cluster-robust standard errors.
  Either a vector of cluster identifiers, one per observation used in
  fitting `model` (e.g. a respondent or unit id); a one-sided formula
  such as `~ cyl`; or a single column name such as `"cyl"`. Formulas and
  column names are evaluated in the `data` the model was fit with, and
  rows dropped from the fit (by `subset` or missing values) are dropped
  from the clusters too. Only one-way clustering is supported. When
  supplied, standard errors are computed using the CR1 cluster-robust
  sandwich estimator (scaled by \\G/(G-1)\\), and the anytime-valid
  p-values and confidence intervals use \\G - 1\\ (clusters minus one)
  as the denominator degrees of freedom, rather than \\n - p\\
  (observations minus coefficients), and the number of clusters \\G\\ as
  the sample size, so `g` is on the cluster scale. The cluster-robust
  e-values and confidence sequences are those of Lindon (2026). Cannot
  be combined with `vcov_estimator`.

- ...:

  Additional arguments passed to or from other methods.

## Value

An enhanced `lm` object of class `avlm` with anytime-valid inference
capabilities.

## References

Lindon, M., Ham, D. W., Tingley, M., and Bojinov, I. (2026).
Anytime-Valid Inference in Linear Models with Applications to
Regression-Adjusted Causal Inference. *Journal of the American
Statistical Association*.
[doi:10.1080/01621459.2026.2692052](https://doi.org/10.1080/01621459.2026.2692052)

Lindon, M. (2026). Anytime-Valid Inference in Linear Models with
Clustered Outcomes: Exact Multivariate Tests and Cluster-Robust
Confidence Sequences. arXiv:2610.25613.
<https://arxiv.org/abs/2610.25613>

## Examples

``` r
# Fit a linear model on the built-in mtcars dataset
fit <- lm(mpg ~ wt + hp, data = mtcars)

# Convert the standard lm object to an anytime-valid avlm object
av_fit <- av(fit, g = 1)

# Print the summary of the anytime-valid model
summary(av_fit)
#> 
#> Call:
#> lm(formula = mpg ~ wt + hp, data = mtcars)
#> 
#> Residuals:
#>    Min     1Q Median     3Q    Max 
#> -3.941 -1.600 -0.182  1.050  5.854 
#> 
#> Coefficients:
#>             Estimate Std. Error t value p value    
#> (Intercept) 37.22727    1.59879  23.285  < 2e-16 ***
#> wt          -3.87783    0.63273  -6.129 3.96e-05 ***
#> hp          -0.03177    0.00903  -3.519   0.0336 *  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> Residual standard error: 2.593 on 29 degrees of freedom
#> Multiple R-squared:  0.8268, Adjusted R-squared:  0.8148
#> F-statistic: 69.21 on 2 and 29 DF,  p-value: 4.226e-10 
#> Anytime-Valid: TRUE (g = 1)

# With cluster-robust standard errors, clustering on number of carburetors
av_fit_cl <- av(fit, g = 1, cluster = ~ carb)
summary(av_fit_cl)
#> 
#> Call:
#> lm(formula = mpg ~ wt + hp, data = mtcars)
#> 
#> Residuals:
#>    Min     1Q Median     3Q    Max 
#> -3.941 -1.600 -0.182  1.050  5.854 
#> 
#> Coefficients:
#>              Estimate Std. Error t value p value  
#> (Intercept) 37.227270   2.334600  15.946   0.0107 *
#> wt          -3.877831   0.804623  -4.819   0.0677 .
#> hp          -0.031773   0.006293  -5.049   0.0602 .
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> Multiple R-squared:  0.8268, Adjusted R-squared:  0.8148
#> F-statistic: 60.46 on 2 and 5 DF,  p-value: 0.03875 
#> Anytime-Valid: TRUE (g = 1)
#> Robust Standard Error Type: CR1 
#> Clusters: 6 

# Equivalent ways of specifying the clusters
av_fit_cl <- av(fit, g = 1, cluster = "carb")
av_fit_cl <- av(fit, g = 1, cluster = mtcars$carb)
```
