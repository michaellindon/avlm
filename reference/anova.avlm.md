# Anytime-Valid ANOVA Table for an avlm Object

Computes the sequential (Type I) ANOVA table of
[`anova.lm`](https://rdrr.io/r/stats/anova.lm.html) and replaces its
p-values with anytime-valid ones, using the same e-values as
[`summary.avaov`](https://michaellindon.github.io/avlm/reference/summary.avaov.md).
Comparing several models (`anova(fit1, fit2)`) is not supported.

## Usage

``` r
# S3 method for class 'avlm'
anova(object, ...)
```

## Arguments

- object:

  An `avlm` object.

- ...:

  Must be empty.

## Value

An `anova` table with anytime-valid p-values.

## Examples

``` r
anova(av(lm(mpg ~ wt + hp, data = mtcars), g = 1))
#> Analysis of Variance Table
#> 
#> Response: mpg
#> Anytime-valid p-values (g = 1)
#> 
#>           Df Sum Sq Mean Sq F value    Pr(>F)    
#> wt         1 847.73  847.73 126.041 4.411e-10 ***
#> hp         1  83.27   83.27  12.381   0.03364 *  
#> Residuals 29 195.05    6.73                      
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
