# Summary Method for Anytime-Valid aov Objects

This method produces a summary for objects of class `avaov`. It first
calls the default `summary.aov` method and then replaces the standard
p-values with anytime-valid p-values calculated using the precision
parameter `g`.

## Usage

``` r
# S3 method for class 'avaov'
summary(object, ...)
```

## Arguments

- object:

  An object of class `avaov` created by `av.aov`.

- ...:

  Additional arguments passed to or from other methods.

## Value

A summary object of class `summary.avaov` that includes the
anytime-valid p-values.

## Examples

``` r

# Fit an ANOVA model to the iris dataset.
# This model tests whether the sepal length differs by species.
fit_aov <- aov(Sepal.Length ~ Species, data = iris)

# Convert the standard aov object to an anytime-valid aov (avaov) with precision parameter g = 1.
av_fit_aov <- av(fit_aov, g = 1)

# Print the summary of the anytime-valid ANOVA model.
# The summary replaces standard p-values with anytime-valid p-values.
summary(av_fit_aov)
#>              Df Sum Sq Mean Sq F value p value    
#> Species       2  63.21  31.606   119.3 <2e-16 ***
#> Residuals   147  38.96   0.265                   
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> Anytime-Valid: TRUE (g = 1)
```
