# Confidence Intervals for Anytime-Valid lm (avlm) Objects

Computes confidence intervals for the coefficients of an `avlm` object
by adjusting the standard errors using the precision parameter `g` and
an anytime-valid approach.

## Usage

``` r
# S3 method for class 'avlm'
confint(object, parm, level = 0.95, ...)
```

## Arguments

- object:

  An `avlm` object.

- parm:

  A specification of which parameters are to be given confidence
  intervals; can be a vector of numbers or names.

- level:

  The confidence level required for the intervals. Defaults to 0.95.

- ...:

  Additional arguments passed to or from other methods.

## Value

A matrix with the lower and upper confidence limits for the specified
parameters.

## Examples

``` r
# Fit a simple linear model using the mtcars dataset
fit <- lm(mpg ~ wt + hp, data = mtcars)

# Convert the standard lm object to an anytime-valid avlm object with precision parameter g = 1
av_fit <- av(fit, g = 1)

# Calculate and print confidence intervals for the coefficients
conf_intervals <- confint(av_fit)
print(conf_intervals)
#>                   2.5 %       97.5 %
#> (Intercept) 31.86338762 42.591152611
#> wt          -6.00063195 -1.755029536
#> hp          -0.06206734 -0.001478552
```
