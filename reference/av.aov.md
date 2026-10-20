# Convert an aov Object to Anytime-Valid aov (avaov)

Converts an object of class `aov` to an anytime-valid version by setting
the precision parameter `g` as an attribute and updating the class.

## Usage

``` r
# S3 method for class 'aov'
av(model, g = 1, ...)
```

## Arguments

- model:

  An `aov` object resulting from an ANOVA analysis.

- g:

  A positive precision parameter for anytime-valid inference (see
  [`av.lm`](https://michaellindon.github.io/avlm/reference/av.lm.md)).
  Default is 1.

- ...:

  Additional arguments passed to or from other methods.

## Value

An object of class `avaov` (inheriting from `avlm`) with anytime-valid
p-values.

## Details

The ANOVA table uses sequential (Type I) sums of squares, as
[`summary.aov`](https://rdrr.io/r/stats/summary.aov.html) does, so in
unbalanced designs each term is tested given the terms listed before it
and the p-values depend on the order of the terms. Each row is the
anytime-valid version of that sequential F-test. The object also
inherits from `avlm`, so
[`confint`](https://rdrr.io/r/stats/confint.html) and
[`vcov`](https://rdrr.io/r/stats/vcov.html) give anytime-valid intervals
for the coefficients.
