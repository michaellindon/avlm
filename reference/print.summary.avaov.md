# Print Method for summary.avaov Objects

This method prints the summary of an `avaov` object. It captures the
output from the default printing method, substitutes the header
"Pr(\>F)" with "p value", and adds a note indicating that anytime-valid
inference is used.

## Usage

``` r
# S3 method for class 'summary.avaov'
print(
  x,
  digits = max(3L, getOption("digits") - 3L),
  signif.stars = getOption("show.signif.stars"),
  ...
)
```

## Arguments

- x:

  An object of class `summary.avaov`.

- digits:

  The number of significant digits to use when printing. Defaults to a
  value based on options.

- signif.stars:

  Logical indicating whether significance stars should be printed.

- ...:

  Additional arguments passed to or from other methods.

## Value

Invisibly returns the summary object.
