# Print Method for summary.avlm Objects

Custom print method for `summary.avlm` objects that mimics the standard
`summary.lm` output, but replaces p-value column headers and includes
anytime-valid annotations.

## Usage

``` r
# S3 method for class 'summary.avlm'
print(
  x,
  digits = max(3L, getOption("digits") - 3L),
  signif.stars = getOption("show.signif.stars"),
  ...
)
```

## Arguments

- x:

  An object of class `summary.avlm` containing the summary of an
  anytime-valid lm object.

- digits:

  The number of significant digits to use when printing. Defaults based
  on system options.

- signif.stars:

  Logical indicating whether significance stars should be printed.

- ...:

  Additional arguments passed to or from other methods.

## Value

Invisibly returns the summary object.
