# Anytime-valid Conversion Generic Function

This generic function converts a fitted model object into an
anytime-valid version. The conversion is performed by the appropriate S3
method based on the class of the input model.

## Usage

``` r
av(model, g = 1, ...)
```

## Arguments

- model:

  A fitted model object (e.g., an object of class `aov` or `lm`).

- g:

  A positive precision parameter for anytime-valid inference (see
  [`av.lm`](https://michaellindon.github.io/avlm/reference/av.lm.md)).
  Defaults to 1.

- ...:

  Additional arguments passed to the method.

## Value

An enhanced version of `model` with anytime-valid inference
capabilities.
