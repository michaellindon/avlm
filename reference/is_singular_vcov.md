# Check Whether a Covariance Matrix Is Numerically Singular

Scale-invariant: checks the smallest eigenvalue of the implied
correlation matrix.

## Usage

``` r
is_singular_vcov(V, tol = 1e-10)
```

## Arguments

- V:

  A symmetric covariance matrix.

- tol:

  Tolerance for the smallest eigenvalue of the correlation matrix.

## Value

`TRUE` if `V` is numerically singular.
