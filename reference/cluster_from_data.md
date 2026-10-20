# Look Up a Cluster Variable in the Data a Model Was Fit With

Look Up a Cluster Variable in the Data a Model Was Fit With

## Usage

``` r
cluster_from_data(model, cluster)
```

## Arguments

- model:

  A fitted `lm` object.

- cluster:

  A one-sided formula such as `~ cyl`, or a single column name.

## Value

The cluster variable, subset to the observations used in the fit.
