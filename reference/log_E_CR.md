# Log Cluster-Robust e-Value for a Joint Wald Test

The cluster-robust e-process of Lindon (2026), eq. (ECR), with \\\nu =
G - 1\\ and exponent \\G/2\\ for any number \\q\\ of tested
coefficients: \$\$E_G = \left(\frac{g}{g+G}\right)^{q/2} \left(\frac{1 +
Q/(G-1)}{1 + \frac{g}{g+G} Q/(G-1)}\right)^{G/2}.\$\$

## Usage

``` r
log_E_CR(Q, q, G, g)
```

## Arguments

- Q:

  The cluster-robust Wald statistic.

- q:

  The number of tested coefficients.

- G:

  The number of clusters.

- g:

  The precision parameter, on the cluster scale.

## Value

The log e-value.

## References

Lindon, M. (2026). Anytime-Valid Inference in Linear Models with
Clustered Outcomes: Exact Multivariate Tests and Cluster-Robust
Confidence Sequences. arXiv:2610.25613.
<https://arxiv.org/abs/2610.25613>
