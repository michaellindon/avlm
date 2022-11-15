# Roshi
An R package for doing fun sequntial/anytime-valid inference in linear models.

[https://arxiv.org/abs/2210.08589](https://arxiv.org/abs/2210.08589).

Contributors:
- Michael Lindon (michael.s.lindon@gmail.com)

## Example Useage
```R
> devtools::install_github("michaellindon/roshi")
> library(roshi)

> set.seed(1)
> n=100
> X = replicate(3, rnorm(n))
> X = scale(X, center=TRUE, scale=FALSE)
> trt = sample(c(0,1), size = n, replace = TRUE)
> y = 1 + 2*X[,1] + 3*X[, 2]+ 4*X[, 3]+ 2.3*trt + 2*X[,1]*trt + 3*X[, 2]*trt+ 0.4*X[, 3]*trt +  rnorm(n)
> df = data.frame(X)
> df$trt = trt
> df$y = y

> lmfit = lm(y~. + trt*.,data=df)
> avsummary(lmfit)


Call:
lm(formula = y ~ . + trt * ., data = df)

Residuals:
     Min       1Q   Median       3Q      Max 
-3.00823 -0.56961  0.07473  0.78458  1.92584 

Coefficients:
            Estimate Std. Error t value Seq. p-value   2.5%  97.5% 
(Intercept)   0.9163     0.1586   5.776    5.632e-06  0.4021  1.430
X1            1.8723     0.1809  10.352    7.565e-15  1.2921  2.453
X2            2.8503     0.1828  15.590    1.838e-24  2.2643  3.436
X3            3.9056     0.1594  24.499    1.051e-37  3.3891  4.422
trt           2.3224     0.2157  10.769    1.842e-15  1.6393  3.006
X1:trt        2.2998     0.2414   9.527    6.110e-13  1.5404  3.059
X2:trt        3.3117     0.2325  14.243    1.523e-21  2.5786  4.045
X3:trt        0.4801     0.2112   2.273    4.471e-01 -0.1899  1.150

Residual standard error: 1.065 on 92 degrees of freedom
Multiple R-squared:  0.9808,	Adjusted R-squared:  0.9793 
F-statistic: 671.1 on 7 and 92 DF,  Seq. p-value: < 2.2e-16

```
