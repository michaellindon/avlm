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
> trt = sample(c(0,1), size = n, replace = TRUE)
> y = 1 + 2*X[,1] + 3*X[, 2]+ 4*X[, 3]+ 2.3*trt + rnorm(n)
> df = data.frame(X)
> df$trt = trt
> df$y = y

> lmfit = lm(y~., data=df)
> avsummary(lmfit)

Call:
lm(formula = y ~ ., data = df)

Residuals:
     Min       1Q   Median       3Q      Max 
-2.90294 -0.58641  0.00453  0.73891  2.02803 

Coefficients:
            Estimate Std. Error t value Seq. p-value  2.5%  97.5% 
(Intercept)    0.924     0.1593    5.80    4.780e-06 0.4082  1.440
X1             2.054     0.1196   17.18    2.747e-28 1.6577  2.451
X2             3.050     0.1125   27.12    5.768e-43 2.6748  3.425
X3             3.946     0.1043   37.82    2.089e-54 3.5953  4.296
trt            2.318     0.2161   10.73    1.526e-15 1.6340  3.002

Residual standard error: 1.067 on 95 degrees of freedom
Multiple R-squared:  0.9651,	Adjusted R-squared:  0.9636 
F-statistic: 656.5 on 4 and 95 DF,  Seq. p-value: < 2.2e-16
```
