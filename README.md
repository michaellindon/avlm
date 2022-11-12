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
    Min      1Q  Median      3Q     Max 
-37.278  -6.875   0.112   6.725  38.622 

Coefficients:
            Estimate Std. Error t value Seq. p-value  2.5%  97.5% 
(Intercept)   0.9612    0.14349   6.699    1.331e-08 0.4148  1.508
X1            1.9552    0.09919  19.712    1.683e-81 1.5680  2.342
X2            3.1247    0.10135  30.831   6.710e-196 2.7296  3.520
X3            3.8641    0.09965  38.777   6.854e-303 3.4753  4.253
trt           2.2120    0.20088  11.012    3.417e-25 1.4650  2.959

Residual standard error: 10.04 on 9995 degrees of freedom
Multiple R-squared:  0.229,	Adjusted R-squared:  0.2287 
F-statistic: 742.4 on 4 and 9995 DF,  Seq. p-value: < 2.2e-16
```
