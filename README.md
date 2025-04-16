# avlm
Anytime-Valid replacements for base R's lm and aov functions (with additional robust standard error functionality).


[https://arxiv.org/abs/2210.08589](https://arxiv.org/abs/2210.08589).

Contributors:
- Michael Lindon (michael.s.lindon@gmail.com)

# Anytime-Valid Linear Regression
## Code
```R
devtools::install_github("michaellindon/avlm")

# Generate Data
set.seed(1)
n = 100
X = scale(replicate(3, rnorm(n)), center = TRUE, scale = FALSE)
trt = sample(c(0, 1), size = n, replace = TRUE)
y = 1 + 0*X[,1] + 0*X[,2] + 1.4*X[,3] + 2.3*trt + 2*X[,1]*trt + 3*X[,2]*trt + 0*X[,3]*trt + 1*rnorm(n)
df = data.frame(X)
df$trt = trt
df$y = y

# Tune g to minimize the width of the confidence interval at n
alpha = 0.05
number_of_coefficients = ncol(model.matrix(y ~ . + trt * ., data = df))
g_star = optimal_g(n, number_of_coefficients, alpha)

# Summary
classic_fit = lm(y ~ . + trt * ., data = df)
av_fit = av(classic_fit, g = g_star)
summary(av_fit)
confint(av_fit)
```

## Output
```
Call:
lm(formula = y ~ . + trt * ., data = df)

Residuals:
     Min       1Q   Median       3Q      Max 
-3.00823 -0.56961  0.07473  0.78458  1.92584 

Coefficients:
            Estimate Std. Error t value p value    
(Intercept)  0.91630    0.15864   5.776 9.77e-06 ***
X1          -0.12768    0.18087  -0.706        1    
X2          -0.14966    0.18283  -0.819        1    
X3           1.30556    0.15941   8.190 8.07e-10 ***
trt          2.32245    0.21566  10.769 3.19e-14 ***
X1:trt       2.29985    0.24141   9.527 3.98e-12 ***
X2:trt       3.31166    0.23251  14.243  < 2e-16 ***
X3:trt       0.08015    0.21122   0.379        1    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
Residual standard error: 1.065 on 92 degrees of freedom
Multiple R-squared:  0.9164,	Adjusted R-squared:   0.91
F-statistic:   144 on 7 and 92 DF,  p-value: < 2.2e-16 
Anytime-Valid: TRUE


                 2.5 %    97.5 %
(Intercept)  0.4238730 1.4087195
X1          -0.6890957 0.4337358
X2          -0.7171420 0.4178268
X3           0.8107460 1.8003705
trt          1.6530529 2.9918408
X1:trt       1.5505190 3.0491788
X2:trt       2.5899720 4.0333555
X3:trt      -0.5754734 0.7357662
```

# Robust Anytime-Valid Linear Regression
## Code
```
av_fit = av(classic_fit, g=g_star, vcov_estimator = "HC0")
summary(av_fit)
confint(av_fit)
```

## Output
```
Call:
lm(formula = y ~ . + trt * ., data = df)

Residuals:
     Min       1Q   Median       3Q      Max 
-3.00823 -0.56961  0.07473  0.78458  1.92584 

Coefficients:
            Estimate Std. Error t value p value    
(Intercept)  0.91630    0.15234   6.015 4.00e-06 ***
X1          -0.12768    0.15327  -0.833        1    
X2          -0.14966    0.17150  -0.873        1    
X3           1.30556    0.13238   9.862 1.07e-12 ***
trt          2.32245    0.20302  11.439 2.54e-15 ***
X1:trt       2.29985    0.21114  10.893 1.99e-14 ***
X2:trt       3.31166    0.20403  16.232  < 2e-16 ***
X3:trt       0.08015    0.18154   0.441        1    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
Multiple R-squared:  0.9164,	Adjusted R-squared:   0.91
F-statistic:  1496 on 7 and 92 DF,  p-value: < 2.2e-16 
Anytime-Valid: TRUE
Robust Standard Error Type: HC0 


                 2.5 %    97.5 %
(Intercept)  0.4434498 1.3891427
X1          -0.6034239 0.3480640
X2          -0.6819698 0.3826546
X3           0.8946528 1.7164637
trt          1.6922727 2.9526209
X1:trt       1.6444954 2.9552024
X2:trt       2.6783800 3.9449474
X3:trt      -0.4833522 0.6436451
```

# Anytime-Valid Anova
## Code
```
classic_aov <- aov(Sepal.Length ~ Species, data = iris)
av_aov = av(classic_aov)
summary(av_aov)
```
## Output
```
             Df Sum Sq Mean Sq F value p value    
Species       2  63.21  31.606   119.3 <2e-16 ***
Residuals   147  38.96   0.265                   
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
Anytime-Valid: TRUE
```

