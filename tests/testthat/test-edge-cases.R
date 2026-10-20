aliased_data <- function() {
  d <- mtcars
  d$wt2 <- 2 * d$wt
  d
}

test_that("aliased coefficients do not shift confidence intervals onto the wrong rows", {
  d <- aliased_data()
  fit <- lm(mpg ~ wt + wt2 + hp, data = d)
  ref <- lm(mpg ~ wt + hp, data = d)
  for (args in list(list(), list(vcov_estimator = "HC1"), list(cluster = ~ carb))) {
    ci <- confint(do.call(av, c(list(fit), args)))
    ci_ref <- confint(do.call(av, c(list(ref), args)))
    expect_true(all(is.na(ci["wt2", ])))
    expect_equal(ci[c("(Intercept)", "wt", "hp"), ], ci_ref, tolerance = 1e-10)
  }
})

test_that("robust and clustered summaries work with aliased coefficients", {
  d <- aliased_data()
  fit <- lm(mpg ~ wt + wt2 + hp, data = d)
  ref <- lm(mpg ~ wt + hp, data = d)
  for (args in list(list(vcov_estimator = "HC3"), list(cluster = ~ carb))) {
    s <- summary(do.call(av, c(list(fit), args)))
    s_ref <- summary(do.call(av, c(list(ref), args)))
    expect_equal(s$coefficients, s_ref$coefficients, tolerance = 1e-10)
    expect_equal(attr(s, "f_pvalue"), attr(s_ref, "f_pvalue"), tolerance = 1e-10)
  }
})

test_that("robust and clustered standard errors account for weights", {
  set.seed(3)
  w <- runif(nrow(mtcars)) + 0.5
  fit <- lm(mpg ~ wt + hp, data = mtcars, weights = w)
  X <- model.matrix(fit) * sqrt(w)
  e <- residuals(fit) * sqrt(w)
  bread <- solve(crossprod(X))

  V_hc0 <- bread %*% crossprod(X * e) %*% bread
  expect_equal(unname(vcov(av(fit, vcov_estimator = "HC0"))), unname(V_hc0), tolerance = 1e-10)

  scores <- rowsum(X * e, mtcars$carb)
  G <- nrow(scores)
  V_cr1 <- G / (G - 1) * bread %*% crossprod(scores) %*% bread
  expect_equal(unname(vcov(av(fit, cluster = ~ carb))), unname(V_cr1), tolerance = 1e-10)
})

test_that("observations with zero weight are dropped from robust inference", {
  w <- rep(1, nrow(mtcars))
  w[1:4] <- 0
  fit_w <- lm(mpg ~ wt + hp, data = mtcars, weights = w)
  fit_sub <- lm(mpg ~ wt + hp, data = mtcars[-(1:4), ])
  for (args in list(list(vcov_estimator = "HC1"), list(vcov_estimator = "HC3"))) {
    expect_equal(summary(do.call(av, c(list(fit_w), args)))$coefficients,
                 summary(do.call(av, c(list(fit_sub), args)))$coefficients, tolerance = 1e-10)
  }
  expect_equal(summary(av(fit_w, cluster = mtcars$carb))$coefficients,
               summary(av(fit_sub, cluster = mtcars$carb[-(1:4)]))$coefficients, tolerance = 1e-10)
})

test_that("intercept-only models work with robust and clustered standard errors", {
  fit <- lm(mpg ~ 1, data = mtcars)
  for (args in list(list(), list(vcov_estimator = "HC0"), list(cluster = ~ carb))) {
    s <- summary(do.call(av, c(list(fit), args)))
    expect_null(s$fstatistic)
    expect_output(print(s), "Anytime-Valid: TRUE")
    expect_equal(nrow(confint(do.call(av, c(list(fit), args)))), 1)
  }
})

test_that("singleton clusters give the HC0 sandwich scaled by n / (n - 1)", {
  fit <- lm(mpg ~ wt + hp, data = mtcars)
  n <- nrow(mtcars)
  V_cl <- vcov(av(fit, cluster = seq_len(n)))
  V_hc0 <- vcov(av(fit, vcov_estimator = "HC0"))
  expect_equal(V_cl, n / (n - 1) * V_hc0, tolerance = 1e-10)
})

test_that("vcov() returns the covariance used for inference", {
  fit <- lm(mpg ~ wt + hp, data = mtcars)
  expect_equal(vcov(av(fit)), vcov(fit))
  s <- summary(av(fit, vcov_estimator = "HC3"))
  expect_equal(sqrt(diag(vcov(av(fit, vcov_estimator = "HC3")))), s$coefficients[, 2])
})

test_that("unsupported inputs give informative errors", {
  fit <- lm(mpg ~ wt, data = mtcars)
  expect_error(av(fit, g = -1), "positive")
  expect_error(av(fit, g = c(1, 2)), "positive")
  expect_error(av(fit, g = "1"), "positive")
  expect_error(av(lm(cbind(mpg, qsec) ~ wt, data = mtcars)), "multivariate")
  expect_error(av(aov(yield ~ N + Error(block), data = npk)), "single-stratum")
  expect_error(av(glm(am ~ wt, data = mtcars, family = binomial)), "glm")
})

test_that("ANOVA e-values use the number of observations as the sample size", {
  fit <- aov(Sepal.Length ~ Species - 1, data = iris)
  s <- summary(av(fit, g = 2))[[1]]
  f <- s$`F value`[1]
  expect_equal(s$`Pr(>F)`[1], p_G_f(log_G_f(f, 3, 147, 150, 2)), tolerance = 1e-12)
})

test_that("the summary reports g", {
  expect_output(print(summary(av(lm(mpg ~ wt, data = mtcars), g = 2.5))), "g = 2.5")
})
