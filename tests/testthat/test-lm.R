test_that("anytime-valid p-values are more conservative than standard p-values for linear models", {
  # Test with different precision parameters
  for (g_val in c(1, 2, 5)) {
    # Fit standard and anytime-valid linear models
    std_fit <- lm(mpg ~ wt + hp, data = mtcars)
    std_summary <- summary(std_fit)
    std_pvals <- std_summary$coefficients[, 4]
    
    av_fit <- av(std_fit, g = g_val)
    av_summary <- summary(av_fit)
    av_pvals <- av_summary$coefficients[, 4]
    
    # Check if all anytime-valid p-values are greater than or equal to standard p-values
    for (i in seq_along(std_pvals)) {
      expect_gte(av_pvals[i], std_pvals[i])
    }
  }
})

test_that("anytime-valid confidence intervals are supersets of standard confidence intervals for linear models", {
  # Test with different precision parameters and confidence levels
  for (g_val in c(1, 2, 5)) {
    for (level_val in c(0.90, 0.95, 0.99)) {
      # Fit standard and anytime-valid linear models
      std_fit <- lm(mpg ~ wt + hp, data = mtcars)
      std_ci <- confint(std_fit, level = level_val)
      
      av_fit <- av(std_fit, g = g_val)
      av_ci <- confint(av_fit, level = level_val)
      
      # Check if each anytime-valid CI encloses the standard CI
      for (i in 1:nrow(std_ci)) {
        coef_name <- rownames(std_ci)[i]
        expect_lte(av_ci[coef_name, 1], std_ci[coef_name, 1])
        expect_gte(av_ci[coef_name, 2], std_ci[coef_name, 2])
      }
    }
  }
})

test_that("sequential p-value computed via log_G_f is greater than the standard F-test p-value", {
  for (g_val in c(1, 2, 5)) {
    # Fit the ordinary linear model and compute its F-test p-value.
    std_fit <- lm(mpg ~ wt + hp, data = mtcars)
    std_f <- summary(std_fit)$fstatistic
    std_f_pvalue <- pf(std_f[1], std_f[2], std_f[3], lower.tail = FALSE)
    n <- length(std_fit$residuals)

    
    av_seq_pvalue <- p_G_f(log_G_f(std_f[1], std_f[2], std_f[3], n, g_val))
    
    # Assert that the standard F p-value is less than the anytime-valid (sequential) p-value.
    expect_lt(std_f_pvalue, av_seq_pvalue)
  }
})


test_that("anytime-valid methods work with increasingly complex linear models", {
  # Test using a more complex linear model with interactions
  complex_lm <- lm(mpg ~ wt * hp + factor(cyl), data = mtcars)
  av_complex_lm <- av(complex_lm)
  expect_no_error(summary(av_complex_lm))
  expect_no_error(confint(av_complex_lm))
})

test_that("heteroskedasticity-robust F-statistic is Q / d, so its e-value uses Q / nu", {
  fit <- lm(mpg ~ wt + hp + qsec, data = mtcars)
  X <- model.matrix(fit)
  e <- residuals(fit)
  n <- nrow(X); k <- ncol(X)
  bread <- solve(crossprod(X))
  V <- (n / (n - k)) * bread %*% crossprod(X * e) %*% bread
  b <- coef(fit)[-1]
  Q <- as.numeric(t(b) %*% solve(V[-1, -1]) %*% b)
  d <- length(b)

  av_summary <- summary(av(fit, g = 1, vcov_estimator = "HC1"))
  expect_equal(unname(av_summary$fstatistic[1]), Q / d, tolerance = 1e-8)

  # Theorem (heteroskedasticity-robust e-process): Q / nu inside, exponent (nu + d) / 2
  nu <- n - k
  r <- 1 / (1 + n)
  expected_log_E <- 0.5 * d * log(r) + 0.5 * (nu + d) * (log(1 + Q / nu) - log(1 + r * Q / nu))
  f <- av_summary$fstatistic
  expect_equal(unname(log_G_f(f[1], f[2], f[3], n, 1)), expected_log_E, tolerance = 1e-8)
})
