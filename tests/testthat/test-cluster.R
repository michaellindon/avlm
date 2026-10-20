manual_vcov_CR1 <- function(fit, cluster) {
  X <- model.matrix(fit)
  e <- residuals(fit)
  k <- ncol(X)
  cluster_f <- factor(cluster)
  G <- nlevels(cluster_f)
  bread <- solve(crossprod(X))
  meat <- matrix(0, k, k)
  for (lvl in levels(cluster_f)) {
    idx <- which(cluster_f == lvl)
    Xeg <- colSums((X * e)[idx, , drop = FALSE])
    meat <- meat + tcrossprod(Xeg)
  }
  (G / (G - 1)) * (bread %*% meat %*% bread)
}

test_that("the CR1 covariance matches a manual sandwich computation", {
  fit <- lm(mpg ~ wt + hp, data = mtcars)
  av_fit <- av(fit, cluster = mtcars$carb)

  expected <- manual_vcov_CR1(fit, mtcars$carb)
  expect_equal(unname(vcov(av_fit)), unname(expected), tolerance = 1e-10)
  expect_equal(av_vcov(av_fit)$G, length(unique(mtcars$carb)))
})

test_that("av.lm with cluster produces cluster-robust standard errors and dfd = G - 1", {
  fit <- lm(mpg ~ wt + hp, data = mtcars)
  G <- length(unique(mtcars$carb))

  av_fit <- av(fit, g = 1, cluster = mtcars$carb)
  av_summary <- summary(av_fit)

  expected_vcov <- manual_vcov_CR1(fit, mtcars$carb)
  expected_se <- sqrt(diag(expected_vcov))
  expect_equal(unname(av_summary$coefficients[, 2]), unname(expected_se), tolerance = 1e-10)

  # coefficient p-values use dfd = G - 1 and sample size G
  t2 <- av_summary$coefficients[, 3] ^ 2
  expected_pvals <- p_G_t(log_G_t(t2, G - 1, G, 1))
  expect_equal(unname(av_summary$coefficients[, 4]), unname(expected_pvals), tolerance = 1e-10)

  # the F-test denominator degrees of freedom is also G - 1
  expect_equal(unname(av_summary$fstatistic[3]), G - 1)
})

test_that("clustered anytime-valid confidence intervals use dfd = G - 1 and sample size G", {
  fit <- lm(mpg ~ wt + hp, data = mtcars)
  G <- length(unique(mtcars$carb))

  av_fit <- av(fit, g = 1, cluster = mtcars$carb)
  se <- summary(av_fit)$coefficients[, 2]

  ci <- confint(av_fit, level = 0.95)
  t_rad <- t_radius(1, G, 0.05, nu = G - 1)

  expected_lo <- coef(fit) - t_rad * se
  expected_hi <- coef(fit) + t_rad * se

  expect_equal(unname(ci[, 1]), unname(expected_lo), tolerance = 1e-10)
  expect_equal(unname(ci[, 2]), unname(expected_hi), tolerance = 1e-10)
})

test_that("cluster and vcov_estimator cannot be combined", {
  fit <- lm(mpg ~ wt + hp, data = mtcars)
  expect_error(
    av(fit, g = 1, cluster = mtcars$carb, vcov_estimator = "HC1"),
    "only one of"
  )
})

test_that("cluster length must match the number of fitted observations", {
  fit <- lm(mpg ~ wt + hp, data = mtcars)
  expect_error(
    av(fit, g = 1, cluster = mtcars$carb[1:10]),
    "length\\(cluster\\)"
  )
})

test_that("cluster must have at least 2 distinct levels", {
  fit <- lm(mpg ~ wt + hp, data = mtcars)
  expect_error(
    av(fit, g = 1, cluster = rep(1, nrow(mtcars))),
    "at least 2 distinct levels"
  )
})

test_that("clustered anytime-valid p-values are more conservative than clustered classical p-values", {
  fit <- lm(mpg ~ wt + hp, data = mtcars)
  for (g_val in c(1, 2, 5)) {
    av_fit <- av(fit, g = g_val, cluster = mtcars$carb)
    av_summary <- summary(av_fit)

    expected_vcov <- manual_vcov_CR1(fit, mtcars$carb)
    se <- sqrt(diag(expected_vcov))
    t <- coef(fit) / se
    G <- length(unique(mtcars$carb))
    classical_pvals <- 2 * pt(-abs(t), df = G - 1)

    for (i in seq_along(classical_pvals)) {
      expect_gte(av_summary$coefficients[i, 4], classical_pvals[i])
    }
  }
})

test_that("cluster can be given as a formula or a column name", {
  fit <- lm(mpg ~ wt + hp, data = mtcars)
  expected <- summary(av(fit, g = 1, cluster = mtcars$carb))$coefficients

  expect_equal(summary(av(fit, g = 1, cluster = ~ carb))$coefficients, expected)
  expect_equal(summary(av(fit, g = 1, cluster = "carb"))$coefficients, expected)
})

test_that("formula and column-name clusters are aligned with rows used in the fit", {
  df <- mtcars
  df$wt[c(3, 10)] <- NA
  fit <- lm(mpg ~ wt + hp, data = df, subset = gear != 5)
  used <- !is.na(df$wt) & df$gear != 5
  expected <- summary(av(fit, g = 1, cluster = df$carb[used]))$coefficients

  expect_equal(summary(av(fit, g = 1, cluster = ~ carb))$coefficients, expected)
  expect_equal(summary(av(fit, g = 1, cluster = "carb"))$coefficients, expected)
})

test_that("cluster formulas can be expressions of columns", {
  fit <- lm(mpg ~ wt + hp, data = mtcars)
  expected <- summary(av(fit, g = 1, cluster = interaction(mtcars$cyl, mtcars$am)))$coefficients
  expect_equal(summary(av(fit, g = 1, cluster = ~ interaction(cyl, am)))$coefficients, expected)
})

test_that("invalid formula and column-name clusters give informative errors", {
  fit <- lm(mpg ~ wt + hp, data = mtcars)
  expect_error(av(fit, cluster = ~ carb + am), "one-way")
  expect_error(av(fit, cluster = gear ~ carb), "one-sided")
  expect_error(av(fit, cluster = "not_a_column"), "not found")

  x <- mtcars$wt
  y <- mtcars$mpg
  fit_no_data <- lm(y ~ x)
  expect_error(av(fit_no_data, cluster = ~ carb), "`data` argument")
})

test_that("cluster must not contain missing values", {
  df <- mtcars
  df$cyl[1] <- NA
  fit <- lm(mpg ~ wt + hp, data = df)
  expect_error(av(fit, cluster = ~ cyl), "missing values")
  expect_error(av(fit, cluster = df$cyl), "missing values")
})

test_that("clustered F-test uses exponent G/2 for any number of tested coefficients", {
  fit <- lm(mpg ~ wt + hp + qsec, data = mtcars)
  G <- length(unique(mtcars$carb))
  for (g_val in c(1, 3)) {
    av_summary <- summary(av(fit, g = g_val, cluster = ~ carb))

    V <- manual_vcov_CR1(fit, mtcars$carb)[-1, -1]
    b <- coef(fit)[-1]
    Q <- as.numeric(t(b) %*% solve(V) %*% b)
    q <- length(b)
    r <- g_val / (g_val + G)
    expected_log_E <- 0.5 * q * log(r) +
      0.5 * G * (log(1 + Q / (G - 1)) - log(1 + r * Q / (G - 1)))

    expect_equal(unname(av_summary$fstatistic[1]), Q / q, tolerance = 1e-10)
    expect_equal(attr(av_summary, "f_pvalue"), min(1, exp(-expected_log_E)), tolerance = 1e-10)
  }
})

test_that("under a common within-cluster design the F-test e-value is the Hotelling Bayes factor", {
  # Proposition (common design): with the G/(G-1) CR1 estimator, Q_G is Hotelling's T^2 of
  # the per-cluster OLS estimates, and the e-value is the exact one-sample Bayes factor
  # (determinant form).
  set.seed(42)
  G <- 8; n0 <- 6
  W0 <- cbind(1, matrix(rnorm(n0 * 2), n0, 2))
  X <- W0[rep(seq_len(n0), G), ]
  cl <- rep(seq_len(G), each = n0)
  y <- as.vector(X %*% c(1, 0.3, -0.2)) + rnorm(G * n0) + rep(rnorm(G), each = n0)
  df <- data.frame(y = y, x1 = X[, 2], x2 = X[, 3], cl = cl)
  fit <- lm(y ~ x1 + x2, data = df)
  q <- 2; g_val <- 2

  b <- t(sapply(split(df, df$cl), function(d) coef(lm(y ~ x1 + x2, data = d))[-1]))
  # determinant form of the Bayes factor for the vectors b_g
  bbar <- colMeans(b)
  S0 <- crossprod(sweep(b, 2, bbar)) + G * tcrossprod(bbar)
  S1 <- S0 - G^2 / (g_val + G) * tcrossprod(bbar)
  expected_log_E <- 0.5 * q * log(g_val / (g_val + G)) + 0.5 * G * log(det(S0) / det(S1))

  av_summary <- summary(av(fit, g = g_val, cluster = ~ cl))
  expect_equal(attr(av_summary, "f_pvalue"), min(1, exp(-expected_log_E)), tolerance = 1e-8)

  # and the CR1 standard errors are those of the per-cluster estimates: sd(b_g) / sqrt(G)
  expect_equal(unname(av_summary$coefficients[-1, 2]), unname(apply(b, 2, sd) / sqrt(G)),
               tolerance = 1e-8)
})

test_that("clustered F-test returns p-value 1 when the cluster-robust covariance is singular", {
  fit <- lm(mpg ~ wt + hp + disp + qsec, data = mtcars)
  av_summary <- suppressWarnings(summary(av(fit, cluster = ~ am)))
  expect_equal(attr(av_summary, "f_pvalue"), 1)
  expect_true(is.na(av_summary$fstatistic[1]))
  expect_output(print(av_summary), "p-value: 1")
})
