manual_vcov_CR1 <- function(fit, cluster) {
  X <- model.matrix(fit)
  e <- residuals(fit)
  n <- nrow(X)
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
  adj <- (G / (G - 1)) * ((n - 1) / (n - k))
  adj * (bread %*% meat %*% bread)
}

test_that("vcov_CR1 matches a manual sandwich computation", {
  fit <- lm(mpg ~ wt + hp, data = mtcars)
  cluster_f <- resolve_cluster_factor(fit, mtcars$cyl)

  computed <- vcov_CR1(fit, cluster_f)
  expected <- manual_vcov_CR1(fit, mtcars$cyl)

  G <- attr(computed, "G")
  attr(computed, "G") <- NULL
  expect_equal(unname(computed), unname(expected), tolerance = 1e-10)
  expect_equal(G, length(unique(mtcars$cyl)))
})

test_that("av.lm with cluster produces cluster-robust standard errors and dfd = G - 1", {
  fit <- lm(mpg ~ wt + hp, data = mtcars)
  G <- length(unique(mtcars$cyl))

  av_fit <- av(fit, g = 1, cluster = mtcars$cyl)
  av_summary <- summary(av_fit)

  expected_vcov <- manual_vcov_CR1(fit, mtcars$cyl)
  expected_se <- sqrt(diag(expected_vcov))
  expect_equal(unname(av_summary$coefficients[, 2]), unname(expected_se), tolerance = 1e-10)

  # denominator degrees of freedom used for the coefficient p-values is G - 1
  t2 <- av_summary$coefficients[, 3] ^ 2
  expected_pvals <- p_G_t(log_G_t(t2, G - 1, length(residuals(fit)), 1))
  expect_equal(unname(av_summary$coefficients[, 4]), unname(expected_pvals), tolerance = 1e-10)

  # the F-test denominator degrees of freedom is also G - 1
  expect_equal(unname(av_summary$fstatistic[3]), G - 1)
})

test_that("clustered anytime-valid confidence intervals use dfd = G - 1", {
  fit <- lm(mpg ~ wt + hp, data = mtcars)
  G <- length(unique(mtcars$cyl))

  av_fit <- av(fit, g = 1, cluster = mtcars$cyl)
  se <- summary(av_fit)$coefficients[, 2]
  n <- length(residuals(fit))

  ci <- confint(av_fit, level = 0.95)
  t_rad <- t_radius(1, n, length(se), 0.05, nu = G - 1)

  expected_lo <- coef(fit) - t_rad * se
  expected_hi <- coef(fit) + t_rad * se

  expect_equal(unname(ci[, 1]), unname(expected_lo), tolerance = 1e-10)
  expect_equal(unname(ci[, 2]), unname(expected_hi), tolerance = 1e-10)
})

test_that("cluster and vcov_estimator cannot be combined", {
  fit <- lm(mpg ~ wt + hp, data = mtcars)
  expect_error(
    av(fit, g = 1, cluster = mtcars$cyl, vcov_estimator = "HC1"),
    "only one of"
  )
})

test_that("cluster length must match the number of fitted observations", {
  fit <- lm(mpg ~ wt + hp, data = mtcars)
  expect_error(
    av(fit, g = 1, cluster = mtcars$cyl[1:10]),
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
    av_fit <- av(fit, g = g_val, cluster = mtcars$cyl)
    av_summary <- summary(av_fit)

    expected_vcov <- manual_vcov_CR1(fit, mtcars$cyl)
    se <- sqrt(diag(expected_vcov))
    t <- coef(fit) / se
    G <- length(unique(mtcars$cyl))
    classical_pvals <- 2 * pt(-abs(t), df = G - 1)

    for (i in seq_along(classical_pvals)) {
      expect_gte(av_summary$coefficients[i, 4], classical_pvals[i])
    }
  }
})
