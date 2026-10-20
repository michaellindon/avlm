#' Coefficient Covariance and Degrees of Freedom for Anytime-Valid Inference
#'
#' Computes, for an \code{avlm} object, the coefficient estimates, their covariance matrix and
#' the two quantities the anytime-valid e-values need: the denominator degrees of freedom
#' \code{nu} and the sample size \code{n} on which the precision \code{g} acts.
#'
#' \itemize{
#'   \item Classical: \eqn{s^2 (X'X)^{-1}}, \eqn{\nu = n - k} (Lindon et al., 2026, Section 4).
#'   \item Heteroskedasticity-robust (\code{HC0}--\code{HC3}): the sandwich
#'   \eqn{(X'X)^{-1} X'\Omega X (X'X)^{-1}}, \eqn{\nu = n - k} (Lindon et al., 2026, Section 5).
#'   \item Cluster-robust (\code{CR1}): the sandwich of the summed cluster scores scaled by
#'   \eqn{G/(G-1)}, \eqn{\nu = G - 1} and sample size \eqn{G} (Lindon, 2026, eq. CR1).
#' }
#'
#' Aliased (\code{NA}) coefficients are dropped. For weighted fits, the rows of \eqn{X} and the
#' residuals are scaled by \eqn{\sqrt{w}} and observations with zero weight are dropped.
#'
#' @param object An \code{avlm} object.
#' @param summ Optionally, \code{summary.lm(object)}, to avoid computing it twice.
#'
#' @return A list with elements \code{coefficients} (non-aliased estimates), \code{vcov},
#' \code{nu}, \code{n}, \code{type} (\code{"classical"}, an \code{HC} type, or \code{"CR1"})
#' and \code{G} (number of clusters, or \code{NULL}).
#' @keywords internal
av_vcov <- function(object, summ = NULL) {
  beta <- coef(object)
  ok <- !is.na(beta)
  beta <- beta[ok]
  k <- length(beta)
  nu_obs <- object$df.residual
  n_obs <- nu_obs + k

  vcov_estimator <- attr(object, "vcov_estimator")
  cluster_f <- attr(object, "cluster")

  if (is.null(vcov_estimator) && is.null(cluster_f)) {
    s <- if (is.null(summ)) summary.lm(object) else summ
    V <- (s$sigma^2 * s$cov.unscaled)[names(beta), names(beta), drop = FALSE]
    return(list(coefficients = beta, vcov = V, nu = nu_obs, n = n_obs,
                type = "classical", G = NULL))
  }

  X <- model.matrix(object)[, ok, drop = FALSE]
  e <- residuals(object)
  w <- object$weights
  if (!is.null(w)) {
    keep <- w != 0
    X <- X[keep, , drop = FALSE] * sqrt(w[keep])
    e <- e[keep] * sqrt(w[keep])
    if (!is.null(cluster_f)) cluster_f <- droplevels(cluster_f[keep])
  }
  bread <- chol2inv(chol(crossprod(X)))

  if (!is.null(cluster_f)) {
    scores <- rowsum(X * e, cluster_f, reorder = FALSE)
    G <- nrow(scores)
    if (G < 2) stop("cluster must have at least 2 distinct levels among the fitted observations.")
    if (G <= ncol(X)) {
      warning("Only ", G, " clusters for ", ncol(X), " coefficients: the cluster-robust ",
              "covariance is singular and the resulting inference is not meaningful.")
    }
    V <- (G / (G - 1)) * (bread %*% crossprod(scores) %*% bread)
    nu <- G - 1
    n <- G
    type <- "CR1"
  } else {
    h <- rowSums((X %*% bread) * X)
    omega <- switch(vcov_estimator,
                    "HC0" = e^2,
                    "HC1" = e^2 * n_obs / nu_obs,
                    "HC2" = e^2 / (1 - h),
                    "HC3" = e^2 / (1 - h)^2)
    V <- bread %*% crossprod(X, X * omega) %*% bread
    nu <- nu_obs
    n <- n_obs
    type <- vcov_estimator
    G <- NULL
  }
  dimnames(V) <- list(names(beta), names(beta))
  list(coefficients = beta, vcov = V, nu = nu, n = n, type = type, G = G)
}

#' Validate the Precision Parameter g
#' @param g Candidate value of \code{g}.
#' @return \code{g}, invisibly, if it is a single finite positive number.
#' @keywords internal
check_g <- function(g) {
  if (!is.numeric(g) || length(g) != 1 || !is.finite(g) || g <= 0) {
    stop("`g` must be a single positive number.")
  }
  invisible(g)
}
