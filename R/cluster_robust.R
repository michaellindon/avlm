#' Resolve a Cluster Argument into a Factor Aligned with the Fitted Model
#'
#' @param model A fitted \code{lm} object.
#' @param cluster A vector of cluster identifiers, one per observation used in the fit,
#' or \code{NULL}.
#'
#' @return A factor of cluster identifiers, or \code{NULL} if \code{cluster} is \code{NULL}.
#' @keywords internal
resolve_cluster_factor <- function(model, cluster) {
  if (is.null(cluster)) return(NULL)

  n_used <- nrow(model.matrix(model))
  if (length(cluster) != n_used) {
    stop(
      "length(cluster) (", length(cluster), ") must equal the number of ",
      "observations used in the fitted model (", n_used, ")."
    )
  }

  cluster_f <- factor(cluster)
  if (nlevels(cluster_f) < 2) {
    stop("cluster must have at least 2 distinct levels.")
  }
  cluster_f
}

#' Compute a CR1 Cluster-Robust Covariance Matrix
#'
#' Computes the cluster-robust ("CRV1") sandwich covariance matrix for the coefficients
#' of a fitted \code{lm} object, following the finite-sample-adjusted estimator of
#' Cameron, Gelbach, and Miller (2011), i.e. the usual sandwich covariance matrix scaled by
#' \eqn{(G/(G-1)) \times ((n-1)/(n-k))}, where \eqn{G} is the number of clusters, \eqn{n} is
#' the number of observations, and \eqn{k} is the number of estimated coefficients.
#'
#' @param model A fitted \code{lm} object.
#' @param cluster_f A factor identifying the cluster membership of each observation used in the fit,
#' as returned by \code{\link{resolve_cluster_factor}}.
#'
#' @return The cluster-robust covariance matrix, with the number of clusters attached as
#' attribute \code{"G"}.
#' @keywords internal
vcov_CR1 <- function(model, cluster_f) {
  X <- model.matrix(model)
  e <- residuals(model)
  n <- nrow(X)
  k <- ncol(X)
  G <- nlevels(cluster_f)

  if (G >= n) {
    stop("Number of clusters must be less than the number of observations.")
  }

  bread <- solve(crossprod(X))
  Xe <- X * e
  meat <- matrix(0, k, k)
  for (lvl in levels(cluster_f)) {
    idx <- which(cluster_f == lvl)
    Xeg <- colSums(Xe[idx, , drop = FALSE])
    meat <- meat + tcrossprod(Xeg)
  }

  adj <- (G / (G - 1)) * ((n - 1) / (n - k))
  vcov <- adj * (bread %*% meat %*% bread)
  dimnames(vcov) <- list(colnames(X), colnames(X))
  attr(vcov, "G") <- G
  vcov
}
