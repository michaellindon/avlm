#' Resolve a Cluster Argument into a Factor Aligned with the Fitted Model
#'
#' @param model A fitted \code{lm} object.
#' @param cluster \code{NULL}, a vector of cluster identifiers with one entry per observation
#' used in the fit, a one-sided formula such as \code{~ cyl}, or a single column name such as
#' \code{"cyl"}. Formulas and column names are looked up in the \code{data} the model was fit
#' with and subset to the observations used in the fit.
#'
#' @return A factor of cluster identifiers, or \code{NULL} if \code{cluster} is \code{NULL}.
#' @keywords internal
resolve_cluster_factor <- function(model, cluster) {
  if (is.null(cluster)) return(NULL)

  if (inherits(cluster, "formula") || (is.character(cluster) && length(cluster) == 1)) {
    cluster <- cluster_from_data(model, cluster)
  }

  n_used <- nrow(model.matrix(model))
  if (length(cluster) != n_used) {
    stop(
      "length(cluster) (", length(cluster), ") must equal the number of ",
      "observations used in the fitted model (", n_used, ")."
    )
  }
  if (anyNA(cluster)) {
    stop("cluster must not contain missing values.")
  }

  cluster_f <- factor(cluster)
  if (nlevels(cluster_f) < 2) {
    stop("cluster must have at least 2 distinct levels.")
  }
  cluster_f
}

#' Look Up a Cluster Variable in the Data a Model Was Fit With
#'
#' @param model A fitted \code{lm} object.
#' @param cluster A one-sided formula such as \code{~ cyl}, or a single column name.
#'
#' @return The cluster variable, subset to the observations used in the fit.
#' @keywords internal
cluster_from_data <- function(model, cluster) {
  if (inherits(cluster, "formula")) {
    if (length(cluster) != 2) {
      stop("cluster formula must be one-sided, e.g. `cluster = ~ id`.")
    }
    rhs <- cluster[[2]]
    if (is.call(rhs) && as.character(rhs[[1]]) %in% c("+", "*", ":", "^", "|", "/")) {
      stop("Only one-way clustering is supported; cluster formula must contain a single variable.")
    }
    label <- deparse(rhs)
  } else {
    label <- cluster
  }

  data_expr <- model$call$data
  if (is.null(data_expr)) {
    stop(
      "cluster = ", if (inherits(cluster, "formula")) "~ " else "", label,
      " requires the model to be fit with a `data` argument; ",
      "otherwise pass the vector of cluster identifiers."
    )
  }
  data <- eval(data_expr, envir = environment(formula(model)))

  if (inherits(cluster, "formula")) {
    values <- eval(rhs, envir = data, enclos = environment(cluster))
  } else {
    if (!label %in% names(data)) {
      stop("cluster variable '", label, "' not found in the model's data.")
    }
    values <- data[[label]]
  }
  if (length(values) != nrow(data)) {
    stop("cluster variable '", label, "' must have one value per row of the model's data.")
  }

  # model.frame() keeps the row names of the rows used in the fit, so this drops
  # rows removed by `subset` or by the `na.action`.
  rows_used <- match(rownames(model.frame(model)), rownames(data))
  if (anyNA(rows_used)) {
    stop("Could not align cluster variable '", label, "' with the observations used in the fit.")
  }
  values[rows_used]
}

#' Log Cluster-Robust e-Value for a Joint Wald Test
#'
#' The cluster-robust e-process of Lindon (2026), eq. (ECR), with \eqn{\nu = G - 1} and exponent
#' \eqn{G/2} for any number \eqn{q} of tested coefficients:
#' \deqn{E_G = \left(\frac{g}{g+G}\right)^{q/2}
#' \left(\frac{1 + Q/(G-1)}{1 + \frac{g}{g+G} Q/(G-1)}\right)^{G/2}.}
#'
#' @param Q The cluster-robust Wald statistic.
#' @param q The number of tested coefficients.
#' @param G The number of clusters.
#' @param g The precision parameter, on the cluster scale.
#'
#' @return The log e-value.
#' @references
#' Lindon, M. (2026). Anytime-Valid Inference in Linear Models with Clustered Outcomes: Exact Multivariate Tests and Cluster-Robust Confidence Sequences. arXiv:2610.25613. \url{https://arxiv.org/abs/2610.25613}
#' @keywords internal
log_E_CR <- function(Q, q, G, g) {
  r <- g / (g + G)
  nu <- G - 1
  0.5 * q * log(r) + 0.5 * G * (log1p(Q / nu) - log1p(r * Q / nu))
}

#' Check Whether a Covariance Matrix Is Numerically Singular
#'
#' Scale-invariant: checks the smallest eigenvalue of the implied correlation matrix.
#'
#' @param V A symmetric covariance matrix.
#' @param tol Tolerance for the smallest eigenvalue of the correlation matrix.
#'
#' @return \code{TRUE} if \code{V} is numerically singular.
#' @keywords internal
is_singular_vcov <- function(V, tol = 1e-10) {
  s <- sqrt(diag(V))
  if (any(!is.finite(s)) || any(s <= 0)) return(TRUE)
  C <- V / tcrossprod(s)
  min(eigen(C, symmetric = TRUE, only.values = TRUE)$values) < tol
}
