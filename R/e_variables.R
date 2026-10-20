#' Computes the value of g such that width of the \eqn{1-\alpha} confidence interval
#' at sample size n is minimized
#'
#' Minimizes the radius \eqn{R(g, n, 1, \alpha)} of the asymptotic t confidence sequence (Lindon et al., 2026, Section 5.1) over \eqn{g}.
#'
#' @param n A positive sample size integer. For cluster-robust inference, the number of clusters.
#' @param number_of_coefficients A positive integer of coefficients in the full model
#' @param alpha A positive numeric scalar in (0,1) for nominal Type I error.
#' @param nu Denominator degrees of freedom. Defaults to \code{n - number_of_coefficients};
#' for cluster-robust inference (\code{av(..., cluster = )}) use \code{nu = n - 1} with
#' \code{n} the number of clusters.
#'
#' @return A positive numeric scalar representing the optimal \eqn{g} that minimizes the CI width.
#'
#' @examples
#' n <- 10000
#' alpha <- 0.05
#' g_star <- optimal_g(n, 5, alpha)
#' cat("The optimal g is:", g_star, "\n")
#'
#' # cluster-robust inference with 200 clusters
#' optimal_g(200, alpha = 0.05, nu = 199)
#'
#' @export
optimal_g <- function(n, number_of_coefficients, alpha, nu = n - number_of_coefficients) {
  if (n <= 0) stop("n must be positive.")
  if (nu <= 0) stop("n must be greater than number_of_coefficients")
  if (alpha <= 0 || alpha >= 1) stop("alpha must be in (0,1).")

  # The radius is finite iff g / (g + n) < alpha^(2 / nu).
  upper_bound <- n * alpha^(2 / nu) / (1 - alpha^(2 / nu))
  # Optimize over log(g): the optimum can be far below 1 for small n.
  opt <- optimize(function(log_g) t_radius(exp(log_g), n, alpha, nu),
                  interval = c(log(upper_bound) - 30, log(upper_bound)))
  exp(opt$minimum)
}

# Radius R(g, n, 1, alpha)^(1/2) of the t confidence sequence on the t scale
# (Lindon et al., 2026, Section 5.1; Lindon, 2026, eq. kappa with q = 1).
t_radius <- function(g, n, alpha, nu) {
  r <- g / (g + n)
  u <- (r * alpha^2)^(1 / (nu + 1))
  sqrt(nu * (1 - u) / max(0, u - r))
}

# log e-value for a single coefficient (d = 1), with the g-prior Phi = (g / n) Z'Z.
log_G_t <- function(t2, nu, n, g) {
  r <- g / (g + n)
  0.5 * log(r) + 0.5 * (nu + 1) * (log1p(t2 / nu) - log1p(r * t2 / nu))
}

p_G_t <- function(log_G_t_values) {
  pmin(1, exp(-log_G_t_values))
}

# log e-value for a d-dimensional F-test with the g-prior (Lindon et al., 2026, Sections 4.3.1 and 5).
log_G_f <- function(f, d, nu, n, g) {
  r <- g / (g + n)
  0.5 * d * log(r) + 0.5 * (nu + d) * (log1p((d / nu) * f) - log1p(r * (d / nu) * f))
}

p_G_f <- function(log_G_f_values) {
  pmin(1, exp(-log_G_f_values))
}
