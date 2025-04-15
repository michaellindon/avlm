#' Anytime-valid Conversion Generic Function
#'
#' This generic function converts a fitted model object into an anytime-valid version.
#' The conversion is performed by the appropriate S3 method based on the class of the input model.
#'
#' @param model A fitted model object (e.g., an object of class \code{aov} or \code{lm}).
#' @param g An integer precision parameter for anytime-valid inference. Defaults to 1.
#' @param ... Additional arguments passed to the method.
#'
#' @return An enhanced version of \code{model} with anytime-valid inference capabilities.
#' @export
av <- function(model, g = 1, ...) {
  UseMethod("av", model)
}


projection_matrix = function(X){
  QRmat <- qr(X)
  Q <- qr.Q(QRmat)
  return(tcrossprod(Q))
}

MPxM = function(M, X){
  QRmat <- qr(X)
  Q <- qr.Q(QRmat)
  QtM = crossprod(Q,M)
  return(crossprod(QtM))
}

MP1M = function(M){
  n=dim(M)[1]
  M = if(is.null(n)) matrix(M, ncol=1) else M
  Mbar = colSums(M)
  return(tcrossprod(Mbar)/n)
}



#' Computes the value of g such that width of the \eqn{1-\alpha} confidence interval 
#' at sample size n is minimized
#'
#' @param n A positive sample size integer.
#' @param number_of_coefficients A positive integer of coefficients in the full model
#' @param alpha A positive numeric scalar in (0,1) for nominal Type I error.
#'
#' @return A positive numeric scalar representing the optimal \( g \) that minimizes the expression.
#'
#' @examples
#' n <- 10000
#' alpha <- 0.05
#' g_star <- optimal_g(n, 5, alpha)
#' cat("The optimal g is:", g_star, "\n")
#'
#' @export
optimal_g <- function(n, number_of_coefficients, alpha) {
  if (n <= 0 ) stop("n must be positive.")
  if (n <= number_of_coefficients ) stop("n must be greater than number_of_coefficients")
  if (alpha < 0 || alpha > 1) stop("alpha must be in (0,1).")
 
  nu = n - number_of_coefficients
  upper_bound <- n * alpha^(2/nu) / (1 - alpha^(2/nu))
  lower_bound <- 1 # A small positive number to avoid division by zero
  
  opt_result <- optimize(t_radius, interval = c(lower_bound, upper_bound), nu = nu, n = n, alpha = alpha)
  g_star <- opt_result$minimum

  return(g_star)
}