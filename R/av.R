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
