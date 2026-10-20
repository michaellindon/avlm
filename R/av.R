#' Anytime-valid Conversion Generic Function
#'
#' This generic function converts a fitted model object into an anytime-valid version.
#' The conversion is performed by the appropriate S3 method based on the class of the input model.
#'
#' @param model A fitted model object (e.g., an object of class \code{aov} or \code{lm}).
#' @param g A positive precision parameter for anytime-valid inference (see \code{\link{av.lm}}). Defaults to 1.
#' @param ... Additional arguments passed to the method.
#'
#' @return An enhanced version of \code{model} with anytime-valid inference capabilities.
#' @export
av <- function(model, g = 1, ...) {
  UseMethod("av", model)
}


#' @export
av.default <- function(model, g = 1, ...) {
  stop(
    "av() supports objects of class \"lm\" and \"aov\" (single-stratum); got class \"",
    paste(class(model), collapse = "\", \""), "\"."
  )
}
