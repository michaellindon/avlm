#' Convert an aov Object to Anytime-Valid aov (avaov)
#'
#' Converts an object of class \code{aov} to an anytime-valid version by setting the
#' precision parameter \code{g} as an attribute and updating the class.
#'
#' @param model An \code{aov} object resulting from an ANOVA analysis.
#' @param g A positive precision parameter for anytime-valid inference (see \code{\link{av.lm}}). Default is 1.
#' @param ... Additional arguments passed to or from other methods.
#'
#' @return An object of class \code{avaov} with anytime-valid p-values.
#' @export
av.aov <- function(model, g = 1, ...) {
  check_g(g)
  attr(model, "g") <- g
  class(model) <- c("avaov", class(model))
  return(model)
}


#' Summary Method for Anytime-Valid aov Objects
#'
#' This method produces a summary for objects of class \code{avaov}. It first calls the
#' default \code{summary.aov} method and then replaces the standard p-values with anytime-valid p-values 
#' calculated using the precision parameter \code{g}.
#'
#' @param object An object of class \code{avaov} created by \code{av.aov}.
#' @param ... Additional arguments passed to or from other methods.
#'
#' @return A summary object of class \code{summary.avaov} that includes the anytime-valid p-values.
#' 
#' @examples
#' 
#' # Fit an ANOVA model to the iris dataset.
#' # This model tests whether the sepal length differs by species.
#' fit_aov <- aov(Sepal.Length ~ Species, data = iris)
#'
#' # Convert the standard aov object to an anytime-valid aov (avaov) with precision parameter g = 1.
#' av_fit_aov <- av(fit_aov, g = 1)
#'
#' # Print the summary of the anytime-valid ANOVA model.
#' # The summary replaces standard p-values with anytime-valid p-values.
#' summary(av_fit_aov)
#' 
#' @export
summary.avaov <- function(object, ...) {
  summ <- NextMethod()
  class(summ) <- c("summary.avaov", class(summ))
  g <- attr(object, "g")
  attr(summ, "g") <- g

  nu <- object$df.residual
  n <- nu + object$rank
  # summary.aov returns a list holding one ANOVA table; each term's F-test is replaced by
  # its anytime-valid counterpart (the residual row has no F value and stays NA).
  summ[] <- lapply(summ, function(tab) {
    f <- tab$`F value`
    is_term <- !is.na(f)
    p <- rep(NA_real_, nrow(tab))
    p[is_term] <- p_G_f(log_G_f(f[is_term], tab$Df[is_term], nu, n, g))
    tab$`Pr(>F)` <- p
    tab
  })
  summ
}

#' Print Method for summary.avaov Objects
#'
#' This method prints the summary of an \code{avaov} object. It captures the output
#' from the default printing method, substitutes the header "Pr(>F)" with "p value", and adds
#' a note indicating that anytime-valid inference is used.
#'
#' @param x An object of class \code{summary.avaov}.
#' @param digits The number of significant digits to use when printing. Defaults to a value based on options.
#' @param signif.stars Logical indicating whether significance stars should be printed.
#' @param ... Additional arguments passed to or from other methods.
#'
#' @return Invisibly returns the summary object.
#' @export
print.summary.avaov <- function(x, digits = max(3L, getOption("digits") - 3L), 
                                signif.stars = getOption("show.signif.stars"), ...) {
  
  output <- capture.output({
    NextMethod()
  })
  
  output <- gsub("Pr\\(>F\\)", "p value", output)
  cat(paste(output, collapse = "\n"))
  cat("\nAnytime-Valid: TRUE\n")

  invisible(x)
}


