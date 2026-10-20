#' Convert an aov Object to Anytime-Valid aov (avaov)
#'
#' Converts an object of class \code{aov} to an anytime-valid version by setting the
#' precision parameter \code{g} as an attribute and updating the class.
#'
#' @param model An \code{aov} object resulting from an ANOVA analysis.
#' @param g A positive precision parameter for anytime-valid inference (see \code{\link{av.lm}}). Default is 1.
#' @param ... Additional arguments passed to or from other methods.
#'
#' @details The ANOVA table uses sequential (Type I) sums of squares, as \code{\link[stats]{summary.aov}}
#' does, so in unbalanced designs each term is tested given the terms listed before it and the
#' p-values depend on the order of the terms. Each row is the anytime-valid version of that
#' sequential F-test. The object also inherits from \code{avlm}, so \code{\link{confint}} and
#' \code{\link{vcov}} give anytime-valid intervals for the coefficients.
#'
#' @return An object of class \code{avaov} (inheriting from \code{avlm}) with anytime-valid p-values.
#' @export
av.aov <- function(model, g = 1, ...) {
  if (inherits(model, "mlm")) {
    stop("av() does not yet support multivariate responses (class \"maov\").")
  }
  check_g(g)
  attr(model, "g") <- g
  class(model) <- c("avaov", "avlm", class(model))
  return(model)
}

# Replace the Pr(>F) column of an ANOVA table (one row per term plus a residual row) with
# anytime-valid p-values; the residual row has no F value and stays NA.
av_anova_pvalues <- function(tab, nu, n, g) {
  f <- tab$`F value`
  is_term <- !is.na(f)
  p <- rep(NA_real_, nrow(tab))
  p[is_term] <- p_G_f(log_G_f(f[is_term], tab$Df[is_term], nu, n, g))
  tab$`Pr(>F)` <- p
  tab
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
  # Not NextMethod(): the next class is avlm, whose summary is not an ANOVA table.
  summ <- summary.aov(object, ...)
  class(summ) <- c("summary.avaov", class(summ))
  g <- attr(object, "g")
  attr(summ, "g") <- g

  nu <- object$df.residual
  n <- nu + object$rank
  # summary.aov returns a list holding one ANOVA table.
  summ[] <- lapply(summ, av_anova_pvalues, nu = nu, n = n, g = g)
  summ
}

#' Anytime-Valid ANOVA Table for an avlm Object
#'
#' Computes the sequential (Type I) ANOVA table of \code{\link[stats]{anova.lm}} and replaces
#' its p-values with anytime-valid ones, using the same e-values as \code{\link{summary.avaov}}.
#' Comparing several models (\code{anova(fit1, fit2)}) is not supported.
#'
#' @param object An \code{avlm} object.
#' @param ... Must be empty.
#'
#' @return An \code{anova} table with anytime-valid p-values.
#' @examples
#' anova(av(lm(mpg ~ wt + hp, data = mtcars), g = 1))
#' @export
anova.avlm <- function(object, ...) {
  if (...length() > 0) {
    stop("anova() on an avlm object takes a single model.")
  }
  tab <- NextMethod()
  attr(tab, "heading") <- c(attr(tab, "heading"),
                            paste0("Anytime-valid p-values (g = ", format(attr(object, "g")), ")\n"))
  av_anova_pvalues(tab, object$df.residual, object$df.residual + object$rank, attr(object, "g"))
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
  cat("\nAnytime-Valid: TRUE (g = ", format(attr(x, "g"), digits = digits), ")\n", sep = "")

  invisible(x)
}


