#' An Anytime-Valid Summary for Linear Models
#'
#' avsummary is similar to summary with some important differences.
#' avsummary provides sequential p-values which permit analysis to be performed 
#' as frequently as desired with no need for multiple testing corrections (coverage 
#' and Type-I guarantees are maintained under continuous monitoring). 
#' For further details see https://arxiv.org/abs/2210.08589
#'
#' @param lmfit an lm object resulting from running lm(y~.,data=)
#' @param phi the precision of the mixture distribution on the standardized coefficients
#' (regression coefficients divided by the residual standard deviation).
#' If you are frequentist, a good rule of thumb is to set phi to 1/MDE^2 where
#' MDE is your estimate of the size of the standardized coefficients.
#' If you are Bayesian, set phi to be the precision of a Gaussian prior over
#' the standardized coefficients.
#' @return an avsummary.lm object, very similar to a summary.lm object
#' @export
avsummary = function(lmfit, phi = 1) {
  mod = summary(lmfit)
  class(mod) = "avsummary.lm"
  p = mod$fstatistic[2] + 1
  n = mod$fstatistic[3] + p
  
  stderrs = mod$coefficients[, 'Std. Error']
  tstats = mod$coefficients[, 't value']
  tstats2 = tstats ^ 2
  estimates = mod$coefficients[, 'Estimate']
  s = mod$sigma
  s2 = s * s
  z2 = (s / stderrs) ^ 2
  nu = n - p - 1
  spvalues = pmin(1, exp(-1 * log_E_t(tstats2, nu, phi, z2)))
  
  # Store parameters for later use in confint method
  mod$cs_params <- list(
    phi = phi,
    nu = nu,
    z2 = z2,
    r = phi / (phi + z2),
    stderrs = stderrs,
    estimates = estimates
  )
  
  mod$f_sequential_p_value = p_F(lmfit, phi, s2)
  
  # Format coefficients table similar to summary.lm (without CIs)
  coefs = cbind(
    mod$coefficients[, -4],  # Keep Estimate, Std. Error, t value
    "Seq. p-value" = spvalues # Add sequential p-values
  )
  
  mod$coefficients = coefs
  
  return(mod)
}

#' @method print avsummary.lm
#' @export
print.avsummary.lm <- function(x, 
                               digits = max(3L, getOption("digits") - 3L), 
                               symbolic.cor = x$symbolic.cor, 
                               ...) {
  # Print the call
  cat("\nCall:\n", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  
  # Print residual summary (similar to summary.lm)
  resid <- x$residuals
  rdf <- x$df[2L]
  cat("Residuals:\n")
  if (rdf > 5L) {
    # Create five-number summary: Min, 1Q, Median, 3Q, Max
    rq <- quantile(resid)
    names(rq) <- c("Min", "1Q", "Median", "3Q", "Max")
    print(rq, digits = digits, ...)
  } else if (rdf > 0L) {
    print(resid, digits = digits, ...)
  } else {
    cat("ALL", x$df[1L], "residuals are 0: no residual degrees of freedom!\n")
  }
  
  # Print coefficients table
  if (length(x$aliased) == 0L) {
    cat("\nNo Coefficients\n")
  } else {
    # Handle possible singularities
    nsingular <- x$df[3L] - x$df[1L]
    if (nsingular)
      cat("\nCoefficients: (", nsingular, " not defined because of singularities)\n", sep = "")
    else
      cat("\nCoefficients:\n")
    
    # Work on the coefficients table:
    # For formatting purposes, rename the sequential p-value column
    coefs <- x$coefficients
    colnames(coefs)[4] <- "Pr(>|t|)"
    
    # Use R’s built-in printCoefmat for standard formatting (alignment, significance stars, etc.)
    printCoefmat(coefs, digits = digits, signif.stars = TRUE, signif.legend = TRUE, ...)
  }
  
  # Print the residual standard error line
  cat("\nResidual standard error:", format(x$sigma, digits), "on", rdf, "degrees of freedom\n")
  
  # Print F-statistic information, if available
  if (!is.null(x$fstatistic)) {
    cat("Multiple R-squared: ", formatC(x$r.squared, digits = digits))
    cat(",\tAdjusted R-squared: ", formatC(x$adj.r.squared, digits = digits), "\n")
    cat("F-statistic:", formatC(x$fstatistic[1L], digits = digits), 
        "on", x$fstatistic[2L], "and", x$fstatistic[3L], "DF,  Seq. p-value:",
        format.pval(x$f_sequential_p_value, digits = digits), "\n")
  }
  
  # Print correlation of coefficients, if available
  if (!is.null(x$correlation)) {
    p <- NCOL(x$correlation)
    if (p > 1L) {
      cat("\nCorrelation of Coefficients:\n")
      if (is.logical(symbolic.cor) && symbolic.cor) {
        print(symnum(x$correlation, abbr.colnames = FALSE))
      } else {
        correl <- format(round(x$correlation, 2), nsmall = 2, digits = digits)
        correl[!lower.tri(correl)] <- ""
        print(correl[-1, -p, drop = FALSE], quote = FALSE)
      }
    }
  }
  cat("\n")
  invisible(x)
}


#' Confidence Intervals for Anytime-Valid Linear Model
#' 
#' Computes confidence sequences for parameters in a fitted anytime-valid linear model.
#' These confidence sequences maintain their coverage guarantees under continuous
#' monitoring, with no need for multiple testing corrections.
#'
#' @param object an avsummary.lm object
#' @param parm a specification of which parameters are to be given confidence intervals,
#'        either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level the confidence level required
#' @param ... additional arguments
#' @return A matrix (or vector) with columns giving lower and upper confidence limits for each parameter.
#' @export
confint.avsummary.lm <- function(object, parm = NULL, level = 0.95, ...) {
  # Extract stored parameters
  params <- object$cs_params
  alpha <- 1 - level  # Convert confidence level to alpha
  
  # Calculate the confidence radii using the stored parameters
  r <- params$r
  radii <- params$stderrs * sqrt(params$nu * ((1 - (r * alpha^2)^(1 / (params$nu + 1))) /
                                                max(0, ((r * alpha^2)^(1 / (params$nu + 1))) - r)))
  
  # Compute confidence limits
  lower <- params$estimates - radii
  upper <- params$estimates + radii
  
  # Create matrix of confidence limits
  cis <- cbind(lower, upper)
  dimnames(cis) <- list(names(params$estimates), c(paste(100 * alpha/2, "%"), paste(100 * (1 - alpha/2), "%")))
  
  # Handle subset of parameters if requested
  if (!is.null(parm)) {
    if (is.numeric(parm)) {
      cis <- cis[parm, , drop = FALSE]
    } else {
      cis <- cis[parm, , drop = FALSE]
    }
  }
  
  return(cis)
}