#' An Anytime-Valid Summary for Linear Models
#'
#' avsummary is similar to summary with some important differences.
#' avsummary provides sequential p-values and confidence sequences which permit
#' analysis to be performed as frequently as desired with no need for multiple 
#' testing corrections (coverage and Type-I guarantees are maintained under
#' continuous monitoring). For further details see
#' https://arxiv.org/abs/2210.08589
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
  spvalues = p_t(tstats2, nu, phi, z2)

  
  mod$f_sequential_p_value = p_F(lmfit, phi, s2)
  coefs = cbind(
    mod$coefficients[, -4],
    "Sp-value" = spvalues
  )
  mod$coefficients = coefs
  return(mod)
}


#' @param lmfit an lm object resulting from running lm(y~.,data=)
#' @param phi the precision of the mixture distribution on the standardized coefficients
#' (regression coefficients divided by the residual standard deviation).
#' If you are frequentist, a good rule of thumb is to set phi to 1/MDE^2 where
#' MDE is your estimate of the size of the standardized coefficients.
#' If you are Bayesian, set phi to be the precision of a Gaussian prior over
#' the standardized coefficients.
#' @return an "avlm" object
#' @export
av = function(lmfit, phi=1){
  class(lmfit) = c("avlm", "lm")
  lmfit$phi = phi
  return(lmfit)
}

#' @export
summary.avlm <- function(object, ...) {
  # First get the regular lm summary
  summ <- NextMethod()
  
  # Change the class to include summary.avlm first 
  class(summ) <- c("summary.avlm", "summary.lm")
  
  # Add any additional elements specific to avlm
  phi <- object$phi
  summ$phi <- phi
  p = summ$fstatistic[2] + 1
  n = summ$fstatistic[3] + p
  
  stderrs = summ$coefficients[, 'Std. Error']
  tstats = summ$coefficients[, 't value']
  tstats2 = tstats ^ 2
  s = summ$sigma
  s2 = s * s
  z2 = (s / stderrs) ^ 2
  nu = n - p - 1
  #summ$f_sequential_p_value = p_F(model, phi, s2)

  summ$coefficients[,4] = p_t(tstats2, nu, phi, z2)
  
  return(summ)
}

#' @export
print.summary.avlm <- function(x, digits = max(3L, getOption("digits") - 3L), 
                               signif.stars = getOption("show.signif.stars"), ...) {
  # Don't call NextMethod() - we'll implement our own printing
  # that closely follows the standard summary.lm printing
  
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  
  # Print residuals
  resid <- x$residuals
  cat("Residuals:\n")
  if (length(resid) > 5) {
    rq <- quantile(resid)
    names(rq) <- c("Min", "1Q", "Median", "3Q", "Max")
    print(rq, digits = digits)
  } else {
    print(resid, digits = digits)
  }
  
  # Print coefficients
  cat("\nCoefficients:\n")
  printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars, ...)
  
  # Print residual standard error
  cat("\nResidual standard error:", format(x$sigma, digits = digits), 
      "on", x$df[2], "degrees of freedom")
  
  # Print NA handling message if needed
  if (nzchar(mess <- naprint(x$na.action))) 
    cat("\n (", mess, ")\n", sep = "")
  else cat("\n")
  
  # Print R-squared and F-statistic with custom p-value
  if (!is.null(x$fstatistic)) {
    cat("Multiple R-squared:  ", formatC(x$r.squared, digits = digits),
        ",\tAdjusted R-squared:  ", formatC(x$adj.r.squared, digits = digits), 
        "\n", sep = "")

    # Calculate custom p-value for F-statistic
    phi <- x$phi
    fstatistic <- x$fstatistic[1]
    numdf <- x$fstatistic[2]
    dendf <- x$fstatistic[3]
    custom_f_pvalue <- 99.9
    
    # Print F-statistic with custom p-value
    cat("F-statistic:", formatC(fstatistic, digits = digits), 
        "on", numdf, "and", dendf, "DF,  p-value:", 
        format.pval(custom_f_pvalue, digits = digits), "\n")
  }
  
  # Print phi parameter
  cat("Precision parameter phi:", format(x$phi, digits = digits), "\n")
  
  invisible(x)
}

#' @export
confint.avlm <- function(object, parm, level = 0.95, ...) {
  # Get the standard confidence intervals from lm
  # We'll use the standard method as a starting point
  orig_confint <- NextMethod()
  
  # Extract needed components for custom CI calculation
  coefs <- coef(object)
  se <- summary(object)$coefficients[, 2]
  phi <- object$phi  # The avlm-specific parameter
  
  # Parameter selection (same as in the standard method)
  if (missing(parm)) {
    parm <- names(coefs)
  } else if (is.numeric(parm)) {
    parm <- names(coefs)[parm]
  }
  
  # Calculate custom confidence intervals that incorporate phi
  alpha <- 1 - level
  critval <- qnorm(1 - alpha/2)
  
  # Example custom calculation (replace with your own method)
  # Here we're adjusting the standard errors using phi
  adjusted_se <- se / sqrt(phi)
  
  # Create the CI matrix
  custom_ci <- matrix(NA, length(parm), 2)
  rownames(custom_ci) <- parm
  colnames(custom_ci) <- c(paste(format(100 * alpha/2, digits = 3), "%"), 
                           paste(format(100 * (1 - alpha/2), digits = 3), "%"))
  
  # Fill in the confidence intervals
  for (i in seq_along(parm)) {
    param_name <- parm[i]
    if (param_name %in% names(coefs)) {
      idx <- which(names(coefs) == param_name)
      custom_ci[i, 1] <- coefs[idx] - 10 * critval * adjusted_se[idx]
      custom_ci[i, 2] <- coefs[idx] + 10 * critval * adjusted_se[idx]
    }
  }
  
  return(custom_ci)
}
