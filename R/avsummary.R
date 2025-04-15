#' @param model an lm object resulting from running lm(y~.,data=)
#' @param g the precision of the mixture distribution on the standardized coefficients
#' (regression coefficients divided by the residual standard deviation).
#' If you are frequentist, a good rule of thumb is to set g to 1/MDE^2 where
#' MDE is your estimate of the size of the standardized coefficients.
#' If you are Bayesian, set g to be the precision of a Gaussian prior over
#' the standardized coefficients.
#' @return an "avlm" object
#' @export
av = function(model, g=1){
  class(model) = c("avlm", "lm")
  model$g = g
  return(model)
}

#' @export
summary.avlm <- function(object, ...) {
  # First get the regular lm summary
  summ <- NextMethod()
  
  # Change the class to include summary.avlm first 
  class(summ) <- c("summary.avlm", "summary.lm")
  
  # Add any additional elements specific to avlm
  g <- object$g
  summ$g <- g
  number_of_coefficients = dim(summ$cov.unscaled)[1]
  n = length(summ$residuals)
  t = summ$coefficients[, 't value']
  t2 = t ^ 2
  
  p = number_of_coefficients - 1
  d = 1
  nu = n - p - d

  log_G_t_values = log_G_t(t2, nu, n, g)
  summ$coefficients[,4] = p_G_t(log_G_t_values)
  
  
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
  regression_table <- capture.output({
    printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars, ...)
  })
  # Replace the p-value column header text
  output <- sub("Pr\\(>\\|t\\|\\)", "p value", regression_table)
  
  # Print the modified output
  cat(paste(output, collapse = "\n"))
  
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
    g <- x$g
    f <- x$fstatistic[1]
    n = length(x$residuals)
    d = x$fstatistic[2]
    nu <- x$fstatistic[3] #n - p - d
    f_spvalue = min(1.0, exp(-1*log_G_f(f, d, nu, n, g)))
    
    # Print F-statistic with custom p-value
    cat("F-statistic:", formatC(f, digits = digits), 
        "on", d, "and", nu, "DF,  p-value:", 
        format.pval(f_spvalue, digits = digits), "\n")
  }
  
  # Print g parameter
  cat("\nAnytime-Valid: TRUE\n")
  # cat("Precision parameter g:", format(x$g, digits = digits), "\n")
  
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
  g <- object$g  # The avlm-specific parameter
  n <- length(object$residuals)
  
  # Parameter selection (same as in the standard method)
  if (missing(parm)) {
    parm <- names(coefs)
  } else if (is.numeric(parm)) {
    parm <- names(coefs)[parm]
  }
  
  # Calculate custom confidence intervals that incorporate g
  alpha <- 1 - level
  critval <- sqrt(((g + n) / n) * log((g + n) / (g * alpha^2)))
  
  # Example custom calculation (replace with your own method)
  # Here we're adjusting the standard errors using g
  adjusted_se <- se / sqrt(g)
  
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
      custom_ci[i, 1] <- coefs[idx] - critval * adjusted_se[idx]
      custom_ci[i, 2] <- coefs[idx] + critval * adjusted_se[idx]
    }
  }
  
  return(custom_ci)
}
