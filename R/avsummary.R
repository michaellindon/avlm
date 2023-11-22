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
#' @param alpha the desired 1-alpha confidence sequences
#' @param phi the precision of the mixture distribution on the standardized coefficients
#' (regression coefficients divided by the residual standard deviation).
#' If you are frequentist, a good rule of thumb is to set phi to 1/MDE^2 where
#' MDE is your estimate of the size of the standardized coefficients.
#' If you are Bayesian, set phi to be the precision of a Gaussian prior over
#' the standardized coefficients.
#' @return an avsummary.lm object, very similar to a summary.lm object
#' @export
avsummary = function(lmfit, alpha = 0.05, phi = 1) {
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
  spvalues = pmin(1, exp(-1 * log_bf(tstats2, nu, phi, z2)))
  r <- phi / (phi + z2)
  radii <- stderrs * sqrt(nu * ((1 - (r * alpha^2)^(1 / (nu + 1))) /
    max(0, ((r * alpha^2)^(1 / (nu + 1))) - r)))
  lowers = estimates - radii
  uppers = estimates + radii
  
  mod$f_sequential_p_value = multivariate_F_sequential_p_value(lmfit, phi, s2)
  coefs = cbind(
    mod$coefficients[, -4],
    "Seq. p-value" = spvalues,
    "CS lower" = lowers,
    "CS upper" = uppers
  )
  colnames(coefs)[5] = sprintf("%s%% ", 100 * alpha / 2)
  colnames(coefs)[6] = sprintf("%s%% ", 100 * (1 - alpha / 2))
  mod$coefficients = coefs
  return(mod)
}


#' @method
#' @export
print.avsummary.lm = function (x, digits = max(3L, getOption("digits") - 3L), symbolic.cor = x$symbolic.cor, 
                               signif.stars = getOption("show.signif.stars"), ...) 
{
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  resid <- x$residuals
  df <- x$df
  rdf <- df[2L]
  cat(if (!is.null(x$weights) && diff(range(x$weights))) 
    "Weighted ", "Residuals:\n", sep = "")
  if (rdf > 5L) {
    nam <- c("Min", "1Q", "Median", "3Q", "Max")
    rq <- if (length(dim(resid)) == 2L) 
      structure(apply(t(resid), 1L, quantile), dimnames = list(nam, 
                                                               dimnames(resid)[[2L]]))
    else {
      zz <- zapsmall(quantile(resid), digits + 1L)
      structure(zz, names = nam)
    }
    print(rq, digits = digits, ...)
  }
  else if (rdf > 0L) {
    print(resid, digits = digits, ...)
  }
  else {
    cat("ALL", df[1L], "residuals are 0: no residual degrees of freedom!")
    cat("\n")
  }
  if (length(x$aliased) == 0L) {
    cat("\nNo Coefficients\n")
  }
  else {
    if (nsingular <- df[3L] - df[1L]) 
      cat("\nCoefficients: (", nsingular, " not defined because of singularities)\n", 
          sep = "")
    else cat("\nCoefficients:\n")
    coefs <- x$coefficients
    if (any(aliased <- x$aliased)) {
      cn <- names(aliased)
      coefs <- matrix(NA, length(aliased), 4, dimnames = list(cn, 
                                                              colnames(coefs)))
      coefs[!aliased, ] <- x$coefficients
    }
    #printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
    #    na.print = "NA", ...)
    print(coefs, digits=digits)
  }
  cat("\nResidual standard error:", format(signif(x$sigma, 
                                                  digits)), "on", rdf, "degrees of freedom")
  cat("\n")
  if (nzchar(mess <- naprint(x$na.action))) 
    cat("  (", mess, ")\n", sep = "")
  if (!is.null(x$fstatistic)) {
    cat("Multiple R-squared: ", formatC(x$r.squared, digits = digits))
    cat(",\tAdjusted R-squared: ", formatC(x$adj.r.squared, 
                                           digits = digits), "\nF-statistic:", formatC(x$fstatistic[1L], 
                                                                                       digits = digits), "on", x$fstatistic[2L], "and", 
        x$fstatistic[3L], "DF,  Seq. p-value:", format.pval(x$f_sequential_p_value, 
                                                            digits = digits))
    cat("\n")
  }
  correl <- x$correlation
  if (!is.null(correl)) {
    p <- NCOL(correl)
    if (p > 1L) {
      cat("\nCorrelation of Coefficients:\n")
      if (is.logical(symbolic.cor) && symbolic.cor) {
        print(symnum(correl, abbr.colnames = NULL))
      }
      else {
        correl <- format(round(correl, 2), nsmall = 2, 
                         digits = digits)
        correl[!lower.tri(correl)] <- ""
        print(correl[-1, -p, drop = FALSE], quote = FALSE)
      }
    }
  }
  cat("\n")
  invisible(x)
}