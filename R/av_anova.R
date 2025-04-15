#' Convert an aov Object to Anytime-Valid aov (avaov)
#'
#' Converts an object of class \code{aov} to an anytime-valid version by setting the
#' precision parameter \code{g} as an attribute and updating the class.
#'
#' @param model An \code{aov} object resulting from an ANOVA analysis.
#' @param g An integer precision parameter for anytime-valid inference. Default is 1.
#' @param ... Additional arguments passed to or from other methods.
#'
#' @return An object of class \code{avaov} with anytime-valid p-values.
#' @export
av.aov <- function(model, g = 1, ...) {
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
#' @export
summary.avaov <- function(object, ...) {
  # First get the regular aov summary
  summ <- NextMethod()
  
  # Change the class to include summary.avaov first 
  class(summ) <- c("summary.avaov", class(summ))
  
  # Add g parameter from the original model
  attr(summ, "g") <- attr(object,"g")
  
  # Process each element in the summary to replace p-values with anytime-valid p-values
  if (is.list(summ) && !is.data.frame(summ)) {
    # For multi-stratum designs, summary returns a list of tables
    for (i in seq_along(summ)) {
      if (is.data.frame(summ[[i]])) {
        # Get the F values from the summary table
        f_values <- summ[[i]]$`F value`
        
        # We need the degrees of freedom and sample size for the anytime-valid calculation
        d <- summ[[i]]$Df  # numerator df
        # The error df is in the last row of the table
        error_row <- nrow(summ[[i]])
        nu <- summ[[i]]$Df[error_row]  # denominator df
        
        # Estimate total sample size (sum of df + 1)
        n <- sum(summ[[i]]$Df) + 1
        
        # Calculate anytime-valid p-values for each F test
        av_pvalues <- numeric(length(f_values))
        for (j in 1:length(f_values)) {
          if (!is.na(f_values[j])) {
            # Skip the error row which doesn't have an F value
            if (j < error_row) {
              av_pvalues[j] <- min(1.0, exp(-1 * log_G_f(f_values[j], d[j], nu, n, attr(object, "g"))))
            } else {
              av_pvalues[j] <- NA
            }
          } else {
            av_pvalues[j] <- NA
          }
        }
        
        # Replace the standard p-values with anytime-valid p-values
        summ[[i]]$`Pr(>F)` <- av_pvalues
      }
    }
  } else if (is.data.frame(summ)) {
    # For single-stratum designs, summary returns a single data frame
    f_values <- summ$`F value`
    
    # Get necessary parameters
    d <- summ$Df  # numerator df
    error_row <- nrow(summ)
    nu <- summ$Df[error_row]  # denominator df
    n <- sum(summ$Df) + 1  # estimate total sample size
    
    # Calculate anytime-valid p-values
    av_pvalues <- numeric(length(f_values))
    for (j in 1:length(f_values)) {
      if (!is.na(f_values[j])) {
        # Skip the error row
        if (j < error_row) {
          log_G_f_value = log_G_f(f_values[j], d[j], nu, n, attr(object,"g"))
          av_pvalues[j] <- p_G_f(log_G_f_values)
        } else {
          av_pvalues[j] <- NA
        }
      } else {
        av_pvalues[j] <- NA
      }
    }
    
    # Replace the standard p-values with anytime-valid p-values
    summ$`Pr(>F)` <- av_pvalues
  }
  
  return(summ)
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


#' Convert an aovlist Object to Anytime-Valid aovlist (avaovlist)
#'
#' Applies the anytime-valid conversion to each component (typically \code{aov} objects)
#' within an \code{aovlist}. The precision parameter \code{g} is stored as an attribute for later use.
#'
#' @param model An object of class \code{aovlist} representing a multi-stratum ANOVA model.
#' @param g An integer precision parameter for anytime-valid inference. Default is 1.
#' @param ... Additional arguments passed to or from other methods.
#'
#' @return An enhanced \code{aovlist} object of class \code{avaovlist} with anytime-valid properties.
#' @export
av.aovlist <- function(model, g = 1, ...) {
  # Apply the 'av' function to each element (which are aov objects)
  new_list <- lapply(model, function(x) {
    if (inherits(x, "aov"))
      av(x, g = g)
    else
      x
  })
  attr(new_list, "g") <- g
  class(new_list) <- c("avaovlist", class(model))
  return(new_list)
}


#' Summary Method for Anytime-Valid aovlist Objects
#'
#' Generates a summary for an \code{avaovlist} object by processing each component's summary
#' to replace the standard p-values with anytime-valid p-values calculated using the precision parameter \code{g}.
#'
#' @param object An object of class \code{avaovlist} created by \code{av.aovlist}.
#' @param ... Additional arguments passed to or from other methods.
#'
#' @return A summary object of class \code{summary.avaovlist} where each table has updated anytime-valid p-values.
#' @export
summary.avaovlist <- function(object, ...) {
  # Call the default summary.aovlist to get a list of summary tables
  summ <- NextMethod()
  
  # Attach the precision parameter
  attr(summ, "g") <- attr(object, "g")
  class(summ) <- c("summary.avaovlist", class(summ))
  
  # The summary for aovlist objects is a list of tables.
  # Iterate over each table and update the p-values.
  if (is.list(summ) && !is.data.frame(summ)) {
    for (i in seq_along(summ)) {
      if (is.data.frame(summ[[i]])) {
        # Extract the F values and degrees of freedom from the table
        f_values <- summ[[i]]$`F value`
        d <- summ[[i]]$Df            # numerator degrees of freedom
        error_row <- nrow(summ[[i]])
        nu <- summ[[i]]$Df[error_row]  # denominator degrees of freedom
        n <- sum(summ[[i]]$Df) + 1     # estimated total sample size
        
        # Calculate anytime-valid p-values for each F test in the table
        av_pvalues <- numeric(length(f_values))
        for (j in 1:length(f_values)) {
          if (!is.na(f_values[j])) {
            # Skip the error row which does not have an F value to compute
            if (j < error_row) {
              log_G_f_value = log_G_f(f_values[j], d[j], nu, n, attr(object,"g"))
              av_pvalues[j] <- p_G_f(log_G_f_value)
            } else {
              av_pvalues[j] <- NA
            }
          } else {
            av_pvalues[j] <- NA
          }
        }
        
        # Replace the standard p-values with the anytime-valid p-values
        summ[[i]]$`Pr(>F)` <- av_pvalues
      }
    }
  }
  
  return(summ)
}


#' Print Method for summary.avaovlist Objects
#'
#' Prints the summary for an \code{avaovlist} object by capturing the default output,
#' modifying the p-value header, and displaying the anytime-valid status.
#'
#' @param x An object of class \code{summary.avaovlist} containing the summary of an anytime-valid aovlist.
#' @param digits The number of significant digits to use when printing numeric values.
#' @param signif.stars Logical indicating whether significance stars should be printed.
#' @param ... Additional arguments passed to or from other methods.
#'
#' @return Invisibly returns the printed summary object.
#' @export
print.summary.avaovlist <- function(x, 
                                    digits = max(3L, getOption("digits") - 3L),
                                    signif.stars = getOption("show.signif.stars"), 
                                    ...) {

  output <- capture.output({
    NextMethod()
  })
  output <- sub("Pr\\(>F\\)", "p value", output)
  cat(paste(output, collapse = "\n"))
  invisible(x)
}
