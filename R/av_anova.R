av <- function(model, g = 1, ...) {
  UseMethod("av", model)
}

#' @export
av.aov <- function(model, g = 1, ...) {
  attr(model, "g") <- g
  class(model) <- c("avaov", class(model))
  return(model)
}


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

#' @export
print.summary.avaov <- function(x, digits = max(3L, getOption("digits") - 3L), 
                                signif.stars = getOption("show.signif.stars"), ...) {
  
  # Use the standard printing method first
  output <- capture.output({
    NextMethod()
  })
  
  output <- gsub("Pr\\(>F\\)", "p value", output)
  cat(paste(output, collapse = "\n"))

  # Then add information about anytime-validity
  cat("\nAnytime-Valid: TRUE\n")

  invisible(x)
}



# Method for objects of class 'aovlist'
av.aovlist <- function(model, g = 1, ...) {
  # Apply the 'av' function to each element (which are aov objects)
  new_list <- lapply(model, function(x) {
    if (inherits(x, "aov"))
      av(x, g = g)
    else
      x
  })
  # Store the precision parameter in the list
  attr(new_list, "g") <- g
  # Set the class to include 'avaovlist'
  class(new_list) <- c("avaovlist", class(model))
  return(new_list)
}

# Summary method for objects of class 'avaovlist'
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

# Print method for summary objects of class 'avaovlist'
print.summary.avaovlist <- function(x, 
                                    digits = max(3L, getOption("digits") - 3L),
                                    signif.stars = getOption("show.signif.stars"), 
                                    ...) {

  # Capture the standard printing output
  output <- capture.output({
    NextMethod()
  })
  
  # Optionally, rename the p-value column for clarity
  output <- sub("Pr\\(>F\\)", "p value", output)
  cat(paste(output, collapse = "\n"))
  


  invisible(x)
}
