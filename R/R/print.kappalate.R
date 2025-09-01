#' @title Print kappalate Results
#' @description Prints a kappalate Object (formatted like Stata)
#' @param x A kappalate object
#' @param ... Additional parameters (not used)
#' @export
print.kappalate <- function(x, ...) {
  # Safety check: Ensure x is a kappalate object
  if (!inherits(x, "kappalate")) {
    stop("Error: Input is not a kappalate object.")
  }
  
  # Extract relevant values from the kappalate object
  coefficients <- x$coefficients
  vcov_matrix <- x$vcov_matrix
  N <- x$N
  metadata <- x$metadata
  
  # Safety check: Ensure matrices exist
  if (is.null(coefficients) || is.null(vcov_matrix)) {
    stop("Error: Coefficients or variance-covariance matrix is missing.")
  }
  
  #-----------------------------------------------------------------------------
  # Preliminaries
  # Estimation method description
  if (metadata$zmodel == "logit") {
    estMeth1 <- "LATE estimation using Logit ML:"
  } else if (metadata$zmodel == "probit") {
    estMeth1 <- "LATE estimation using Probit ML:"
  } else if (metadata$zmodel == "cbps") {
    estMeth1 <- "LATE estimation using Logit Covariate Balancing:"
  } else {
    estMeth1 <- "LATE estimation using an unknown method:"
  }
  
  #-----------------------------------------------------------------------------
  # Output header
  cat("\n Call:\n")
  print(metadata$formula)
  
  cat("----------------------------------------------------------------------\n")
  cat("\n", estMeth1, "\n")
  cat("\nWeighting estimation of the LATE\n\n")
  
  cat("Outcome      : ", metadata$depvar, "\n")
  cat("Treatment    : ", metadata$tvar, "\n")
  cat("Instrument   : ", metadata$zvar, "\n")
  cat("Covariates   : ", metadata$xvarsips, "\n")
  cat("Number of obs = ", N, "\n\n")
  
  #-----------------------------------------------------------------------------
  # Extract coefficient names and values
  coef_names <- colnames(coefficients)
  coef_values <- as.numeric(coefficients)
  
  # Ensure that standard errors are computable
  if (ncol(vcov_matrix) != length(coef_values)) {
    stop("Error: Variance-covariance matrix dimensions do not match coefficients.")
  }
  
  # Extract standard errors
  std_errors <- sqrt(diag(vcov_matrix))
  
  # Compute z-values
  z_values <- coef_values / std_errors
  
  # Compute p-values
  p_values <- 2 * (1 - pnorm(abs(z_values)))
  
  # Compute 95% confidence intervals
  conf_low <- coef_values - 1.96 * std_errors
  conf_high <- coef_values + 1.96 * std_errors
  
  #-----------------------------------------------------------------------------
  # Print Stata-style table with `z` moved slightly left
  cat("----------------------------------------------------------------------\n")
  cat(sprintf("%-10s %10s %10s %6s %8s %22s\n",
              " ", "Coefficient", "Std. Err.", "z", "P>|z|", "[95% conf. interval]"))
  cat("----------------------------------------------------------------------\n")
  
  for (i in seq_along(coef_names)) {
    cat(sprintf("%-10s %10.7f %10.7f %6.2f %8.3f %10.7f %10.7f\n",
                coef_names[i], coef_values[i], std_errors[i],
                z_values[i], p_values[i], conf_low[i], conf_high[i]))
  }
  
  cat("----------------------------------------------------------------------\n")
  cat("\nAnalytical standard error.\n")
}
