delta_method_se <- function(outcome, exposure) {
  # Save useful constants
  n <- length(exposure) ## total number of observations

  # Calculate rankings based on exposure
  R <- (rank(exposure) - 1) / n + 1 / (2 * n)
  varR <- var(R) ## and their variance 
  
  # Fit linear regression of outcome ~ rank(exposure)
  mu_hat <- mean(outcome) ## mean health outcome
  fit_ci <- lm(outcome ~ R) 
  
  ## Separate pieces needed for standard errors
  premult <- 2 * varR / (mu_hat ^ 2) 
  beta0_hat <- fit_ci$coefficients[1]
  beta1_hat <- fit_ci$coefficients[2]
  cov_hat <- vcov(fit_ci)
  
  ## Calculate standard errors from delta method formula
  var_ci <- premult ^ 2 * 
    (beta1_hat ^ 2 * cov_hat[1, 1] - 
       2 * beta0_hat * beta1_hat * cov_hat[1, 2] + 
       beta0_hat ^ 2 * cov_hat[2, 2])
  
  ## Return 
  return(sqrt(var_ci))
}