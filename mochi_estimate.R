mochi_estimate = function(health, unval_exposure, val_exposure, return_naive) {
  # Save useful constants
  n = length(unval_exposure) ## total number of observations
  nv = sum(!is.na(val_exposure)) ## size of the validation sample
  
  # Define validation indicator
  V = !is.na(val_exposure)
  
  # Calculate rankings
  ## Error-prone exposures' fractional ranks (all observations)
  Rstar = (rank(unval_exposure) - 1) / n + 1 / (2 * n)
  ## Error-free exposures' fractional ranks (validation subsample)
  Rval = rep(x = NA, length = n) ## initialize as NA
  Rval[V] = (rank(val_exposure[V]) - 1) / nv + 1 / (2 * nv)
  
  # Calculate error magnitude (validation subsample)
  Wval = Rstar - Rval ## W = R* - R
  
  # Use validation subset to estimate quantities in the bias factor
  varRval = var(Rval, na.rm = TRUE) ## Var(R)
  varWval = var(Wval, na.rm = TRUE) ## Var(W)
  covRWval = cov(Rval, Wval, use = "complete.obs") ## Cov(R,W)
  lambdahat_varRval = (varRval + covRWval) /
    (varRval + varWval + 2 * covRWval) ## Estimated bias factor
  
  # Estimate the concentratioin index
  ## Naive: Using ranks based on error-prone for all observations
  mu_hat <- mean(health)
  varRstar <- var(Rstar) ## Var(R*) = Var(R)
  fit_ci_naive <- lm(health ~ Rstar)
  beta1star_hat <- fit_ci_naive$coefficients[2]
  ci_premult = 2 * varRstar / mu_hat
  ci_xstar <- ci_premult * beta1star_hat # Error-prone CI
  
  ## Moment-based: Divide naive by estimate of bias factor
  ci_xmb <-  ci_xstar / lambdahat_varRval
  
  # Return estimates
  if (return_naive) {
    return(list(
      ci.moment = as.numeric(ci_xmb),
      ci.naive = as.numeric(ci_xstar))
    )
  } else {
    return(list(
      ci.moment = as.numeric(ci_xmb)
    )
    )
  }
}
