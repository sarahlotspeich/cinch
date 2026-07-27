mochi_estimate = function(outcome, unval_exposure, val_exposure, return_naive, return_cc) {
  # Save useful constants
  n <- length(unval_exposure) ## total number of observations
  nv <- sum(!is.na(val_exposure)) ## size of the validation sample
  
  # Define validation indicator
  V <- !is.na(val_exposure)
  
  # Calculate rankings
  ## Error-prone exposures' fractional ranks (all observations)
  Rstar <- (rank(unval_exposure) - 1) / n + 1 / (2 * n)
  ## Error-free exposures' fractional ranks (validation subsample)
  Rval <- rep(x = NA, length = n) ## initialize as NA
  Rval[V] <- (rank(val_exposure[V]) - 1) / nv + 1 / (2 * nv) ## replace for validated
  
  # Calculate error magnitude (validation subsample)
  Wval <- Rstar - Rval ## W = R* - R
  
  # Use validation subset to estimate quantities in the bias factor
  varRval <- var(Rval, na.rm = TRUE) ## Var(R)
  varWval <- var(Wval, na.rm = TRUE) ## Var(W)
  covRWval <- cov(Rval, Wval, use = "complete.obs") ## Cov(R,W)
  lambdahat_varRval <- (varRval + covRWval) /
    (varRval + varWval + 2 * covRWval) ## Estimated bias factor
  
  # Estimate the concentration index
  ## Naive: Using ranks based on error-prone for all observations
  mu_hat <- mean(outcome)
  varRstar <- var(Rstar) ## Var(R*) = Var(R)
  fit_ci_naive <- lm(outcome ~ Rstar)
  beta1star_hat <- fit_ci_naive$coefficients[2]
  ci_premult <- 2 * varRstar / mu_hat
  ci_xstar <- as.numeric(ci_premult * beta1star_hat) ## Error-prone CI
  
  ## Moment-based: Divide naive by estimate of bias factor
  ci_xmb <- as.numeric(ci_xstar / lambdahat_varRval)
  
  ## Complete-case: Keep only observations with non-missing error-free 
  outcome_cc <- outcome[V]
  val_exposure_cc <- val_exposure[V]
  rank_cc <- (rank(val_exposure_cc) - 1) / nv + 1 / (2 * nv)
  mu_hat_cc <- mean(outcome_cc) 
  varR_cc <- var(rank_cc)
  fit_ci_cc <- lm(outcome_cc ~ rank_cc)
  beta1_hat_cc <- fit_ci_cc$coefficients[2]
  ci_premult <- 2 * varR_cc / mu_hat_cc
  ci_cc <- as.numeric(ci_premult * beta1_hat_cc) ## Complete-case CI
  
  # Return estimates
  if (return_naive) {
    ### Estimate the standard error for the complete-case concentration index 
    ### using jacknknife (sample generally too small for delta method)
    se_ci_naive <- delta_method_se(outcome = outcome, 
                                   exposure = unval_exposure)
    if (return_cc) {
      ### Estimate the standard error for the naive concentration index
      se_ci_cc <- cc_jacknife_se(outcome = outcome_cc, 
                                 val_exposure = val_exposure_cc, 
                                 conf_level = conf_level)
      return(
        list(
          ci_moment = ci_xmb,
          ci_naive = ci_xstar, 
          se_ci_naive = se_ci_naive, 
          ci_cc = ci_cc, 
          se_ci_cc = se_ci_cc[["se"]], 
          lb_ci_cc = se_ci_cc[["lb"]], 
          ub_ci_cc = se_ci_cc[["ub"]]
        )
      )
    } else if (return_cc) {
      ### Estimate the standard error for the naive concentration index
      se_ci_cc <- cc_jacknife_se(outcome = outcome_cc, 
                                 val_exposure = val_exposure_cc, 
                                 conf_level = conf_level)
      return(
        list(
          ci_moment = ci_xmb,
          ci_cc = ci_cc, 
          se_ci_cc = se_ci_cc[["se"]], 
          lb_ci_cc = se_ci_cc[["lb"]], 
          ub_ci_cc = se_ci_cc[["ub"]]
        )
      )
    } else {
      return(
        list(
          ci_moment = as.numeric(ci_xmb),
          ci_naive = as.numeric(ci_xstar), 
          se_ci_naive = as.numeric(se_ci_naive)
        )
      )
    }
  } else {
    return(
      list(
        ci_moment = as.numeric(ci_xmb)
        )
      )
  }
}
