cc_jacknife_se <- function(outcome, val_exposure, conf_level) {
  n <- length(outcome) ## total number of observations
  jack_ci <- rep(NA, n) ## initialize empty vector
  
  for(i in 1:n) { ## loop over sample, leaving out each row one-by-one
    ## Define vectors excluding row i 
    Y <- outcome[-i] 
    X <- val_exposure[-i]
    
    ## Compute complete-case concentration index excluding row i
    outcome_cc <- Y[V]
    val_exposure_cc <- X[V]
    rank_cc <- (rank(val_exposure_cc) - 1) / (n - 1) + 1 / (2 * (n - 1))
    mu_hat_cc <- mean(outcome_cc) 
    varR_cc <- var(rank_cc)
    fit_ci_cc <- lm(outcome_cc ~ rank_cc)
    beta1_hat_cc <- fit_ci_cc$coefficients[2]
    ci_premult <- 2 * varR_cc / mu_hat_cc
    jack_ci[i] <-  ci_premult * beta1_hat_cc ## Complete-case CI
  }
  ## Calculate standard error from vector using jacknife formula
  se <- sqrt((n - 1) / n * sum((jack_ci - mean(jack_ci)) ^ 2))
  lb <- quantile(x = jack_ci, probs = (1 - conf_level) / 2)
  ub <- quantile(x = jack_ci, probs = (1 - (1 - conf_level) / 2))
  
  ## Return 
  to_return <- c(se, lb, ub)
  names(to_return) <- c("se", "lb", "ub")
  return(to_return)
}