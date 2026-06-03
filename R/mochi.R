#' Moment-based measurement error correction for the concentration index
#' This function returns the moment-based concentration index estimates, correcting for measurement error in the ranking variable by incorporating partial validation data.
#'
#' @param outcome vector containing the outcome outcomes for all observations
#' @param unval_exposure vector (of the same length as \code{outcome}) containing the error-prone, unvalidated exposure on which all observations will be ranked
#' @param val_exposure vector (of the same length as \code{outcome}) containing the error-free, validated exposure on which all observations will be ranked. For observations that were not validated, this vector should contain \code{NA}.
#' @param return_naive logical for whether the naive estimate (based only on \code{unval_exposure}) should be returned (if \code{TRUE}). The default is \code{FALSE}, in which case only the moment-based estimate is returned.
#' @param include_se logical for whether jackknife standard error estimated should be returned
#' @return a list with the concentration index estimates requested (and standard errors, where applicable)
#' @export

mochi <- function(outcome, unval_exposure, val_exposure, return_naive = FALSE, include_se = FALSE) {
  # Estimate concentration index 
  c_hat <- mochi_estimate(outcome = outcome, 
                          unval_exposure = unval_exposure, 
                          val_exposure = val_exposure, 
                          return_naive = return_naive) 
  
  # If requested, calculate jacknife standard error
  if(include_se) {
    n <- length(outcome) ## total number of observations
    jack_ci_mb <- rep(NA, n) ## initialize empty vector
    for(i in 1:n) { ## loop over sample, leaving out each row one-by-one
      ## Define vectors excluding row i 
      Y <- outcome[-i] 
      Xstar <- unval_exposure[-i]
      X <- val_exposure[-i]
      
      ## Compute moment-based concentration index excluding row i
      jack_ci_mb[i] <- mochi_estimate(
        outcome = Y, 
        unval_exposure = Xstar,
        val_exposure = X,
        return_naive = FALSE)[[1]]
    }
    ## Calculate standard error from vector using jacknife formula
    se_jack <- sqrt((n - 1) / n * sum((jack_ci_mb - mean(jack_ci_mb)) ^ 2))
  } else {
    se_jack <- NA
  }
  
  # Return estimates
  if (return_naive) {
    return(
      list(
        ci_moment = c_hat[["ci_moment"]],
        se_ci_moment = se_jack,
        ci_naive = c_hat[["ci_naive"]], 
        se_ci_naive = c_hat[["se_ci_naive"]]
        )
      )
  } else {
    return(
      list(
        ci_moment = c_hat[["ci_moment"]],
        se_ci_moment = se_jack
        )
    )
  }
}
