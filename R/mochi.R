#' Moment-based measurement error correction for the concentration index
#' This function returns the moment-based concentration index estimates, correcting for measurement error in the rankingn variable by incorporating partial validation data.
#'
#' @param health vector containing the health outcomes for all observations
#' @param unval_exposure vector (of the same length as \code{health}) containing the error-prone, unvalidated exposure on which all observations will be ranked
#' @param val_exposure vector (of the same length as \code{health}) containing the error-free, validated exposure on which all observations will be ranked. For observations that were not validated, this vector should contain \code{NA}.
#' @param return_naive logical for whether the naive estimate (based only on \code{unval_exposure}) should be returned (if \code{TRUE}). The default is \code{FALSE}, in which case only the moment-based estimate is returned.
#' @param include_se logical for whether jackknife standard error estimated should be returned
#' @return a list with the concentration index estimates requested
#' @export

mochi = function(health, unval_exposure, val_exposure, return_naive = FALSE, include_se = FALSE) {
  c_hat <- mochi_estimate(health, unval_exposure, val_exposure, return_naive) 
  if(include_se) {
    n <- length(health)
    jack_ci_mb <- rep(NA, n)
    for(i in 1:n) {
      Y <- health[-i]
      Xstar <- unval_exposure[-i]
      X <- val_exposure[-i]
      
      # Compute moment-based concentration index
      jack_ci_mb[i] <- mochi_estimate(health = Y, 
                             unval_exposure = Xstar, 
                             val_exposure = X,
                             return_naive = FALSE)[[1]]
    }
    se_jack <- sqrt((n - 1)/n * sum((jack_ci_mb - mean(jack_ci_mb))^2))
  } else {
    se_jack = NA
  }
  
  # Return estimates
  if (return_naive) {
    return(list(
      ci.moment = c_hat[[1]],
      ci.moment.se = se_jack,
      ci.naive = c_hat[[2]]
    ))
  } else {
    return(list(
      ci.moment = c_hat[[1]],
      ci.moment.se = se_jack
    )
    )
  }
}
