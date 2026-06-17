#' Simulate data to test moment-based measurement error correction for the concentration index
#' 
#' This function returns a dataframe with all columns needed to apply the \code{mochi} function.
#'
#' @param var_exposure scalar for the variance of the exposure. Default is \code{var_exposure = 0.6}. 
#' @param var_error scalar for the variance of the errors in the exposure. Default is \code{var_error = 1}. 
#' @param corr_exposure_error scalar for the correlation between the exposure and its measurement errors. Default is \code{corr_exposure_error = 0.9}.
#' @param approx_disparity scalar for the approximate concentration index underlying the error-free exposure and outcome. Default is \code{approx_disparity = -0.5}.
#' @param n scalar for total sample size generated. Default is \code{n = 1000}. 
#' @param val_prop scalar for proportion of \code{n} to be sampled for the validation study (via simple random sampling). Default is \code{val_prop = 0.1}. 
#' @return a dataframe
#' @export

sim_mochi_data <- function(var_error = 1, var_exposure = 0.6, corr_exposure_error = 0.9, approx_disparity = -0.5, n = 1000, val_prop = 0.1) {
  # Get beta0/beta1 params from approx_disparity 
  if (approx_disparity == -0.5) {
    beta0 <- 2.5
    beta1 <- -3
  } else if (approx_disparity == 0) {
    beta0 <- 3
    beta1 <- 0
  } else if (approx_disparity == 0.5) {
    beta0 <- -0.5
    beta1 <- 3
  }
  
  # Define the covariance matrix between error-free disadvantage and error 
  error_sd <- sqrt(var_error)
  x_sd <- sqrt(var_exposure)
  Sigma_xu <- matrix(data = c(x_sd ^ 2, corr_exposure_error * x_sd * error_sd, 
                              corr_exposure_error * x_sd * error_sd, error_sd ^ 2), 
                     nrow = 2)
  
  # Simulate error-free disadvantage and error from bivariate normal 
  XU <- MASS::mvrnorm(n = n, 
                      mu = c(1.8, 0), 
                      Sigma = Sigma_xu)
  X <- XU[, 1] ## Error-free disadvantage
  U <- XU[, 2] ## Errors in disadvantage

  # Error-free fractional rank
  R <- (rank(X) - 1) / n + 1 / (2 * n)

  # Health outcome | Fractional rank of true proximity
  eps <- rnorm(n = n, mean = 0, sd = 1)
  Y <- beta0 + beta1 * R + eps
  Ybin <- as.numeric(Y >= median(Y)) ## Ybin = 1 if Y > median and = 0 otherwise

  # Error-prone disadvantage
  Xstar <- X + U

  # Error-prone fractional rank
  Rstar <- (rank(Xstar) - 1) / n + 1 / (2 * n)
  W <- Rstar - R

  dat <- data.frame(Y, X, R, U, Xstar, Rstar, W)
  
  # Partially validate (SRS)
  V <- as.numeric(seq(1, n) %in% 
                    sample(x = seq(1, n), size = val_prop * n, replace = F))
  V <- as.logical(V) ## coerce from 0/1 --> FALSE/TRUE 
  Xval <- dat$X ## initialize Xval = X 
  Xval[!V] <- NA ## but then redact Xval if V = FALSE (unvalidated)
  Wval <- dat$W ## initialize Wval = w 
  Wval[!V] <- NA ## but then redact Xval if V = FALSE (unvalidated)
  
  # Calculate ranks based on X in the validation subsample
  nv <- sum(V) ## sample size validated
  Rval <- NA ## initialize as NA
  Rval[V] <- (rank(Xval[V]) - 1) / nv + 1 / (2 * nv)
  
  # Return complete dataset, including mean of Y
  return(data.frame(Y, X, R, U, Xstar, Rstar, W, Xval, Wval, Rval))
}
