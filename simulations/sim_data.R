# Function to simulate data 
## Choose standard deviation of errors, total sample size, approx CI, validation proportion, and correlation between disadvantage and errors
sim_data <- function(error_sd, n, approx_ci, pv = 0.1, rho_XU = 0) {
  # Get beta0/beta1 params from approx_ci 
  if (approx_ci == -0.5) {
    beta0 <- 2.5
    beta1 <- -3
  } else if (approx_ci == 0) {
    beta0 <- 3
    beta1 <- 0
  } else if (approx_ci == 0.5) {
    beta0 <- -0.5
    beta1 <- 3
  }
  
  # Define the covariance matrix between error-free disadvantage and error 
  Sigma_XU <- matrix(data = c(1, rho_XU * error_sd, rho_XU * error_sd, error_sd ^ 2), 
                     nrow = 2)
  
  # Simulate error-free disadvantage and error from bivariate normal 
  XU <- MASS::mvrnorm(n = n, 
                      mu = c(1.8, 0), 
                      Sigma = Sigma_XU)
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
  Xstarbin <- as.numeric(Xstar >= median(Xstar)) ## Xstarbin = 1 if Xstar > median and = 0 otherwise

  # Error-prone fractional rank
  Rstar <- (rank(Xstar) - 1) / n + 1 / (2 * n)
  W <- Rstar - R

  dat <- data.frame(Y, X, R, U, Xstar, Rstar, W, Ybin, Xstarbin)
  
  # Partially validate (SRS)
  V <- sample_srs(phI = n, 
                  phII = pv * n)
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
  return(data.frame(Y, X, R, U, Xstar, Rstar, W, Xval, Wval, Rval, Ybin, Xstarbin))
}
