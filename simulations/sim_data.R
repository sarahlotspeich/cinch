# Function to simulate data 
## Choose standard deviation of errors, total sample size, approx CI, validation proportion, and validation design
sim_data <- function(error_sd, n, approx_ci, pv = 0.1, design = "SRS") {
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
  
  # Error-free exposure
  X <- rnorm(n = n, mean = 1.8, sd = 1)

  # Error-free fractional rank
  R <- (rank(X) - 1) / n + 1 / (2 * n)

  # Health outcome | Fractional rank of true proximity
  eps <- rnorm(n = n, mean = 0, sd = 1)
  Y <- beta0 + beta1 * R + eps
  Ybin <- as.numeric(Y >= median(Y)) ## Ybin = 1 if Y > median and = 0 otherwise

  # Errors and error-prone exposure
  U <- rnorm(n = n, mean = 0, sd = error_sd)
  Xstar <- X + U
  Xstarbin <- as.numeric(Xstar >= median(Xstar)) ## Xstarbin = 1 if Xstar > median and = 0 otherwise

  # Error-prone fractional rank
  Rstar <- (rank(Xstar) - 1) / n + 1 / (2 * n)
  W <- Rstar - R

  dat <- data.frame(Y, X, R, U, Xstar, Rstar, W, Ybin, Xstarbin)
  
  # Partially validate 
  if (design == "SRS") {
    V <- sample_srs(phI = n, 
                    phII = pv * n)
  } else if (design == "CC") {
    V <- sample_cc(dat = dat, 
                   phI = n, 
                   phII = pv * n, 
                   sample_on = "Ybin")
  } else if (design == "BCC") {
    V <- sample_bcc(dat = dat, 
                    phI = n, 
                    phII = pv * n, 
                    sample_on = c("Xstarbin"))
  } else if (design == "ETS") {
    V <- sample_ets(ets_dat = dat$Xstar, 
                    phI = n, 
                    phII = pv * n)
  } else if (design == "RS") {
    V <- sample_resid(formula = Y ~ Rstar,
                      family = "gaussian",
                      dat = dat, 
                      phI = n, 
                      phII = pv * n)
  }
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
