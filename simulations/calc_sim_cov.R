N = 1000
set.seed(1)

# Function to simulate data (normal errors)
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
covRW = varW = vector()
for (error_var in c(0.001, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 10)) {
  temp = sim_data(error_sd = sqrt(error_var), 
                  n = N, 
                  approx_ci = ci, 
                  pv = pv, 
                  design = "SRS")
  covRW = append(covRW, with(temp, cov(R, W)))
  varW = append(varW, with(temp, var(W)))
}
data.frame(sigma2U = c(0.001, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 10), 
           varW, covRW) |> 
  ggplot(aes(x = sigma2U)) + 
  geom_line(aes(y = varW), color = "blue") + 
  geom_line(aes(y = covRW), color = "red")

# Function to simulate data (uniform errors)
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
  U <- runif(n = n, min = 0, max = error_sd) 
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
covRW = varW = vector()
for (error_var in c(0.001, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 10)) {
  temp = sim_data(error_sd = sqrt(error_var), 
                  n = N, 
                  approx_ci = ci, 
                  pv = pv, 
                  design = "SRS")
  covRW = append(covRW, with(temp, cov(R, W)))
  varW = append(varW, with(temp, var(W)))
}
data.frame(sigma2U = c(0.001, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 10), 
           varW, covRW) |> 
  ggplot(aes(x = sigma2U)) + 
  geom_line(aes(y = varW), color = "blue") + 
  geom_line(aes(y = covRW), color = "red")

# Function to simulate data (gamma errors)
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
  U <- rgamma(n = n, shape = 1, scale = error_sd) 
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
covRW = varW = vector()
for (error_var in c(0.001, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 10)) {
  temp = sim_data(error_sd = sqrt(error_var), 
                  n = N, 
                  approx_ci = ci, 
                  pv = pv, 
                  design = "SRS")
  covRW = append(covRW, with(temp, cov(R, W)))
  varW = append(varW, with(temp, var(W)))
}
data.frame(sigma2U = c(0.001, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 10), 
           varW, covRW) |> 
  ggplot(aes(x = sigma2U)) + 
  geom_line(aes(y = varW), color = "blue") + 
  geom_line(aes(y = covRW), color = "red")

# Function to simulate data (normal errors)
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
covRW = varW = vector()
for (error_var in c(0.001, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 10)) {
  temp = sim_data(error_sd = sqrt(error_var), 
                  n = N, 
                  approx_ci = ci, 
                  pv = pv, 
                  design = "SRS")
  covRW = append(covRW, with(temp, cov(R, W)))
  varW = append(varW, with(temp, var(W)))
}
data.frame(sigma2U = c(0.001, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 10), 
           varW, covRW) |> 
  ggplot(aes(x = sigma2U)) + 
  geom_line(aes(y = varW), color = "blue") + 
  geom_line(aes(y = covRW), color = "red")

temp = sim_data(error_sd = sqrt(0.1), 
         n = N, 
         approx_ci = ci, 
         pv = pv, 
         design = "SRS")
with(temp, cov(R, W)) ## -0.004146908
with(temp, -cov(R, W) < var(W)) ## TRUE (attenuation)

temp = sim_data(error_sd = sqrt(0.5), 
                n = N, 
                approx_ci = ci, 
                pv = pv, 
                design = "SRS")
with(temp, cov(R, W)) ## -0.01629967
with(temp, -cov(R, W) < var(W)) ## TRUE (attenuation)

temp = sim_data(error_sd = sqrt(1), 
                n = N, 
                approx_ci = ci, 
                pv = pv, 
                design = "SRS")
with(temp, cov(R, W)) ## -0.02891362
with(temp, -cov(R, W) < var(W)) ## TRUE (attenuation)

temp = sim_data(error_sd = sqrt(2), 
                n = N, 
                approx_ci = ci, 
                pv = pv, 
                design = "SRS")
with(temp, cov(R, W)) ## -0.02891362
with(temp, -cov(R, W) < var(W)) ## TRUE (attenuation)

results <- expand.grid(error_sd = c(0.01, 0.05, 0.1),
                       n        = c(50, 100, 200),
                       x_sd     = c(2, 5, 10))

results$met <- mapply(function(esd, nn, xsd) {
  # You'd need to parameterize sd of X in sim_data for this
  set.seed(42)
  X     <- rnorm(nn, mean = 1.8, sd = xsd)
  U     <- rnorm(nn, mean = 0,   sd = esd)
  Xstar <- X + U
  R     <- (rank(X)     - 1) / nn + 1 / (2 * nn)
  Rstar <- (rank(Xstar) - 1) / nn + 1 / (2 * nn)
  W     <- Rstar - R
  -cov(R, W) >= var(W)
}, results$error_sd, results$n, results$x_sd)

results[results$met, ]