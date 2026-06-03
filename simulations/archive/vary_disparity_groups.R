# Load packages
## RUN ONCE: devtools::install_github("sarahlotspeich/mochi")
library(mochi) ## for moment-based correction
## RUN ONCE: devtools::install_github("sarahlotspeich/auditDesignR")
library(auditDesignR) ## for validation study designs 

# Source script for simulate_data() function
#devtools::source_url("https://raw.githubusercontent.com/sarahlotspeich/mochi/refs/heads/main/simulations/sim_data.R")
source("~/Documents/cinch/simulations/sim_data.R")

# Set useful constants 
sim_seed = 11422 ## For reproducibility, seed to start each setting
num_rep = 1000 ## Number of replications per setting

# Set constants that are not varied between simulation settings 
N = 1000 ## total sample size
pv = 0.1 ## validation study size
sigma2U = 0.5 ## error variance (less error-prone group)

# Loop over different error variances and concentration indices 
df_res = data.frame( ## initialize empty dataframe to hold results
  true_ci = NA, 
  less_error_var = NA,
  more_error_var = NA,
  gs_ci = NA, 
  se_gs_ci = NA, 
  mb_ci = NA, 
  se_mb_ci = NA, 
  nv_ci = NA, 
  se_nv_ci = NA
)  
for (ci in c(-0.5, 0, 0.5)) {
  for (mult_ci in c(1.25, 1.1, 1)) {
    set.seed(sim_seed) ## for reproducibility 
    for (r in 1:num_rep) {
      ## Simulate data for group 1 (larger disparity)
      dat1 = sim_data(error_sd = sqrt(sigma2U), 
                      n = N / 2, 
                      approx_ci = mult_ci * ci, 
                      pv = pv, 
                      design = "SRS")
      
      ## Calculate oracle/fully validated CI
      mu_hat <- mean(dat1$Y)
      varR <- var(dat1$R)
      fit_ci_oracle <- lm(Y ~ R, data = dat1)
      beta1_hat <- fit_ci_oracle$coefficients[2] 
      oracle_ci1 <- 2 * varR / mu_hat * beta1_hat 
      se_oracle_ci1 <- delta_method_se(
        outcome = dat1$Y, 
        exposure = dat1$X)
      
      ## Calculate moment-based and naive CI using mochi()
      mochi_res1 <- mochi(outcome = dat1$Y, 
                          unval_exposure = dat1$Xstar, 
                          val_exposure = dat1$Xval, 
                          include_se = TRUE, 
                          return_naive = TRUE)
      
      ## Simulate data for group 2 (smaller disparity)
      dat2 = sim_data(error_sd = sqrt(sigma2U), 
                      n = N / 2, 
                      approx_ci = ci, 
                      pv = pv, 
                      design = "SRS")
      
      ## Calculate oracle/fully validated CI
      mu_hat <- mean(dat2$Y)
      varR <- var(dat2$R)
      fit_ci_oracle <- lm(Y ~ R, data = dat2)
      beta1_hat <- fit_ci_oracle$coefficients[2] 
      oracle_ci2 <- 2 * varR / mu_hat * beta1_hat 
      se_oracle_ci2 <- delta_method_se(
        outcome = dat2$Y, 
        exposure = dat2$X)
      
      ## Calculate moment-based and naive CI using mochi()
      mochi_res2 <- mochi(outcome = dat2$Y, 
                          unval_exposure = dat2$Xstar, 
                          val_exposure = dat2$Xval, 
                          include_se = TRUE, 
                          return_naive = TRUE)
      
      ## Combine, row stack, and save 
      df_res = data.frame(
        true_ci = ci, 
        gr1_ci = mult_ci * ci,
        gr1_gs_ci = oracle_ci1, 
        gr1_se_gs_ci = se_oracle_ci1, 
        gr1_mb_ci = mochi_res1$ci_moment, 
        gr1_se_mb_ci = mochi_res1$se_ci_moment, 
        gr1_nv_ci = mochi_res1$ci_naive, 
        gr1_se_nv_ci = mochi_res1$se_ci_naive, 
        gr2_ci = ci,
        gr2_gs_ci = oracle_ci2, 
        gr2_se_gs_ci = se_oracle_ci2, 
        gr2_mb_ci = mochi_res2$ci_moment, 
        gr2_se_mb_ci = mochi_res2$se_ci_moment, 
        gr2_nv_ci = mochi_res2$ci_naive, 
        gr2_se_nv_ci = mochi_res2$se_ci_naive) |> 
        dplyr::bind_rows(df_res)
      df_res |> 
        write.csv(file = paste0("simulations/vary_disparity_groups.csv"), 
                row.names = FALSE)  
    }
  }
}
