# Load packages
## RUN ONCE: devtools::install_github("sarahlotspeich/mochi")
library(mochi) ## for moment-based correction

# Set useful constants 
sim_seed = 11422 ## For reproducibility, seed to start each setting
num_rep = 1000 ## Number of replications per setting

# Loop over different error variances and concentration indices 
df_res = data.frame( ## initialize empty dataframe to hold results
  true_ci = NA, 
  corr_disadv_error = NA,
  gs_ci = NA, 
  se_gs_ci = NA, 
  mb_ci = NA, 
  se_mb_ci = NA, 
  nv_ci = NA, 
  se_nv_ci = NA
)  
for (ci in c(-0.5, 0, 0.5)) {
  for (rho in c(-0.9, -0.6, -0.3, 0.3, 0.6, 0.9)) {
    set.seed(sim_seed) ## for reproducibility 
    for (r in 1:num_rep) {
      ## Simulate data 
      dat <- sim_mochi_data(
        var_error = 1, ## fixed: Var(U) = 1 
        var_exposure = 0.6, ## fixed: Var(X) = 0.6
        corr_exposure_error = rho,
        approx_disparity = ci,
        n = 1000, ## fixed: N = 1000
        val_prop = 0.1 ## fixed: n/N = 0.1
      )
      ## Calculate oracle/fully validated CI
      mu_hat <- mean(dat$Y)
      varR <- var(dat$R)
      fit_ci_oracle <- lm(Y ~ R, data = dat)
      beta1_hat <- fit_ci_oracle$coefficients[2] 
      oracle_ci <- 2 * varR / mu_hat * beta1_hat 
      se_oracle_ci <- delta_method_se(
        outcome = dat$Y, 
        exposure = dat$X)
      
      ## Calculate moment-based and naive CI using mochi()
      mochi_res <- mochi(outcome = dat$Y, 
                         unval_exposure = dat$Xstar, 
                         val_exposure = dat$Xval, 
                         include_se = TRUE, 
                         return_naive = TRUE)
      
      ## Combine, row stack, and save 
      df_res = data.frame(
        true_ci = ci, 
        corr_disadv_error = rho,
        gs_ci = oracle_ci, 
        se_gs_ci = se_oracle_ci, 
        mb_ci = mochi_res$ci_moment, 
        se_mb_ci = mochi_res$se_ci_moment, 
        nv_ci = mochi_res$ci_naive, 
        se_nv_ci = mochi_res$se_ci_naive) |> 
        dplyr::bind_rows(df_res)
      df_res |> 
        write.csv(file = paste0("simulations/vary_disadvantage_error_corr.csv"), 
                row.names = FALSE)  
    }
  }
}