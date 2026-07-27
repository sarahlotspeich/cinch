# Load packages
## RUN ONCE: devtools::install_github("sarahlotspeich/mochi")
library(mochi) ## for moment-based correction

# Set useful constants 
sim_seed = 11422 ## For reproducibility, seed to start each setting
num_rep = 1000 ## Number of replications per setting

# Loop over different error variances and concentration indices 
df_res = data.frame( ## initialize empty dataframe to hold results
  true_ci = NA, 
  val_prop = NA,
  gs_lambda = NA,
  gs_ci = NA, 
  se_gs_ci = NA, 
  cc_ci = NA, 
  se_cc_ci = NA,
  mb_ci = NA, 
  se_mb_ci = NA, 
  nv_ci = NA, 
  se_nv_ci = NA
)  
for (ci in c(-0.5, 0, 0.5)) {
  for (val_prop in c(0.05, 0.1, 0.25)) {
    set.seed(sim_seed) ## for reproducibility 
    for (r in 1:num_rep) {
      ## Simulate data 
      dat <- sim_mochi_data(
        var_error = 1, ## fixed: Var(U) = 1 
        var_exposure = 0.6, ## fixed: Var(X) = 0.6
        diff_exposure_error = FALSE, ## fixed: nondifferential errors
        approx_disparity = ci,
        n = 1000, ## fixed: N = 1000
        val_prop = val_prop
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

      ## Calculate oracle/fully validated bias factor 
      varRstar <- var(dat$Rstar)
      covR_Rstar <- cov(dat$R, dat$Rstar)
      lambda <- covR_Rstar / varRstar

      ## Calculate moment-based, complete-case, and naive CI using mochi()
      mochi_res <- mochi(outcome = dat$Y,
                         unval_exposure = dat$Xstar,
                         val_exposure = dat$Xval,
                         include_se = TRUE,
                         return_naive = TRUE,
                         return_cc = TRUE,
                         bootstraps = 0) ### use jackknife
      
      ## Combine, row stack, and save 
      df_res = data.frame(
        true_ci = ci, 
        val_prop = val_prop,
        gs_lambda = lambda,
        gs_ci = oracle_ci, 
        se_gs_ci = se_oracle_ci,
        cc_ci = mochi_res$ci_cc, 
        se_cc_ci = mochi_res$se_ci_cc,
        mb_ci = mochi_res$ci_moment, 
        se_mb_ci = mochi_res$se_ci_moment, 
        nv_ci = mochi_res$ci_naive, 
        se_nv_ci = mochi_res$se_ci_naive) |> 
        dplyr::bind_rows(df_res)
      df_res |> 
        write.csv(file = paste0("simulations/vary_validation_size.csv"), 
                row.names = FALSE)  
    }
  }
}
