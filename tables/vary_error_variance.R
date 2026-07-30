# Load packages 
library(dplyr) ## for data manipulation
library(kableExtra) ## for pretty tables 

# Read in simulation results
res = read.csv("~/Documents/mochi/simulations/vary_error_variance.csv") |> 
  filter(!is.na(true_ci)) ## get rid of all NA row used to set dataframe structure

# Check that 1000 replications ran per setting 
table(res$true_ci, res$error_var) ## Yes! 

# Summarize results per setting
n = 1000 # total sample size 
R = (1:n - 1) / n + 1 / (2 * n)
meanR = mean(R)
varR = var(R)
res_summ = res |> 
  group_by(true_ci, error_var) |> 
  mutate(
    approx_ci = true_ci,
    beta0 = case_when(
      true_ci == -0.5 ~ 2.5, 
      true_ci == 0 ~ 3, 
      true_ci == 0.5 ~ -0.5 
    ), 
    beta1 = case_when(
      true_ci == -0.5 ~ -3, 
      true_ci == 0 ~ 0, 
      true_ci == 0.5 ~ 3
    ), 
    true_ci = 2 * varR / (beta0 + beta1 / 2) * beta1 ## only different by 0.0005 from "approx"
  ) |> 
  summarize(
    ## Gold standard lambda 
    mean_gs_lambda = mean(gs_lambda),
    ## Absolute bias 
    bias_gs_ci = mean(gs_ci - true_ci),
    bias_cc_ci = mean(cc_ci - true_ci), 
    bias_mb_ci = mean(mb_ci - true_ci),
    bias_nv_ci = mean(nv_ci - true_ci),
    ## Relative (%) bias 
    perc_bias_gs_ci = 100 * mean((gs_ci - true_ci) / true_ci),
    perc_bias_cc_ci = 100 * mean((cc_ci - true_ci) / true_ci),
    perc_bias_mb_ci = 100 * mean((mb_ci - true_ci) / true_ci),
    perc_bias_nv_ci = 100 * mean((nv_ci - true_ci) / true_ci),
    ## Empirical standard error
    ese_gs_ci = sd(gs_ci),
    ese_cc_ci = sd(cc_ci),
    ese_mb_ci = sd(mb_ci),
    ese_nv_ci = sd(nv_ci),
    ## Average standard error estimator
    ase_gs_ci = mean(se_gs_ci),
    ase_cc_ci = mean(se_cc_ci),
    ase_mb_ci = mean(se_mb_ci),
    ase_nv_ci = mean(se_nv_ci),
    ## 95% confidence interval coverage 
    cp_gs_ci = mean((gs_ci - 1.96 * se_gs_ci) <= true_ci & 
                      true_ci <= (gs_ci + 1.96 * se_gs_ci)),
    cp_cc_ci = mean((cc_ci - 1.96 * se_cc_ci) <= true_ci & 
                      true_ci <= (cc_ci + 1.96 * se_cc_ci)),
    cp_mb_ci = mean((mb_ci - 1.96 * se_mb_ci) <= true_ci & 
                      true_ci <= (mb_ci + 1.96 * se_mb_ci)),
    cp_nv_ci = mean((nv_ci - 1.96 * se_nv_ci) <= true_ci & 
                      true_ci <= (nv_ci + 1.96 * se_nv_ci))
  )

# Format for LaTex 
latex_num = function(x, digits = 3) {
  paste0("$", round(x = x, digits = digits), "$") 
}
latex_num_perc = function(x, digits = 0) {
  paste0("$(", abs(round(x = x, digits = digits)), "\\%)$") 
}
res_summ |> 
  ungroup() |>
  arrange(true_ci, error_var) |> 
  select(error_var, mean_gs_lambda,
                contains(c("gs", "nv", "mb", "cc")), 
                -ase_gs_ci, -cp_gs_ci) |> 
  mutate(
    across(starts_with("perc_bias"), latex_num_perc),
    across(-starts_with("perc_bias"), latex_num)
  ) |> 
  kbl(
    format = "latex", 
    booktabs = TRUE, 
    escape = FALSE,
    col.names = c(
      "ErrorVar", "Lambda",
      "Bias", "(%)", "ESE",  
      "Bias", "(%)", "ESE", "ASE", "CP", 
      "Bias", "(%)", "ESE", "ASE", "CP",
      "Bias", "(%)", "ESE", "ASE", "CP")
    ) |>
  add_header_above(
    c(" " = 2,
      "Gold Standard" = 3,
      "Naive" = 5, 
      "Moment-Based" = 5, 
      "Complete-Case" = 5)) |> 
  row_spec(0, bold = TRUE)
