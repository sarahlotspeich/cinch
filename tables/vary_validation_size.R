# Read in simulation results
res = read.csv("~/Documents/cinch/simulations/vary_validation_size.csv")

# Check that 1000 replications ran per setting 
table(res$true_ci, res$val_prop) ## Yes! 

# Summarize results per setting
res_summ = res |> 
  dplyr::group_by(true_ci, val_prop) |> 
  dplyr::summarize(
    ## Absolute bias 
    bias_gs_ci = mean(gs_ci - true_ci),
    bias_mb_ci = mean(mb_ci - true_ci),
    bias_nv_ci = mean(nv_ci - true_ci),
    ## Relative (%) bias 
    perc_bias_gs_ci = mean((gs_ci - true_ci) / true_ci) * 100,
    perc_bias_mb_ci = mean((mb_ci - true_ci) / true_ci) * 100,
    perc_bias_nv_ci = mean((nv_ci - true_ci) / true_ci) * 100,
    ## Empirical standard error
    ese_gs_ci = sd(gs_ci),
    ese_mb_ci = sd(mb_ci),
    ese_nv_ci = sd(nv_ci),
    ## Average standard error estimator
    ase_gs_ci = mean(se_gs_ci),
    ase_mb_ci = mean(se_mb_ci),
    ase_nv_ci = mean(se_nv_ci),
    ## 95% confidence interval coverage 
    cp_gs_ci = mean((gs_ci - 1.96 * se_gs_ci) <= true_ci & 
                      true_ci <= (gs_ci + 1.96 * se_gs_ci)),
    cp_mb_ci = mean((mb_ci - 1.96 * se_mb_ci) <= true_ci & 
                      true_ci <= (mb_ci + 1.96 * se_mb_ci)),
    cp_nv_ci = mean((nv_ci - 1.96 * se_nv_ci) <= true_ci & 
                      true_ci <= (nv_ci + 1.96 * se_nv_ci))
  )

# Format for LaTex 
latex_num = function(x, digits = 3) {
  paste0("$", round(x = x, digits = digits), "$") 
}
res_summ |> 
  dplyr::arrange(true_ci, val_prop) |> 
  dplyr::select(true_ci, val_prop, 
                dplyr::contains(c("gs", "nv", "mb"))) |> 
  dplyr::mutate_all(.funs = latex_num) |>
  kableExtra::kbl(
    format = "latex", 
    booktabs = TRUE, 
    escape = FALSE,
    col.names = c(
      "CI", "ValProp", 
      "Bias", "(%)", "ESE", "ASE", "CP", 
      "Bias", "(%)", "ESE", "ASE", "CP",
      "Bias", "(%)", "ESE", "ASE", "CP")
    ) |>
  kableExtra::add_header_above(
    c(" " = 2,
      "Gold Standard" = 5,
      "Naive" = 5, 
      "Moment-Based" = 5)) |> 
  kableExtra::row_spec(0, bold = TRUE)
