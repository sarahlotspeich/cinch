# Read in simulation results
res = read.csv("~/Documents/mochi/simulations/vary_error_variance.csv") |> 
  dplyr::filter(!is.na(true_ci)) ## exclude one all-NA row (didn't run)

# Check that 1000 replications ran per setting 
table(res$true_ci, res$error_var) ## Yes! 

# Read in simulation results
res = read.csv("~/Documents/mochi/simulations/diff_error.csv") |> 
  dplyr::filter(!is.na(true_ci)) ## exclude one all-NA row (didn't run)

# Check that 1000 replications ran per setting 
table(res$true_ci, res$diff_mult) ## Yes! 

# Read in simulation results
res = read.csv("~/Documents/mochi/simulations/vary_validation_size.csv") |> 
  dplyr::filter(!is.na(true_ci)) ## exclude one all-NA row (didn't run)

# Check that 1000 replications ran per setting 
table(res$true_ci, res$val_prop) ## Yes! 

# Read in simulation results
res = read.csv("~/Documents/mochi/simulations/vary_total_size.csv")
# Check that 1000 replications ran per setting 
table(res$true_ci, res$total_size) ## Yes! 

# Read in simulation results
res = read.csv("~/Documents/mochi/simulations/try_bootstrap.csv") |> 
  dplyr::filter(!is.na(true_ci)) ## exclude one all-NA row (didn't run)

# Check that 1000 replications ran per setting 
table(res$true_ci, res$total_size) ## Yes! 