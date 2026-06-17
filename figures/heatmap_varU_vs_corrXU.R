# Build a grid of SD(X), Corr(X,U) combinations to tru
grid <- expand.grid(
  error_var = seq(0.25, 3, by = 0.25),
  rho_xu   = seq(-1, 1, by = 0.1)
)

## Calculate attenuation factor based on N = 1000 people 
grid$attenuation <- mapply(function(esd, rho) {
  dat <- sim_data(x_sd = sqrt(0.6), error_sd = sqrt(esd), n = 1000, approx_ci = 0.5, rho_xu = rho)
  cov(dat$R, dat$Rstar) / var(dat$Rstar)
}, grid$error_var, grid$rho_xu)

# Plot
grid |>
  filter(rho_xu != 0, abs(rho_xu) < 1) |> 
  mutate(error_sd = factor(round(error_var, 2)), ### variance!
         rho_xu   = factor(round(rho_xu, 2))) |>
  ggplot() +
  geom_tile(aes(x = error_sd, y = rho_xu, fill = attenuation), 
            color = "white", linewidth = 0.3) +
  geom_text(aes(x = error_sd, y = rho_xu, label = formatC(attenuation, digits = 2, format = "f")),
            size = 3, color = "black") +
  scale_fill_gradientn(
    colors   = c("#fcd061", "white", "#378ADD"),
    values   = scales::rescale(c(-1, 0, 1)),
    limits   = c(-1, 1),
    guide = "none"
  ) + 
  labs(x = "Error Variance",
       y = "Correlation Between Disadvantage and Error") +
  theme_minimal(base_size = 12) +
  theme(panel.grid  = element_blank(), 
        axis.title = element_text(face = "bold")) 
ggsave("~/Documents/mochi/figures/heatmap_varU_vs_corrXU.png",
       device = "png",
       width = 7,
       height = 5,
       units = "in")
