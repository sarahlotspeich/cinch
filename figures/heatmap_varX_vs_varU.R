# Build a grid of SD(X), Corr(X,U) combinations to try
grid <- expand.grid(
  error_var = seq(0.25, 3, by = 0.25),
  x_var = seq(0.25, 3, by = 0.25)
)

## Calculate attenuation factor based on N = 1000 people 
grid$attenuation <- mapply(function(esd, esx) {
  dat <- sim_data(x_sd = sqrt(esx), error_sd = sqrt(esd), n = 1000, approx_ci = 0.5, rho_xu = -0.9)
  cov(dat$R, dat$Rstar) / var(dat$Rstar)
}, grid$error_var, grid$x_var)

# Plot
grid |>
  mutate(error_var = factor(round(error_var, 2)), 
         x_var = factor(round(x_var, 2))) |>
  ggplot() +
  geom_tile(aes(x = x_var, y = error_var, fill = attenuation), 
            color = "white", linewidth = 0.3) +
  geom_text(aes(x = x_var, y = error_var, 
                label = formatC(attenuation, digits = 2, format = "f")),
            size = 3, color = "black") +
  scale_fill_gradientn(
    colors   = c("#fcd061", "white", "#378ADD"),
    values   = scales::rescale(c(-1, 0, 1)),
    limits   = c(-1, 1),
    guide = "none"
  ) + 
  labs(x = "Disadvantage Variance",
       y = "Error Variance") +
  theme_minimal(base_size = 12) +
  theme(panel.grid  = element_blank(), 
        axis.title = element_text(face = "bold")) 
ggsave("~/Documents/mochi/figures/heatmap_varX_vs_varU.png",
       device = "png",
       width = 7,
       height = 5,
       units = "in")
