## Sensitivity analyses (refit from final model) ----------------------------
source("R/00_setup.R")
ms_final <- readRDS("results/ms_final.rds")

## 1) Size preference (gam_a, gam_b)
refit_gam <- function(gam_a, gam_b) {
  inits <- ms_final$estimated_params
  inits$log_gam_a <- c(0, log(gam_a))
  inits$log_gam_b <- c(0, log(gam_b))
  fit_ms(ms_final$data_list, inits, ms_final$map)
}

sens_gam_median <- refit_gam(gam_a = c(3.35, 2.94), gam_b = c(1.83, 1.12))  # median ratio
sens_gam_broad  <- refit_gam(gam_a = c(3.7, 3.1),   gam_b = c(1.1, 0.5))    # broader preference

## 2) Other food x2
dat_food <- ms_final$data_list
dat_food$other_food <- 2 * dat_food$other_food
sens_food <- fit_ms(dat_food, ms_final$estimated_params, ms_final$map)

saveRDS(sens_gam_median, "results/sens_gam_median.rds")
saveRDS(sens_gam_broad,  "results/sens_gam_broad.rds")
saveRDS(sens_food,       "results/sens_other_food.rds")

## Figures: size preference
mod_list  <- list(ms_final, sens_gam_median, sens_gam_broad)
mod_names <- c("ms model", "gam median", "gam_b broad")

plot_biomass(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE, species = 1)
plot_ssb(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE, species = 1)
plot_recruitment(Rceattle = mod_list, model_names = mod_names, species = 1,
                 width = 5, height = 4, lwd = 2, mod_cex = 1, alpha = 0.3,
                 file = "figures/sensitivity/recruitment_gam")
plot_b_eaten_prop(Rceattle = mod_list, model_names = mod_names, species =1,
                  width = 7, height = 7, lwd = 2)
plot_diet_fit(sens_gam_broad, file = "figures/sensitivity/diet_fit_gam_b_broad")

## Figures: other food
mod_list  <- list(ms_final, sens_food)
mod_names <- c("ms model", "ms model (x2 other food)")

plot_biomass(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE, species = 1)
plot_ssb(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE, species = 1)
plot_recruitment(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE, species = 1)
plot_b_eaten_prop(Rceattle = mod_list, model_names = mod_names,species =1,
                  width = 7, height = 7, lwd = 2,
                  file = "figures/sensitivity/biomass_consumed_other_food")
plot_m_at_age(Rceattle = mod_list, model_names = mod_names, age = 6, species = 1,
              file = "figures/sensitivity/M_age6_other_food")
plot_m2_at_age_prop(Rceattle = mod_list, model_names = mod_names, species = 1, age = 1)
