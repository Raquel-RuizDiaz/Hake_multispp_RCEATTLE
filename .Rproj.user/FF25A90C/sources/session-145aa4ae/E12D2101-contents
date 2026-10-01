## Multispecies models ------------------------------------------------------
source("R/00_setup.R")
ss_fix <- readRDS("results/ss_fix.rds")

## 1) Empirical suitability, M1 estimated (starting values for step 2)
ms_emp <- fit_mod(data_list = dat, inits = ss_fix$estimated_params, file = NULL,
                  compFun = compFun,
                  M1Fun = build_M1(M1_model = 1, updateM1 = FALSE,
                                   M1_use_prior = FALSE, M2_use_prior = FALSE),
                  estimateMode = 0, msmMode = 1, niter = 3, suitMode = 0,
                  suit_styr = c(1980, 1980, 1980), suit_endyr = c(2019, 2019, 2019),
                  initMode = 2, random_rec = FALSE, verbose = 1)

## 2) Lognormal size-based suitability
inits <- ms_emp$estimated_params
map   <- ms_emp$map

# Predator:prey log weight ratio, mean (gam_a) and sd (gam_b); fixed
inits$log_gam_a <- c(0, log(3.7), log(3.1))
inits$log_gam_b <- c(0, log(1.83), log(1.12))

# Vulnerability log_phi[predator, prey]: only ATF->hake and SBF->hake estimated
no_pred <- rbind(c(1, 2), c(2, 2), c(1, 3), c(3, 3), c(2, 3), c(3, 2))
inits$log_phi[no_pred] <- -999
map$mapList$log_phi[] <- seq_along(map$mapList$log_phi)
map$mapList$log_phi[rbind(c(1, 1), no_pred)] <- NA
map$mapFactor$log_phi <- factor(map$mapList$log_phi)

ms_ln       <- fit_ms(dat, inits, map)
ms_ln_refit <- fit_ms(ms_ln$data_list, ms_ln$estimated_params, ms_ln$map)  # reduce M1 gradient

## 3) Estimate diet DM weights for ATF and SBF (hake fixed)
map_w <- ms_ln_refit$map
map_w$mapList$diet_comp_weights <- c(NA, 2, 3)
map_w$mapFactor$diet_comp_weights <- factor(map_w$mapList$diet_comp_weights)

ms_ln_w  <- fit_ms(ms_ln_refit$data_list, ms_ln_refit$estimated_params, map_w)
ms_final <- fit_ms(ms_ln_w$data_list, ms_ln_w$estimated_params, ms_ln_w$map)  # final refit

## Checks
ms_final$quantities$jnll
ms_final$quantities$jnll_comp
ms_final$quantities$M1
ms_final$quantities$vulnerability
ms_final$sdrep
summary(ms_final)
head(grad_table(ms_final))

saveRDS(ms_emp,   "results/ms_emp.rds")
saveRDS(ms_ln,    "results/ms_ln.rds")
saveRDS(ms_final, "results/ms_final.rds")
