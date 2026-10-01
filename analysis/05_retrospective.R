## Retrospective analysis (5 peels) and Mohn's rho ---------------------------
source("R/00_setup.R")
ms_final <- readRDS("results/ms_final.rds")

retro <- retrospective(Rceattle = ms_final, peels = 5)
retro$mohns
saveRDS(retro, "results/retrospective.rds")

# Use peel end year so each trajectory is clipped correctly
retro_list <- lapply(retro$Rceattle_list, function(x) {
  if (!is.null(x$data_list$endyr_peel)) x$data_list$endyr <- x$data_list$endyr_peel
  x
})

plot_biomass(retro_list, species = 1, legend.pos = "topright", width = 5, height = 3,
             file = "figures/diagnostics/retrospective")
