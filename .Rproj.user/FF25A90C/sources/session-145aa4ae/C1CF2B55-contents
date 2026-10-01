## Setup: packages, data, shared model settings -----------------------------
library(remotes)
remotes::install_github("afsc-assessments/Rceattle@5.45.3")

library(Rceattle)
library(ggplot2)
library(dplyr)
library(patchwork)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
dir.create("figures/diagnostics", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/sensitivity", recursive = TRUE, showWarnings = FALSE)

# Species: 1 = hake, 2 = arrowtooth flounder (ATF), 3 = sablefish (SBF)
dat <- read_data(file = "data/300426_SBF_ATF_Hake_Final.xlsx")
dat$index_data <- dat$index_data %>%
  dplyr::select(-Q_block)

# Dirichlet-multinomial (DM) likelihood for age and diet composition
dat$fleet_control$Comp_distribution <- "DirichletMultinomial"
dat$Diet_distribution <- rep(1, dat$nspp)

# DM weights: one per fleet / predator, lognormal(0, 2) prior
compFun <- build_composition(linkages = list(
  theta_comp = linkage_spec(formula = ~ 1, by = ~ fleet,
                            fleet = dat$fleet_control$Fleet_code,
                            priors = list(`(Intercept)` = prior_lognormal(0, 2))),
  theta_diet = linkage_spec(formula = ~ 1, by = ~ species,
                            species = seq_len(dat$nspp),
                            priors = list(`(Intercept)` = prior_lognormal(0, 2)))))

# Residual natural mortality (M1): estimated, prior N(0.2, 0.1)
M1_prior <- build_M1(M1_model = 1, M1_use_prior = TRUE, M_prior = 0.2, M_prior_sd = 0.1)

# Multispecies fit with lognormal size-based suitability (ATF, SBF)
fit_ms <- function(data_list, inits, map) {
  fit_mod(data_list = data_list, inits = inits, map = map, file = NULL,
          M1Fun = M1_prior, compFun = compFun,
          estimateMode = 0,
          msmMode = 1, niter = 3,                 # set multispecies model
          suitMode = c(0, 4, 4),                  # empirical (hake), lognormal (ATF, SBF)
          suit_styr = c(1980, 2013, 2005),
          suit_endyr = c(2019, 2018, 2008),
          initMode = 2,                           # fished start with initial deviations
          random_rec = FALSE, loopnum = 5, phase = TRUE, verbose = 1)
}

# Absolute gradients, largest first
grad_table <- function(fit) {
  gr <- data.frame(par = names(fit$obj$par), est = fit$obj$par,
                   gr = abs(fit$obj$gr()[1, ]))
  gr[order(-gr$gr), ]
}

