## Single-species models (no predation) --------------------------------------
source("R/00_setup.R")

# M fixed at input value
ss_fix <- fit_mod(data_list = dat, inits = NULL, file = NULL, compFun = compFun,
                  estimateMode = 0, msmMode = 0, random_rec = FALSE,
                  phase = TRUE, verbose = 1)

# M estimated with prior
ss_Mprior <- fit_mod(data_list = dat, inits = NULL, file = NULL, compFun = compFun,
                     M1Fun = M1_prior,
                     estimateMode = 0, msmMode = 0, random_rec = FALSE,
                     phase = TRUE, verbose = 1)

# M estimated without prior (free)
ss_Mfree <- fit_mod(data_list = dat, inits = NULL, file = NULL, compFun = compFun,
                    M1Fun = build_M1(M1_model = 1, M1_use_prior = FALSE),
                    estimateMode = 0, msmMode = 0, random_rec = FALSE,
                    phase = TRUE, verbose = 1)

ss_models <- list(ss_fix = ss_fix, ss_Mprior = ss_Mprior, ss_Mfree = ss_Mfree)
lapply(ss_models, function(m) list(jnll = m$quantities$jnll, M1 = m$quantities$M1))

for (nm in names(ss_models)) saveRDS(ss_models[[nm]], file.path("results", paste0(nm, ".rds")))
