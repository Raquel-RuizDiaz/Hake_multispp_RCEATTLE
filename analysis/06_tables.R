################# Outputs extraction #############
## Model fit and likelihood component tables ---------------------------------
models <- list("SS: M fixed"  = ss_fix,
               "SS: M prior"  = ss_Mprior,
               "SS: M free"   = ss_Mfree,
               "Multispecies" = ms_final)

tabs <- report_tables(models, model_names = names(models))

## Fit summary (AIC only comparable among models fit to the same data: SS models)
fit_tab <- tabs$model[, c("model", "n_parameters", "joint_nll", "AIC",
                          "max_gradient", "pdHess", "converged")]
fit_tab$M1_hake <- sapply(models, function(m) m$quantities$M1[1])
fit_tab

## Likelihood components (weighted NLL summed over fleets/species), one column per model
lik <- tabs$likelihood
lik$component <- factor(lik$component, levels = c(setdiff(unique(lik$component), "Total"), "Total"))
lik$model     <- factor(lik$model, levels = names(models))
nll_tab <- as.data.frame.matrix(round(tapply(lik$weighted, lik[, c("component", "model")], sum), 2))
nll_tab

write.csv(fit_tab, "results/table_model_fit.csv", row.names = FALSE)
write.csv(nll_tab, "results/table_likelihood_components.csv")


## Diet effective sample size (linear DM: ESS = (1 + theta * N) / (1 + theta), theta = exp(weight))
ms    <- ms_final
theta <- exp(ms$estimated_params$diet_comp_weights)
diet  <- unique(ms$data_list$diet_data[ms$data_list$diet_data$Pred %in% c(2, 3),
                                       c("Pred", "Pred_age", "Sample_size")])
diet$Predator <- ms$data_list$spnames[diet$Pred]
diet$theta    <- theta[diet$Pred]
diet$ESS      <- (1 + diet$theta * diet$Sample_size) / (1 + diet$theta)
diet$ESS_N    <- diet$ESS / diet$Sample_size
diet

