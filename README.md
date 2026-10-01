# Hake–ATF–sablefish multispecies model (CEATTLE)

Code to reproduce the models and figures in *XXX*.

## Requirements
- R ≥ 4.x
- Rceattle v5.45.3: `remotes::install_github("afsc-assessments/Rceattle@5.45.3")`

## Structure
```
data/   300426_SBF_ATF_Hake_Final.xlsx   model input
R/      00_setup.R          data, DM likelihood, shared settings
        01_single_species.R M fixed / M with prior / M free
        02_multispecies.R   empirical → lognormal suitability → final model
        03_figures.R        comparison and fit diagnostics
        04_sensitivity.R    size preference and other-food sensitivity
        05_retrospective.R  retrospective peels, Mohn's rho
        06_tables.R         model fit and likelihood component tables
results/  fitted models (.rds, created by scripts)
figures/  output figures (created by scripts)
```

## Run
Open R in the repository root and run the scripts in order (`01` → `06`).

Species order: 1 = Pacific hake, 2 = arrowtooth flounder, 3 = sablefish.
