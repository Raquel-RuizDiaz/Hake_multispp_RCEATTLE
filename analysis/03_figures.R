## Model comparison and diagnostics (hake = species 1) -----------------------
source("R/00_setup.R")

ss_fix    <- readRDS("results/ss_fix.rds")
ss_Mprior <- readRDS("results/ss_Mprior.rds")
ss_Mfree  <- readRDS("results/ss_Mfree.rds")
ms_final  <- readRDS("results/ms_final.rds")

mod_list  <- list(ss_fix, ss_Mprior, ss_Mfree, ms_final)
mod_names <- c("ss_fix", "ss_Mprior", "ss_Mfree", "ms_Mprior")

## Single- vs multispecies
plot_ssb(Rceattle = mod_list, model_names = mod_names, species = 1)
plot_biomass(Rceattle = mod_list, model_names = mod_names, species = 1)
plot_recruitment(Rceattle = mod_list, model_names = mod_names, species = 1)
plot_depletionSSB(Rceattle = mod_list, model_names = mod_names, species = 1)

## Predation
plot_m2_at_age_prop(Rceattle = list(ms_final), model_names = "ms_Mprior", species = 1)
plot_b_eaten_prop(Rceattle = list(ms_final), model_names = "ms_Mprior", species = 1)
plot_mortality(ms_final, type = 3) #, width = 8, height = 8,
              # file = "figures/M1_M2_plot")

## Fit diagnostics
plot_diet_comp2(ms_final, file = "figures/diagnostics/diet_fit")
plot_catch(Rceattle = mod_list, model_names = mod_names, width = 5, height = 4,
           lwd = 2, alpha = 0.3, file = "figures/diagnostics/catch")
plot_index(Rceattle = mod_list, model_names = mod_names, width = 5, height = 4,
           file = "figures/diagnostics/survey_index")
plot_logindex(Rceattle = mod_list, model_names = mod_names, width = 5, height = 3,
              file = "figures/diagnostics/survey_logindex")
plot_indexresidual(Rceattle = mod_list, model_names = mod_names, width = 5.5,
                   height = 4, top_adj = -0.8,
                   file = "figures/diagnostics/survey_index_residual")
plot_selectivity(ms_final, width = 4, height = 4,
                 file = "figures/diagnostics/selectivity")

plot_comp(ms_final, cex = 3, lwd = 3, file = "figures/diagnostics/comps")


#### FIGURES #########################################################################################
## 3-D perspective plot of predation mortality (M2) at age by year, with M1 label
## (replaces plot_mortality(type = 3), removed in newer Rceattle versions)
plot_m_persp <- function(fit, species = 1, sex = 1, incl_proj = FALSE, maxage = NULL,
                         zlim = NULL, theta = 155, phi = 20, m1 = NULL,
                         file = NULL, width = 8, height = 8) {
  q  <- fit$quantities
  dl <- fit$data_list

  get_m <- function(nm) {
    a <- q[[paste0(nm, "_at_age")]]
    if (is.null(a)) stop(nm, "_at_age not found in fit$quantities")
    if (length(dim(a)) == 4) a[species, sex, , , drop = TRUE] else a[species, , , drop = TRUE]
  }

  # Ages and years
  nages <- dl$nages[species]
  ages  <- dl$minage[species] + seq_len(nages) - 1
  yrs   <- dl$styr:(if (incl_proj) dl$projyr else dl$endyr)

  z  <- get_m("M2")[seq_len(nages), seq_along(yrs), drop = FALSE]
  # M1 label: age-invariant M1 for this species (arrays are zero-padded to the max ages)
  if (is.null(m1)) {
    m1 <- if (!is.null(q[["M1"]])) {
      a <- q[["M1"]]
      if (is.null(dim(a))) a[species] else if (length(dim(a)) == 2) a[species, 1] else a[species, sex, 1]
    } else {
      a <- get_m("M1")
      if (is.null(dim(a))) a[1] else a[1, 1]
    }
  }
  if (!is.null(maxage)) { keep <- ages <= maxage; z <- z[keep, , drop = FALSE]; ages <- ages[keep] }
  if (is.null(zlim)) zlim <- c(0, max(z))

  if (!is.null(file)) {
    png(paste0(file, ".png"), width = width, height = height, units = "in", res = 300)
    on.exit(dev.off())
  }
  par(mar = c(1, 1, 1, 1))
  persp(x = ages, y = yrs, z = z, zlim = zlim, theta = theta, phi = phi,
        col = "white", border = "black", ticktype = "detailed",
        xlab = "Age", ylab = "Year", zlab = "M2")
  legend("topleft", legend = sprintf("M1 = %.3f", m1), bty = "n", inset = c(0.08, 0.1))
  invisible(z)
}

plot_m_persp(ms_final, species = 1, file = "figures/M1_M2_plot")  # 3-D M2 at age, M1 label

## Composition Pearson residuals without the +/-6 cap used by plot_comp()
## Same residuals as plot_comp() (DM variance included); bubble area ~ |residual|
plot_pearson_full <- function(fit, species = NULL, file = NULL, width = 9, height = 6.5) {
  d <- Rceattle:::.comp_resid_long(fit, species = species)
  d <- d[is.finite(d$pearson), , drop = FALSE]
  d$sign <- ifelse(d$pearson >= 0, "positive", "negative")

  p <- ggplot2::ggplot(d, ggplot2::aes(Year, bin)) +
    ggplot2::geom_point(ggplot2::aes(size = abs(pearson), colour = sign), alpha = 0.7) +
    ggplot2::scale_colour_manual(values = c(positive = "#d7301f", negative = "#2c7fb8"),
                                 guide = "none") +
    ggplot2::scale_size_area(max_size = 8, breaks = c(1, 2, 4, 8, 12), name = "|Pearson|") +
    ggplot2::facet_wrap(~source, scales = "free_y") +
    ggplot2::labs(x = "Year", y = "Age bin") +
    ggplot2::theme_bw(base_size = 10)

  if (!is.null(file)) ggplot2::ggsave(paste0(file, ".png"), p, width = width, height = height, dpi = 300)
  print(p)
  invisible(d[order(-abs(d$pearson)), ])
}

big_resid <- plot_pearson_full(ms_final, file = "figures/diagnostics/comps_pearson_uncapped")
head(big_resid, 5)   # largest residuals

## Observed vs predicted proportion of hake in ATF and SBF diets, by predator age group
## Predictions: predator-N-weighted mean over ages in each group (mean N over suitability years)
library(ggplot2)
library(patchwork)

# Predator age groups
atf_groups <- list(breaks = c(0, 2, 7, 16, Inf),  labels = c("1-2", "3-7", "8-16", "17+"))
sbf_groups <- list(breaks = c(0, 10, 20, 31, Inf), labels = c("1-10", "11-20", "21-31", "32+"))

# Mean predator N-at-age (sexes summed) over the diet-fitting years
pred_N <- function(fit, sp) {
  dl  <- fit$data_list
  yrs <- (dl$suit_styr[sp]:dl$suit_endyr[sp]) - dl$styr + 1
  N   <- fit$quantities$avgN_at_age[sp, seq_len(dl$nsex[sp]), seq_len(dl$nages[sp]), yrs, drop = FALSE]
  apply(N, 3, function(x) sum(x) / length(yrs))
}

# SD of predicted diet proportion (diet_hat[, 2]) from sdrep; NA if not reported
diet_hat_sd <- function(fit) {
  nr  <- nrow(fit$data_list$diet_data)
  sdv <- fit$sdrep$sd[names(fit$sdrep$value) == "diet_hat"]
  nc  <- ncol(as.matrix(fit$quantities$diet_hat))
  if (length(sdv) == nr * nc && nc >= 2) return(sdv[(nr + 1):(2 * nr)])  # column 2 of diet_hat
  if (length(sdv) == nr) return(sdv)
  warning("diet_hat SDs not found in sdrep (", length(sdv), " values for ", nr,
          " diet rows); predicted intervals omitted.", call. = FALSE)
  rep(NA_real_, nr)
}

diet_fit_table <- function(fit, pred_sp, groups, prey_sp = 1) {
  dd <- fit$data_list$diet_data
  d  <- Rceattle:::.diet_plot_data(fit)            # obs, Est, row-aligned with diet_data
  d$Est_sd <- diet_hat_sd(fit)
  if (all(c("lower_95", "upper_95") %in% names(dd))) d[, c("lower_95", "upper_95")] <- dd[, c("lower_95", "upper_95")]
  d  <- d[d$Pred == pred_sp & d$Prey == prey_sp, ]

  N  <- pred_N(fit, pred_sp)
  d$age_group <- cut(d$Pred_age, groups$breaks, groups$labels)
  d$w <- N[d$Pred_age - fit$data_list$minage[pred_sp] + 1]

  do.call(rbind, lapply(split(d, d$age_group, drop = TRUE), function(g) data.frame(
    age_group  = g$age_group[1],
    obs        = g$Observed[1],                      # observed diet is reported per age group
    obs_lower  = g$lower_95[1],
    obs_upper  = g$upper_95[1],
    pred       = weighted.mean(g$Est, g$w),
    pred_sd    = sqrt(weighted.mean(g$Est_sd^2, g$w)))))
}

make_diet_plot <- function(tab, title) {
  long <- rbind(
    data.frame(age_group = tab$age_group, Source = "Observed",  Mean = tab$obs,
               Lower = tab$obs_lower, Upper = tab$obs_upper),
    data.frame(age_group = tab$age_group, Source = "Predicted", Mean = tab$pred,
               Lower = pmax(0, tab$pred - 1.96 * tab$pred_sd),
               Upper = tab$pred + 1.96 * tab$pred_sd))
  y_max <- max(long$Mean, long$Upper, na.rm = TRUE) * 1.25
  pd <- position_dodge(width = 0.4)

  ggplot(long, aes(age_group, Mean, colour = Source, shape = Source, group = Source)) +
    geom_errorbar(aes(ymin = Lower, ymax = pmin(y_max, Upper)), width = 0.2, linewidth = 0.8, position = pd) +
    geom_point(size = 3.5, position = pd) +
    scale_colour_manual(name = NULL, values = c(Observed = "black", Predicted = "#C0392B")) +
    scale_shape_manual(name = NULL, values = c(Observed = 16, Predicted = 17)) +
    coord_cartesian(ylim = c(0, y_max)) +
    labs(title = title, x = "Predator age group", y = "Proportion of hake in diet") +
    theme_bw(base_size = 13) +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5),
          panel.grid = element_blank())
}


plot_diet_fit <- function(fit, file = NULL, width = 8, height = 5) {
  p <- (make_diet_plot(diet_fit_table(fit, 2, atf_groups), "Arrowtooth flounder") +
          make_diet_plot(diet_fit_table(fit, 3, sbf_groups), "Sablefish")) +
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  if (!is.null(file)) ggsave(paste0(file, ".png"), p, width = width, height = height, dpi = 300)
  print(p)
  invisible(p)
}

plot_diet_fit(ms_final)

plot_suitability <- function(fit, prey = 1, year = NULL, max_age = 20, normalise = FALSE,
                             file = NULL, width = 11, height = 5) {
  d <- suitability_table(fit, prey, year, max_age)
  if (normalise) {  # scale to max within each predator x sex
    mx <- ave(d$Suitability, d$Predator, d$Sex, FUN = max)
    d$Suitability <- ifelse(mx > 0, d$Suitability / mx, 0)
  }
  prey_name <- fit$data_list$spnames[prey]

  # One plot per predator, sexes stacked; single-sex predators fill the full height
  one <- function(dp) {
    ggplot(droplevels(dp), aes(Prey_Age, Predator_Age, fill = Suitability)) +
      geom_tile() +
      facet_wrap(~ Sex, ncol = 1) +
      scale_x_continuous(breaks = seq(0, max_age, 2), expand = c(0, 0)) +
      scale_y_continuous(breaks = seq(0, max_age, 2), expand = c(0, 0)) +
      scale_fill_viridis_c(option = "plasma", limits = c(0, max(d$Suitability)),
                           name = if (normalise) "Relative\nsuitability" else "Suitability") +
      labs(title = dp$Predator[1], x = paste(prey_name, "age (prey)"), y = "Predator age") +
      theme_bw(base_size = 12) +
      theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))
  }
  p <- wrap_plots(lapply(split(d, d$Predator, drop = TRUE), one), nrow = 1) +
    plot_layout(guides = "collect", axis_titles = "collect")

  if (!is.null(file)) ggsave(paste0(file, ".png"), p, width = width, height = height, dpi = 300)
  print(p)
  invisible(d)
}

plot_suitability(ms_final) # file = "figures/suitability")

# =============================================================================
# plot M2 by predator
# =============================================================================

plot_M2_by_predator <- function(Rceattle,
                                species = NULL,
                                ages = 1:5,
                                file = NULL,
                                minyr = NULL,
                                incl_proj = FALSE,
                                width = 8,
                                height = 6,
                                ...) {

  # --- 1. Data Extraction & Setup ---
  if(class(Rceattle) != "Rceattle"){
    stop("This function supports a single Rceattle object.")
  }

  data_list <- Rceattle$data_list
  spnames <- data_list$spnames
  nspp <- data_list$nspp

  if(is.null(minyr)){ minyr <- data_list$styr }
  endyr <- if(incl_proj) data_list$projyr else data_list$endyr
  years <- minyr:endyr

  if(is.null(species)){ species <- 1:nspp }
  target_spnames <- spnames[species]

  # --- 2. Helper Function for Tidying M2_prop ---
  tidy_m2_prop <- function(m2_prop_array, m1_dimnames, all_spnames, target_spnames, target_ages, target_years) {

    species_names <- m1_dimnames[[1]]
    sex_names     <- m1_dimnames[[2]]
    age_labels    <- m1_dimnames[[3]]
    year_labels   <- m1_dimnames[[4]]

    sp_sex_labels <- paste(
      rep(species_names, times = length(sex_names)),
      rep(sex_names, each = length(species_names)),
      sep = "_"
    )

    dn_m2_prop <- list(
      Pred_SpSex = sp_sex_labels,
      Prey_SpSex = sp_sex_labels,
      Pred_Age   = age_labels,
      Prey_Age   = age_labels,
      Year       = year_labels
    )

    tryCatch({
      dimnames(m2_prop_array) <- dn_m2_prop
    }, error = function(e) {
      print("Error applying dimnames. Array dimensions and names do not match.")
      stop(e)
    })

    df <- as.data.frame.table(m2_prop_array, responseName = "Mortality", stringsAsFactors = FALSE)
    if(nrow(df) == 0) return(data.frame())

    spname_pattern <- paste(all_spnames, collapse = "|")

    df_processed <- df %>%
      mutate(
        Year         = as.numeric(Year),
        Prey_Age     = as.numeric(gsub("\\D", "", Prey_Age)),
        Pred_Species = stringr::str_extract(Pred_SpSex, spname_pattern),
        Prey_Species = stringr::str_extract(Prey_SpSex, spname_pattern)
      ) %>%
      filter(
        Prey_Species %in% target_spnames,
        Prey_Age     %in% target_ages,
        Year         %in% target_years,
        grepl("Sex combined or females", Prey_SpSex)
      )

    if(nrow(df_processed) == 0) return(data.frame())

    df_agg <- df_processed %>%
      group_by(Prey_Species, Prey_Age, Year, Pred_Species) %>%
      summarise(Mortality = sum(Mortality, na.rm = TRUE), .groups = 'drop') %>%
      filter(Mortality > 1e-10) %>%
      rename(Species = Prey_Species, Age = Prey_Age) %>%
      mutate(Type = paste0("M2 (from ", Pred_Species, ")")) %>%
      select(Species, Age, Year, Type, Mortality)

    return(df_agg)
  }

  # --- 3. Extract M2 data ---
  m1_dims <- dimnames(Rceattle$quantities$M1_at_age)

  if(is.null(Rceattle$quantities$M2_prop)) {
    stop("`Rceattle$quantities$M2_prop` not found.")
  }

  combined_mortality_df <- tidy_m2_prop(
    m2_prop_array  = Rceattle$quantities$M2_prop,
    m1_dimnames    = m1_dims,
    all_spnames    = spnames,
    target_spnames = target_spnames,
    target_ages    = ages,
    target_years   = years
  )

  if(nrow(combined_mortality_df) == 0) {
    stop("Data frame is empty after filtering. Check species, ages, and years.")
  }

  # --- 4. Extract SD from sdrep for M2_at_age ---
  combined_mortality_df <- combined_mortality_df %>%
    mutate(Lower_95 = NA_real_, Upper_95 = NA_real_)

  if ("sdrep" %in% names(Rceattle) && !is.null(Rceattle$sdrep)) {

    m2_indices <- which(names(Rceattle$sdrep$value) == "M2_at_age")

    if (length(m2_indices) > 0) {
      m2_sds       <- Rceattle$sdrep$sd[m2_indices]
      m2_array     <- Rceattle$quantities$M2_at_age
      m2_sd_array  <- array(m2_sds, dim = dim(m2_array), dimnames = dimnames(m2_array))

      tidy_total_m2_sd <- function(mort_array, sd_array, target_spnames, target_ages, target_years) {
        df_val <- as.data.frame.table(mort_array, responseName = "M2_total",    stringsAsFactors = FALSE)
        df_sd  <- as.data.frame.table(sd_array,   responseName = "M2_total_sd", stringsAsFactors = FALSE)
        colnames(df_val) <- c("Species_Name_Raw", "Sex", "Age_Str", "Year", "M2_total")
        colnames(df_sd)  <- c("Species_Name_Raw", "Sex", "Age_Str", "Year", "M2_total_sd")

        left_join(df_val, df_sd, by = c("Species_Name_Raw", "Sex", "Age_Str", "Year")) %>%
          mutate(
            Year    = as.numeric(Year),
            Age     = as.numeric(gsub("\\D", "", Age_Str)),
            Species = case_when(
              grepl("Hake", Species_Name_Raw) ~ "Hake",
              grepl("ATF",  Species_Name_Raw) ~ "ATF",
              TRUE ~ "Unknown"
            )
          ) %>%
          filter(
            grepl("Sex combined or females", Sex),
            Species %in% target_spnames,
            Age     %in% target_ages,
            Year    %in% target_years
          ) %>%
          select(Species, Age, Year, M2_total, M2_total_sd)
      }

      total_m2_sd_df <- tidy_total_m2_sd(m2_array, m2_sd_array, target_spnames, ages, years)

      combined_mortality_df <- combined_mortality_df %>%
        left_join(total_m2_sd_df, by = c("Species", "Age", "Year")) %>%
        mutate(
          prop_contribution = ifelse(M2_total > 1e-10, Mortality / M2_total, 0),
          M2_sd_approx      = prop_contribution * M2_total_sd,
          Lower_95          = pmax(0, Mortality - 1.96 * M2_sd_approx),
          Upper_95          = Mortality + 1.96 * M2_sd_approx
        ) %>%
        select(-M2_total, -M2_total_sd, -prop_contribution, -M2_sd_approx)

      message("SD ribbon extracted from sdrep and apportioned per predator.")
    } else {
      message("'M2_at_age' not found in sdrep$value. Plotting without confidence intervals.")
    }
  } else {
    message("No sdrep found in model object. Plotting without confidence intervals.")
  }

  has_ci <- !all(is.na(combined_mortality_df$Lower_95))

  # --- 5. Factor levels, colors, and labels ---
  type_levels <- sort(unique(combined_mortality_df$Type))
  combined_mortality_df <- combined_mortality_df %>%
    mutate(
      Type      = factor(Type, levels = type_levels),
      Age_label = factor(paste("Age", Age),
                         levels = paste("Age", sort(unique(Age))))
    )

  n_types    <- length(type_levels)
  # In section 5, replace:
  all_colors <- RColorBrewer::brewer.pal(max(3, n_types), "Dark2")
  names(all_colors) <- type_levels

  # With:
  all_colors <- c(
    "M2 (from ATF)"      = "#009E73",  # teal green
    "M2 (from Hake)"     = "#21295C",  # deep navy/purple
    "M2 (from Sablefish)"= "#D55E00"   # warm orange
  )

  # --- 6. Plot ---
  p <- ggplot(combined_mortality_df, aes(x = Year, color = Type, fill = Type)) +

    {if(has_ci)
      geom_ribbon(aes(ymin = Lower_95, ymax = Upper_95), alpha = 0.2, color = NA)
    } +

    geom_line(aes(y = Mortality), linewidth = 1.1) +

    # facet on Age_label only — removes the redundant "Hake" strip
    facet_wrap(
      ~ Age_label,
      scales   = "free_y",
      labeller = label_value
    ) +

    labs(
      # title = "M2 (Predation Mortality) by Predator and Prey Age",  # commented out
      x    = "Year",
      y    = "M2 Rate",
      color = "Predator",
      fill  = "Predator"
    ) +
    expand_limits(y = 0) +
    scale_color_manual(values = all_colors) +
    scale_fill_manual(values  = all_colors) +
    theme_bw(base_size = 14) +
    theme(
      legend.position  = "bottom",         # legend at the bottom
      legend.title     = element_text(face = "bold"),
      strip.text       = element_text(face = "bold"),
      axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid = element_blank(),
      panel.spacing    = unit(0.5, "lines")
    )

  # --- 7. Mean summary tables ---
  m2_mean_by_age <- combined_mortality_df %>%
    group_by(Species, Age, Type) %>%
    summarise(Mean_Mortality = mean(Mortality, na.rm = TRUE), .groups = 'drop')

  m2_mean_overall <- combined_mortality_df %>%
    group_by(Species, Type) %>%
    summarise(Mean_Mortality = mean(Mortality, na.rm = TRUE), .groups = 'drop')

  # --- 8. Save & return ---
  if (!is.null(file)) {
    ggsave(
      filename = paste0(file, "_M2_by_predator.png"),
      plot     = p,
      width    = width,
      height   = height,
      units    = "in",
      dpi      = 300
    )
  }

  print(p)

  invisible(list(
    plot            = p,
    m2_timeseries   = combined_mortality_df,
    m2_mean_by_age  = m2_mean_by_age,
    m2_mean_overall = m2_mean_overall
  ))
}


# =============================================================================
# plot_mortality_components
# =============================================================================

plot_mortality_components <- function(Rceattle,
                                      species = NULL,
                                      ages = 1:5,
                                      file = NULL,
                                      minyr = NULL,
                                      incl_proj = FALSE,
                                      width = 10,
                                      height = 5,
                                      ...) {

  # --- 1. Data Extraction & Setup ---
  if(class(Rceattle) != "Rceattle"){
    stop("This function supports a single Rceattle object.")
  }

  data_list <- Rceattle$data_list
  spnames   <- data_list$spnames
  nspp      <- data_list$nspp

  if(is.null(minyr)){ minyr <- data_list$styr }
  endyr <- if(incl_proj) data_list$projyr else data_list$endyr
  years <- minyr:endyr

  if(is.null(species)){ species <- 1:nspp }

  # --- 2. Tidy helper ---
  tidy_mortality <- function(mort_array, mort_name) {
    df <- as.data.frame.table(mort_array, responseName = "Mortality", stringsAsFactors = FALSE)
    colnames(df) <- c("Species_Name_Raw", "Sex", "Age_Str", "Year", "Mortality")

    df %>%
      mutate(
        Year    = as.numeric(Year),
        Age     = as.numeric(gsub("\\D", "", Age_Str)),
        Species = case_when(
          grepl("Hake", Species_Name_Raw) ~ "Hake",
          grepl("ATF",  Species_Name_Raw) ~ "ATF",
          TRUE ~ "Unknown"
        ),
        Type = mort_name
      ) %>%
      filter(grepl("Sex combined or females", Sex)) %>%
      filter(Species %in% spnames[species], Age %in% ages, Year %in% years) %>%
      select(Species, Age, Year, Type, Mortality)
  }

  m1_df      <- tidy_mortality(Rceattle$quantities$M1_at_age, "M1 (Residual)")
  m2_df      <- tidy_mortality(Rceattle$quantities$M2_at_age, "M2 (Predation)")
  total_m_df <- tidy_mortality(Rceattle$quantities$M_at_age,  "Total M")

  combined_mortality_df <- bind_rows(m1_df, m2_df, total_m_df) %>%
    mutate(Type = factor(Type, levels = c("M1 (Residual)", "M2 (Predation)", "Total M")))

  if(nrow(combined_mortality_df) == 0) {
    stop("Data frame is empty after filtering.")
  }

  # --- 3. Extract SD from sdrep ---
  combined_mortality_df <- combined_mortality_df %>%
    mutate(Lower_95 = NA_real_, Upper_95 = NA_real_)

  if ("sdrep" %in% names(Rceattle) && !is.null(Rceattle$sdrep)) {

    # Helper to attach SD for a given reported quantity name
    attach_sd <- function(quantity_name, type_label, target_array) {
      idx <- which(names(Rceattle$sdrep$value) == quantity_name)
      if (length(idx) == 0) return(NULL)

      sds      <- Rceattle$sdrep$sd[idx]
      sd_array <- array(sds, dim = dim(target_array), dimnames = dimnames(target_array))

      df_val <- as.data.frame.table(target_array, responseName = "Mortality",   stringsAsFactors = FALSE)
      df_sd  <- as.data.frame.table(sd_array,     responseName = "Mortality_sd", stringsAsFactors = FALSE)
      colnames(df_val) <- c("Species_Name_Raw", "Sex", "Age_Str", "Year", "Mortality")
      colnames(df_sd)  <- c("Species_Name_Raw", "Sex", "Age_Str", "Year", "Mortality_sd")

      left_join(df_val, df_sd, by = c("Species_Name_Raw", "Sex", "Age_Str", "Year")) %>%
        mutate(
          Year    = as.numeric(Year),
          Age     = as.numeric(gsub("\\D", "", Age_Str)),
          Species = case_when(
            grepl("Hake", Species_Name_Raw) ~ "Hake",
            grepl("ATF",  Species_Name_Raw) ~ "ATF",
            TRUE ~ "Unknown"
          ),
          Type     = type_label,
          Lower_95 = pmax(0, Mortality - 1.96 * Mortality_sd),
          Upper_95 = Mortality + 1.96 * Mortality_sd
        ) %>%
        filter(
          grepl("Sex combined or females", Sex),
          Species %in% spnames[species],
          Age     %in% ages,
          Year    %in% years
        ) %>%
        select(Species, Age, Year, Type, Lower_95, Upper_95)
    }

    sd_m1    <- attach_sd("M1_at_age", "M1 (Residual)", Rceattle$quantities$M1_at_age)
    sd_m2    <- attach_sd("M2_at_age", "M2 (Predation)", Rceattle$quantities$M2_at_age)
    sd_total <- attach_sd("M_at_age",  "Total M",        Rceattle$quantities$M_at_age)

    sd_df <- bind_rows(sd_m1, sd_m2, sd_total)

    if (!is.null(sd_df) && nrow(sd_df) > 0) {
      combined_mortality_df <- combined_mortality_df %>%
        select(-Lower_95, -Upper_95) %>%
        left_join(sd_df, by = c("Species", "Age", "Year", "Type"))
      message("SD ribbons attached for M1, M2, and Total M.")
    } else {
      message("Mortality quantities not found in sdrep. Plotting without confidence intervals.")
    }

  } else {
    message("No sdrep found. Plotting without confidence intervals.")
  }

  has_ci <- !all(is.na(combined_mortality_df$Lower_95))

  # --- 4. Labels ---
  combined_mortality_df <- combined_mortality_df %>%
    mutate(
      Type      = factor(Type, levels = c("M1 (Residual)", "M2 (Predation)", "Total M")),
      Age_label = paste("Age", Age)
    )

  # --- 5. Plot ---
  p <- ggplot(combined_mortality_df, aes(x = Year, color = Type, fill = Type)) +

    {if(has_ci)
      geom_ribbon(aes(ymin = Lower_95, ymax = Upper_95), alpha = 0.2, color = NA)
    } +

    geom_line(aes(y = Mortality), linewidth = 1.1) +

    # facet on Age_label only — removes the redundant species strip
    facet_grid(
      ~ Age_label,
      scales   = "free_y",
      labeller = label_value
    ) +

    labs(
      # title = "Mortality Components by Species and Age",  # commented out
      x     = "Year",
      y     = "Mortality Rate",
      color = "Mortality Component",
      fill  = "Mortality Component"
    ) +
    expand_limits(y = 0) +
    scale_color_manual(values = c(
      "M1 (Residual)"  = "cornflowerblue",
      "M2 (Predation)" = "yellow4",
      "Total M"        = "black"
    )) +
    scale_fill_manual(values = c(
      "M1 (Residual)"  = "cornflowerblue",
      "M2 (Predation)" = "yellow4",
      "Total M"        = "black"
    )) +
    theme_bw(base_size = 14) +
    theme(
      legend.position  = "bottom",         # legend at the bottom
      legend.title     = element_text(face = "bold"),
      strip.text       = element_text(face = "bold"),
      axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid = element_blank(),
      panel.spacing    = unit(0.5, "lines")
    )

  # --- 6. Save & return ---
  if (!is.null(file)) {
    ggsave(
      filename = paste0(file, "_mortality_components.png"),
      plot     = p,
      width    = width,
      height   = height,
      units    = "in",
      dpi      = 300
    )
  }

  print(p)
  invisible(p)
}

m2_results <- plot_M2_by_predator(
  Rceattle = ms_final,
  species = 1, # Hake
  ages = 1:20) #,
#file = "figures/Manuscript1/M2_by_pred_test2")

# This will now produce the correct plot
plot_mortality_components(
  Rceattle = ms_final,
  species = 1, # Hake
  ages = 1:9)
#file = "figures/Manuscript2/M_components")


## Hake population estimates: single-species configurations vs final multispecies model
## Shaded bands: +/- 2 SD from sdrep
library(ggplot2)
library(patchwork)

# Hake (species 1) time series + sdrep SD over the hindcast years
series_sd <- function(fit, quantity, label, sp = 1) {
  dl  <- fit$data_list
  yrs <- dl$styr:dl$endyr
  val <- fit$quantities[[quantity]][sp, seq_along(yrs)]
  sdv <- fit$sdrep$sd[names(fit$sdrep$value) == quantity]          # [species x year], column-major
  sd  <- if (length(sdv) %% dl$nspp == 0 && length(sdv) >= dl$nspp * length(yrs))
    matrix(sdv, nrow = dl$nspp)[sp, seq_along(yrs)] else NA
  data.frame(Year = yrs, Value = as.numeric(val), SD = as.numeric(sd), Model = label)
}

ts_panel <- function(models, quantity, ylab, scale, cols) {
  d <- do.call(rbind, Map(series_sd, models, quantity, names(models)))
  d$Model <- factor(d$Model, levels = names(models))
  d[, c("Value", "SD")] <- d[, c("Value", "SD")] / scale
  d$lo <- pmax(0, d$Value - 2 * d$SD)
  d$hi <- d$Value + 2 * d$SD

  ggplot(d, aes(Year, Value, colour = Model, fill = Model)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.12, colour = NA) +
    geom_line(linewidth = 0.9) +
    scale_colour_manual(values = cols) +
    scale_fill_manual(values = cols) +
    expand_limits(y = 0) +
    labs(x = NULL, y = ylab) +
    theme_bw(base_size = 13) +
    theme(legend.title = element_blank(), panel.grid = element_blank())
}

plot_population <- function(models, consumption = TRUE, file = NULL, width = 10, height = 8,
                            cols = c("#56B4E9", "#E69F00", "#009E73", "#440154")) {
  cols <- setNames(cols[seq_along(models)], names(models))
  p_ssb <- ts_panel(models, "ssb",     "Spawning biomass (million t)", 1e6, cols)
  p_bio <- ts_panel(models, "biomass", "Age-1+ biomass (million t)",   1e6, cols)
  p_rec <- ts_panel(models, "R",       "Recruitment (millions)",       1e6, cols)  # R in thousands

  if (consumption) {
    # Hake eaten by all predators (last model = multispecies)
    ms  <- models[[length(models)]]
    dl  <- ms$data_list
    yrs <- dl$styr:dl$endyr
    con <- data.frame(Year = yrs, Value = apply(
      ms$quantities$B_eaten_as_prey[1, , , seq_along(yrs), drop = FALSE], 4, sum, na.rm = TRUE) / 1e6)
    p_con <- ggplot(con, aes(Year, Value)) +
      geom_line(linewidth = 0.9, colour = cols[[length(cols)]]) +
      expand_limits(y = 0) +
      labs(x = NULL, y = "Hake consumed (million t)") +
      theme_bw(base_size = 13) +
      theme(panel.grid = element_blank())
    p <- (p_ssb + p_bio) / (p_rec + p_con)
  } else {
    p <- p_ssb + p_bio + p_rec
  }
  p <- p + plot_layout(guides = "collect") & theme(legend.position = "bottom")

  if (!is.null(file)) ggsave(paste0(file, ".png"), p, width = width, height = height, dpi = 300)
  print(p)
  invisible(p)
}

## Mean habitat temperature by species (env_data columns named by species)
plot_temperature <- function(fit, file = NULL, width = 8, height = 7,
                             cols = c(Hake = "#21295C", ATF = "#009E73", Sablefish = "#D55E00")) {
  env <- fit$data_list$env_data
  spp <- intersect(names(cols), names(env))
  d <- do.call(rbind, lapply(spp, function(s) data.frame(Year = env$Year, Temperature = env[[s]], Species = s)))
  d$Species <- factor(d$Species, levels = spp)

  p <- ggplot(d, aes(Year, Temperature, colour = Species)) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.5) +
    facet_wrap(~ Species, ncol = 1, scales = "free_y") +
    scale_colour_manual(values = cols) +
    labs(x = NULL, y = "Mean habitat temperature (°C)") +
    theme_bw(base_size = 14) +
    theme(legend.position = "none", strip.text = element_text(face = "bold"),
          panel.grid = element_blank())

  if (!is.null(file)) ggsave(paste0(file, ".png"), p, width = width, height = height, dpi = 300)
  print(p)
  invisible(p)
}

# Run
plot_population(list("SS: M fixed" = ss_fix, "SS: M prior" = ss_Mprior,
                     "SS: M free"  = ss_Mfree, "Multispecies" = ms_final),
                file = "figures/population_estimates")
plot_temperature(ms_final, file = "figures/temperature_timeseries")

