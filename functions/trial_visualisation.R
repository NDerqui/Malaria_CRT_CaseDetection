
# INPUT: df(s) from the different trial analyses' outcomes.
# OUTPUT: plot(s) tracking each outcome and effect size.


# DESCRIPTION:

# Functions to do some basic plots for verbose trail analyses' outputs.



# GENERAL -----------------------------------------------------------------


# First, a function to save plots with
# (default) width and height, and a resolution of 1200 dpi.

save_plot <- function(plot, filename, width = 12, height = 8) {
  
  png(filename = filename, width = width, height = height, units = "in", res = 1200)
  print(plot)
  dev.off()

  }

# Shared label-tidying helper: clean population names, order facets consistently.
# Apply to any df before plotting to avoid repeating this block in every function.

tidy_plot_labels <- function(df) {
  
  df %>%
    mutate(
      # Get rid of ugly age labels
      analysis_population = dplyr::case_when(
        grepl("^all", as.character(analysis_population)) ~ "All ages",
        TRUE ~ gsub("_", "-", gsub("age_", "Ages ", as.character(analysis_population)))
      ),
      # Re-order analysis cohorts so they come in specific orders
      analysis_population = forcats::fct_reorder(
        analysis_population,
        dplyr::case_when(
          grepl("All", analysis_population) ~ 1L,
          grepl("0-5", analysis_population) ~ 2L,
          grepl("5-15", analysis_population) ~ 3L,
          grepl("15-", analysis_population) ~ 4L,
          TRUE ~ 5L
        ),
        .fun = min
      ),
      # Re-order measures so they come in specific order if several in same plot
      effect_order = dplyr::case_when(
        (grepl("Prevalence", measure) & grepl("True", type_measure))  ~ 1L,
        (grepl("Prevalence", measure) & grepl("survey", type_measure)) ~ 2L,
        (grepl("Incidence",  measure) & grepl("True", type_measure))  ~ 3L,
        (grepl("Incidence",  measure) & grepl("ACD",  type_measure))  ~ 4L,
        (grepl("Hazard",     measure) & grepl("True", type_measure))  ~ 5L,
        (grepl("Hazard",     measure) & grepl("ACD",  type_measure))  ~ 6L,
        TRUE ~ 7L
      ),
      type_measure = forcats::fct_reorder(as.character(type_measure), effect_order, .fun = min),
      measure      = forcats::fct_reorder(as.character(measure),      effect_order, .fun = min)
    ) %>%
    select(-effect_order)
}



# DISEASE LEVELS ------------------------------------------------------


## Following functions are designed to compare disease levels (prev and incidence)
## across age populations and clusters.


# Across clustes (just as check - only with true instantaneous measures)

plot_cluster_variability <- function(estimates_df, key_intervention_time, sim_length, trial_title,
                                     measures_to_show = c("Infection Prevalence", "Case Incidence p.p.day"),
                                     type_to_show     = "True Instantaneous") {
  
  # Requirements
  
  year <- 365
  require(ggplot2)
  require(dplyr)
  require(ggh4x)
  require(rcartocolor)
  
  # Data preparations
  
  df <- estimates_df %>%
    filter(!is.na(mean),
           as.character(measure)      %in% measures_to_show,
           as.character(type_measure) %in% type_to_show) %>%
    tidy_plot_labels()
  
  df_mean <- df %>%
    group_by(run, timestep, measure, type_measure, analysis_population) %>%
    summarise(mean = mean(mean, na.rm = TRUE), .groups = "drop")
  
  # Plot
  
  ggplot(df, aes(x = timestep, y = mean,
                 group = interaction(cluster_id, run), color = run)) +
    geom_line(alpha = 0.3, linewidth = 0.4) +
    geom_line(data = df_mean, aes(group = run), linewidth = 1.1) +
    geom_vline(xintercept = key_intervention_time * year, color = "firebrick", linetype = "dashed") +
    scale_color_manual(values = carto_pal(name = "Safe")[c(11, 10)]) +
    scale_x_continuous(breaks = seq(0, sim_length * year, by = year), labels = 0:sim_length) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "Year", y = "Estimate",
         title = paste0(trial_title, ": Cluster variability")) +
    theme_bw() +
    theme(legend.position = "bottom", legend.title = element_blank()) +
    ggh4x::facet_grid2(analysis_population ~ measure, scales = "free_y", independent = "y")
}


# Shared internals for cohort-comparison plots (prev, true inc, ACD inc).

.cohort_comparison_base <- function(df, key_intervention_time, sim_length,
                                    title, errorbar = FALSE) {
  
  # Requirements
  
  year      <- 365
  
  # Getting an idea on number of age cohorts and assigning colors
  
  n_pops    <- length(unique(df$analysis_population))
  pop_colors <- rcartocolor::carto_pal(name = "Safe")[seq_len(n_pops)]
  
  # Plot
  
  p <- ggplot(df, aes(x = timestep, y = mean,
                      color = analysis_population, fill = analysis_population,
                      linetype = run,
                      group = interaction(analysis_population, run))) +
    geom_vline(xintercept = key_intervention_time * year, color = "firebrick", linetype = "dashed") +
    scale_color_manual(values = pop_colors) +
    scale_fill_manual(values  = pop_colors) +
    scale_x_continuous(breaks = seq(0, sim_length * year, by = year), labels = 0:sim_length) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "Year", y = "Estimate", title = title,
         color = "Population", fill = "Population", linetype = "Arm") +
    theme_bw() +
    theme(legend.position = "bottom") +
    ggh4x::facet_grid2(type_measure + measure ~ ., scales = "free_y", independent = "y")
  
  # Error bar or ribbon depending on whether continuous or point estimate
  
  if (errorbar) {
    p <- p +
      geom_errorbar(aes(ymin = lower_95quant, ymax = upper_95quant),
                    width = 0.2, position = position_dodge(width = 50)) +
      geom_point(position = position_dodge(width = 50)) +
      geom_line()
  } else {
    p <- p +
      geom_ribbon(aes(ymin = lower_95quant, ymax = upper_95quant), alpha = 0.15, color = NA) +
      geom_line()
  }
  p
}

# Bul of functions to get different true or observed prevalence / incidence

plot_cohort_comparison_true_prev <- function(estimates_summary, key_intervention_time, sim_length, trial_title,
                                             measures_to_show = c("Infection Prevalence", "Case Prevalence")) {
  require(ggplot2); require(dplyr); require(ggh4x); require(rcartocolor)
  
  df <- estimates_summary %>%
    filter(!is.na(mean),
           as.character(measure) %in% measures_to_show,
           type_measure == "True Instantaneous") %>%
    tidy_plot_labels()
  
  .cohort_comparison_base(df, key_intervention_time, sim_length,
                          title = paste0(trial_title, ": True instantaneous estimates by age cohort"),
                          errorbar = FALSE)
}

plot_cohort_comparison_survey_prev <- function(estimates_summary, key_intervention_time, sim_length, trial_title,
                                               measures_to_show = c("Infection Prevalence", "Case Prevalence")) {
  require(ggplot2); require(dplyr); require(ggh4x); require(rcartocolor)
  
  df <- estimates_summary %>%
    filter(!is.na(mean),
           as.character(measure) %in% measures_to_show,
           grepl("survey", type_measure)) %>%
    tidy_plot_labels()
  
  .cohort_comparison_base(df, key_intervention_time, sim_length,
                          title = paste0(trial_title, ": Estimates observed by cross-sectional surveys"),
                          errorbar = TRUE)
}

plot_cohort_comparison_true_inc <- function(estimates_summary, key_intervention_time, sim_length, trial_title) {
  require(ggplot2); require(dplyr); require(ggh4x); require(rcartocolor)
  
  df <- estimates_summary %>%
    filter(!is.na(mean),
           grepl("p.p.y", as.character(measure)),
           grepl("True.*[Aa]g", type_measure)) %>%
    tidy_plot_labels()
  
  .cohort_comparison_base(df, key_intervention_time, sim_length,
                          title = paste0(trial_title, ": True aggregate incidence by age cohort"),
                          errorbar = TRUE)
}


plot_cohort_comparison_acd_inc <- function(estimates_summary, key_intervention_time, sim_length, trial_title) {
  require(ggplot2); require(dplyr); require(ggh4x); require(rcartocolor)
  
  df <- estimates_summary %>%
    filter(!is.na(mean),
           grepl("p.p.y", as.character(measure)),
           grepl("ACD", type_measure)) %>%
    tidy_plot_labels()
  
  .cohort_comparison_base(df, key_intervention_time, sim_length,
                          title = paste0(trial_title, ": Incidence observed by ACD by age cohort"),
                          errorbar = TRUE)
}



# PE / TIME-TO-EVENT ------------------------------------------------------


## These plots' section are related to the estimated effect size or HR,
## not disease levels (incidence and prevalence).

# Get general functions for protective effect

plot_protective_effect <- function(protective_effect, key_intervention_time,
                                   sim_length, trial_title) {
  
  # Requirements
  
  year <- 365
  
  require(ggplot2)
  require(ggh4x)
  require(rcartocolor)
  
  # Cleaning

  protective_effect <- tidy_plot_labels(protective_effect)
  
  # Plot (simple, overall effect)

  ggplot(data = protective_effect,
         aes(x = timestep, y = mean)) +
    geom_ribbon(data = filter(protective_effect, grepl("Ins", type_measure)),
      aes(ymin = lower_95quant, ymax = upper_95quant), alpha = 0.3) +
    geom_errorbar(data = filter(protective_effect, !grepl("Ins", type_measure)),
              aes(ymin = lower_95quant, ymax = upper_95quant), width = 0.2) +
    geom_point() + geom_line() +
    geom_vline(xintercept = key_intervention_time * year, color = "firebrick", linetype = "dashed") +
    scale_x_continuous(breaks = seq(0, sim_length * year, by = year), labels = 0:sim_length) +
    scale_y_continuous(labels = scales::percent, limits = 0:1) +
    labs(x = "Year", y = "Effect size",
         title = paste0(trial_title, ": Intervention Protective Effect")) +
    theme_bw() +
    theme(legend.position = "bottom", legend.title = element_blank()) +
    ggh4x::facet_grid2(type_measure + measure ~ ., scales = "free_y", independent = "y")
}

# EFFECT BY COHORT
# Effect sizes for all populations overlaid; colour = population.


plot_effect_by_cohort <- function(effects_df, key_intervention_time, sim_length, trial_title) {
  
  # Requirements
  
  year <- 365
  
  require(ggplot2)
  require(ggh4x)
  require(dplyr)
  require(rcartocolor)
  
  # Cleaning

  df <- effects_df %>%
    filter(!is.na(mean)) %>%
    tidy_plot_labels()

  # Get the number of populations and assign colours to each population

  n_pops     <- length(unique(df$analysis_population))
  pop_colors <- carto_pal(name = "Safe")[seq_len(n_pops)]

  # Plot
  
  ggplot(df, aes(x = timestep, y = mean,
                 color = analysis_population, fill = analysis_population,
                 group = analysis_population)) +
    geom_ribbon(data = filter(df, grepl("Ins", type_measure)),
                aes(ymin = lower_95quant, ymax = upper_95quant), alpha = 0.2, color = NA) +
    geom_errorbar(data = filter(df, !grepl("Ins", type_measure)),
                  aes(ymin = lower_95quant, ymax = upper_95quant),
                  width = 0.2, position = position_dodge(width = 50)) +
    geom_point(position = position_dodge(width = 50)) +
    geom_line() +
    geom_vline(xintercept = key_intervention_time * year, color = "firebrick", linetype = "dashed") +
    scale_color_manual(values = pop_colors) +
    scale_fill_manual(values = pop_colors) +
    scale_x_continuous(breaks = seq(0, sim_length * year, by = year), labels = 0:sim_length) +
    scale_y_continuous(labels = scales::percent, limits = 0:1) +
    labs(x = "Year", y = "Effect size",
         title = paste0(trial_title, ": Intervention Protective Effect by age cohort"),
         color = "Population", fill = "Population") +
    theme_bw() +
    theme(legend.position = "bottom") +
    ggh4x::facet_grid2(type_measure + measure ~ ., scales = "free_y", independent = "y")
}

# TIME-TO-EVENT
# Kaplan Meier curves

plot_time_to_event_pair <- function(tte_df, sim_length, trial_title, x_label) {
  
  # Requirements
  
  require(dplyr)
  require(ggplot2)
  require(ggpubr)
  require(rcartocolor)
  require(survival)
  require(tidycmprsk)

  year <- 365
  
  # Run survival model and plot the K-M curve

  surv_fit_infection <- survfit(Surv(time_to_infection, ever_infected) ~ run, data = tte_df) %>%
    tidy() %>% mutate(strata = gsub("run=", "", strata))

  surv_fit_case <- survfit(Surv(time_to_case, ever_case) ~ run, data = tte_df) %>%
    tidy() %>% mutate(strata = gsub("run=", "", strata))

  make_km_plot <- function(surv_fit, y_label) {
    ggplot(surv_fit, aes(x = time, y = estimate, color = strata, fill = strata)) +
      geom_line(linewidth = 1) +
      geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
      geom_point(data = filter(surv_fit, n.censor != 0),
                 shape = 4, size = 4) +
      scale_color_manual(values = carto_pal(name = "Safe")[c(11, 10)]) +
      scale_fill_manual(values  = carto_pal(name = "Safe")[c(11, 10)]) +
      scale_x_continuous(breaks = seq(0, sim_length * year, by = year), labels = 0:sim_length) +
      scale_y_continuous(labels = scales::percent) +
      labs(x = x_label, y = y_label) +
      theme_bw() +
      theme(legend.position = "bottom", legend.title = element_blank())
  }

  annotate_figure(
    ggarrange(
      make_km_plot(surv_fit_infection, "Proportion without infection"),
      make_km_plot(surv_fit_case,      "Proportion without clinical case"),
      nrow = 2
    ),
    top = trial_title
  )
}



# WRAPPER -----------------------------------------------------------------


save_two_arm_trial_plots <- function(trial_results, trial_slug,
                                     key_intervention_time, sim_length,
                                     trial_title) {
  
  # Requirements
  
  require(dplyr)
  require(ggplot2)
  require(ggh4x)
  require(rcartocolor)

  source("functions/trial_tidy_outputs.R")
  make_output_dirs()

  year <- 365


  ## Cross-population: cluster variability

  save_plot(
    plot_cluster_variability(
      estimates_df        = trial_results$estimates_true,
      key_intervention_time = key_intervention_time,
      sim_length          = sim_length,
      trial_title         = trial_title
    ),
    paste0("outputs/plots/prevalence_incidence/", trial_slug, "_cluster_variability.png"),
    height = 6
  )


  ## Cross-population: cohort comparisons (three separate plots)

  save_plot(
    plot_cohort_comparison_true_prev(
      estimates_summary     = trial_results$estimates_summary,
      key_intervention_time = key_intervention_time,
      sim_length            = sim_length,
      trial_title           = trial_title
    ),
    paste0("outputs/plots/prevalence_incidence/", trial_slug, "_cohort_comparison_true_prev.png"),
    height = 6
  )
  
  save_plot(
    plot_cohort_comparison_survey_prev(
      estimates_summary     = trial_results$estimates_summary,
      key_intervention_time = key_intervention_time,
      sim_length            = sim_length,
      trial_title           = trial_title
    ),
    paste0("outputs/plots/prevalence_incidence/", trial_slug, "_cohort_comparison_survey_prev.png"),
    height = 6
  )

  save_plot(
    plot_cohort_comparison_true_inc(
      estimates_summary     = trial_results$estimates_summary,
      key_intervention_time = key_intervention_time,
      sim_length            = sim_length,
      trial_title           = trial_title
    ),
    paste0("outputs/plots/prevalence_incidence/", trial_slug, "_cohort_comparison_true_inc.png"),
    height = 6
  )

  save_plot(
    plot_cohort_comparison_acd_inc(
      estimates_summary     = trial_results$estimates_summary,
      key_intervention_time = key_intervention_time,
      sim_length            = sim_length,
      trial_title           = trial_title
    ),
    paste0("outputs/plots/prevalence_incidence/", trial_slug, "_cohort_comparison_acd_inc.png"),
    height = 6
  )


  ## Effect sizes across populations

  save_plot(
    plot_effect_by_cohort(
      effects_df = trial_results$all_effects %>%
        filter(
          !is.na(mean),
          !(type_measure == "True Instantaneous" & grepl("Incidence", measure)),
          grepl("Infection Prev", measure) |
            grepl("Case Incidence p.p.y", measure) |
            (grepl("ime-to", type_measure) & grepl("Hazard", measure))
        ),
      key_intervention_time = key_intervention_time,
      sim_length            = sim_length,
      trial_title           = trial_title
    ),
    paste0("outputs/plots/effect_size/", trial_slug, "_effect_by_cohort.png"),
    height = 10
  )


  ## Kaplan-Meier curves (one file per analysis population)

  tte_configs <- list(
    list(key = "tte_true_1", label = "True time-to-event",          x = "Year after trial start",         suffix = "_true_1_intervention"),
    list(key = "tte_true_2", label = "True time-to-event",          x = "Year after second intervention",  suffix = "_true_2_intervention"),
    list(key = "tte_acd_1",  label = "Time-to-event w/ ACD visits", x = "Year after trial start",         suffix = "_acd_1_intervention"),
    list(key = "tte_acd_2",  label = "Time-to-event w/ ACD visits", x = "Year after second intervention",  suffix = "_acd_2_intervention")
  )

  for (cfg in tte_configs) {
    tte_df <- trial_results[[cfg$key]]
    if (is.null(tte_df)) next
    for (ap in unique(tte_df$analysis_population)) {
      save_plot(
        plot_time_to_event_pair(
          tte_df      = filter(tte_df, analysis_population == ap),
          sim_length  = sim_length,
          trial_title = paste0(trial_title, " [", ap, "]: ", cfg$label),
          x_label     = cfg$x
        ),
        paste0("outputs/plots/time_to_event/", trial_slug, cfg$suffix, "_", ap, ".png"),
        width = 8, height = 8
      )
    }
  }
}
