
# INPUT: df(s) from the different trial analyses' outcomes.
# OUTPUT: plot(s) tracking each outcome and effect size.


# DESCRIPTION:

# Functions to do some basic plots for verbose trail analyses' outputs.


# First, a function to save plots with
# (default) width and height, and a resolution of 1200 dpi.

save_plot <- function(plot, filename, width = 12, height = 8) {
  png(filename = filename, width = width, height = height, units = "in", res = 1200)
  print(plot)
  dev.off()
}


# Get general functions for protective effect and time-to-event

plot_protective_effect <- function(protective_effect, key_intervention_time,
                                   sim_length, trial_title) {
  year <- 365
  
  require(ggplot2)
  require(rcartocolor)
  
  # So that facets come somehow ordered...
  
  protective_effect <- protective_effect %>%
    mutate(
      effect_order = case_when(
        grepl("Prevalence", measure) ~ 1,
        grepl("Incidence", measure) ~ 2,
        grepl("Hazard Ratio", measure) ~ 3,
        TRUE ~ 4
      ),
      type_measure = forcats::fct_reorder(type_measure, effect_order, .fun = min),
      measure = forcats::fct_reorder(measure, effect_order, .fun = min)
    )
  
  # Actual plot - general for effect sizes
  
  ggplot(data = protective_effect,
         aes(x = timestep, y = mean)) +
    geom_ribbon(data = filter(protective_effect, grepl("Ins", type_measure)),
      aes(ymin = lower_95quant, ymax = upper_95quant), alpha = 0.3) +
    geom_errorbar(data = filter(protective_effect, !grepl("Ins", type_measure)),
              aes(ymin = lower_95quant, ymax = upper_95quant), width = 0.2) +
    geom_point() + geom_line() +
    geom_vline(xintercept = key_intervention_time*year, color = "firebrick", linetype = "dashed") +
    scale_x_continuous(breaks = seq(0, sim_length * year, by = year),
                       labels = (0:sim_length)) +
    scale_y_continuous(labels = scales::percent, limits = 0:1) +
    labs(x = "Year", y = "Intervention Protective Effect",
         title = trial_title) +
    theme_bw() +
    theme(legend.position = "bottom", legend.title = element_blank()) +
    facet_nested(type_measure + measure ~ ., scales = "free")
}

plot_time_to_event_pair <- function(tte_df, sim_length, trial_title, x_label) {
  
  require(dplyr)
  require(ggplot2)
  require(ggpubr)
  require(rcartocolor)
  require(survival)
  require(tidycmprsk)

  year <- 365

  surv_fit_infection <- survfit(Surv(time_to_infection, ever_infected) ~ run, data = tte_df) %>%
    tidy() %>%
    mutate(strata = gsub("run=", "", strata))

  surv_fit_case <- survfit(Surv(time_to_case, ever_case) ~ run, data = tte_df) %>%
    tidy() %>%
    mutate(strata = gsub("run=", "", strata))

  p_infection <- ggplot(data = surv_fit_infection,
                        aes(x = time, y = estimate, color = strata, fill = strata)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
    geom_point(data = surv_fit_infection[surv_fit_infection$n.censor != 0,],
               aes(x = time, y = estimate, color = strata, fill = strata),
               shape = 4, size = 4) +
    scale_color_manual(values = carto_pal(name = "Safe")[c(11, 10)]) +
    scale_fill_manual(values = carto_pal(name = "Safe")[c(11, 10)]) +
    scale_x_continuous(breaks = seq(0, sim_length * year, by = year),
                       labels = 0:sim_length) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = x_label, y = "Proportion without infection") +
    theme_bw() +
    theme(legend.position = "bottom", legend.title = element_blank())

  p_case <- ggplot(data = surv_fit_case,
                   aes(x = time, y = estimate, color = strata, fill = strata)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
    geom_point(data = surv_fit_case[surv_fit_case$n.censor != 0,],
               aes(x = time, y = estimate, color = strata, fill = strata),
               shape = 4, size = 4) +
    scale_color_manual(values = carto_pal(name = "Safe")[c(11, 10)]) +
    scale_fill_manual(values = carto_pal(name = "Safe")[c(11, 10)]) +
    scale_x_continuous(breaks = seq(0, sim_length * year, by = year),
                       labels = 0:sim_length) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = x_label, y = "Proportion without clinical case") +
    theme_bw() +
    theme(legend.position = "bottom", legend.title = element_blank())

  annotate_figure(ggarrange(p_infection, p_case, nrow = 2), top = trial_title)
}


# CLUSTER VARIABILITY
# Spaghetti of sim-averaged prevalence per cluster (thin) + overall mean (bold), by arm.

plot_cluster_variability <- function(estimates_df, key_intervention_time, sim_length, trial_title,
                                     measure_filter = "Infection Prevalence",
                                     type_filter = "True Instantaneous") {
  year <- 365
  require(ggplot2)
  require(dplyr)
  require(rcartocolor)

  df <- estimates_df %>%
    filter(as.character(measure) == measure_filter,
           type_measure == type_filter,
           !is.na(mean))

  df_mean <- df %>%
    group_by(run, timestep) %>%
    summarise(mean = mean(mean, na.rm = TRUE), .groups = "drop")

  ggplot(df, aes(x = timestep, y = mean,
                 group = interaction(cluster_id, run), color = run)) +
    geom_line(alpha = 0.3, linewidth = 0.4) +
    geom_line(data = df_mean, aes(group = run), linewidth = 1.2) +
    geom_vline(xintercept = key_intervention_time * year, color = "firebrick", linetype = "dashed") +
    scale_color_manual(values = carto_pal(name = "Safe")[c(11, 10)]) +
    scale_x_continuous(breaks = seq(0, sim_length * year, by = year), labels = 0:sim_length) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "Year", y = measure_filter,
         title = paste0(trial_title, ": Cluster variability (", measure_filter, ")")) +
    theme_bw() +
    theme(legend.position = "bottom", legend.title = element_blank()) +
    facet_wrap(~ run)
}


# COHORT COMPARISON
# All analysis populations on one plot; colour = population, linetype = arm.

plot_cohort_comparison <- function(estimates_summary, key_intervention_time, sim_length, trial_title,
                                   measures_to_show = c("Infection Prevalence", "Case Incidence p.p.year"),
                                   type_to_show = c("True Instantaneous", "True Aggregate Incidence")) {
  year <- 365
  require(ggplot2)
  require(dplyr)
  require(rcartocolor)

  n_pops <- length(unique(estimates_summary$analysis_population))
  pop_colors <- carto_pal(name = "Safe")[seq_len(n_pops)]

  df <- estimates_summary %>%
    filter(as.character(measure) %in% measures_to_show,
           type_measure %in% type_to_show,
           !is.na(mean))

  ggplot(df, aes(x = timestep, y = mean,
                 color = analysis_population, fill = analysis_population,
                 linetype = run,
                 group = interaction(analysis_population, run))) +
    geom_ribbon(aes(ymin = lower_95quant, ymax = upper_95quant), alpha = 0.15, color = NA) +
    geom_line() +
    geom_vline(xintercept = key_intervention_time * year, color = "firebrick", linetype = "dashed") +
    scale_color_manual(values = pop_colors) +
    scale_fill_manual(values = pop_colors) +
    scale_x_continuous(breaks = seq(0, sim_length * year, by = year), labels = 0:sim_length) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "Year", y = NULL,
         title = paste0(trial_title, ": Estimates by age cohort"),
         color = "Population", fill = "Population", linetype = "Arm") +
    theme_bw() +
    theme(legend.position = "bottom") +
    facet_wrap(~ measure, scales = "free_y")
}


# EFFECT BY COHORT
# Effect sizes for all populations overlaid; colour = population.

plot_effect_by_cohort <- function(effects_df, key_intervention_time, sim_length, trial_title) {
  year <- 365
  require(ggplot2)
  require(ggh4x)
  require(dplyr)
  require(rcartocolor)

  n_pops <- length(unique(effects_df$analysis_population))
  pop_colors <- carto_pal(name = "Safe")[seq_len(n_pops)]

  effects_df <- effects_df %>%
    filter(!is.na(mean)) %>%
    mutate(
      effect_order = case_when(
        grepl("Prevalence", measure) ~ 1,
        grepl("Incidence", measure) ~ 2,
        grepl("Hazard Ratio", measure) ~ 3,
        TRUE ~ 4
      ),
      type_measure = forcats::fct_reorder(type_measure, effect_order, .fun = min),
      measure = forcats::fct_reorder(measure, effect_order, .fun = min)
    )

  ggplot(effects_df, aes(x = timestep, y = mean,
                          color = analysis_population, fill = analysis_population,
                          group = analysis_population)) +
    geom_ribbon(data = filter(effects_df, grepl("Ins", type_measure)),
                aes(ymin = lower_95quant, ymax = upper_95quant), alpha = 0.2, color = NA) +
    geom_errorbar(data = filter(effects_df, !grepl("Ins", type_measure)),
                  aes(ymin = lower_95quant, ymax = upper_95quant),
                  width = 0.2, position = position_dodge(width = 50)) +
    geom_point(position = position_dodge(width = 50)) +
    geom_line() +
    geom_vline(xintercept = key_intervention_time * year, color = "firebrick", linetype = "dashed") +
    scale_color_manual(values = pop_colors) +
    scale_fill_manual(values = pop_colors) +
    scale_x_continuous(breaks = seq(0, sim_length * year, by = year), labels = 0:sim_length) +
    scale_y_continuous(labels = scales::percent, limits = 0:1) +
    labs(x = "Year", y = "Protective Effect",
         title = paste0(trial_title, ": Effect size by age cohort"),
         color = "Population", fill = "Population") +
    theme_bw() +
    theme(legend.position = "bottom") +
    facet_nested(type_measure + measure ~ ., scales = "free")
}


# A wrapper to save all the basic plots for a two-arm trial analysis,
# given the outputs of the trial analyses and the survey/ACD protocol used.

save_two_arm_trial_plots <- function(trial_results, trial_slug,
                                     key_intervention_time, sim_length,
                                     trial_title) {
  # General options
  
  require(dplyr)
  require(ggplot2)
  require(ggh4x)
  require(rcartocolor)

  source("functions/trial_tidy_outputs.R")
  make_output_dirs()

  year <- 365
  
  
  analysis_populations <- unique(trial_results$estimates_summary$analysis_population)

  for (ap in analysis_populations) {

    ## Prev/Inc estimates plots

    # One with the true estimates

    plot_true_estimates <- ggplot(
      data = filter(
        trial_results$estimates_summary,
        analysis_population == ap,
        (type_measure == "True Instantaneous" & grepl("Prev", measure)) |
          (grepl("aggregate", type_measure) & grepl("p.p.y", measure) & !grepl("ACD", type_measure))),
      aes(x = timestep, y = mean, group = run, color = run)) +
      geom_ribbon(aes(ymin = lower_95quant, ymax = upper_95quant, fill = run),
                  alpha = 0.3, color = NA) +
      geom_point() + geom_line() +
      geom_vline(xintercept = key_intervention_time * year, color = "firebrick", linetype = "dashed") +
      scale_color_manual(values = carto_pal(name = "Safe")[c(11, 10)]) +
      scale_fill_manual(values = carto_pal(name = "Safe")[c(11, 10)]) +
      scale_x_continuous(breaks = seq(0, sim_length * year, by = year),
                         labels = 0:sim_length) +
      labs(x = "Year", y = NULL, title = paste0(trial_title, " [", ap, "]: True estimates")) +
      theme_bw() +
      theme(legend.position = "bottom", legend.title = element_blank()) +
      facet_nested(type_measure + measure ~ ., scales = "free", drop = TRUE)

    save_plot(plot_true_estimates,
              paste0("outputs/plots/prevalence_incidence/", trial_slug, "_true_estimates_", ap, ".png"))

    # One for incidence estimates only

    plot_incidence_estimates <- ggplot(
      data = filter(
        trial_results$estimates_summary,
        analysis_population == ap, !is.na(mean),
        grepl("p.p.y", measure) & !grepl("Ins", type_measure)),
      aes(x = timestep, y = mean, group = run, color = run)) +
      geom_ribbon(aes(ymin = lower_95quant, ymax = upper_95quant, fill = run),
                  alpha = 0.3, color = NA) +
      geom_point() + geom_line() +
      geom_vline(xintercept = key_intervention_time * year, color = "firebrick", linetype = "dashed") +
      scale_color_manual(values = carto_pal(name = "Safe")[c(11, 10)]) +
      scale_fill_manual(values = carto_pal(name = "Safe")[c(11, 10)]) +
      scale_x_continuous(breaks = seq(0, sim_length * year, by = year),
                         labels = 0:sim_length) +
      labs(x = "Year", y = NULL, title = paste0(trial_title, " [", ap, "]: Incidence estimates")) +
      theme_bw() +
      theme(legend.position = "bottom", legend.title = element_blank()) +
      facet_nested(type_measure + measure ~ ., scales = "free", drop = TRUE)

    save_plot(plot_incidence_estimates,
              paste0("outputs/plots/prevalence_incidence/", trial_slug, "_incidence_estimates_", ap, ".png"))

    # One for prevalence estimates only

    plot_prevalence_estimates <- ggplot(
      data = filter(
        trial_results$estimates_summary,
        analysis_population == ap, !is.na(mean),
        grepl("Prevalence", measure)),
      aes(x = timestep, y = mean, group = run, color = run)) +
      geom_ribbon(aes(ymin = lower_95quant, ymax = upper_95quant, fill = run),
                  alpha = 0.3, color = NA) +
      geom_point(aes(shape = type_measure, size = type_measure)) +
      geom_line() +
      geom_vline(xintercept = key_intervention_time * year, color = "firebrick", linetype = "dashed") +
      scale_shape_manual(breaks = c("True Instantaneous", "Cross-sectional surveys"),
                         values = c(16:17)) +
      scale_size_manual(breaks = c("True Instantaneous", "Cross-sectional surveys"),
                        values = c(1, 4)) +
      scale_color_manual(values = carto_pal(name = "Safe")[c(11, 10)]) +
      scale_fill_manual(values = carto_pal(name = "Safe")[c(11, 10)]) +
      scale_x_continuous(breaks = seq(0, sim_length * year, by = year),
                         labels = 0:sim_length) +
      labs(x = "Year", y = NULL, title = paste0(trial_title, " [", ap, "]: Prevalence estimates")) +
      theme_bw() +
      theme(legend.position = "bottom", legend.title = element_blank()) +
      facet_grid(measure ~ ., scales = "free")

    save_plot(plot_prevalence_estimates,
              paste0("outputs/plots/prevalence_incidence/", trial_slug, "_prevalence_estimates_", ap, ".png"),
              height = 5)


    ## Intervention protective effect plots

    # One for the prev/inc-based relative protective effect

    plot_relative_effect <- plot_protective_effect(
      protective_effect = trial_results$relative_effect %>%
        filter(
          analysis_population == ap,
          !is.na(mean),
          !(type_measure == "True Instantaneous" & grepl("Incidence", measure)),
          grepl("Infection Prev", measure) | grepl("Case Incidence p.p.y", measure)
        ),
      key_intervention_time = key_intervention_time,
      sim_length = sim_length,
      trial_title = paste0(trial_title, " [", ap, "]: Relative protective effect")
    )

    save_plot(plot_relative_effect,
              paste0("outputs/plots/effect_size/", trial_slug, "_relative_effect_", ap, ".png"))

    # One with HRs

    plot_all_effects <- plot_protective_effect(
      protective_effect = trial_results$all_effects %>%
        filter(
          analysis_population == ap,
          !is.na(mean),
          !(type_measure == "True Instantaneous" & grepl("Incidence", measure)),
          grepl("Infection Prev", measure) |
            grepl("Case Incidence p.p.y", measure) |
            (grepl("ime-to", type_measure) & grepl("Case", measure))
        ),
      key_intervention_time = key_intervention_time,
      sim_length = sim_length,
      trial_title = paste0(trial_title, " [", ap, "]: Protective effect with hazard ratios")
    )

    save_plot(plot_all_effects,
              paste0("outputs/plots/effect_size/", trial_slug, "_all_effects_with_hr_", ap, ".png"),
              height = 10)

  }
  
  
  ## Kaplan Maier curves (one file per analysis population)

  tte_configs <- list(
    list(data = "tte_true_1",  label = "True time-to-event",        x = "Year after trial start",        suffix = "_true_1_intervention"),
    list(data = "tte_true_2",  label = "True time-to-event",        x = "Year after second intervention", suffix = "_true_2_intervention"),
    list(data = "tte_acd_1",   label = "Time-to-event w/ ACD visits", x = "Year after trial start",       suffix = "_acd_1_intervention"),
    list(data = "tte_acd_2",   label = "Time-to-event w/ ACD visits", x = "Year after second intervention", suffix = "_acd_2_intervention")
  )

  for (cfg in tte_configs) {
    tte_df <- trial_results[[cfg$data]]
    for (ap in unique(tte_df$analysis_population)) {
      save_plot(
        plot_time_to_event_pair(
          tte_df = filter(tte_df, analysis_population == ap),
          sim_length = sim_length,
          trial_title = paste0(trial_title, " [", ap, "]: ", cfg$label),
          x_label = cfg$x
        ),
        paste0("outputs/plots/time_to_event/", trial_slug, cfg$suffix, "_", ap, ".png"),
        width = 8, height = 8
      )
    }
  }


  ## Cross-population plots

  # Cluster variability per population
  for (ap in analysis_populations) {
    save_plot(
      plot_cluster_variability(
        estimates_df = trial_results$estimates_true %>% filter(analysis_population == ap),
        key_intervention_time = key_intervention_time,
        sim_length = sim_length,
        trial_title = paste0(trial_title, " [", ap, "]")
      ),
      paste0("outputs/plots/prevalence_incidence/", trial_slug, "_cluster_variability_", ap, ".png"),
      height = 5
    )
  }

  # All populations on one plot
  save_plot(
    plot_cohort_comparison(
      estimates_summary = trial_results$estimates_summary,
      key_intervention_time = key_intervention_time,
      sim_length = sim_length,
      trial_title = trial_title
    ),
    paste0("outputs/plots/prevalence_incidence/", trial_slug, "_cohort_comparison.png"),
    height = 6
  )

  # Effect sizes across populations
  save_plot(
    plot_effect_by_cohort(
      effects_df = trial_results$all_effects %>%
        filter(
          !is.na(mean),
          !(type_measure == "True Instantaneous" & grepl("Incidence", measure)),
          grepl("Infection Prev", measure) |
            grepl("Case Incidence p.p.y", measure) |
            (grepl("ime-to", type_measure) & grepl("Case", measure))
        ),
      key_intervention_time = key_intervention_time,
      sim_length = sim_length,
      trial_title = trial_title
    ),
    paste0("outputs/plots/effect_size/", trial_slug, "_effect_by_cohort.png"),
    height = 10
  )
}
