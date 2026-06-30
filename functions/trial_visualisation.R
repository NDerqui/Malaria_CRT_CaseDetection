# functions/trial_visualisation.R

save_plot <- function(plot, filename, width = 12, height = 8) {
  png(filename = filename, width = width, height = height, units = "in", res = 1200)
  print(plot)
  dev.off()
}

plot_protective_effect <- function(protective_effect, key_intervention_time,
                                   sim_length, trial_title) {
  year <- 365
  
  ggplot(data = protective_effect,
         aes(x = timestep, y = effect)) +
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

save_two_arm_trial_plots <- function(trial_results, trial_slug,
                                     key_intervention_time, sim_length,
                                     trial_title) {
  require(dplyr)
  require(ggplot2)
  require(ggh4x)
  require(rcartocolor)

  source("functions/trial_tidy_outputs.R")
  source("functions/trial_effects.R")
  make_output_dirs()

  year <- 365

  plot_all_estimates <- trial_results$estimates_all %>%
    tidy_outcomes_for_effect() %>%
    filter(!is.na(value))

  plot_true_estimates <- ggplot(
    data = filter(
      plot_all_estimates,
      (type_measure == "True Instantaneous" & grepl("Prev", measure)) |
        (grepl("aggregate", type_measure) & grepl("p.p.y", measure) & !grepl("ACD", type_measure))
    ),
    aes(x = timestep, y = value, group = run, color = run)
  ) +
    geom_point() + geom_line() +
    geom_vline(xintercept = key_intervention_time * year, color = "firebrick", linetype = "dashed") +
    scale_color_manual(values = carto_pal(name = "Safe")[c(11, 10)]) +
    scale_x_continuous(breaks = seq(0, sim_length * year, by = year),
                       labels = 0:sim_length) +
    labs(x = "Year", y = NULL, title = paste0(trial_title, ": True estimates")) +
    theme_bw() +
    theme(legend.position = "bottom", legend.title = element_blank()) +
    facet_nested(type_measure + measure ~ ., scales = "free", drop = TRUE)

  save_plot(plot_true_estimates,
            paste0("outputs/plots/prevalence_incidence/", trial_slug, "_true_estimates.png"))

  plot_incidence_estimates <- ggplot(
    data = filter(plot_all_estimates, grepl("p.p.y", measure) & !grepl("Ins", type_measure)),
    aes(x = timestep, y = value, group = run, color = run)
  ) +
    geom_point() + geom_line() +
    geom_vline(xintercept = key_intervention_time * year, color = "firebrick", linetype = "dashed") +
    scale_color_manual(values = carto_pal(name = "Safe")[c(11, 10)]) +
    scale_x_continuous(breaks = seq(0, sim_length * year, by = year),
                       labels = 0:sim_length) +
    labs(x = "Year", y = NULL, title = paste0(trial_title, ": Incidence estimates")) +
    theme_bw() +
    theme(legend.position = "bottom", legend.title = element_blank()) +
    facet_nested(type_measure + measure ~ ., scales = "free", drop = TRUE)

  save_plot(plot_incidence_estimates,
            paste0("outputs/plots/prevalence_incidence/", trial_slug, "_incidence_estimates.png"))

  plot_prevalence_estimates <- ggplot(
    data = filter(plot_all_estimates, grepl("Prevalence", measure)),
    aes(x = timestep, y = value, group = run, color = run)
  ) +
    geom_point(aes(shape = type_measure, size = type_measure)) +
    geom_line() +
    geom_vline(xintercept = key_intervention_time * year, color = "firebrick", linetype = "dashed") +
    scale_color_manual(values = carto_pal(name = "Safe")[c(11, 10)]) +
    scale_x_continuous(breaks = seq(0, sim_length * year, by = year),
                       labels = 0:sim_length) +
    labs(x = "Year", y = NULL, title = paste0(trial_title, ": Prevalence estimates")) +
    theme_bw() +
    theme(legend.position = "bottom", legend.title = element_blank()) +
    facet_grid(measure ~ ., scales = "free")

  save_plot(plot_prevalence_estimates,
            paste0("outputs/plots/prevalence_incidence/", trial_slug, "_prevalence_estimates.png"),
            height = 5)

  plot_relative_effect <- plot_protective_effect(
    protective_effect = trial_results$relative_effect %>%
      filter(
        !is.na(effect),
        !(type_measure == "True Instantaneous" & grepl("Incidence", measure)),
        grepl("Infection Prev", measure) | grepl("Case Incidence p.p.y", measure)
      ),
    key_intervention_time = key_intervention_time,
    sim_length = sim_length,
    trial_title = paste0(trial_title, ": Relative protective effect")
  )

  save_plot(plot_relative_effect,
            paste0("outputs/plots/effect_size/", trial_slug, "_relative_effect.png"))

  plot_all_effects <- plot_protective_effect(
    protective_effect = trial_results$all_effects %>%
      filter(
        !is.na(effect),
        !(type_measure == "True Instantaneous" & grepl("Incidence", measure)),
        grepl("Infection Prev", measure) |
          grepl("Case Incidence p.p.y", measure) |
          (grepl("ime-to", type_measure) & grepl("Case", measure))
      ),
    key_intervention_time = key_intervention_time,
    sim_length = sim_length,
    trial_title = paste0(trial_title, ": Protective effect with hazard ratios")
  )

  save_plot(plot_all_effects,
            paste0("outputs/plots/effect_size/", trial_slug, "_all_effects_with_hr.png"),
            height = 10)

  save_plot(
    plot_time_to_event_pair(
      trial_results$tte_true_1,
      sim_length = sim_length,
      trial_title = paste0(trial_title, ": True time-to-event"),
      x_label = "Year after trial start"
    ),
    paste0("outputs/plots/time_to_event/", trial_slug, "_true_1_intervention.png"),
    width = 8, height = 8
  )

  save_plot(
    plot_time_to_event_pair(
      trial_results$tte_true_2,
      sim_length = sim_length,
      trial_title = paste0(trial_title, ": True time-to-event"),
      x_label = "Year after second intervention"
    ),
    paste0("outputs/plots/time_to_event/", trial_slug, "_true_2_intervention.png"),
    width = 8, height = 8
  )

  save_plot(
    plot_time_to_event_pair(
      trial_results$tte_acd_1,
      sim_length = sim_length,
      trial_title = paste0(trial_title, ": Time-to-event w/ ACD visits"),
      x_label = "Year after trial start"
    ),
    paste0("outputs/plots/time_to_event/", trial_slug, "_acd_1_intervention.png"),
    width = 8, height = 8
  )

  save_plot(
    plot_time_to_event_pair(
      trial_results$tte_acd_2,
      sim_length = sim_length,
      trial_title = paste0(trial_title, ": Time-to-event w/ ACD visits"),
      x_label = "Year after second intervention"
    ),
    paste0("outputs/plots/time_to_event/", trial_slug, "_acd_2_intervention.png"),
    width = 8, height = 8
  )
}
