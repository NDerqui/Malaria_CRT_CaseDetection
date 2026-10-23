
# INPUT: filtered verbose output w/ infections signaled (for two arms)
# OUTPUT: list with main outcomes and effect sizes computed


# DESCRIPTION:

# Wrapper function to run all functions that estimate all outcomes
# (prevalence, incidence, time-to-event), true or with PCD/ACD protocols,
# and protective effect size of intervention.

analyse_two_arm_trial <- function(trial_id,
                                  # Trial and sim characteristics
                                  trial_start, trial_second_intervention,
                                  sim_length,
                                  # Protocols for PCD and ACD
                                  survey_protocol, acd_protocol,
                                  # Hazard ratio window (in years; mirrors aggregate incidence periods)
                                  window_years_hr = 1) {
  
  year <- 365
  require(rlang)
  require(dplyr)
  
  # Fetch all necessary functions
  
  source("functions/trial_outcomes.R")
  source("functions/trial_outcomes_event_time.R")
  source("functions/trial_effects.R")
  
  ## DATA
  
  # Read the data just with the trial name
  
  file_path <- paste0("outputs/cohort_data/", trial_id, ".csv")
  
  if (file.exists(file_path)) {
    
    infections <- read.csv(file_path)
  
  } else {
    
    infections <- readRDS(paste0("outputs/cohort_data/", trial_id, ".rds"))
  }
  
  infections_control <- infections %>% filter(run == "Control")
  infections_intervention <- infections %>% filter(run == "Intervention")
  
  rm(infections)
  
  ## OUTCOMES (prev and inc)
  
  # Get all true estimates:
  # Instantaneous (prevalence and incidence) and cumulative (only incidence)
  
  estimates_true <- bind_rows(
    infections_control %>% estimate_true_realtime_outcomes() %>% mutate(run = "Control"),
    infections_intervention %>% estimate_true_realtime_outcomes() %>% mutate(run = "Intervention"),
    infections_control %>% estimate_true_aggregate_incidence() %>% mutate(run = "Control"),
    infections_intervention %>% estimate_true_aggregate_incidence() %>% mutate(run = "Intervention")
  )
  
  # Get estimates as if measured with PCD:
  # Use a cross-sectional survey to estimate prevalence (protocol can be modified)
  
  estimates_survey <- bind_rows(
    inject(infections_control %>% estimate_survey_prevalence(trial_start = trial_start, !!!survey_protocol)) %>% mutate(run = "Control"),
    inject(infections_intervention %>% estimate_survey_prevalence(trial_start = trial_start, !!!survey_protocol)) %>% mutate(run = "Intervention")
  )
  
  # Get estimates as if measured with ACD:
  # Use routine visits to estimate incidence (protocol can be modified)
  
  estimates_acd <- bind_rows(
    inject(infections_control %>% estimate_acd_incidence(trial_start = trial_start, !!!acd_protocol)) %>% mutate(run = "Control"),
    inject(infections_intervention %>% estimate_acd_incidence(trial_start = trial_start, !!!acd_protocol)) %>% mutate(run = "Intervention")
  )
  
  # Put all estimates together
  
  estimates_all <- bind_rows(estimates_true, estimates_survey, estimates_acd)
  
  ## PROTECTIVE EFFECT (prev and inc)
  
  # Prepare our estimates to get protective effect sizes
  
  relative_effect <- estimates_all %>%
    tidy_outcomes_for_effect() %>%
    estimate_relative_effect() %>%
    filter(!is.na(mean))
  
  ## SUMARIES (prev and inc across simulations)
  
  # Now that we have calculated effect for each simulation using each sim's prev and inc,
  # we can now get summaries (mean and 95%CI) across simulations for each measure and timestep
  
  possible_grouping_vars <- c(
    "run", "analysis_population", "type_measure", "cluster_id",
    "sim", "timestep", "period", "period_label")
  measures <- c("prevalence_infection", "prevalence_case",
                "incidence_ppd_infection", "incidence_ppd_case",
                "incidence_ppy_infection", "incidence_ppy_case")
  measures_labels <- c("Infection Prevalence", "Case Prevalence",
                       "Infection Incidence p.p.day", "Case Incidence p.p.day",
                       "Infection Incidence p.p.year", "Case Incidence p.p.year")
  
  estimates_all <- estimates_all %>%
    select(-any_of(c("n", "person_days_at_risk", "infections", "cases",
                     "new_infections", "new_cases"))) %>%
    # Pivot to get one row per timestep/sim and per measure
    pivot_longer(-any_of(possible_grouping_vars),
                 names_to = "measure", values_to = "value") %>%
    # Refactor measure to have a more readable label for plotting
    mutate(measure = factor(measure, levels = measures, labels = measures_labels)) %>%
    # Create mean and 95%CI across simulations
    group_by(across(c(all_of(possible_grouping_vars[!(possible_grouping_vars == "sim")]), "measure"))) %>%
    mutate(mean = mean(value, na.rm = TRUE),
           lower_95quant = quantile(value, probs = 0.025, na.rm = TRUE),
           upper_95quant = quantile(value, probs = 0.975, na.rm = TRUE)) %>%
    filter(row_number() == 1) %>% ungroup() %>%
    select(all_of(possible_grouping_vars[!(possible_grouping_vars == "sim")]), measure, mean, lower_95quant, upper_95quant)
  
  estimates_true <- estimates_all %>% filter(grepl("True", type_measure)) %>% select_if(function(x){!all(is.na(x))})
  estimates_acd <- estimates_all %>% filter(grepl("ACD", type_measure)) %>% select_if(function(x){!all(is.na(x))})
  estimates_survey <- estimates_all %>% filter(grepl("Cross", type_measure)) %>% select_if(function(x){!all(is.na(x))})

  # Cluster-averaged summary for plotting: one row per (run, population, type, timestep, measure)
  estimates_summary <- estimates_all %>%
    group_by(across(any_of(c("run", "analysis_population", "type_measure",
                              "timestep", "period", "period_label", "measure")))) %>%
    summarise(mean = mean(mean, na.rm = TRUE),
              lower_95quant = mean(lower_95quant, na.rm = TRUE),
              upper_95quant = mean(upper_95quant, na.rm = TRUE),
              .groups = "drop")

  rm(estimates_all); gc()

  ## OUTCOMES (time to event)

  has_second_intervention <- !is.null(trial_second_intervention) && !is.na(trial_second_intervention)

  # Get true time-to-event estimates for first infection/case from the first intervention

  tte_true_1 <- bind_rows(
    infections_control %>% estimate_true_time_to_event(time_inter = trial_start*year) %>% mutate(run = "Control"),
    infections_intervention %>% estimate_true_time_to_event(time_inter = trial_start*year) %>% mutate(run = "Intervention")
  ) %>%
    prepare_time_to_event_survival(time_inter = trial_start*year, sim_length = sim_length*year)

  # Get true time-to-event estimates for first infection/case from the second intervention (if any)

  tte_true_2 <- if (has_second_intervention) {
    t2 <- (trial_start + trial_second_intervention) * year
    bind_rows(
      infections_control %>% estimate_true_time_to_event(time_inter = t2) %>% mutate(run = "Control"),
      infections_intervention %>% estimate_true_time_to_event(time_inter = t2) %>% mutate(run = "Intervention")
    ) %>%
      prepare_time_to_event_survival(time_inter = t2, sim_length = sim_length*year)
  } else NULL

  # Get time-to-event estimates for first infection/case from the first intervention as if doing ACD

  tte_acd_1 <- bind_rows(
    inject(infections_control %>% estimate_acd_time_to_event(time_inter = trial_start*year, !!!acd_protocol)) %>% mutate(run = "Control"),
    inject(infections_intervention %>% estimate_acd_time_to_event(time_inter = trial_start*year, !!!acd_protocol)) %>% mutate(run = "Intervention")
  ) %>%
    prepare_time_to_event_survival(time_inter = trial_start*year, sim_length = sim_length*year)

  # Get time-to-event estimates for first infection/case from the second intervention as if doing ACD (if any)

  tte_acd_2 <- if (has_second_intervention) {
    t2 <- (trial_start + trial_second_intervention) * year
    bind_rows(
      inject(infections_control %>% estimate_acd_time_to_event(time_inter = t2, !!!acd_protocol)) %>% mutate(run = "Control"),
      inject(infections_intervention %>% estimate_acd_time_to_event(time_inter = t2, !!!acd_protocol)) %>% mutate(run = "Intervention")
    ) %>%
      prepare_time_to_event_survival(time_inter = t2, sim_length = sim_length*year)
  } else NULL

  rm(infections_control, infections_intervention); gc()

  # Windowed (landmark) hazard ratios: one per year window from trial start.
  # Uses tte_true_1 / tte_acd_1 (full follow-up from first intervention).
  # window_years_hr controls window size, default 1 year (mirrors aggregate incidence).

  hazard_ratio <- bind_rows(
    estimate_hazard_ratio_by_sim_windowed(tte_true_1, "time_to_case",      "ever_case",
      time_inter = trial_start*year, sim_length = sim_length*year, window_years = window_years_hr) %>%
      mutate(type_measure = "True Time-to-event",          measure = "Case Hazard Ratio"),
    estimate_hazard_ratio_by_sim_windowed(tte_true_1, "time_to_infection", "ever_infected",
      time_inter = trial_start*year, sim_length = sim_length*year, window_years = window_years_hr) %>%
      mutate(type_measure = "True Time-to-event",          measure = "Infection Hazard Ratio"),
    estimate_hazard_ratio_by_sim_windowed(tte_acd_1,  "time_to_case",      "ever_case",
      time_inter = trial_start*year, sim_length = sim_length*year, window_years = window_years_hr) %>%
      mutate(type_measure = "Time-to-event w/ ACD visits", measure = "Case Hazard Ratio"),
    estimate_hazard_ratio_by_sim_windowed(tte_acd_1,  "time_to_infection", "ever_infected",
      time_inter = trial_start*year, sim_length = sim_length*year, window_years = window_years_hr) %>%
      mutate(type_measure = "Time-to-event w/ ACD visits", measure = "Infection Hazard Ratio")
  ) %>%
    mutate(effect = 1 - hazard_ratio) %>%
    group_by(across(c("analysis_population", "type_measure", "measure", "timestep"))) %>%
    mutate(mean = mean(effect, na.rm = TRUE),
           lower_95quant = quantile(effect, probs = 0.025, na.rm = TRUE),
           upper_95quant = quantile(effect, probs = 0.975, na.rm = TRUE),
           power_any = mean(significant_any, na.rm = TRUE),
           power_benefit = mean(significant_benefit, na.rm = TRUE),
           n_sim = n_distinct(sim)) %>%
    filter(row_number() == 1) %>% ungroup() %>%
    select(analysis_population, type_measure, measure, timestep,
           mean, lower_95quant, upper_95quant, power_any, power_benefit, n_sim)
  
  ## SAVE
  
  # Get a list with all the estimates and effect sizes for further uses
  
  list(
    estimates_true = estimates_true,
    estimates_survey = estimates_survey,
    estimates_acd = estimates_acd,
    estimates_summary = estimates_summary,
    tte_true_1 = tte_true_1,
    tte_true_2 = tte_true_2,
    tte_acd_1 = tte_acd_1,
    tte_acd_2 = tte_acd_2,
    relative_effect = relative_effect,
    hazard_ratio = hazard_ratio,
    all_effects = bind_rows(relative_effect, hazard_ratio)
  )
}


# DESCRIPTION:

# A function to save all the estimates and effect sizes to corresponding folder/file

save_two_arm_trial <- function(trial_results, trial_slug) {
  
  
  # Create folder structure if not present already
  source("functions/trial_tidy_outputs.R")
  make_output_dirs()
  
  # Save all estimates and effect sizes to corresponding folder/file
  write.csv(trial_results$estimates_true, file = paste0("outputs/estimates/prevalence_incidence/", trial_slug, "_true.csv"), row.names = FALSE)
  write.csv(trial_results$estimates_survey, file = paste0("outputs/estimates/prevalence_incidence/", trial_slug, "_survey.csv"), row.names = FALSE)
  write.csv(trial_results$estimates_acd, file = paste0("outputs/estimates/prevalence_incidence/", trial_slug, "_acd.csv"), row.names = FALSE)
  write.csv(trial_results$estimates_summary, file = paste0("outputs/estimates/prevalence_incidence/", trial_slug, "_summary.csv"), row.names = FALSE)
  write.csv(trial_results$tte_true_1, file = paste0("outputs/estimates/time_to_event/", trial_slug, "_true_1_intervention.csv"), row.names = FALSE)
  if (!is.null(trial_results$tte_true_2))
    write.csv(trial_results$tte_true_2, file = paste0("outputs/estimates/time_to_event/", trial_slug, "_true_2_intervention.csv"), row.names = FALSE)
  write.csv(trial_results$tte_acd_1, file = paste0("outputs/estimates/time_to_event/", trial_slug, "_acd_1_intervention.csv"), row.names = FALSE)
  if (!is.null(trial_results$tte_acd_2))
    write.csv(trial_results$tte_acd_2, file = paste0("outputs/estimates/time_to_event/", trial_slug, "_acd_2_intervention.csv"), row.names = FALSE)
  write.csv(trial_results$relative_effect, file = paste0("outputs/estimates/effect_size/", trial_slug, "_relative_effect.csv"), row.names = FALSE)
  write.csv(trial_results$hazard_ratio, file = paste0("outputs/estimates/effect_size/", trial_slug, "_hazard_ratio.csv"), row.names = FALSE)
  write.csv(trial_results$all_effects, file = paste0("outputs/estimates/effect_size/", trial_slug, "_all_effects.csv"), row.names = FALSE)
}