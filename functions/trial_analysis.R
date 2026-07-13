
# INPUT: filtered verbose output w/ infections signaled (for two arms)
# OUTPUT: list with main outcomes and effect sizes computed


# DESCRIPTION:

# Wrapper function to run all functions that estimate all outcomes
# (prevalence, incidence, time-to-event), true or with PCD/ACD protocols,
# and protective effect size of intervention.

analyse_two_arm_trial <- function(trial_slug,
                                  # Trial and sim characteristics
                                  trial_start, trial_second_intervention,
                                  sim_length,
                                  # Protocols for PCD and ACD
                                  survey_protocol, acd_protocol) {
  
  year <- 365
  require(rlang)
  require(dplyr)
  
  # Fetch all necessary functions
  
  source("functions/trial_outcomes.R")
  source("functions/trial_outcomes_event_time.R")
  source("functions/trial_effects.R")
  
  ## DATA
  
  # Read the data just with the trial name
  
  infections_control <- read.csv(paste0("outputs/cohort_data/", trial_slug, "_control.csv"))
  infections_intervention <- read.csv(paste0("outputs/cohort_data/", trial_slug, "_intervention.csv"))
  
  ## OUTCOMES (prev and inc)
  
  # Get all true estimates:
  # Instantaneous (prevalence and incidence) and cumulative (only incidence)
  
  estimates_true <- bind_rows(
    infections_control %>% estimate_true_realtime_outcomes() %>% mutate(run = "Control"),
    infections_intervention %>% estimate_true_realtime_outcomes() %>% mutate(run = "Intervention"),
    infections_control %>% estimate_true_aggregate_incidence(trial_start = trial_start) %>% mutate(run = "Control"),
    infections_intervention %>% estimate_true_aggregate_incidence(trial_start = trial_start) %>% mutate(run = "Intervention")
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
    mutate(mean = ifelse(is.nan(mean), NA, ifelse(is.numeric(mean), mean, NA))) %>%
    filter(!is.na(mean) & is.numeric(mean))
  
  ## SUMARIES (prev and inc across simulations)
  
  # Now that we have calculated effect for each simulation using each sim's prev and inc,
  # we can now get summaries (mean and 95%CI) across simulations for each measure and timestep
  
  possible_grouping_vars <- c("run", "sim", "timestep", "type_measure", "period", "period_label")
  
  estimates_all <- estimates_all %>%
    # Pivot to get one row per timestep/sim and per measure
    pivot_longer(-any_of(possible_grouping_vars),
                 names_to = "measure", values_to = "value") %>%
    # Create mean and 95%CI across simulations
    group_by(across(c(all_of(possible_grouping_vars[!(possible_grouping_vars == "sim")]), "measure"))) %>%
    mutate(mean = mean(value, na.rm = TRUE),
           lower_95quant = quantile(value, probs = 0.025, na.rm = TRUE),
           upper_95quant = quantile(value, probs = 0.975, na.rm = TRUE)) %>%
    filter(row_number() == 1) %>% ungroup() %>%
    select(all_of(possible_grouping_vars[!(possible_grouping_vars == "sim")]), mean, lower_95quant, upper_95quant)
  
  estimates_true <- estimates_all %>% filter(grepl("True", type_measure)) %>% select_if(function(x){!all(is.na(x))})
  estimates_acd <- estimates_all %>% filter(grepl("ACD", type_measure)) %>% select_if(function(x){!all(is.na(x))})
  estimates_survey <- estimates_all %>% filter(grepl("Cross", type_measure)) %>% select_if(function(x){!all(is.na(x))})
  
  ## OUTCOMES (time to event)
  
  # Get true time-to-event estimates for first infection/case from the first intervention
  
  tte_true_1 <- bind_rows(
    infections_control %>% estimate_true_time_to_event(time_inter = trial_start*year) %>% mutate(run = "Control"),
    infections_intervention %>% estimate_true_time_to_event(time_inter = trial_start*year) %>% mutate(run = "Intervention")
  ) %>%
    prepare_time_to_event_survival(time_inter = trial_start*year, sim_length = sim_length*year)
  
  # Get true time-to-event estimates for first infection/case from the second intervention
  
  tte_true_2 <- bind_rows(
    infections_control %>% estimate_true_time_to_event(time_inter = (trial_start + trial_second_intervention)*year) %>% mutate(run = "Control"),
    infections_intervention %>% estimate_true_time_to_event(time_inter = (trial_start + trial_second_intervention)*year) %>% mutate(run = "Intervention")
  ) %>%
    prepare_time_to_event_survival(time_inter = (trial_start + trial_second_intervention)*year, sim_length = sim_length*year)
  
  # Get time-to-event estimates for first infection/case from the first intervention as if doing ACD
  
  tte_acd_1 <- bind_rows(
    inject(infections_control %>% estimate_acd_time_to_event(time_inter = trial_start*year, !!!acd_protocol)) %>% mutate(run = "Control"),
    inject(infections_intervention %>% estimate_acd_time_to_event(time_inter = trial_start*year, !!!acd_protocol)) %>% mutate(run = "Intervention")
  ) %>%
    prepare_time_to_event_survival(time_inter = trial_start*year, sim_length = sim_length*year)
  
  # Get time-to-event estimates for first infection/case from the second intervention as if doing ACD
  
  tte_acd_2 <- bind_rows(
    inject(infections_control %>% estimate_acd_time_to_event(time_inter = (trial_start + trial_second_intervention)*year, !!!acd_protocol)) %>% mutate(run = "Control"),
    inject(infections_intervention %>% estimate_acd_time_to_event(time_inter = (trial_start + trial_second_intervention)*year, !!!acd_protocol)) %>% mutate(run = "Intervention")
  ) %>%
    prepare_time_to_event_survival(time_inter = (trial_start + trial_second_intervention)*year, sim_length = sim_length*year)
  
  # Put all estimates together and calculate hazard ratios
  # Once we have all HR for all sim, estimate summaries across sims
  
  hazard_ratio <- bind_rows(
    estimate_hazard_ratio_by_sim(tte_true_1, "time_to_case", "ever_case") %>% mutate(type_measure = "True Time-to-event", timestep = (trial_start + trial_second_intervention)*year),
    estimate_hazard_ratio_by_sim(tte_true_2, "time_to_case", "ever_case") %>% mutate(type_measure = "True Time-to-event", timestep = sim_length*year),
    estimate_hazard_ratio_by_sim(tte_acd_1, "time_to_case", "ever_case") %>% mutate(type_measure = "Time-to-event w/ ACD visits", timestep = (trial_start + trial_second_intervention)*year),
    estimate_hazard_ratio_by_sim(tte_acd_2, "time_to_case", "ever_case") %>% mutate(type_measure = "Time-to-event w/ ACD visits", timestep = sim_length*year)
  ) %>%
    mutate(measure = "Case Hazard Ratio",
           effect = 1 - hazard_ratio) %>%
    # Create mean and 95%CI across simulations
    group_by(across(c("type_measure", "measure", "timestep"))) %>%
    mutate(mean = mean(effect, na.rm = TRUE),
           lower_95quant = quantile(effect, probs = 0.025, na.rm = TRUE),
           upper_95quant = quantile(effect, probs = 0.975, na.rm = TRUE),
           power_any = mean(significant_any, na.rm = TRUE),
           power_benefit = mean(significant_benefit, na.rm = TRUE),
           n_sim = n_distinct(sim)) %>%
    filter(row_number() == 1) %>% ungroup() %>%
    select(type_measure, measure, timestep, mean, lower_95quant, upper_95quant, power_any, power_benefit, n_sim)
  
  ## SAVE
  
  # Get a list with all the estimates and effect sizes for further uses
  
  list(
    estimates_true = estimates_true,
    estimates_survey = estimates_survey,
    estimates_acd = estimates_acd,
    estimates_all = estimates_all,
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
  write.csv(trial_results$estimates_all, file = paste0("outputs/estimates/prevalence_incidence/", trial_slug, "_all.csv"), row.names = FALSE)
  write.csv(trial_results$tte_true_1, file = paste0("outputs/estimates/time_to_event/", trial_slug, "_true_1_intervention.csv"), row.names = FALSE)
  write.csv(trial_results$tte_true_2, file = paste0("outputs/estimates/time_to_event/", trial_slug, "_true_2_intervention.csv"), row.names = FALSE)
  write.csv(trial_results$tte_acd_1, file = paste0("outputs/estimates/time_to_event/", trial_slug, "_acd_1_intervention.csv"), row.names = FALSE)
  write.csv(trial_results$tte_acd_2, file = paste0("outputs/estimates/time_to_event/", trial_slug, "_acd_2_intervention.csv"), row.names = FALSE)
  write.csv(trial_results$relative_effect, file = paste0("outputs/estimates/effect_size/", trial_slug, "_relative_effect.csv"), row.names = FALSE)
  write.csv(trial_results$hazard_ratio, file = paste0("outputs/estimates/effect_size/", trial_slug, "_hazard_ratio.csv"), row.names = FALSE)
  write.csv(trial_results$all_effects, file = paste0("outputs/estimates/effect_size/", trial_slug, "_all_effects.csv"), row.names = FALSE)
}