
# INPUT: filtered verbose output w/ infections signaled (for two arms)
# OUTPUT: list with main outcomes and effect sizes computed


# DESCRIPTION:

# Wrapper function to run all functions that estimate all outcomes
# (prevalence, incidence, time-to-event), true or with PCD/ACD protocols,
# and protective effect size of intervention.

analyse_two_arm_trial <- function(infections_control, infections_intervention,
                                  # Trial and sim characteristics
                                  trial_start, trial_second_intervention,
                                  sim_length,
                                  # Protocols for PCD and ACD
                                  survey_protocol, acd_protocol) {
  
  year <- 356
  
  # Fetch all necessary functions
  
  source("functions/trial_outcomes.R")
  source("functions/trial_outcomes_event_time.R")
  source("functions/trial_effects.R")
  
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
    infections_control %>% estimate_survey_prevalence(trial_start = trial_start, !!!survey_protocol) %>% mutate(run = "Control"),
    infections_intervention %>% estimate_survey_prevalence(trial_start = trial_start, !!!survey_protocol) %>% mutate(run = "Intervention")
  )
  
  # Get estimates as if measured with ACD:
  # Use routine visits to estimate incidence (protocol can be modified)
  
  estimates_acd <- bind_rows(
    infections_control %>% estimate_acd_incidence(trial_start = trial_start, !!!acd_protocol) %>% mutate(run = "Control"),
    infections_intervention %>% estimate_acd_incidence(trial_start = trial_start, !!!acd_protocol) %>% mutate(run = "Intervention")
  )
  
  # Put all estimates together
  
  estimates_all <- bind_rows(estimates_true, estimates_survey, estimates_acd)
  
  ## PROTECTIVE EFFECT (prev and inc)
  
  # Prepare our estimates to get protective effect sizes
  
  relative_effect <- estimates_all %>%
    tidy_outcomes_for_effect() %>%
    estimate_relative_effect() %>%
    filter(!is.na(effect) & is.numeric(effect))
  
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
    infections_control %>% estimate_acd_time_to_event(time_inter = trial_start*year, !!!acd_protocol) %>% mutate(run = "Control"),
    infections_intervention %>% estimate_acd_time_to_event(time_inter = trial_start*year, !!!acd_protocol) %>% mutate(run = "Intervention")
  ) %>%
    prepare_time_to_event_survival(time_inter = trial_start*year, sim_length = sim_length*year)
  
  # Get time-to-event estimates for first infection/case from the second intervention as if doing ACD
  
  tte_acd_2 <- bind_rows(
    infections_control %>% estimate_acd_time_to_event(time_inter = (trial_start + trial_second_intervention)*year, !!!acd_protocol) %>% mutate(run = "Control"),
    infections_intervention %>% estimate_acd_time_to_event(time_inter = (trial_start + trial_second_intervention)*year, !!!acd_protocol) %>% mutate(run = "Intervention")
  ) %>%
    prepare_time_to_event_survival(time_inter = (trial_start + trial_second_intervention)*year, sim_length = sim_length*year)
  
  # Put all estimates together and calculate hazard ratios
  
  hazard_ratio <- bind_rows(
    estimate_hazard_ratio(tte_true_1, "time_to_case", "ever_case") %>% mutate(type_measure = "True Time-to-event", timestep = (trial_start + trial_second_intervention)*year),
    estimate_hazard_ratio(tte_true_2, "time_to_case", "ever_case") %>% mutate(type_measure = "True Time-to-event", timestep = sim_length*year),
    estimate_hazard_ratio(tte_acd_1, "time_to_case", "ever_case") %>% mutate(type_measure = "Time-to-event w/ ACD visits", timestep = (trial_start + trial_second_intervention)*year),
    estimate_hazard_ratio(tte_acd_2, "time_to_case", "ever_case") %>% mutate(type_measure = "Time-to-event w/ ACD visits", timestep = sim_length*year)
  ) %>%
    mutate(measure = "Case Hazard Ratio",
           effect = 1 - hazard_ratio)
  
  ## SAVE
  
  # Get a list with all the estimates and effect sizes for further uses
  
  list(
    estimates_true = estimates_true,
    estimates_survey = estimates_survey,
    estimates_acd = estimates_acd,
    estimates_all = estimates_all,
    tte_true_1 = tte_true_1,
    tte_true_2 = tte_true_2,
    tte_acd_2 = tte_acd_2,
    relative_effect = relative_effect,
    hazard_ratio = hazard_ratio,
    all_effects = bind_rows(relative_effect, hazard_ratio)
  )
}