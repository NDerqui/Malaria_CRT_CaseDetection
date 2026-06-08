#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#-#-#-#-#-#-#-#-# Modelling Malaria CRT - Case Detection #-#-#-#-#-#-#-#-#-#-#-
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#


## Modelling case detection in malaria CRTs using verbose simulations:
## statistical analyses from analyses cohort (read from previous step).



# SET-UP ------------------------------------------------------------------


rm(list = ls())

# Packages

library(tidyverse)
library(ggpubr)


#### trial conditions ####

year <- 365
month <- 30

# General inputs: init_EIR, and
# total sim size versus pop followed (analyses performed) at trial.

init_EIR <- 25

human_population <- 10000
trial_size <- 100

# Length of sim, length of trial, time of start interventions.

sim_length <- 9

trial_start <- 3
trial_second_intervention <- 3

key_intervention_time <- c(trial_start, trial_start+trial_second_intervention)

# Trial name of the simulation we want to analyse

trial_name <- paste0("Seasonal Init EIR ", init_EIR)


#### data ####

## Read the analyses cohort data (age cohort) from previous step

analyses_cohort_control <- read.csv(paste0("outputs_agecohort_data/", gsub(" ", "_", tolower(trial_name)), "_control.csv"))
analyses_cohort_bednet <- read.csv(paste0("outputs_agecohort_data/", gsub(" ", "_", tolower(trial_name)), "_bednet.csv"))


#### infections/cases ####

## Functions to signal infection/case, at time and overall

source("functions/verbose_infection_state.R")

# Run for each

infections_control <- analyses_cohort_control %>%
  detect_ever_malaria() %>% detect_infection()

infections_bednet <- analyses_cohort_bednet %>%
  detect_ever_malaria() %>% detect_infection()



# INCIDENCE / PREVALENCE ----------------------------------------------------------------


## Functions to get incidence/prevalence at each timestep

source("functions/verbose_prevalence_incidence.R")

## Factor levels for outcomes

measures <- c("prevalence_infection", "prevalence_case",
              "incidence_ppd_infection", "incidence_ppd_case",
              "incidence_ppy_infection", "incidence_ppy_case")
measures_labels <- c("Infection Prevalence", "Case Prevalence",
                     "Infection Incidence p.p.day", "Case Incidence p.p.day",
                     "Infection Incidence p.p.yea", "Case Incidence p.p.yea")


#### all estimates ####

## true estimates

estimates_control <- infections_control %>%
  true_realtime_measures() %>% mutate(run = "Control")

estimates_bednet <- infections_bednet %>%
  true_realtime_measures() %>% mutate(run = "Intervention")

## aggregate incidence

estimates_aggr_control <- infections_control %>%
  aggregate_incidence_period(trial_start = trial_start) %>% mutate(run = "Control")

estimates_aggr_bednet <- infections_bednet %>%
  aggregate_incidence_period(trial_start = trial_start) %>% mutate(run = "Intervention")

## cross-sectional survey prevalence

estimates_survey_control <- infections_control %>%
  survey_prevalence(trial_start = trial_start,
                    cross_surveys_in_years = seq(0.5, 6, 0.5)) %>% # Surveys every 6 months
  mutate(run = "Control")
estimates_survey_bednet <- infections_bednet %>%
  survey_prevalence(trial_start = trial_start,
                    cross_surveys_in_years = seq(0.5, 6, 0.5)) %>% # Surveys every 6 months
  mutate(run = "Intervention")

## routine ACD visits

estimates_visits_control <- infections_control %>%
  visits_incidence(trial_start = trial_start,
                   routine_visits_in_weeks = seq(4, 6*52, 4), # Surveys every 2 weeks
                   days_catchment = 2) %>%           # Cases appearing in last 48 h
  mutate(run = "Control")
estimates_visits_bednet <- infections_bednet %>%
  visits_incidence(trial_start = trial_start,
                   routine_visits_in_weeks = seq(4, 6*52, 4), # Surveys every 2 weeks
                   days_catchment = 2) %>%           # Cases appearing in last 48 h
  mutate(run = "Intervention")

## merge

estimates_all <- bind_rows(estimates_control, estimates_bednet,
                           estimates_aggr_control, estimates_aggr_bednet,
                           estimates_survey_control, estimates_survey_bednet,
                           estimates_visits_control, estimates_visits_bednet)

dir.create("outputs_effect_size/", showWarnings = FALSE)
write.csv(estimates_all, row.names = FALSE,
          file = paste0("outputs_effect_size/all_outputs_", gsub(" ", "_", tolower(trial_name)), ".csv"))


#### quick vis ####

source("functions/verbose_effect_size.R")
dir.create("outputs_effect_plots/", showWarnings = FALSE)

# Wrap up function (useful for effect estimate later)

plot_all_estimates <- wrap_for_plot_effect(estimates_all) %>%
  filter(!is.na(value))

require(rcartocolor)
require(ggh4x)

# all true estimates

png(filename = paste0("outputs_effect_plots/all_true_estimates_", gsub(" ", "_", tolower(trial_name)), ".png"),
    width = 12, height = 8, units = "in", res = 1200)
ggplot(data = filter(plot_all_estimates,
                     (type_measure == "True Instantaneous" & grepl("Prev", measure)) |
                       (grepl("aggr", type_measure) & grepl("p.p.y", measure) & !(grepl("ACD", type_measure)))),
       aes(x = timestep, y = value, group = run, color = run)) +
  geom_point() + geom_line() +
  geom_vline(xintercept = key_intervention_time*year, color = "firebrick", linetype = "dashed") +
  scale_color_manual(values = carto_pal(name = "Safe")[c(11, 10)]) +
  scale_x_continuous(breaks = seq(0, sim_length * year, by = year),
                     labels = (0:sim_length)) +
  labs(x = "Year", y = NULL,
       title = paste0("Simulated a ", human_population, " population, Sampled ", trial_size, " for trial, ", trial_name, ": True estimates")) +
  theme_bw() + theme(legend.position = "bottom", legend.title = element_blank()) +
  facet_nested(type_measure + measure ~ ., scales = "free", drop = TRUE)
dev.off()

# comparing incidence measurements

png(filename = paste0("outputs_effect_plots/incidence_estimates_", gsub(" ", "_", tolower(trial_name)), ".png"),
    width = 12, height = 8, units = "in", res = 1200)
ggplot(data = filter(plot_all_estimates,
                     grepl("p.p.yea", measure) & !grepl("Ins", type_measure)),
       aes(x = timestep, y = value, group = run, color = run)) +
  geom_point() + geom_line() +
  geom_vline(xintercept = key_intervention_time*year, color = "firebrick", linetype = "dashed") +
  scale_color_manual(values = carto_pal(name = "Safe")[c(11, 10)]) +
  scale_x_continuous(breaks = seq(0, sim_length * year, by = year),
                     labels = (0:sim_length)) +
  labs(x = "Year", y = NULL,
       title = paste0("Simulated a ", human_population, " population, Sampled ", trial_size, " for trial, ", trial_name, ": Incidence estimates")) +
  theme_bw() + theme(legend.position = "bottom", legend.title = element_blank()) +
  facet_nested(type_measure + measure ~ ., scales = "free", drop = TRUE)
dev.off()

# Overlay cross-sectional

png(filename = paste0("outputs_effect_plots/prevalence_estimates_", gsub(" ", "_", tolower(trial_name)), ".png"),
    width = 12, height = 5, units = "in", res = 1200)
ggplot(data = filter(plot_all_estimates, grepl("Prevalence", measure)),
       aes(x = timestep, y = value, group = run, color = run)) +
  geom_point(aes(shape = type_measure, size = type_measure)) + geom_line() +
  scale_shape_manual(breaks = c("True Instantaneous", "Cross-sectional surveys"),
                     values = c(16:17)) +
  scale_size_manual(breaks = c("True Instantaneous", "Cross-sectional surveys"),
                    values = c(1, 4)) +
  geom_vline(xintercept = key_intervention_time*year, color = "firebrick", linetype = "dashed") +
  scale_color_manual(values = carto_pal(name = "Safe")[c(11, 10)]) +
  scale_x_continuous(breaks = seq(0, sim_length * year, by = year),
                     labels = (0:sim_length)) +
  labs(x = "Year", y = NULL,
       title = paste0("Simulated a ", human_population, " population, Sampled ", trial_size, " for trial, ", trial_name, ": Prevalence estimates")) +
  theme_bw() + theme(legend.position = "bottom", legend.title = element_blank()) +
  facet_grid(measure ~ ., scales = "free")
dev.off()


#### effect size ####

source("functions/verbose_effect_size.R")

# Get the protective effect of our intervention,
# comparing incidence and prevalence

plot_all_estimates <- wrap_for_plot_effect(estimates_all)

protective_effect <- get_relative_effect(df = plot_all_estimates) %>%
  filter(!is.na(effect) & is.numeric(effect))

write.csv(protective_effect, row.names = FALSE,
          file = paste0("outputs_effect_size/protective_effect_", gsub(" ", "_", tolower(trial_name)), ".csv"))

png(filename = paste0("outputs_effect_plots/protective_effect_", gsub(" ", "_", tolower(trial_name)), ".png"),
    width = 12, height = 8, units = "in", res = 1200)
ggplot(data = filter(protective_effect, !is.na(effect) &
                       !(type_measure == "True Instantaneous" & grepl("Incidence", measure)) &
                       (grepl("Infection Prev", measure) | grepl("Case Incidence p.p.y", measure))),
       aes(x = timestep, y = effect)) +
  geom_point() + geom_line() +
  geom_vline(xintercept = key_intervention_time*year, color = "firebrick", linetype = "dashed") +
  scale_x_continuous(breaks = seq(0, sim_length * year, by = year),
                     labels = (0:sim_length)) +
  scale_y_continuous(labels = scales::percent, limits = 0:1) +
  labs(x = "Year", y = "Intervention Protective Effect",
       title = paste0("Simulated a ", human_population, " population, Sampled ", trial_size, " for trial, ", trial_name, ": Intervention Effect Size")) +
  theme_bw() + theme(legend.position = "bottom", legend.title = element_blank()) +
  facet_nested(type_measure + measure ~ ., scales = "free")
dev.off()


#### toying ####

# Modify parameters in function

results <- data.frame()
cross_surveys <- c(0.25,0.5, 0.75, 1)

for (test in 1:length(cross_surveys)) {
  
  # Iter over prevalence surveys
  
  control <- infections_control %>%
    survey_prevalence(trial_start = trial_start,
                      cross_surveys_in_years = seq(cross_surveys[test], 6, cross_surveys[test])) %>%
    mutate(run = "Control")
  bednet <- infections_bednet %>%
    survey_prevalence(trial_start = trial_start,
                      cross_surveys_in_years = seq(cross_surveys[test], 6, cross_surveys[test])) %>%
    mutate(run = "Intervention")
  
  prevalence <- rbind(control, bednet) %>%
    wrap_for_plot_effect() %>% get_relative_effect() %>%
    mutate(var = paste0(cross_surveys[test]*12, " months"))
  
  # Add to previous results
  
  results <- rbind(results, prevalence)
  rm(control, bednet, prevalence)
  
}
rm(test)

p1 <- ggplot(data = filter(results, grepl("Infection Prev", measure)),
       aes(x = timestep, y = effect, color = var)) +
  geom_jitter(size = 3) +
  geom_vline(xintercept = key_intervention_time*year, color = "firebrick", linetype = "dashed") +
  scale_color_manual(breaks = paste0(cross_surveys*12, " months"),
                     values = carto_pal(name = "Safe")[c(2, 3, 8, 11)]) +
  scale_x_continuous(breaks = seq(0, sim_length * year, by = year),
                     labels = (0:sim_length)) +
  scale_y_continuous(labels = scales::percent, limits = 0:1) +
  labs(x = "Year", y = "Intervention Protective Effect",
       title = paste0("Simulated a ", human_population, " population, Sampled ", trial_size, " for trial, ", trial_name)) +
  theme_bw() + theme(legend.position = "bottom", legend.title = element_blank()) +
  facet_nested(type_measure + measure ~ ., scales = "free")

# Modify parameters in function

results2 <- data.frame()
visits <- c(2, 4, 8)
visits_period <- c(2, 4, 7)

for (test in 1:length(visits)) {
  
  for (test2 in 1:length(visits_period)) {
    
    # Iter over incidence visits
    
    control <- infections_control %>%
      visits_incidence(trial_start = trial_start,
                       routine_visits_in_weeks = seq(visits[test], 6*52, visits[test]),
                       days_catchment = visits_period[test2]) %>%
      mutate(run = "Control")
    bednet <- infections_bednet %>%
      visits_incidence(trial_start = trial_start,
                       routine_visits_in_weeks = seq(visits[test], 6*52, visits[test]),
                       days_catchment = visits_period[test2]) %>%
      mutate(run = "Intervention")
    
    incidence <- rbind(control, bednet) %>%
      wrap_for_plot_effect() %>% get_relative_effect() %>%
      mutate(var = paste0(visits[test], " weeks - ", visits_period[test2], " days window"))
    
    # Add to previous results
    
    results2 <- rbind(results2, incidence)
    rm(control, bednet, incidence)
  }
  
}
rm(test, test2)

results2$type_measure <- "ACD visits, aggr. over 12 mos."

p2 <- ggplot(data = filter(results2, grepl("Case Incidence p.p.ye", measure)),
       aes(x = timestep, y = effect, color = var)) +
  geom_jitter(size = 3) +
  geom_vline(xintercept = key_intervention_time*year, color = "firebrick", linetype = "dashed") +
  scale_color_manual(values = c(carto_pal(name = "Safe")[c(9, 10, 2, 3, 8, 4, 7, 1, 11, 5, 12)], "black")) +
  scale_x_continuous(breaks = seq(0, sim_length * year, by = year),
                     labels = (0:sim_length)) +
  scale_y_continuous(labels = scales::percent, limits = 0:1) +
  labs(x = "Year", y = "Intervention Protective Effect") +
  theme_bw() + theme(legend.position = "bottom", legend.title = element_blank()) +
  facet_nested(type_measure + measure ~ ., scales = "free")

png(filename = paste0("outputs_effect_plots/prot_eff_pars_test_", gsub(" ", "_", tolower(trial_name)), ".png"),
    width = 12, height = 8, units = "in", res = 1200)
ggarrange(p1, p2, ncol = 1, common.legend = FALSE)
dev.off()
rm(p1, p2)



# TIME-TO-EVENT -----------------------------------------------------------


## Function to get time to first infection / case

source("functions/verbose_infection_time.R")


#### time-to-infection ####

## For time since (first intervention)

event_time_control <- infections_control %>%
  get_time_to_event(time_inter = trial_start*year) %>% mutate(run = "Control")

event_time_bednet <- infections_bednet %>%
  get_time_to_event(time_inter = trial_start*year) %>% mutate(run = "Intervention")

# Survival analysis

library(survival)
library(tidycmprsk)

# Prepare for survival analysis (with censoring time)

plot_survival_data <- rbind(event_time_control, event_time_bednet) %>%
  prepare_survival(time_inter = trial_start*year, sim_length = sim_length*year)

# Run models

surv_fit_1infection <-  survfit(Surv(time_to_infection, ever_infected) ~ run, data = plot_survival_data)
surv_fit_1case <-  survfit(Surv(time_to_case, ever_case) ~ run, data = plot_survival_data)

surv_fit_1infection <- tidy(surv_fit_1infection) %>%
  mutate(strata = gsub("run=", "", strata))
surv_fit_1case <- tidy(surv_fit_1case) %>%
  mutate(strata = gsub("run=", "", strata))

# Plot

plot_survival_a <- ggplot(data = surv_fit_1infection,
                          aes(x = time, y = estimate, color = strata, fill = strata)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  geom_point(data = surv_fit_1infection[surv_fit_1infection$n.censor != 0,],
             aes(x = time, y = estimate,color = strata, fill = strata),
             shape = 4, size = 4) +
  scale_color_manual(values = carto_pal(name = "Safe")[c(11, 10)]) +
  scale_fill_manual(values = carto_pal(name = "Safe")[c(11, 10)]) +
  scale_x_continuous(breaks = seq(0, sim_length * year, by = year),
                     labels = (0:sim_length)) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Year after trial start", y = "Proportion without infection") +
  theme_bw() + theme(legend.position = "bottom", legend.title = element_blank()) 

plot_survival_b <- ggplot(data = surv_fit_1case,
                          aes(x = time, y = estimate, color = strata, fill = strata)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  geom_point(data = surv_fit_1case[surv_fit_1case$n.censor != 0,],
             aes(x = time, y = estimate,color = strata, fill = strata),
             shape = 4, size = 4) +
  scale_color_manual(values = carto_pal(name = "Safe")[c(11, 10)]) +
  scale_fill_manual(values = carto_pal(name = "Safe")[c(11, 10)]) +
  scale_x_continuous(breaks = seq(0, sim_length * year, by = year),
                     labels = (0:sim_length)) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Year after trial start", y = "Proportion without clinical case") +
  theme_bw() + theme(legend.position = "bottom", legend.title = element_blank())

png(filename = paste0("outputs_effect_plots/true_time_to_event_", gsub(" ", "_", tolower(trial_name)), ".png"),
    width = 8, height = 8, units = "in", res = 1200)
annotate_figure(ggarrange(plot_survival_a, plot_survival_b, nrow = 2),
                top = paste0("Simulated a ", human_population, " population, Sampled ", trial_size, " for trial, ", trial_name, ": Time-to-event"))
dev.off()

## For time since second intervention

event_time_control <- infections_control %>%
  get_time_to_event(time_inter = (trial_start + trial_second_intervention)*year) %>% mutate(run = "Control")

event_time_bednet <- infections_bednet %>%
  get_time_to_event(time_inter = (trial_start + trial_second_intervention)*year) %>% mutate(run = "Intervention")

# Survival analysis

# Prepare for survival analysis (with censoring time)

plot_survival_data_2 <- rbind(event_time_control, event_time_bednet) %>%
  prepare_survival(time_inter = (trial_start + trial_second_intervention)*year,
                   sim_length = sim_length*year)

# Run models

surv_fit_2infection <-  survfit(Surv(time_to_infection, ever_infected) ~ run, data = plot_survival_data_2)
surv_fit_2case <-  survfit(Surv(time_to_case, ever_case) ~ run, data = plot_survival_data_2)

surv_fit_2infection <- tidy(surv_fit_2infection) %>%
  mutate(strata = gsub("run=", "", strata))
surv_fit_2case <- tidy(surv_fit_2case) %>%
  mutate(strata = gsub("run=", "", strata))

# Plot

plot_survival_2a <- ggplot(data = surv_fit_2infection,
                           aes(x = time, y = estimate, color = strata, fill = strata)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  geom_point(data = surv_fit_2infection[surv_fit_2infection$n.censor != 0,],
             aes(x = time, y = estimate,color = strata, fill = strata),
             shape = 4, size = 4) +
  scale_color_manual(values = carto_pal(name = "Safe")[c(11, 10)]) +
  scale_fill_manual(values = carto_pal(name = "Safe")[c(11, 10)]) +
  scale_x_continuous(breaks = seq(0, sim_length * year, by = year),
                     labels = (0:sim_length)) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Year after second intervention", y = "Proportion without infection") +
  theme_bw() + theme(legend.position = "bottom", legend.title = element_blank()) 

plot_survival_2b <- ggplot(data = surv_fit_2case,
                           aes(x = time, y = estimate, color = strata, fill = strata)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  geom_point(data = surv_fit_2case[surv_fit_2case$n.censor != 0,],
             aes(x = time, y = estimate,color = strata, fill = strata),
             shape = 4, size = 4) +
  scale_color_manual(values = carto_pal(name = "Safe")[c(11, 10)]) +
  scale_fill_manual(values = carto_pal(name = "Safe")[c(11, 10)]) +
  scale_x_continuous(breaks = seq(0, sim_length * year, by = year),
                     labels = (0:sim_length)) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Year after second intervention", y = "Proportion without clinical case") +
  theme_bw() + theme(legend.position = "bottom", legend.title = element_blank())

png(filename = paste0("outputs_effect_plots/true_time_to_event2_", gsub(" ", "_", tolower(trial_name)), ".png"),
    width = 8, height = 8, units = "in", res = 1200)
annotate_figure(ggarrange(plot_survival_2a, plot_survival_2b, nrow = 2),
                top = paste0("Simulated a ", human_population, " population, Sampled ", trial_size, " for trial, ", trial_name, ": Time-to-event"))
dev.off()


#### hazard ratio effect ####

source("functions/verbose_effect_size.R")

# HR for first event from first intervention

hr_1infection <- get_hazard_ratio(df = plot_survival_data,
                                  time_col = "time_to_infection",
                                  event_col = "ever_infected") %>%
  mutate(timestep = (trial_start+trial_second_intervention) * year,
         type_measure = "True Time-to-event",
         measure = "Infection Hazard Ratio",
         effect = 1 - hazard_ratio)

hr_1case <- get_hazard_ratio(df = plot_survival_data,
                             time_col = "time_to_case",
                             event_col = "ever_case") %>%
  mutate(timestep = (trial_start+trial_second_intervention) * year,
         type_measure = "True Time-to-event",
         measure = "Case Hazard Ratio",
         effect = 1 - hazard_ratio)

# HR for first event from second intervention

hr_2infection <- get_hazard_ratio(df = plot_survival_data_2,
                                  time_col = "time_to_infection",
                                  event_col = "ever_infected") %>%
  mutate(timestep = sim_length * year,
         type_measure = "True Time-to-event",
         measure = "Infection Hazard Ratio",
         effect = 1 - hazard_ratio)

hr_2case <- get_hazard_ratio(df = plot_survival_data_2,
                             time_col = "time_to_case",
                             event_col = "ever_case") %>%
  mutate(timestep = sim_length * year,
         type_measure = "True Time-to-event",
         measure = "Case Hazard Ratio",
         effect = 1 - hazard_ratio)

hr_effect <- bind_rows(hr_1infection, hr_1case,
                       hr_2infection, hr_2case)
rm(hr_1infection, hr_1case,
   hr_2infection, hr_2case)

# Save as before

write.csv(hr_effect, row.names = FALSE,
          file = paste0("outputs_effect_size/timetoevent_hr_", gsub(" ", "_", tolower(trial_name)), ".csv"))

# Quiz viz

# png(filename = paste0("outputs_effect_plots/protective_effect_w_hr_", gsub(" ", "_", tolower(trial_name)), ".png"),
#     width = 12, height = 8, units = "in", res = 1200)
ggplot(data = filter(bind_rows(protective_effect, hr_effect), !is.na(effect) &
                     !(type_measure == "True Instantaneous" & grepl("Incidence", measure)) &
                     (grepl("Infection Prev", measure) | grepl("Case Incidence p.p.y", measure) | grepl("Time", type_measure))),
       aes(x = timestep, y = effect)) +
  geom_point() + geom_line() +
  geom_vline(xintercept = key_intervention_time*year, color = "firebrick", linetype = "dashed") +
  scale_x_continuous(breaks = seq(0, sim_length * year, by = year),
                     labels = (0:sim_length)) +
  scale_y_continuous(labels = scales::percent, limits = 0:1) +
  labs(x = "Year", y = "Intervention Protective Effect",
       title = paste0("Simulated a ", human_population, " population, Sampled ", trial_size, " for trial, ", trial_name, ": Intervention Effect Size")) +
  theme_bw() + theme(legend.position = "bottom", legend.title = element_blank()) +
  facet_nested(type_measure + measure ~ ., scales = "free")
# dev.off()


#### time-to-infection with ACDs ####

## For time since second intervention

event_time_visit_control <- infections_control %>%
  get_time_to_visit_case(time_inter = (trial_start + trial_second_intervention)*year,
                         routine_visits_in_weeks = seq(4, 6*52, 4),
                         days_catchment = 2) %>% mutate(run = "Control")

event_time_visit_bednet <- infections_bednet %>%
  get_time_to_visit_case(time_inter = (trial_start + trial_second_intervention)*year,
                         routine_visits_in_weeks = seq(4, 6*52, 4),
                         days_catchment = 2) %>% mutate(run = "Intervention")

# Survival analysis

# Prepare for survival analysis (with censoring time)

plot_survival_data_3 <- rbind(event_time_visit_control, event_time_visit_bednet) %>%
  prepare_survival(time_inter = (trial_start + trial_second_intervention)*year,
                   sim_length = sim_length*year)

# Run models

surv_fit_3infection <-  survfit(Surv(time_to_infection, ever_infected) ~ run, data = plot_survival_data_3)
surv_fit_3case <-  survfit(Surv(time_to_case, ever_case) ~ run, data = plot_survival_data_3)

surv_fit_3infection <- tidy(surv_fit_3infection) %>%
  mutate(strata = gsub("run=", "", strata))
surv_fit_3case <- tidy(surv_fit_3case) %>%
  mutate(strata = gsub("run=", "", strata))

# Plot

plot_survival_3a <- ggplot(data = surv_fit_3infection,
                           aes(x = time, y = estimate, color = strata, fill = strata)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  geom_point(data = surv_fit_3infection[surv_fit_3infection$n.censor != 0,],
             aes(x = time, y = estimate,color = strata, fill = strata),
             shape = 4, size = 4) +
  scale_color_manual(values = carto_pal(name = "Safe")[c(11, 10)]) +
  scale_fill_manual(values = carto_pal(name = "Safe")[c(11, 10)]) +
  scale_x_continuous(breaks = seq(0, sim_length * year, by = year),
                     labels = (0:sim_length)) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Year after second intervention", y = "Proportion without infection") +
  theme_bw() + theme(legend.position = "bottom", legend.title = element_blank()) 

plot_survival_3b <- ggplot(data = surv_fit_3case,
                           aes(x = time, y = estimate, color = strata, fill = strata)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  geom_point(data = surv_fit_3case[surv_fit_3case$n.censor != 0,],
             aes(x = time, y = estimate,color = strata, fill = strata),
             shape = 4, size = 4) +
  scale_color_manual(values = carto_pal(name = "Safe")[c(11, 10)]) +
  scale_fill_manual(values = carto_pal(name = "Safe")[c(11, 10)]) +
  scale_x_continuous(breaks = seq(0, sim_length * year, by = year),
                     labels = (0:sim_length)) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Year after second intervention", y = "Proportion without clinical case") +
  theme_bw() + theme(legend.position = "bottom", legend.title = element_blank())

png(filename = paste0("outputs_effect_plots/time_to_event2_w_acd_", gsub(" ", "_", tolower(trial_name)), ".png"),
    width = 8, height = 8, units = "in", res = 1200)
annotate_figure(ggarrange(plot_survival_3a, plot_survival_3b, nrow = 2),
                top = paste0("Simulated a ", human_population, " population, Sampled ", trial_size, " for trial, ", trial_name, ": Time-to-event w/ ACD visits"))
dev.off()

# HR for first event from second intervention

hr_3infection <- get_hazard_ratio(df = plot_survival_data_3,
                                  time_col = "time_to_infection",
                                  event_col = "ever_infected") %>%
  mutate(timestep = sim_length * year,
         type_measure = "Time-to-event w/ ACD visits",
         measure = "Infection Hazard Ratio",
         effect = 1 - hazard_ratio)

hr_3case <- get_hazard_ratio(df = plot_survival_data_3,
                             time_col = "time_to_case",
                             event_col = "ever_case") %>%
  mutate(timestep = sim_length * year,
         type_measure = "Time-to-event w/ ACD visits",
         measure = "Case Hazard Ratio",
         effect = 1 - hazard_ratio)

# Save as before

write.csv(bind_rows(hr_effect, hr_3infection, hr_3case), row.names = FALSE,
          file = paste0("outputs_effect_size/timetoevent_hr_", gsub(" ", "_", tolower(trial_name)), ".csv"))

# Quiz viz

png(filename = paste0("outputs_effect_plots/protective_effect_w_hr_", gsub(" ", "_", tolower(trial_name)), ".png"),
    width = 12, height = 10, units = "in", res = 1200)
ggplot(data = filter(bind_rows(protective_effect, hr_effect, hr_3case), !is.na(effect) &
                       !(type_measure == "True Instantaneous" & grepl("Incidence", measure)) &
                       (grepl("Infection Prev", measure) | grepl("Case Incidence p.p.y", measure) |
                          (grepl("ime-to", type_measure) & grepl("Case",measure)))),
       aes(x = timestep, y = effect)) +
  geom_point() + geom_line() +
  geom_vline(xintercept = key_intervention_time*year, color = "firebrick", linetype = "dashed") +
  scale_x_continuous(breaks = seq(0, sim_length * year, by = year),
                     labels = (0:sim_length)) +
  scale_y_continuous(labels = scales::percent, limits = 0:1) +
  labs(x = "Year", y = "Intervention Protective Effect",
       title = paste0("Simulated a ", human_population, " population, Sampled ", trial_size, " for trial, ", trial_name, ": Intervention Effect Size")) +
  theme_bw() + theme(legend.position = "bottom", legend.title = element_blank()) +
  facet_nested(type_measure + measure ~ ., scales = "free")
dev.off()

rm(hr_3infection, hr_3case)
