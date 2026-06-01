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

init_EIR <- 1

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



measures <- c("prevalence_infection", "prevalence_case", "incidence_infection", "incidence_case")
measures_labels <- c("Infection Prevalence", "Case Prevalence", "Infection Incidence", "Case Incidence")

## Functions to get incidence/prevalence at each timestep

source("functions/verbose_prevalence_incidence.R")


#### true estimates ####

# Apply

estimates_control <- infections_control %>%
  true_realtime_measures() %>% mutate(run = "Control")

estimates_bednet <- infections_bednet %>%
  true_realtime_measures() %>% mutate(run = "Intervention")

# Quick vis

plot_prev_inc <- rbind(estimates_control, estimates_bednet) %>%
  select(-c(n, person_days_at_risk, infections, cases, new_infections, new_cases)) %>%
  pivot_longer(-c(timestep, type_measure, run), names_to = "measure", values_to = "value") %>%
  mutate(measure = factor(measure,
                          levels = measures,
                          labels = measures_labels))

require(rcartocolor)

png(filename = paste0("outputs_plots/outcomes_prev_inc_", gsub(" ", "_", tolower(trial_name)), ".png"),
    width = 12, height = 8, units = "in", res = 1200)
ggplot(data = plot_prev_inc,
       aes(x = timestep, y = value, group = run, color = run)) +
  geom_point() + geom_line() +
  geom_vline(xintercept = key_intervention_time*year, color = "firebrick", linetype = "dashed") +
  scale_color_manual(values = carto_pal(name = "Safe")[c(11, 10)]) +
  scale_x_continuous(breaks = seq(0, sim_length * year, by = year),
                     labels = (0:sim_length)) +
  labs(x = "Year", y = NULL,
       title = paste0("Simulated a ", human_population, " population, Sampled ", trial_size, " for trial, ", trial_name)) +
  theme_bw() + theme(legend.position = "bottom", legend.title = element_blank()) +
  facet_grid(type_measure*measure ~ ., scales = "free")
dev.off()


#### with cross-sectional surveys ####

# Get prevalence estimates as if by a cross sectional survey 
# (imp, cross surveys passed as year measure, i.e. 6 month = 0.5 year)

estimates_control_prev_survey <- infections_control %>%
  get_prev_survey(trial_start = trial_start, cross_surveys = seq(0.5, 6, 0.5)) %>% # Surveys every 6 months
  mutate(run = "Prevalence Survey - Control")
estimates_bednet_prev_survey <- infections_bednet %>%
  get_prev_survey(trial_start = trial_start, cross_surveys = seq(0.5, 6, 0.5)) %>% # Surveys every 6 months
  mutate(run = "Prevalence Survey - Intervention")

plot_prev_survey <- rbind(estimates_control_prev_survey, estimates_bednet_prev_survey) %>%
  select(-c(n, infections, cases)) %>%
  pivot_longer(-c(timestep, run), names_to = "measure", values_to = "value") %>%
  mutate(measure = factor(measure,
                          levels = c("prevalence_infec", "prevalence_case", "incidence_infec", "incidence_case"),
                          labels = c("Infection Prevalence", "Case Prevalence", "Infection Incidence", "Case Incidence")))

# Get incidence estimates as if by visits

estimates_control_inc_survey <- infections_control %>%
  get_inc_survey(trial_start = trial_start,
                 routine_visits = seq(4, 6*52, 4), # Surveys every 2 weeks
                 days_catchment = 2) %>%           # Cases appearing in last 48 h
  mutate(run = "Incidence Routine Visits - Control")
estimates_bednet_inc_survey <- infections_bednet %>%
  get_inc_survey(trial_start = trial_start,
                 routine_visits = seq(4, 6*52, 4), # Surveys every 2 weeks
                 days_catchment = 2) %>%           # Cases appearing in last 48 h
  mutate(run = "Incidence Routine Visits - Intervention")

plot_inc_survey <- rbind(estimates_control_inc_survey, estimates_bednet_inc_survey) %>%
  select(-c(n, at_risk, new_infections, new_cases)) %>%
  pivot_longer(-c(timestep, run), names_to = "measure", values_to = "value") %>%
  mutate(measure = factor(measure,
                          levels = c("prevalence_infec", "prevalence_case", "incidence_infec", "incidence_case"),
                          labels = c("Infection Prevalence", "Case Prevalence", "Infection Incidence", "Case Incidence")))

# Overlay to old plot to see if matches

png(filename = paste0("outputs_plots/outcomes_survey_prev_inc_", gsub(" ", "_", tolower(trial_name)), ".png"),
    width = 12, height = 8, units = "in", res = 1200)
ggplot(data = plot_prev_inc,
       aes(x = timestep, y = value, group = run, color = run)) +
  geom_point(alpha = 0.75) + geom_line(alpha = 0.75) +
  geom_vline(xintercept = key_intervention_time*year, color = "firebrick", linetype = "dashed") +
  geom_point(data = plot_prev_survey %>%
               mutate(run = gsub(" - Intervention", "", gsub(" - Control", "", run))),
             aes(x = timestep, y = value, group = run, color = run), size = 3) +
  geom_point(data = plot_inc_survey %>%
               mutate(run = gsub(" - Intervention", "", gsub(" - Control", "", run))),
             aes(x = timestep, y = value, group = run, color = run), size = 1.5) +
  scale_color_manual(values = carto_pal(name = "Safe")[c(11, 10, 3, 8)],
                     breaks = c("Control", "Intervention", "Prevalence Survey", "Incidence Routine Visits")) +
  scale_x_continuous(breaks = seq(0, sim_length * year, by = year),
                     labels = (0:sim_length)) +
  labs(x = "Year", y = NULL,
       title = paste0("Simulated a ", human_population, " population, Sampled ", trial_size, " for trial, ", trial_name)) +
  theme_bw() + theme(legend.position = "bottom", legend.title = element_blank()) +
  facet_grid(measure ~ ., scales = "free")
dev.off()

# Save estimates with effect size

prev_ic_effect_survey <- rbind(plot_prev_survey, plot_inc_survey) %>%
  mutate(run = gsub("Incidence Routine Visits - ", "", gsub("Prevalence Survey - ", "", run))) %>%
  arrange(timestep, run) %>%
  group_by(timestep, measure) %>%
  mutate(
    est_control = value[run == "Control"],
    est_intervention = value[run == "Intervention"],
    effect = (1 - est_intervention/est_control)
  ) %>% filter(row_number() == 1) %>%
  ungroup() %>% select(-c(est_control, est_intervention, value, run)) %>%
  mutate(measure = gsub(" ", "", measure)) %>%
  mutate(effect = as.numeric(gsub("NaN", NA, effect))) %>%
  pivot_wider(names_from = "measure", values_from = "effect") 

dir.create("outputs_effect_size/", showWarnings = FALSE)
write.csv(prev_ic_effect_survey, row.names = FALSE,
          file = paste0("outputs_effect_size/survey_prev_inc_", gsub(" ", "_", tolower(trial_name)), ".csv"))

png(filename = paste0("outputs_plots/effect_survey_prev_inc_", gsub(" ", "_", tolower(trial_name)), ".png"),
    width = 12, height = 8, units = "in", res = 1200)
ggplot(data = prev_ic_effect_survey %>% pivot_longer(-timestep, names_to = "measure", values_to = "value"),
       aes(x = timestep, y = value)) +
  geom_point() + geom_line() +
  geom_vline(xintercept = key_intervention_time*year, color = "firebrick", linetype = "dashed") +
  scale_x_continuous(breaks = seq(0, sim_length * year, by = year),
                     labels = (0:sim_length)) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Year", y = NULL,
       title = paste0("Simulated a ", human_population, " population, Sampled ", trial_size, " for trial, ", trial_name)) +
  theme_bw() + theme(legend.position = "bottom", legend.title = element_blank()) +
  facet_grid(measure ~ ., scales = "free")
dev.off()


#### effect size ####

# Save estimates with effect size

prev_ic_effect <- plot_prev_inc %>%
  arrange(timestep, run) %>%
  group_by(timestep, measure) %>%
  mutate(
    est_control = value[run == "Control"],
    est_intervention = value[run == "Intervention"],
    effect = (1 - est_intervention/est_control)
  ) %>% filter(row_number() == 1) %>%
  ungroup() %>% select(-c(est_control, est_intervention, value, run)) %>%
  mutate(measure = gsub(" ", "", measure)) %>%
  mutate(effect = as.numeric(gsub("NaN", NA, effect))) %>%
  pivot_wider(names_from = "measure", values_from = "effect") 

dir.create("outputs_effect_size/", showWarnings = FALSE)
write.csv(prev_ic_effect, row.names = FALSE,
          file = paste0("outputs_effect_size/prev_inc_", gsub(" ", "_", tolower(trial_name)), ".csv"))

png(filename = paste0("outputs_plots/effect_prev_inc_", gsub(" ", "_", tolower(trial_name)), ".png"),
    width = 12, height = 8, units = "in", res = 1200)
ggplot(data = prev_ic_effect %>% pivot_longer(-timestep, names_to = "measure", values_to = "value"),
       aes(x = timestep, y = value)) +
  geom_point() + geom_line() +
  geom_vline(xintercept = key_intervention_time*year, color = "firebrick", linetype = "dashed") +
  scale_x_continuous(breaks = seq(0, sim_length * year, by = year),
                     labels = (0:sim_length)) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Year", y = NULL,
       title = paste0("Simulated a ", human_population, " population, Sampled ", trial_size, " for trial, ", trial_name)) +
  theme_bw() + theme(legend.position = "bottom", legend.title = element_blank()) +
  facet_grid(measure ~ ., scales = "free")
dev.off()

# With new functions

source("functions/verbose_effect_size.R")

get_relative_effect(df = plot_prev_inc) %>%
  ggplot(aes(x = timestep, y = effect)) +
  geom_point() + geom_line() +
  geom_vline(xintercept = key_intervention_time*year, color = "firebrick", linetype = "dashed") +
  scale_x_continuous(breaks = seq(0, sim_length * year, by = year),
                     labels = (0:sim_length)) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Year", y = NULL,
       title = paste0("Simulated a ", human_population, " population, Sampled ", trial_size, " for trial, ", trial_name)) +
  theme_bw() + theme(legend.position = "bottom", legend.title = element_blank()) +
  facet_grid(measure ~ ., scales = "free")






# TIME-TO-EVENT -----------------------------------------------------------


## Function to get time to first infection / case

source("functions/verbose_infection_time.R")

#### time-to-first-infection ####

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

png(filename = paste0("outputs_plots/outcomes_time_to_", gsub(" ", "_", tolower(trial_name)), ".png"),
    width = 8, height = 8, units = "in", res = 1200)
annotate_figure(ggarrange(plot_survival_a, plot_survival_b, nrow = 2),
                top = paste0("Simulated a ", human_population, " population, Sampled ", trial_size, " for trial, ", trial_name))
dev.off()

# Cox and HZ

# For infection (with and without age)
summary(coxph(Surv(time_to_infection, ever_infected) ~ run, data = plot_survival_data))
summary(coxph(Surv(time_to_infection, ever_infected) ~ run + age_at_first_infection, data = plot_survival_data))

# For clinical case (with and without age)
summary(coxph(Surv(time_to_case, ever_case) ~ run, data = plot_survival_data))
summary(coxph(Surv(time_to_case, ever_case) ~ run + age_at_first_case, data = plot_survival_data))


#### time-to-second-itn ####

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
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.5) +
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
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.5) +
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

png(filename = paste0("outputs_plots/outcomes_time_to_2_", gsub(" ", "_", tolower(trial_name)), ".png"),
    width = 8, height = 8, units = "in", res = 1200)
annotate_figure(ggarrange(plot_survival_2a, plot_survival_2b, nrow = 2),
                top = paste0("Simulated a ", human_population, " population, Sampled ", trial_size, " for trial, ", trial_name))
dev.off()

# Cox and HZ

# For infection (with and without age)
summary(coxph(Surv(time_to_infection, ever_infected) ~ run, data = plot_survival_data_2))
summary(coxph(Surv(time_to_infection, ever_infected) ~ run + age_at_first_infection, data = plot_survival_data_2))

# For clinical case (with and without age)
summary(coxph(Surv(time_to_case, ever_case) ~ run, data = plot_survival_data_2))
summary(coxph(Surv(time_to_case, ever_case) ~ run + age_at_first_case, data = plot_survival_data_2))