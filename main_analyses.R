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



# INCIDENCE / PREVALENCE ----------------------------------------------------------------


#### infections/cases ####

## Functions to signal infection/case, at time and overall

source("functions/verbose_infection_state.R")

# Run for each

infections_control <- analyses_cohort_control %>%
  detect_ever_malaria() %>% detect_infection()

infections_bednet <- analyses_cohort_bednet %>%
  detect_ever_malaria() %>% detect_infection()


#### incidence/prevalence ####

## Functions to get incidence/prevalence at each timestep

source("functions/verbose_prevalence_incidence.R")

# Apply

estimates_control <- infections_control %>%
  get_prev_inc() %>% mutate(run = "Control")

estimates_bednet <- infections_bednet %>%
  get_prev_inc() %>% mutate(run = "Intervention")

# Quick vis

plot_prev_inc <- rbind(estimates_control, estimates_bednet) %>%
  select(-c(n, at_risk, infections, cases, new_infections, new_cases)) %>%
  pivot_longer(-c(timestep, run), names_to = "measure", values_to = "value") %>%
  mutate(measure = factor(measure,
                          levels = c("prevalence_infec", "prevalence_case", "incidence_infec", "incidence_case"),
                          labels = c("Infection Prevalence", "Case Prevalence", "Infection Incidence", "Case Incidence")))

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
  facet_grid(measure ~ ., scales = "free")
dev.off()

# Apply to age estimates

estimates_control_by_age <- infections_control %>%
  get_prev_inc_by_age() %>% mutate(run = "Control")

estimates_bednet_by_age <- infections_bednet %>%
  get_prev_inc_by_age() %>% mutate(run = "Intervention")

# Quick vis

plot_prev_inc_by_age <- rbind(estimates_control_by_age, estimates_bednet_by_age) %>%
  select(-c(n, at_risk, infections, cases, new_infections, new_cases)) %>%
  pivot_longer(-c(timestep, age_at_time_year, run), names_to = "measure", values_to = "value") %>%
  mutate(measure = factor(measure,
                          levels = c("prevalence_infec", "prevalence_case", "incidence_infec", "incidence_case"),
                          labels = c("Infection Prevalence", "Case Prevalence", "Infection Incidence", "Case Incidence")))

png(filename = paste0("outputs_plots/outcomes_prev_by_age_", gsub(" ", "_", tolower(trial_name)), ".png"),
    width = 12, height = 8, units = "in", res = 1200)
ggplot(data = filter(plot_prev_inc_by_age, age_at_time_year %in% c(1, 2, 5, 10)),
       aes(x = timestep, y = value, group = run, color = run)) +
  geom_point() + geom_line() +
  geom_vline(xintercept = key_intervention_time*year, color = "firebrick", linetype = "dashed") +
  scale_color_manual(values = carto_pal(name = "Safe")[c(11, 10)]) +
  scale_x_continuous(breaks = seq(0, sim_length * year, by = year),
                     labels = (0:sim_length)) +
  labs(x = "Year", y = NULL,
       title = paste0("Simulated a ", human_population, " population, Sampled ", trial_size, " for trial, ", trial_name)) +
  theme_bw() + theme(legend.position = "bottom", legend.title = element_blank()) +
  facet_grid(age_at_time_year ~ measure, scales = "free")

dev.off()



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