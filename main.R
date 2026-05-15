#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#-#-#-#-#-#-#-#-# Modelling Malaria CRT - Case Detection #-#-#-#-#-#-#-#-#-#-#-
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#


## Modelling case detection in malaria CRTs using verbose simulations.



# SET-UP ------------------------------------------------------------------


rm(list = ls())

# Packages

# remotes::install_github("mrc-ide/individual@dev")
# remotes::install_github("JasonRWood/malariasimulation@verbose_simulations")
library(individual)
# library(devtools)
# setwd("~/Code_base/malariasimulation/")
# load_all(".")
library(malariasimulation)

library(tidyverse)
library(ggpubr)


## Set up a couple of basic general options for simulations

year <- 365
month <- 30



# TRIAL SIM -----------------------------------------------------------------


#### trial conditions ####

# General inputs: init_EIR, and
# total sim size versus pop followed (analyses performed) at trial.

init_EIR <- 15

human_population <- 10000
trial_size <- 100

# Set some basics, like length of sim vs length of trial,
# and when do our interventions start.

# Example, run the sim for total 9 years, with 3 years "without intervention"
# (i.e., start the trial on the third year of the sim),
# and doing two rounds of intervention separated by 3 years.

sim_length <- 9

trial_start <- 3
trial_second_intervention <- 3

key_intervention_time <- c(trial_start, trial_start+trial_second_intervention)

# Control when we get our age snapshot (best at start of trial, timestep = 1)

snapshot_time <- 1

# Add a "trial name" to keep track of results

trial_name <- paste0("Basic Init EIR ", init_EIR)


#### run sim ####

## Functions to set up parameters and run the verbose simulation

source("functions/verbose_set_parameters.R")
source("functions/verbose_runsim.R")

## Basic parameters

# With this function, we set baseline parameters
# (including seasonality, treatment use, base bednet use, etc.)
# to be used in intervention and control runs.

baseline_parameters <- set_baseline_pars(sim_length = sim_length,
                                         init_EIR = init_EIR,
                                         human_population = human_population,
                                         seasonality = TRUE,
                                         treatment = TRUE,
                                         bednets = TRUE)

## Run sim

# Control arm

out <- run_verbose_sim(simparams = baseline_parameters, sim_length = sim_length,
                       snapshot_time = snapshot_time,
                       key_bednet = FALSE, run_note = "control")

df_control <- read.csv("outputs_verbose_sims/verbose_dumping_control.csv")

df_control$process <- out$process_vector[df_control$process_index]
df_control$state <- out$state_list[df_control$state_index]

write.csv(df_control, "outputs_verbose_sims/verbose_dumping_control.csv", row.names = FALSE)

# Intervention arm

out <- run_verbose_sim(simparams = baseline_parameters, sim_length = sim_length,
                       snapshot_time = snapshot_time,
                       key_bednet = TRUE, key_intervention_time = key_intervention_time,
                       bed_coverage = 0.95, run_note = "bednet")

df_bednet <- read.csv("outputs_verbose_sims/verbose_dumping_bednet.csv")

df_bednet$process <- out$process_vector[df_bednet$process_index]
df_bednet$state <- out$state_list[df_bednet$state_index]

write.csv(df_bednet, "outputs_verbose_sims/verbose_dumping_bednet.csv", row.names = FALSE)

## Read the verbose files

rm(out)

df_control <- read.csv("outputs_verbose_sims/verbose_dumping_control.csv")
df_control_age <- read.csv("outputs_verbose_sims/verbose_dumping_snapshot_control.csv")

df_bednet <- read.csv("outputs_verbose_sims/verbose_dumping_bednet.csv")
df_bednet_age <- read.csv("outputs_verbose_sims/verbose_dumping_snapshot_bednet.csv")

gc()


#### analyses cohort ####

## Simple clean to subtract to the cohort we can follow with age.

source("functions/verbose_dataclean_agecohort.R")
source("functions/verbose_dataclean_trial_sample.R")
source("functions/verbose_vis.R")

# Filter individuals born / with age from snapshot,
# estimate their age at each timestep and final age (at death or sim end),
# and sample (say 100).

analyses_cohort_control <- df_control %>%
  get_birth_death() %>%
  get_age_cohort(age_snapshot = df_control_age, snapshot_time = snapshot_time) %>%
  trial_sample(alive_by = trial_start * year, trial_size = trial_size)

analyses_cohort_bednet <- df_bednet %>%
  get_birth_death() %>%
  get_age_cohort(age_snapshot = df_bednet_age, snapshot_time = snapshot_time) %>%
  trial_sample(alive_by = trial_start * year, trial_size = trial_size)

# Clean space

rm(df_control, df_bednet, df_control_age, df_bednet_age)
gc()

# Save for future

dir.create("outputs_agecohort_data", showWarnings = FALSE)
write.csv(analyses_cohort_control, row.names = FALSE,
          file = paste0("outputs_agecohort_data/", gsub(" ", "_", tolower(trial_name)), "_control.csv"))
write.csv(analyses_cohort_bednet, row.names = FALSE,
          file = paste0("outputs_agecohort_data/", gsub(" ", "_", tolower(trial_name)), "_bednet.csv"))


## Check

# Plot the control (no bednet)

png(filename = "outputs_plots/agecohort_overtime_control.png",
    width = 8, height = 5, units = "in", res = 1200)
plot_verbose_itn(df = analyses_cohort_control,
                 note = "Control", sim_length = sim_length,
                 human_population = human_population, trial_size = trial_size,
                 bednetstimesteps = seq(0, sim_length, 3)*year)
dev.off()

# Plot the bednet

png(filename = "outputs_plots/agecohort_overtime_bednet.png",
    width = 8, height = 5, units = "in", res = 1200)
plot_verbose_itn(df = analyses_cohort_bednet,
                 note = "Intervention", sim_length = sim_length,
                 human_population = human_population, trial_size = trial_size,
                 bednetstimesteps = seq(0, sim_length, 3)*year) +
  geom_vline(xintercept = key_intervention_time*year, color = "firebrick", linetype = "dashed")
dev.off()



# ANALYSES ----------------------------------------------------------------


#### infections/cases ####

## Functions to signal infection/case, at time and overall

source("functions/verbose_infection_state.R")

# Run for each

infections_control <- analyses_cohort_control %>%
  ever_malaria() %>% detect_infection()

infections_bednet <- analyses_cohort_bednet %>%
  ever_malaria() %>% detect_infection()


#### incidence/prevalence ####

## Functions to get incidence/prevalence at each timestep

source("functions/verbose_prevalence_incidence.R")

# Apply

estimates_control <- infections_control %>%
  get_prev_inc() %>% mutate(run = "Control")

estimates_bednet <- infections_bednet %>%
  get_prev_inc() %>% mutate(run = "Intervention")

# Quick vis

plot <- rbind(estimates_control, estimates_bednet) %>%
  select(-c(n, at_risk, infections, cases, new_infections, new_cases)) %>%
  pivot_longer(-c(timestep, run), names_to = "measure", values_to = "value") %>%
  mutate(measure = factor(measure,
                          levels = c("prevalence_infec", "prevalence_case", "incidence_infec", "incidence_case"),
                          labels = c("Infection Prevalence", "Case Prevalence", "Infection Incidence", "Case Incidence")))

require(rcartocolor)

png(filename = "outputs_plots/outcomes.png",
    width = 12, height = 8, units = "in", res = 1200)
ggplot(data = plot,
       aes(x = timestep, y = value, group = run, color = run)) +
  geom_point() + geom_line() +
  geom_vline(xintercept = key_intervention_time*year, color = "firebrick", linetype = "dashed") +
  scale_color_manual(values = carto_pal(name = "Safe")[c(11, 10)]) +
  scale_x_continuous(breaks = seq(0, sim_length * year, by = year),
                     labels = (0:sim_length)) +
  labs(x = "Year", y = NULL,
       title = paste0("Simulated a ", human_population, " population, Sampled ", trial_size, " for trial")) +
  theme_bw() + theme(legend.position = "bottom", legend.title = element_blank()) +
  facet_grid(measure ~ ., scales = "free")
dev.off()

# Apply to age estimates

estimates_control_by_age <- infections_control %>%
  get_prev_inc_by_age() %>% mutate(run = "Control")

estimates_bednet_by_age <- infections_bednet %>%
  get_prev_inc_by_age() %>% mutate(run = "Intervention")

# Quick vis

plot_by_age <- rbind(estimates_control_by_age, estimates_bednet_by_age) %>%
  select(-c(n, at_risk, infections, cases, new_infections, new_cases)) %>%
  pivot_longer(-c(timestep, age_at_time_year, run), names_to = "measure", values_to = "value") %>%
  mutate(measure = factor(measure,
                          levels = c("prevalence_infec", "prevalence_case", "incidence_infec", "incidence_case"),
                          labels = c("Infection Prevalence", "Case Prevalence", "Infection Incidence", "Case Incidence")))

png(filename = "outputs_plots/outcomes_by_age.png",
    width = 12, height = 8, units = "in", res = 1200)
ggplot(data = filter(plot_by_age, age_at_time_year %in% c(1, 2, 5, 10)),
       aes(x = timestep, y = value, group = run, color = run)) +
  geom_point() + geom_line() +
  geom_vline(xintercept = key_intervention_time*year, color = "firebrick", linetype = "dashed") +
  scale_color_manual(values = carto_pal(name = "Safe")[c(11, 10)]) +
  scale_x_continuous(breaks = seq(0, sim_length * year, by = year),
                     labels = (0:sim_length)) +
  labs(x = "Year", y = NULL,
       title = paste0("Simulated a ", human_population, " population, Sampled ", trial_size, " for trial")) +
  theme_bw() + theme(legend.position = "bottom", legend.title = element_blank()) +
  facet_grid(age_at_time_year ~ measure, scales = "free")

dev.off()


#### time-to-event ####

## Function to get time to first infection / case

source("functions/verbose_infection_time.R")

# Apply

event_time_control <- infections_control %>%
  get_time_to_event(time_inter = trial_start*year) %>% mutate(run = "Control")

event_time_bednet <- infections_bednet %>%
  get_time_to_event(time_inter = trial_start*year) %>% mutate(run = "Intervention")

# Quick vis

library(survival)
library(tidycmprsk)

plot2 <- rbind(event_time_control, event_time_bednet)

surv_fit_a <-  survfit(Surv(time_to_infection, ever_infected) ~ run, data = plot2)
surv_fit_b <-  survfit(Surv(time_to_case, ever_case) ~ run, data = plot2)

surv_fit_a <- tidy(surv_fit_a) %>%
  mutate(strata = gsub("run=", "", strata))
surv_fit_b <- tidy(surv_fit_b) %>%
  mutate(strata = gsub("run=", "", strata))

plot2a <- ggplot(data = surv_fit_a, aes(x = time, y = estimate, color = strata, fill = strata)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.5) +
  scale_color_manual(values = carto_pal(name = "Safe")[c(11, 10)]) +
  scale_fill_manual(values = carto_pal(name = "Safe")[c(11, 10)]) +
  scale_x_continuous(breaks = seq(0, sim_length * year, by = year),
                     labels = (0:sim_length)) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Year after trial start", y = "Proportion without infection") +
  theme_bw() + theme(legend.position = "bottom", legend.title = element_blank()) 

plot2b <- ggplot(data = surv_fit_b, aes(x = time, y = estimate, color = strata, fill = strata)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.5) +
  scale_color_manual(values = carto_pal(name = "Safe")[c(11, 10)]) +
  scale_fill_manual(values = carto_pal(name = "Safe")[c(11, 10)]) +
  scale_x_continuous(breaks = seq(0, sim_length * year, by = year),
                     labels = (0:sim_length)) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Year after trial start", y = "Proportion without clinical case") +
  theme_bw() + theme(legend.position = "bottom", legend.title = element_blank())

png(filename = "outputs_plots/outcomes2.png",
    width = 8, height = 8, units = "in", res = 1200)
annotate_figure(ggarrange(plot2a, plot2b, nrow = 2),
                top = paste0("Simulated a ", human_population, " population, Sampled ", trial_size, " for trial"))
dev.off()

# Cox and HZ

# For infection (with and without age)
summary(coxph(Surv(time_to_infection, ever_infected) ~ run, data = plot2))
summary(coxph(Surv(time_to_infection, ever_infected) ~ run + age_at_first_infection, data = plot2))

# For clinical case (with and without age)
summary(coxph(Surv(time_to_case, ever_case) ~ run, data = plot2))
summary(coxph(Surv(time_to_case, ever_case) ~ run + age_at_first_case, data = plot2))