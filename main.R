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


#### run sim ####

## Functions to set up parameters and run the verbose simulation

source("functions/verbose_par_set.R")
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
source("functions/verbose_plots_basic.R")

# Filter individuals born / with age from snapshot,
# estimate their age at each timestep and final age (at death or sim end),
# and sample (say 100).

analyses_cohort_control <- df_control %>%
  birth_death() %>%
  get_age_cohort(age_snapshot = df_control_age, snapshot_time = snapshot_time) %>%
  trial_sample(alive_by = trial_start * year, trial_size = trial_size)

analyses_cohort_bednet <- df_bednet %>%
  birth_death() %>%
  get_age_cohort(age_snapshot = df_bednet_age, snapshot_time = snapshot_time) %>%
  trial_sample(alive_by = trial_start * year, trial_size = trial_size)

# Clean space

rm(df_control, df_bednet, df_control_age, df_bednet_age)
gc()

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


#### true no infections / cases ####

## Signal all new infections and new clinical-infections.

source("functions/verbose_dataclean_cases_detect.R")

infection_cases_control <- age_cohort_control %>%
  # Filter to only get new infections AFTER trial starts
  filter(timestep >= key_intervention_time[1]) %>%
  detect_new_infection()

infection_cases_bednet <- age_cohort_bednet %>%
  # Filter to only get new infections AFTER trial starts
  filter(timestep >= key_intervention_time[1]) %>%
  detect_new_infection()

clinical_cases_control <- age_cohort_control %>%
  # Filter to only get new cases AFTER trial starts
  filter(timestep >= key_intervention_time[1]) %>%
  detect_new_clinical()

clinical_cases_bednet <- age_cohort_bednet %>%
  # Filter to only get new cases AFTER trial starts
  filter(timestep >= key_intervention_time[1]) %>%
  detect_new_clinical()

## Count per individual the true number new infections / cases

source("functions/verbose_dataclean_cases_sum.R")

# Infections

# All
count_infections(infection_cases_control) %>% select(total_infection) %>%
  summarise(total = sum(total_infection),
            min = min(total_infection), max = max(total_infection),
            mean = mean(total_infection), sd = sd(total_infection),
            median = median(total_infection)) %>% t() %>% round(digits = 1)
count_infections(infection_cases_bednet) %>% select(total_infection) %>%
  summarise(total = sum(total_infection),
            min = min(total_infection), max = max(total_infection), 
            mean = mean(total_infection), sd = sd(total_infection),
            median = median(total_infection)) %>% t() %>% round(digits = 1)

# By age
count_infection_by_age(infection_cases_control) %>%
  select(total_infection_year, age_year) %>% group_by(age_year) %>%
  summarise(total = sum(total_infection_year),
            min = min(total_infection_year), max = max(total_infection_year),
            mean = mean(total_infection_year), sd = sd(total_infection_year),
            median = median(total_infection_year)) %>% print(n = nrow(.))
count_infection_by_age(infection_cases_bednet) %>%
  select(total_infection_year, age_year) %>% group_by(age_year) %>%
  summarise(total = sum(total_infection_year),
            min = min(total_infection_year), max = max(total_infection_year),
            mean = mean(total_infection_year), sd = sd(total_infection_year),
            median = median(total_infection_year)) %>% print(n = nrow(.))

# Clinical cases

# All
count_clinical(clinical_cases_control) %>% select(total_clinical) %>%
  summarise(total = sum(total_clinical),
            min = min(total_clinical), max = max(total_clinical),
            mean = mean(total_clinical), sd = sd(total_clinical),
            median = median(total_clinical)) %>% t() %>% round(digits = 1)
count_clinical(clinical_cases_bednet) %>% select(total_clinical) %>%
  summarise(total = sum(total_clinical),
            min = min(total_clinical), max = max(total_clinical), 
            mean = mean(total_clinical), sd = sd(total_clinical),
            median = median(total_clinical))  %>% t() %>% round(digits = 1)

# By age
count_clinical_by_age(clinical_cases_control) %>%
  select(total_clin_year, age_year) %>% group_by(age_year) %>%
  summarise(total = sum(total_clin_year),
            min = min(total_clin_year), max = max(total_clin_year),
            mean = mean(total_clin_year), sd = sd(total_clin_year),
            median = median(total_clin_year)) %>% print(n = nrow(.))
count_clinical_by_age(clinical_cases_bednet) %>%
  select(total_clin_year, age_year) %>% group_by(age_year) %>%
  summarise(total = sum(total_clin_year),
            min = min(total_clin_year), max = max(total_clin_year),
            mean = mean(total_clin_year), sd = sd(total_clin_year),
            median = median(total_clin_year)) %>% print(n = nrow(.))

## Plot

# Result by age of infection / case

result_inf_age_control <- count_infection_by_age(infection_cases_control) %>% mutate(run = "Control")
result_inf_age_bednet <- count_infection_by_age(infection_cases_bednet) %>% mutate(run = "ITNs")

result_clin_age_control <- count_clinical_by_age(clinical_cases_control) %>% mutate(run = "Control")
result_clin_age_bednet <- count_clinical_by_age(clinical_cases_bednet) %>% mutate(run = "ITNs")

# Put together and clean

result_age_inf <- rbind(result_inf_age_control, result_inf_age_bednet)
rm(result_inf_age_control, result_inf_age_bednet)

result_age_clin <- rbind(result_clin_age_control, result_clin_age_bednet)
rm(result_clin_age_control, result_clin_age_bednet)

result_age_clin_long <- result_age_clin %>%
  select(-c(timestep_born, timestep_died, max_age, age_gr_final,
            D_by_time, T_by_time)) %>%
  pivot_longer(-c(individual_index, age_year, recieved_net, removed_net, run),
               names_to = "measure", values_to = "count") %>%
  mutate(measure = gsub("total_D_year", "D", measure)) %>%
  mutate(measure = gsub("total_T_year", "T", measure)) %>%
  mutate(measure = gsub("total_clin_year", "Clinical cases", measure))

# Plots to see infection / clinical incidence by age

png(filename = "outputs_plots/clinical_incidence_age.png",
    width = 12, height = 5, units = "in", res = 1200)

ggplot(filter(result_age_clin_long, measure == "Clinical cases"),
       aes(x = as.factor(age_year), y = count, color = run, fill = run)) +
  geom_boxplot(alpha = 0.3) +
  scale_color_manual(breaks = c("Control", "ITNs"),
                     values = carto_pal(name = "Safe")[c(2, 8)]) +
  scale_fill_manual(breaks = c("Control", "ITNs"),
                    values = carto_pal(name = "Safe")[c(2, 8)]) +
  labs(x = "Age (year)", y = "Total no. of new clinical infection (D and T) per year of age") +
  theme_bw() + theme(legend.title = element_blank(), legend.position = "bottom")

dev.off()

png(filename = "outputs_plots/infection_incidence_age.png",
    width = 12, height = 5, units = "in", res = 1200)

ggplot(filter(result_age_inf),
       aes(x = as.factor(age_year), y = total_infection_year, color = run, fill = run)) +
  geom_boxplot(alpha = 0.3) +
  scale_color_manual(breaks = c("Control", "ITNs"),
                     values = carto_pal(name = "Safe")[c(2, 8)]) +
  scale_fill_manual(breaks = c("Control", "ITNs"),
                    values = carto_pal(name = "Safe")[c(2, 8)]) +
  labs(x = "Age (year)", y = "Total no. of new infection per year of age") +
  theme_bw() + theme(legend.title = element_blank(), legend.position = "bottom")

dev.off()


#### cross-sectional ####

## Count per individual the number of infections detected with a cross sectional

source("functions/verbose_dataclean_cross_sectional.R")

# Use the already detected clinical and new infections dataset (subset to time)
# and then use the cross sectional command

# Infections

# All
infection_cases_control %>%
  cross_survey(survey_time = 25*year, period = 15) %>%
  detect_all_infection() %>%
  group_by(individual_index) %>% filter(row_number() == 1) %>% ungroup() %>%
  summarise(total = sum(infection))
infection_cases_bednet %>%
  cross_survey(survey_time = 25*year, period = 15) %>%
  detect_all_infection() %>%
  group_by(individual_index) %>% filter(row_number() == 1) %>% ungroup() %>%
  summarise(total = sum(infection))

# By age
infection_cases_control %>%
  cross_survey(survey_time = 25*year, period = 15) %>%
  detect_all_infection() %>%
  group_by(individual_index) %>% filter(row_number() == 1) %>% ungroup() %>%
  group_by(age_year) %>% summarise(total = sum(infection))
infection_cases_bednet %>%
  cross_survey(survey_time = 25*year, period = 15) %>%
  detect_all_infection() %>%
  group_by(individual_index) %>% filter(row_number() == 1) %>% ungroup() %>%
  group_by(age_year) %>% summarise(total = sum(infection))

# Clinical cases

# All
clinical_cases_control %>%
  cross_survey(survey_time = 25*year, period = 15) %>%
  detect_all_clinical() %>%
  group_by(individual_index) %>% filter(row_number() == 1) %>% ungroup() %>%
  summarise(total = sum(case))
clinical_cases_bednet %>%
  cross_survey(survey_time = 25*year, period = 15) %>%
  detect_all_clinical() %>%
  group_by(individual_index) %>% filter(row_number() == 1) %>% ungroup() %>%
  summarise(total = sum(case))

# By age

clinical_cases_control %>%
  cross_survey(survey_time = 25*year, period = 15) %>%
  count_clinical_by_age() %>%
  select(total_clin_year, age_year) %>% group_by(age_year) %>%
  summarise(total = sum(total_clin_year),
            min = min(total_clin_year), max = max(total_clin_year),
            mean = mean(total_clin_year), sd = sd(total_clin_year),
            median = median(total_clin_year))
clinical_cases_bednet %>%
  cross_survey(survey_time = 25*year, period = 15) %>%
  count_clinical_by_age() %>%
  select(total_clin_year, age_year) %>% group_by(age_year) %>%
  summarise(total = sum(total_clin_year),
            min = min(total_clin_year), max = max(total_clin_year),
            mean = mean(total_clin_year), sd = sd(total_clin_year),
            median = median(total_clin_year))


#### time to infection ####

clinical_cases_control %>%
  # Only interested in time to first clinical infection,
  # so filter to that timepoint for each individual
  filter(new_D | new_T) %>%
  arrange(individual_index, timestep) %>%
  group_by(individual_index) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  # Calculate time to infection
  mutate(time_to_infection = timestep - key_intervention_time[1]) %>%
  summarise(min = min(time_to_infection), max = max(time_to_infection),
            mean = mean(time_to_infection), sd = sd(time_to_infection),
            median = median(time_to_infection)) %>% t() %>% round(digits = 1)
clinical_cases_bednet %>%
  # Only interested in time to first clinical infection,
  # so filter to that timepoint for each individual
  filter(new_D | new_T) %>%
  arrange(individual_index, timestep) %>%
  group_by(individual_index) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  # Calculate time to infection
  mutate(time_to_infection = timestep - key_intervention_time[1]) %>%
  summarise(min = min(time_to_infection), max = max(time_to_infection),
            mean = mean(time_to_infection), sd = sd(time_to_infection),
            median = median(time_to_infection)) %>% t() %>% round(digits = 1)