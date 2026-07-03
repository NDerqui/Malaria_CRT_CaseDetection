#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#-#-#-#-#-#-#-#-# Modelling Malaria CRT - Case Detection #-#-#-#-#-#-#-#-#-#-#-
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#


## Modelling case detection in malaria CRTs using verbose simulations:
## set up parameters, run sim, clean and save analyses cohort.



# SET-UP ------------------------------------------------------------------


rm(list = ls())

# Packages

# remotes::install_github("mrc-ide/individual@dev")
# remotes::install_github("JasonRWood/malariasimulation@verbose_simulations")
library(individual)
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

init_EIR <- 25

human_population <- 10000
trial_size <- 500

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

trial_name <- paste0("Seasonal Init EIR ", init_EIR)


#### run sim ####

## Functions to set up parameters and run the verbose simulation

source("functions/verbose_set_parameters.R")
source("functions/verbose_simulation.R")

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
                       bed_coverage = 0.95, run_note = "intervention")

df_intervention <- read.csv("outputs_verbose_sims/verbose_dumping_intervention.csv")

df_intervention$process <- out$process_vector[df_intervention$process_index]
df_intervention$state <- out$state_list[df_intervention$state_index]

write.csv(df_intervention, "outputs_verbose_sims/verbose_dumping_intervention.csv", row.names = FALSE)

## Read the verbose files only

rm(out)

df_control <- read.csv("outputs_verbose_sims/verbose_dumping_control.csv")
df_control_age <- read.csv("outputs_verbose_sims/verbose_dumping_snapshot_control.csv")

df_intervention <- read.csv("outputs_verbose_sims/verbose_dumping_intervention.csv")
df_intervention_age <- read.csv("outputs_verbose_sims/verbose_dumping_snapshot_intervention.csv")

gc()



# CLEAN SIM OUTPUT --------------------------------------------------------


#### analyses cohort ####

## Simple clean to subtract to the cohort we can follow with age.

source("functions/verbose_analysis_cohort.R")

# Filter individuals born / with age from snapshot,
# estimate their age at each timestep and final age (at death or sim end),
# and sample (say 100).

analyses_cohort_control <- df_control %>%
  get_birth_death() %>%
  get_age_cohort(age_snapshot = df_control_age, snapshot_time = snapshot_time) %>%
  get_enrol_sample(alive_by = trial_start * year, trial_size = trial_size,
                   age_min = 0, age_max = 10)

analyses_cohort_intervention <- df_intervention %>%
  get_birth_death() %>%
  get_age_cohort(age_snapshot = df_intervention_age, snapshot_time = snapshot_time) %>%
  get_enrol_sample(alive_by = trial_start * year, trial_size = trial_size,
                   age_min = 0, age_max = 10)

# Clean space

rm(df_control, df_intervention, df_control_age, df_intervention_age)
gc()

# Save for future

dir.create("outputs_agecohort_data", showWarnings = FALSE)
write.csv(analyses_cohort_control, row.names = FALSE,
          file = paste0("outputs_agecohort_data/", gsub(" ", "_", tolower(trial_name)), "_control.csv"))
write.csv(analyses_cohort_intervention, row.names = FALSE,
          file = paste0("outputs_agecohort_data/", gsub(" ", "_", tolower(trial_name)), "_intervention.csv"))


#### vis ####

source("functions/verbose_visualisation.R")

# Plot the control (no intervention)

png(filename = paste0("outputs_plots/agecohort_overtime_control_", gsub(" ", "_", tolower(trial_name)), ".png"),
    width = 8, height = 5, units = "in", res = 1200)
plot_verbose_itn(df = analyses_cohort_control,
                 note = paste0("Control: ", trial_name), sim_length = sim_length,
                 human_population = human_population, trial_size = trial_size,
                 bednetstimesteps = seq(0, sim_length, 3)*year)
dev.off()

# Plot the intervention

png(filename = paste0("outputs_plots/agecohort_overtime_intervention_", gsub(" ", "_", tolower(trial_name)), ".png"),
    width = 8, height = 5, units = "in", res = 1200)
plot_verbose_itn(df = analyses_cohort_intervention,
                 note = paste0("Intervention: ", trial_name), sim_length = sim_length,
                 human_population = human_population, trial_size = trial_size,
                 bednetstimesteps = seq(0, sim_length, 3)*year) +
  geom_vline(xintercept = key_intervention_time*year, color = "firebrick", linetype = "dashed")
dev.off()