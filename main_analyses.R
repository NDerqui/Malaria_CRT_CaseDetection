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

# Functions

source("functions/trial_tidy_outputs.R")
source("functions/trial_analysis.R")
source("functions/trial_visualisation.R")

# General options

year <- 365
month <- 30



# TRIAL CONDITIONS --------------------------------------------------------

# General inputs: init_EIR, total sim size, and pop followed in trial.

init_EIR <- 25

human_population <- 10000
trial_size <- 200

# Length of sim, length of trial, time of start interventions.

sim_length <- 6

trial_start <- 2
trial_second_intervention <- 2

key_intervention_time <- c(trial_start, trial_start + trial_second_intervention)

# Trial name of the simulation we want to analyse

trial_name <- paste0("Seasonal Init EIR ", init_EIR)
trial_slug <- make_trial_slug(trial_name)
trial_title <- paste0("Simulated a ", human_population,
                      " population, Sampled ", trial_size,
                      " for trial, ", trial_name)

make_output_dirs()



# TRIAL OUTCOMES ----------------------------------------------------------

# Define a PCD cross-sectional survey and an ACD routine visit protocol

survey_protocol <- list(
  cross_surveys_in_years = seq(0.5, 6, 0.5)
)

acd_protocol <- list(
  routine_visits_in_weeks = seq(4, 6 * 52, 4),
  days_catchment = 2
)

# Run function for all estimates just from the trial name and protocols defined above

trial_results <- analyse_two_arm_trial(
  trial_slug = trial_slug,
  trial_start = trial_start,
  trial_second_intervention = trial_second_intervention,
  sim_length = sim_length,
  survey_protocol = survey_protocol,
  acd_protocol = acd_protocol
)



# SAVE OUTPUTS ------------------------------------------------------------

# Save all estimates and effect sizes to corresponding folder/file.

save_two_arm_trial(
  trial_results = trial_results,
  trial_slug = trial_slug
)

# Save all plots to corresponding folder/file.

save_two_arm_trial_plots(
  trial_results = trial_results,
  trial_slug = trial_slug,
  key_intervention_time = key_intervention_time,
  sim_length = sim_length,
  trial_title = trial_title
)
