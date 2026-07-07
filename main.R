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

human_population <- 100
trial_size <- 10

# Set some basics, like length of sim vs length of trial,
# and when do our interventions start.

# Example, run the sim for total 9 years, with 3 years "without intervention"
# (i.e., start the trial on the third year of the sim),
# and doing two rounds of intervention separated by 3 years.

sim_length <- 6

trial_start <- 2
trial_second_intervention <- 2

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

## Run sim and clean the data

# Put the options for the verbose run, the intervention and the analysis cohort in lists to pass to the wrapper function.

verbose_protocol <- list(
  simparams = baseline_parameters,
  sim_length = sim_length,
  snapshot_time = snapshot_time)

intervention_protocol <- list(
  key_intervention_time = key_intervention_time,
  bed_coverage = 0.95)

analysis_cohort_protocol <- list(
  alive_by = trial_start*year,
  trial_size = trial_size,
  age_min = 0, age_max = 10)

# Running with the wrapper function(s)

source("functions/trial_tidy_outputs.R")
source("functions/verbose_tidy_outputs.R")

trial_slug <- make_trial_slug(trial_name = trial_name)
make_output_dirs()

sim_two_arm_trial(trial_slug = trial_slug,
                  n_power = 10,
                  verbose_protocol = verbose_protocol,
                  intervention_protocol = intervention_protocol,
                  analysis_cohort_protocol = analysis_cohort_protocol)


#### vis ####

source("functions/verbose_visualisation.R")

# Plot the control (no intervention)

png(filename = paste0("outputs/plots/cohort/agecohort_", trial_slug, "_control.png"),
    width = 8, height = 5, units = "in", res = 1200)
read.csv(paste0("outputs/cohort_data/", trial_slug, "_control.csv")) %>%
  plot_verbose_itn(note = paste0("Control: ", trial_name), sim_length = sim_length,
                   human_population = human_population, trial_size = trial_size,
                   bednetstimesteps = seq(0, sim_length, 3)*year)
dev.off()

# Plot the intervention

png(filename = paste0("outputs/plots/cohort/agecohort_", trial_slug, "_intervention.png"),
    width = 8, height = 5, units = "in", res = 1200)
read.csv(paste0("outputs/cohort_data/", trial_slug, "_intervention.csv")) %>%
plot_verbose_itn(note = paste0("Intervention: ", trial_name), sim_length = sim_length,
                 human_population = human_population, trial_size = trial_size,
                 bednetstimesteps = seq(0, sim_length, 3)*year) +
  geom_vline(xintercept = key_intervention_time*year, color = "firebrick", linetype = "dashed")
dev.off()