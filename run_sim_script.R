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


#### set trial conditions ####

# General inputs: init_EIR and total sim size (different to pop w/ analyses performed).
# Then, general baseline conditions (shared across control and intervention).

init_EIR <- 25

human_population <- 10000

boolean_seasonality <- TRUE
boolean_treatment <- TRUE
boolean_bednets <- TRUE

# Set length of sim and when do our interventions start.
# We can do two interventions, timing of second calculated with respect to first.

sim_length <- 6

trial_start <- 2
trial_second_intervention <- 2

key_intervention_time <- c(trial_start, trial_start+trial_second_intervention)

# Our follow-up cohort size (no individuals followed from total)

trial_size <- 200

# Control when we get our age snapshot (best at start of trial, timestep = 1)

snapshot_time <- 1

# Number of simulations over which we repeat analysis

n_power <- 10

# Add a "trial name" to keep track of results

trial_name <- "High transmission test for workflow"


#### sim pars and metadata ####

## Functions to set up parameters

source("functions/verbose_set_parameters.R")

## Basic parameters

# With this function, we set baseline parameters
# (including seasonality, treatment use, base bednet use, etc.)
# to be used in intervention and control runs.

baseline_parameters <- set_baseline_pars(sim_length = sim_length,
                                         init_EIR = init_EIR,
                                         human_population = human_population,
                                         seasonality = boolean_seasonality,
                                         treatment = boolean_treatment,
                                         bednets = boolean_bednets)

## Create protocols for our simulation and intervention

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

## Functions to save the protocols and all trial metadata

source("functions/trial_metadata.R")
source("functions/trial_tidy_outputs.R")

## Trial metadata

# With this function, we save all the above information so we can retrieve after the simulation is run.

metadata <- create_trial_metadata(
  
  trial_name = trial_name,
  
  simulation = list(
    init_EIR = init_EIR,
    seasonality = boolean_seasonality,
    treatment = boolean_treatment,
    bednets = boolean_bednets,
    human_population = human_population,
    sim_length = sim_length,
    snapshot_time = snapshot_time
  ),
  
  trial = list(
    trial_start = trial_start,
    trial_second_intervention = trial_second_intervention
  ),
  intervention = intervention_protocol,
  
  analysis_cohort = analysis_cohort_protocol,
  n_power = n_power
)

# Create our trial id and slug to use in file names and directories

trial_id <- metadata$trial_id
trial_slug <- make_trial_slug(trial_name = trial_name)

make_output_dirs()
save_trial_metadata(metadata)


#### run sim ####

# Running with the wrapper function(s)

source("functions/verbose_tidy_outputs.R")

sim_two_arm_trial(trial_id = trial_id,
                  n_power = n_power,
                  verbose_protocol = verbose_protocol,
                  intervention_protocol = intervention_protocol,
                  analysis_cohort_protocol = analysis_cohort_protocol)


#### vis ####

source("functions/verbose_visualisation.R")

# Plot the control (no intervention)

png(filename = paste0("outputs/plots/cohort/agecohort_", trial_slug, "_control.png"),
    width = 8, height = 5, units = "in", res = 1200)
read.csv(paste0("outputs/cohort_data/", trial_id, "_control.csv")) %>%
  plot_verbose_itn(note = paste0("Control: ", trial_name), sim_length = sim_length,
                   human_population = human_population, trial_size = trial_size,
                   bednetstimesteps = seq(0, sim_length, 3)*year)
dev.off()

# Plot the intervention

png(filename = paste0("outputs/plots/cohort/agecohort_", trial_slug, "_intervention.png"),
    width = 8, height = 5, units = "in", res = 1200)
read.csv(paste0("outputs/cohort_data/", trial_id, "_intervention.csv")) %>%
plot_verbose_itn(note = paste0("Intervention: ", trial_name), sim_length = sim_length,
                 human_population = human_population, trial_size = trial_size,
                 bednetstimesteps = seq(0, sim_length, 3)*year) +
  geom_vline(xintercept = key_intervention_time*year, color = "firebrick", linetype = "dashed")
dev.off()