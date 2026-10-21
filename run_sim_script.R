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


# Seed for reproducible simulations

simulation_seed <- 1234


#### set trial conditions ####

# General inputs: init_EIR and general baseline conditions (shared across control and intervention).

init_EIR <- 25

boolean_seasonality <- TRUE
boolean_treatment <- TRUE
boolean_mass_clearance <- TRUE
boolean_bednets <- TRUE
boolean_vaccine <- FALSE

# Set length of sim

sim_length <- 5

# Set the intervention and when do our interventions start.
# We can do two interventions, timing of second calculated with respect to first.

key_intervention <- "bednet"

trial_start <- 3
trial_second_intervention <- NULL

key_intervention_time <- c(trial_start, trial_start+trial_second_intervention)

# Set a time for mass clearance of infection before (2 weeks) the start of the trial
# and the target population of such mass clearance (particular age group?)

mass_clearance_timestep <- (trial_start * year) - 14
mass_clearance_min_age <- 0
mass_clearance_min_age <- 15

# Number of clusters

n_clusters <- 20

# Human population size (per cluster)

human_population <- 1000

# Our follow-up cohort size (no individuals followed from total, per cluster)

trial_size <- 30

# Control when we get our age snapshot (best at start of trial, timestep = 1)

snapshot_time <- 1

# Number of simulations over which we repeat analysis

n_power <- 10

# Add a "trial name" to keep track of results

trial_name <- "Test of new treatment"


#### sim pars ####

## Functions to set up parameters

source("functions/verbose_set_parameters.R")

## Basic parameters

# With this function, we set baseline parameters
# (including seasonality, treatment use, base bednet use, etc.)
# to be used in intervention and control runs.
# Also, include as mass drug administration at the start of the trial to clear infection if desired.

baseline_parameters <- set_baseline_pars(sim_length = sim_length,
                                         init_EIR = init_EIR,
                                         human_population = human_population,
                                         seasonality = boolean_seasonality,
                                         treatment = boolean_treatment,
                                         mass_clearance = boolean_mass_clearance,
                                         mass_clearance_timestep = mass_clearance_timestep,
                                         mass_clearance_min_age = mass_clearance_min_age,
                                         mass_clearance_max_age = mass_clearance_max_age,
                                         bednets = boolean_bednets,
                                         baseline_bed_coverage = 0.6,
                                         baseline_bed_retention = 4,
                                         baseline_bed_dn0 = 0.2218,
                                         baseline_bed_rn = 0.7208,
                                         baseline_bed_rnm = 0.24,
                                         baseline_bed_gamman = (2.05/log(2)),
                                         vaccine = boolean_vaccine)

## Create protocols for our simulation and intervention

# Put the options for the verbose run,
# the intervention and the analysis cohort in lists to pass to the wrapper function.
# Note the intervention protocol is only passed to the intervention arm,
# so any modification to default pars for baseline conditions need to be passed in verbose_protocol.
# This includes the population(s) age boundary(ies) for the follow-up cohort.

verbose_protocol <- list(
  simparams = baseline_parameters,
  sim_length = sim_length,
  snapshot_time = snapshot_time,
  # Add to both arms the same baseline bednet coverage and retention parameters
  baseline_bed_coverage = 0.6,
  baseline_bed_retention = 4,
  baseline_bed_dn0 = 0.2218,
  baseline_bed_rn = 0.7208,
  baseline_bed_rnm = 0.24,
  baseline_bed_gamman = (2.05/log(2)))

intervention_protocol <- list(
  key_intervention_time = key_intervention_time,
  intervention_bed_coverage = 0.8,  
  intervention_bed_dn0 = 0.4616,     
  intervention_bed_rn = 0.5319)

analysis_populations <- tibble::tribble(
  ~analysis_population, ~age_min, ~age_max,  ~all_cleared,
  "age_0_5",                    0,        5,         TRUE,  
  "age_5_15",                   5,       15,         TRUE,
  "age_15_plus",               15,      100,        FALSE,
  "all_ages",                   0,      100,        FALSE
)

analysis_cohort_protocol <- list(
  alive_by = trial_start * year,
  trial_size = trial_size,
  analysis_populations = analysis_populations
)


#### metadata ####

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
    vaccine = boolean_vaccine,
    human_population = human_population,
    sim_length = sim_length,
    snapshot_time = snapshot_time,
    simulation_seed = simulation_seed
  ),
  
  trial = list(
    trial_start = trial_start,
    trial_second_intervention = trial_second_intervention,
    n_clusters = n_clusters,
    allocation = "1:1"
  ),
  intervention = c(
    list(key_intervention = key_intervention),
    intervention_protocol
  ),
  
  analysis_cohort = analysis_cohort_protocol,
  n_power = n_power
)

# Create our trial id and slug to use in file names and directories

trial_id <- metadata$trial_id
trial_slug <- make_trial_slug(trial_name = trial_name)

# Create output file structure

make_output_dirs()


#### run sim ####

# Configure local parallel execution

future::plan(future::multisession, workers = 2)

set.seed(simulation_seed)

# Running with the wrapper function(s)

source("functions/verbose_tidy_outputs.R")

tryCatch(
  {sim_two_arm_trial(
      trial_id = trial_id,
      n_power = n_power,
      n_clusters = n_clusters,
      verbose_protocol = verbose_protocol,
      key_intervention = key_intervention,
      intervention_protocol = intervention_protocol,
      analysis_cohort_protocol = analysis_cohort_protocol)
  },
  finally = {future::plan(future::sequential)}
)

# Only save the metadata after a successful sim runs

save_trial_metadata(metadata)



# CHECKS ------------------------------------------------------------------

# Not necessary for the workflow, but useful to check that the sim ran as expected.


#### vis ####

source("functions/verbose_visualisation.R")

# Plot the control (no intervention)

png(filename = paste0("outputs/plots/cohort/agecohort_", trial_slug, "_control.png"),
    width = 12, height = 10, units = "in", res = 1200)
read.csv(paste0("outputs/cohort_data/", trial_id, ".csv")) %>%
  filter(run == "Control") %>%
  filter(sim == sample(x = unique(sim), size = 1)) %>%
  filter(cluster_id == sample(x = unique(cluster_id), size = 5)) %>%
  plot_verbose_itn(note = paste0("Control: ", trial_name), sim_length = sim_length,
                   human_population = human_population, trial_size = trial_size,
                   bednetstimesteps = seq(0, sim_length, 3)*year) +
  facet_wrap(analysis_population ~ cluster_id, scales = "free")
dev.off()

# Plot the intervention

png(filename = paste0("outputs/plots/cohort/agecohort_", trial_slug, "_intervention.png"),
    width = 12, height = 10, units = "in", res = 1200)
read.csv(paste0("outputs/cohort_data/", trial_id, ".csv")) %>%
  filter(run == "Intervention") %>%
  filter(sim == sample(x = unique(sim), size = 1)) %>%
  filter(cluster_id == sample(x = unique(cluster_id), size = 5)) %>%
  plot_verbose_itn(note = paste0("Intervention: ", trial_name), sim_length = sim_length,
                 human_population = human_population, trial_size = trial_size,
                 bednetstimesteps = seq(0, sim_length, 3)*year) +
  geom_vline(xintercept = key_intervention_time*year, color = "firebrick", linetype = "dashed") +
  facet_wrap(analysis_population ~ cluster_id, scales = "free")
dev.off()