
# INPUT: Basic simparams and other key trial info (sim_length, intervention, age to follow, etc).
# OUTPUT: Age-cohort cleaned data for both control and intervention arms.


# DESCRIPTION:

# This function runs verbose sims and cleans the basic output to extract state,
# then gets the analysis age-cohort to minimise what we need for each simulation.

run_and_clean_verbose <- function(run_note,
                                  ## Shared verbose sim parameters
                                  verbose_protocol,
                                  ## Key intervention parameters
                                  key_bednet, intervention_protocol = NULL,
                                  ## Parameters to follow only analysis cohort
                                  analysis_cohort_protocol) {
  
  ## Basic inputs
  
  # Packages and year/month options
  
  require(malariasimulation)
  require(dplyr)
  
  month <- 30
  year <- 365
  
  # Source our functions to run and clean the verbose simulation
  
  source("functions/verbose_simulation.R")
  source("functions/verbose_analysis_cohort.R")
  
  
  ## Run sim
  
  # Run and extract state
  
  sim_args <- c(verbose_protocol,
                list(run_note = run_note, key_bednet = key_bednet))
  
  if (!is.null(intervention_protocol)) {
    sim_args <- c(sim_args, intervention_protocol)
  }
  
  out <- do.call(run_verbose_sim, sim_args)
  
  df <- read.csv(paste0("verbose_dump/", run_note, "_full_output.csv"))
  
  df$process <- out$process_vector[df$process_index]
  df$state <- out$state_list[df$state_index]
  
  rm(out)
  
  # Read the age snapshot
  
  df_age <- read.csv(paste0("verbose_dump/", run_note, "_snapshot_age.csv"))
  
  
  ## Simple clean to subtract to the cohort we can follow with age.
  
  # Filter individuals born / with age from snapshot,
  # estimate their age at each timestep and final age (at death or sim end),
  # and sample trial size.
  
  analyses_cohort <- df %>%
    get_birth_death() %>%
    get_age_cohort(age_snapshot = df_age,
                   snapshot_time = verbose_protocol$snapshot_time)
  
  analyses_cohort <- do.call(
    get_enrol_sample,
    c(list(df = analyses_cohort), analysis_cohort_protocol)
  )
  
  
  ## Return strictly necessary
  
  return(analyses_cohort)
  
}

# DESCRIPTION:

# This function runs one verbose sim for the control arm and one for the intervention,
# over n simulations and cleaning each to extract only the analysis age cohort data.

# Finally, it appends the data of each sim and write a .csv for each control and intervention.

sim_two_arm_trial <- function(trial_slug, n_power,
                              ## Verbose sim parameters
                              verbose_protocol,
                              ## Key intervention parameters
                              intervention_protocol,
                              ## Parameters to follow only analysis cohort
                              analysis_cohort_protocol) {
  
  ## Basic inputs
  
  # Packages and year/month options
  
  require(dplyr)
  require(purrr)
  
  month <- 30
  year <- 365
  
  
  ## Run sim over n simulations
  
  # Control
  
  purrr::map_dfr(seq_len(n_power), \(i) {
    run_and_clean_verbose(
      verbose_protocol = verbose_protocol,
      run_note = "control",
      key_bednet = FALSE,
      analysis_cohort_protocol = analysis_cohort_protocol
    ) %>%
      dplyr::mutate(sim = i)
  }) %>%
    write.csv(paste0("outputs/cohort_data/", trial_slug, "_control.csv"), row.names = FALSE)
  
  # Intervention
  
  purrr::map_dfr(seq_len(n_power), \(i) {
    run_and_clean_verbose(
      verbose_protocol = verbose_protocol,
      run_note = "intervention",
      key_bednet = TRUE,
      intervention_protocol = intervention_protocol,
      analysis_cohort_protocol = analysis_cohort_protocol
    ) %>%
      dplyr::mutate(sim = i)
  }) %>%
    write.csv(paste0("outputs/cohort_data/", trial_slug, "_intervention.csv"), row.names = FALSE)
  
}