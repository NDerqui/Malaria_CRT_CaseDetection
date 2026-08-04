
# INPUT: Basic simparams and other key trial info (sim_length, intervention, age to follow, etc).
# OUTPUT: Age-cohort cleaned data for both control and intervention arms.


# DESCRIPTION:

# This function runs verbose sims and cleans the basic output to extract state,
# then gets the analysis age-cohort to minimise what we need for each simulation.
# Finally, detect infection states to run all verbose-functions in one step.

run_and_clean_verbose <- function(run_note,
                                  ## Shared verbose sim parameters
                                  verbose_protocol,
                                  ## Key intervention parameters
                                  key_intervention, intervention_protocol = NULL,
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
  source("functions/verbose_detect_event.R")
  
  
  ## Run sim
  
  # Run and extract state
  
  sim_args <- c(verbose_protocol,
                list(run_note = run_note, key_intervention = key_intervention))
  
  if (!is.null(intervention_protocol)) {
    sim_args <- c(sim_args, intervention_protocol)
  }
  
  out <- do.call(run_verbose_sim, sim_args)
  
  # Option with local .csv
  # df <- read.csv(paste0("verbose_dump/", run_note, "_full_output.csv"))
  
  # Option without csv dependency
  df <- out$verbose_data
  
  df$process <- out$process_vector[df$process_index]
  df$state <- out$state_list[df$state_index]
  
  # Some final clean
  df <- df %>%
    select(-process_index, -state_index) %>%
    arrange(individual_index, timestep)
  
  rm(out)
  
  # Read the age snapshot
  
  df_age <- read.csv(paste0("verbose_dump/", run_note, "_snapshot_age.csv"))
  
  
  ## Simple clean to subtract to the cohort with age we want to follow.
  
  # Filter individuals born / with age from snapshot,
  # estimate their age at each timestep and final age (at death or sim end),
  # and sample trial size.
  
  analyses_cohort <- df %>%
    get_birth_death() %>%
    get_age_cohort(age_snapshot = df_age,
                   snapshot_time = verbose_protocol$snapshot_time)
  
  analyses_cohort <- do.call(
    get_enrol_samples,
    c(list(df = analyses_cohort), analysis_cohort_protocol)
  )
  
  
  ## Detect the infection events
  
  # Run all the verbose-related functions in the early stages,
  # so future analyses are only for trial outcomes and effects.
  
  analyses_cohort <- analyses_cohort %>%
    detect_ever_malaria() %>%
    detect_infection()
  
  
  ## Return strictly necessary
  
  return(analyses_cohort)
  
}

# DESCRIPTION:

# This function runs one verbose sim for the control arm and one for the intervention,
# over n simulations and cleaning each to extract only the analysis age cohort data.

# Finally, it appends the data of each sim and write a .csv for each control and intervention.

sim_two_arm_trial <- function(trial_id, n_power, n_clusters,
                              ## Verbose sim parameters
                              verbose_protocol,
                              ## Key intervention parameters
                              key_intervention, intervention_protocol,
                              ## Parameters to follow only analysis cohort
                              analysis_cohort_protocol) {
  
  ## Basic inputs
  
  # Packages and year/month options
  
  require(dplyr)
  require(purrr)
  
  month <- 30
  year <- 365
  
  
  # 1. Define the two arms
  
  arms <- c("Control", "Intervention")
  
  # 2. Run every sim and cluster within each arm
  
  # Simple, sequential, lapply
  # arm_results <- lapply(arms, function(run) {
  
  # Local, parallel, lapply
  arm_results <- future.apply::future_lapply(arms, function(run) {
    
    intervention_arm <- run == "Intervention"
    
    purrr::map_dfr(seq_len(n_power), function(sim) {
      
      purrr::map_dfr(seq_len(n_clusters), function(cluster_within_arm) {
        
        # This assigns _id to be unique across both arms,
        # i.e. sequential starting with control then intervention.
        cluster_id <- if (intervention_arm) {
          n_clusters + cluster_within_arm
        } else {
          cluster_within_arm
        }
        
        run_and_clean_verbose(
          # Same verbose protocol for both arms
          run_note = paste(tolower(run), sim, cluster_id, sep = "_"),
          verbose_protocol = verbose_protocol,
          # Assign intervention parameters only to the intervention arm
          key_intervention = if (intervention_arm) {
            key_intervention
          } else {
            "none"
          },
          intervention_protocol = if (intervention_arm) {
            intervention_protocol
          } else {
            NULL
          },
          # Same analysis cohort protocol for both arms
          analysis_cohort_protocol = analysis_cohort_protocol
        ) %>%
          mutate(
            sim = sim,
            cluster_id = cluster_id,
            run = run
          )
      })
    })
  },
  future.seed = TRUE)
  
  # 3. Combine both arms
  
  cohort_data <- bind_rows(arm_results)
  
  # 4. Save one cohort file
  
  write.csv( cohort_data,
    paste0("outputs/cohort_data/", trial_id, ".csv"),
    row.names = FALSE)

  # 5. Clean up the verbose_dump folder to save space
  
  unlink("verbose_dump/*")

}