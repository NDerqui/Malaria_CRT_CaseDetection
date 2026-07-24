
# INPUT: Basic simparams and other key information (sim_length).
# OUTPUT: Verbose dump file.


# DESCRIPTION:

# This function runs verbose sims modulating one/several key ITN or other pars.
# It uses parameters already fixed (baseline), so that key intervention
# can be analysed on a control and an intervention run.

# Verbose parameters are set here, so modify in advance.

run_verbose_sim <- function(simparams, sim_length,
                            ## Note for verbose dump file,
                            # and timing for the age snapshot (time at which we get age for all alive)
                            run_note = "", snapshot_time,
                            # IMP: if we modified any of the default following pars,
                            # we need to set them here too, so that the verbose run is consistent with the baseline run.
                            baseline_bednets_timesteps = seq(0, sim_length, 3), # By default, bednet rounds every 3 years
                            baseline_bed_coverage = 0.5,  # Each round is distributed to 50% of the population.
                            baseline_bed_retention = 5,   # Nets are kept on average 5 years
                            baseline_bed_dn0 = 0.352,     # Death probabilities for each mosquito species 
                            baseline_bed_rn = 0.568,      # Repelling probabilities for each mosquito species 
                            baseline_bed_rnm = 0.24,      # Minimum repelling probabilities for each mosquito species
                            baseline_bed_gamman = 2.64,   # Bed net half-lives
                            ## Key intervention
                            # Here, we only have one intervention (vaccine or bednet)
                            key_intervention, key_intervention_time = NA,
                            # IMP: we need to have an options to modify the following pars,
                            # but by default are kept the same as the set_baseline_par
                            intervention_bed_coverage = baseline_bed_coverage,  
                            intervention_bed_retention = baseline_bed_retention,   
                            intervention_bed_dn0 = baseline_bed_dn0,     
                            intervention_bed_rn = baseline_bed_rn,      
                            intervention_bed_rnm = baseline_bed_rnm,      
                            intervention_bed_gamman = baseline_bed_gamman    
                            ) {
  
  require(malariasimulation)
  
  month <- 30
  year <- 365
  
 
  ## Key intervention
  
  if (key_intervention == "bednet") {
    key_bednet <- TRUE
    key_vaccine <- FALSE
  } else {
    if (key_intervention == "vaccine") {
      key_bednet <- FALSE
      key_vaccine <- TRUE
    } else {
      stop("Key intervention must be either 'bednet' or 'vaccine'.")
    }
  }
  
  # Re-set bednets pars if this is the key intervention of interest
  
  if (key_bednet) {
  
    # Add our key intervention timepoint
    # (if already there, no adding; if not, add)
    bednets_timesteps <- sort(unique(c(baseline_bednets_timesteps, key_intervention_time)))
    
    # Override the bednet pars with the extra (or not) timepoint
    # The pars here should match those used in the baseline function
    simparams <- set_bednets(
      simparams,
      timesteps = bednets_timesteps * year,
      coverages = rep(baseline_bed_coverage, times = length(bednets_timesteps)),
      retention = baseline_bed_retention * year, 
      dn0 = matrix(rep(baseline_bed_dn0, times = length(bednets_timesteps)), nrow = length(bednets_timesteps), ncol = 1), # Matrix of death probabilities
      rn = matrix(rep(baseline_bed_rn, times = length(bednets_timesteps)), nrow = length(bednets_timesteps), ncol = 1), # Matrix of repelling probabilities 
      rnm = matrix(rep(baseline_bed_rnm, times = length(bednets_timesteps)), nrow = length(bednets_timesteps), ncol = 1), # Matrix of minimum repelling probabilities
      gamman = rep(baseline_bed_gamman * year, times = length(bednets_timesteps)) # Vector of bed net half-lives for each distribution timestep
    )
    
    # Get the index for that intervention timepoint(s)
    index_key_intervention <- which(bednets_timesteps %in% key_intervention_time)
    
    # And subscribe that parameter in that timepoint(s)
    simparams[["bednet_coverages"]][index_key_intervention] <- intervention_bed_coverage
    simparams[["bednet_retention"]] <- intervention_bed_retention * year
    simparams[["bednet_dn0"]][index_key_intervention] <- intervention_bed_dn0
    simparams[["bednet_rn"]][index_key_intervention] <- intervention_bed_rn
    simparams[["bednet_rnm"]][index_key_intervention] <- intervention_bed_rnm
    simparams[["bednet_gamman"]][index_key_intervention] <- intervention_bed_gamman * year

    }
  
  
  ## Verbose sims options
  
  # Keep as basic:
  # i.e. retrieve all states, process and ITNs pars for people <80 years.
  
  simparams$progress_bar <- TRUE
  
  simparams$infection_verbose <- TRUE
  simparams$biting_verbose <- FALSE
  simparams$mortality_verbose <- TRUE
  simparams$progression_verbose <- FALSE
  simparams$spraying_verbose <- FALSE
  simparams$nets_verbose <- FALSE
  simparams$pev_verbose <- FALSE
  simparams$states_verbose <- TRUE
  simparams$snapshot_verbose <- TRUE
  simparams$snapshot_times <- snapshot_time
  simparams$start_time <- 0
  simparams$lower_age_bound <- 0
  simparams$upper_age_bound <- 1000*year
  simparams$state_recording_freq <- 1
  
  # Set a directory to dump the verbose file
  
  folder <- "verbose_dump"
  
  dir.create(paste0(folder, "/"), showWarnings = FALSE)
  
  # Set the verbose file name 
  
  simparams$file_name <- paste0(folder, "/", run_note, "_full_output.csv")
  simparams$snapshot_file_name <- paste0(folder, "/", run_note, "_snapshot_age.csv")
  
  
  ## Run simulation
  
  output <- malariasimulation:::run_verbose_simulation(timesteps = sim_length * year,
                                                       parameters = simparams)
  
  return(output)
}
