
# INPUT: Basic simparams and other key information (sim_length).
# OUTPUT: Verbose dump file.

# DESCRIPTION:

# This function runs verbose sims modulating one/several key ITN or other pars.
# It uses parameters already fixed (baseline), so that key intervention
# can be analysed on a control and an intervention run.

# Verbose parameters are set here, so modify in advance.

run_verbose_sim <- function(simparams, sim_length,
                            ## Note for verbose dump file
                            run_note = "", snapshot_time,
                            ## Key intervention
                            # Here, only set so that key intervention is bednet
                            key_bednet, key_intervention_time = NA,
                            # IMP: we need to have an options to modify the following pars,
                            # but by default are kept the same as the set_baseline_par
                            bed_coverage = 0.5,  # Each round is distributed to 50% of the population.
                            bed_retention = 5,   # Nets are kept on average 5 years
                            bed_dn0 = 0.352,     # Death probabilities for each mosquito species 
                            bed_rn = 0.568,      # Repelling probabilities for each mosquito species 
                            bed_rnm = 0.24,      # Minimum repelling probabilities for each mosquito species
                            bed_gamman = 2.64    # Bed net half-lives
                            ) {
  
  require(malariasimulation)
  
  month <- 30
  year <- 365
  
 
  ## Key intervention 
  
  # Re-set bednets pars if this is the key intervention of interest
  
  if (key_bednet) {
  
    # These are our bednets distribution rounds as per baseline par
    bednets_timesteps <- seq(0, sim_length, 3)*year
    
    # Add our key intervention timepoint
    # (if already there, no adding, if not, add)
    bednets_timesteps <- sort(unique(c(bednets_timesteps, key_intervention_time)))
    
    # Override the bednet pars with the extra (or not) timepoint
    # The pars here should match those used in the baseline function
    simparams <- set_bednets(
      simparams,
      timesteps = bednets_timesteps,
      coverages = rep(0.5, times = length(bednets_timesteps)),
      retention = 5 * year, 
      dn0 = matrix(rep(0.352, times = length(bednets_timesteps)), nrow = length(bednets_timesteps), ncol = 1), # Matrix of death probabilities
      rn = matrix(rep(0.568, times = length(bednets_timesteps)), nrow = length(bednets_timesteps), ncol = 1), # Matrix of repelling probabilities 
      rnm = matrix(rep(0.24, times = length(bednets_timesteps)), nrow = length(bednets_timesteps), ncol = 1), # Matrix of minimum repelling probabilities
      gamman = rep(2.64 * year, times = length(bednets_timesteps)) # Vector of bed net half-lives for each distribution timestep
    )
    
    # Get the index for that intervention timepoint(s)
    index_key_intervention <- which(bednets_timesteps %in% key_intervention_time)
    
    # And subscribe that parameter in that timepoint(s)
    simparams[["bednet_coverages"]][index_key_intervention] <- bed_coverage
    simparams[["bednet_retention"]] <- bed_retention * year
    simparams[["bednet_dn0"]][index_key_intervention] <- bed_dn0
    simparams[["bednet_rn"]][index_key_intervention] <- bed_rn
    simparams[["bednet_rnm"]][index_key_intervention] <- bed_rnm
    simparams[["bednet_gamman"]][index_key_intervention] <- bed_gamman * year

    }
  
  
  ## Verbose sims options
  
  # Keep as basic:
  # i.e. retrieve all states, process and ITNs pars for people <80 years.
  
  simparams$progress_bar <- TRUE
  
  simparams$infection_verbose <- FALSE
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
  
  folder <- "outputs_verbose_sims"
  
  dir.create(paste0(folder, "/"), showWarnings = FALSE)
  
  # Set the verbose file name 
  
  simparams$file_name <- paste0(folder, "/verbose_dumping_", run_note, ".csv")
  simparams$snapshot_file_name <- paste0(folder, "/verbose_dumping_snapshot_", run_note, ".csv")
  
  
  ## Run simulation
  
  output <- malariasimulation:::run_verbose_simulation(timesteps = sim_length*year,
                                                       parameters = simparams)
}
