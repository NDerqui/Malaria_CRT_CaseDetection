
# INPUT: Basic info, like sim_length, init_EIR, etc.
# OUTPUT: Setting background parameters for running malariasim.


# DESCRIPTION:

# This function sets all basic parameters for malariasim:
# bednet usage, treatment, seasonality, etc.
# Some default parameters have been already written, but these can be modified too.

# On future functions, possible to modify one/several parameters,
# like analysing the effect of a key intervention (better bednet or vaccine)
# in a run while keeping all basic parameters the same in a control run. 

set_baseline_pars <- function(sim_length, init_EIR, human_population,
                              ## Some basic options
                              seasonality,
                              season_g0 = 0,
                              season_g = c(1, 0, 0),
                              season_h = c(0, 0, 0),
                              ## Treatment pars
                              # Using default SP-AQ as per JC m/s
                              treatment, treat_timesteps = 0, # Treatment introduced at start (and not remove)
                              treat_coverage = 0.75,          # Treatment administered to 75%
                              ## Bednet pars (coverage, etc. default)
                              # By default in these sims, parameters constant over time
                              bednets,
                              baseline_bednets_timesteps = seq(1, sim_length, 3), # By default, bednet rounds every 3 years
                              baseline_bed_coverage = 0.5,  # Each round is distributed to 50% of the population.
                              baseline_bed_retention = 5,   # Nets are kept on average 5 years
                              baseline_bed_dn0 = 0.352,     # Death probabilities for each mosquito species 
                              baseline_bed_rn = 0.568,      # Repelling probabilities for each mosquito species 
                              baseline_bed_rnm = 0.24,      # Minimum repelling probabilities for each mosquito species
                              baseline_bed_gamman = 2.64,   # Bed net half-lives
                              ## Vaccine routine admin pars (coverage, etc. default)
                              # By default in these sims, parameters constant over time
                              vaccine,
                              baseline_vax_timesteps = 1,   # Round of vaccination is at 1st year into the simulation.
                              baseline_vax_profile = rtss_profile, # Implementation of the RTSS vaccine.
                              baseline_vax_coverage = 1,    # The vaccine is given to 100% of the population between the specified ages.
                              baseline_vax_min_wait = 0,    # The minimum acceptable time since the last vaccination is 0 because in our case we are only implementing one round of vaccination.
                              baseline_vax_min_age_month = 5,   # The minimum age for the target population to be vaccinated. 
                              baseline_vax_max_age_year = 5,    # The maximum age for the target population to be vaccinated.
                              baseline_vax_booster_spacing_month = 12, # The booster is given at 12 months after the primary series. 
                              baseline_vax_booster_coverage = 0.95,    # Coverage of the booster dose is 95%.
                              baseline_vax_booster_profile = rtss_booster_profile # We will model implementation of the RTSS booster.
                              ) {
  
  require(malariasimulation)
  
  month <- 30
  year <- 365
  
  
  ## Set basic sim parameters
  
  simparams <- get_parameters(list(
    human_population = human_population))
  
  # Set seasonality
  
  if (seasonality) {
    
    simparams$model_seasonality <- TRUE
    # Define the seasonality parameters as per malariaverse
    simparams$g0 <- 0
    simparams$g <- c(1, 0, 0)
    simparams$h <- c(0, 0, 0)
  }
  
  # Run the equilibrium after setting the seasonality
  
  simparams <- simparams  %>%
    set_equilibrium(init_EIR = init_EIR)
  
  # Set treatment pars if adding that option
  
  if (treatment) {
    
    simparams <- simparams %>%
      set_drugs(list(SP_AQ_params)) %>%
      # Initial coverage (before our first introduction of treatment) is default 0%,
      # then SP-AQ is introduced at coverage (which in this case is at time 0)
      set_clinical_treatment(drug = 1,
                             timesteps = treat_timesteps * year,
                             coverages = treat_coverage) %>% 
      set_equilibrium(init_EIR = init_EIR)
  }
  
  # Set bednets pars
  
  if (bednets) {
  
  simparams <- set_bednets(
    simparams,
    timesteps = baseline_bednets_timesteps * year,
    coverages = rep(baseline_bed_coverage, times = length(baseline_bednets_timesteps)),
    retention = baseline_bed_retention * year, 
    dn0 = matrix(rep(baseline_bed_dn0, times = length(baseline_bednets_timesteps)), nrow = length(baseline_bednets_timesteps), ncol = 1), # Matrix of death probabilities
    rn = matrix(rep(baseline_bed_rn, times = length(baseline_bednets_timesteps)), nrow = length(baseline_bednets_timesteps), ncol = 1), # Matrix of repelling probabilities 
    rnm = matrix(rep(baseline_bed_rnm, times = length(baseline_bednets_timesteps)), nrow = length(baseline_bednets_timesteps), ncol = 1), # Matrix of minimum repelling probabilities
    gamman = rep(baseline_bed_gamman * year, times = length(baseline_bednets_timesteps)) # Vector of bed net half-lives for each distribution timestep
  )
  
  }
  
  # Set vaccine pars
  
  if (vaccine) {
    
    simparams <- set_mass_pev(
      simparams,
      profile = baseline_vax_profile,
      timesteps = baseline_vax_timesteps * year,
      coverages = baseline_vax_coverage,
      min_wait = baseline_vax_min_wait,
      min_ages = baseline_vax_min_age_month * month,
      max_ages = baseline_vax_max_age_year * year,
      booster_spacing = baseline_vax_booster_spacing_month * month,
      booster_coverage = matrix(baseline_vax_booster_coverage),
      booster_profile = list(baseline_vax_booster_profile)
    )
    
  }
  
  
  ## Output is the sim pars
  
  return(simparams)

}
