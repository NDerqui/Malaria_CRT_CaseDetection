
# INPUT: a df with new infections or cases flagged (still one obs per time per person).
# OUTPUT: a df with one observation per timestep, with information on prevcalence and incidence.


# DESCRIPTION:

# Functions to get:
# infection prevalence at each timestep (infected at timestep/population),
# infection incidence at each timestep (newly infected at timestep/population at risk),
# case prevalence at each timestep (case at timestep/population),
# case incidence at each timestep (newly clinical case at timestep/population at risk).

get_prev_inc <- function(df) {
  
  require(dplyr)
  
  # Prepare to iter over all timesteps
  
  timesteps <- unique(df$timestep)
  result <- data.frame()
  
  # Loop over all timesteps to get prevalence and incidence at each timstep
  
  for (time in 1:length(timesteps)) {
    
    #new_infection_at_time
    #new_case_at_time
    #infected_at_time
    #case_at_time
    
    add <- df %>%
      ## Get no. at risk in this timestep considering prev timstep:
      # To later get no. at risk, get the state at our prior timestep.
      # No at risk of incident infection won't include infected in previous timestep.
      filter(timestep %in% (timesteps[time] - 1):timesteps[time]) %>% # No need to get lag for all the df
      group_by(individual_index) %>%
      mutate(prev_state = lag(state)) %>%
      ungroup() %>%
      ## Ready for the each timestep calcs
      # Filter to each timestep and remove everyone dead by then (for denom)
      filter(timestep == timesteps[time]) %>%
      filter(!is.na(timestep_died) | timestep_died > timesteps[time]) %>%
      # Get some basic counts
      mutate(n = n()) %>%
      mutate(at_risk = sum(!(prev_state %in% c("U", "A", "D", "Tr")))) %>%
      mutate(infections = sum(infected_at_time)) %>%
      mutate(cases = sum(case_at_time)) %>%
      mutate(new_infections = sum(new_infection_at_time)) %>%
      mutate(new_cases = sum(new_case_at_time)) %>%
      # Incidence and Prevalence calcs
      mutate(prevalence_infec = infections/n) %>%
      mutate(prevalence_case = cases/n) %>%
      mutate(incidence_infec = new_infections/at_risk) %>%
      mutate(incidence_case = new_cases/at_risk) %>%
      # Cleaning
      group_by(timestep) %>% filter(row_number() == 1) %>% ungroup() %>%
      select(timestep, n, at_risk,
             infections, cases, new_infections, new_cases,
             prevalence_infec, prevalence_case, incidence_infec, incidence_case)
    
    result <- rbind(result, add)
    rm(add)
  }
  
  return(result)
}

