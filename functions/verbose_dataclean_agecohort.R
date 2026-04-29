
# INPUT: a df result of a verbose simulation.
# OUTPUT: a df with verbose output with age calcs.

# DESCRIPTION:

# From the total output of a verbose sim, function that gets age of each individual.
# Age is estimated as timestep - timestep_born.
# The function notes age at each timestep and age at final step.
# Note the final step can be death or end of sim.

get_ages <- function(df, age_snapshot) {
  
  require(dplyr)
  
  # First, signal if individuals are born during the sim
  
  df <- df %>%
    # Flag if individuals were born during the sim
    group_by(individual_index) %>%
    mutate(ever_born = case_when("born" %in% process ~ 1,
                                 !("born" %in% process) ~ 0)) %>%
    ungroup()
  
  # Check if we can have age for everyone!
  # (Either by being born or being in the age_snapshot df)
  
  those_born <- unique(df$individual_index[df$ever_born == 1])
  those_snapshot <- unique(age_snapshot$individual_index)
  
  if (!(FALSE %in% (unique(df$individual_index) %in% c(those_born, those_snapshot)))) {
    
    print("Great! We can get age for everyone!")
    
    # Add age at start to our dataset (ensure timstep is 1)
    
    age_snapshot <- age_snapshot %>%
      filter(timestep == 1) %>%
      select(individual_index, age) %>%
      rename(age_at_start = age)
    
    age_cohort <- merge(df, age_snapshot, all = TRUE)
    
    # To calculate age, get timestep at which individuals are born
    # (Which is <0 for those who are not born in the sim)
    
    age_cohort <- age_cohort %>%
      # For those that were not born, timestep born is -age_at_start
      mutate(timestep_born = case_when(ever_born == 0 ~ (-age_at_start),
                                       ever_born == 1  ~ NA)) %>%
      # For those born, get the timestep for that process
      mutate(timestep_born = case_when(ever_born == 1 & process == "born" ~ timestep,
                                       ever_born == 0 | process != "born" ~ timestep_born)) %>%
      group_by(individual_index) %>% mutate(timestep_born = max(timestep_born, na.rm = TRUE)) %>% ungroup()
    
    # Now that we have timestep born, get age at each timestep and final
    
    age_cohort <- age_cohort %>%
      # Flag if individuals were born during the sim
      group_by(individual_index) %>%
      mutate(ever_died = case_when("died" %in% process ~ 1,
                                   !("died" %in% process) ~ 0)) %>%
      ungroup() %>%
      # Because we want to stop aging after death, get timestep of death
      mutate(timestep_died = case_when(process == "died" ~ timestep,
                                       process != "died" ~ 0)) %>%
      group_by(individual_index) %>% mutate(timestep_died = max(timestep_died)) %>% ungroup() %>%
      # If timestep_died = 0 (ie, individual not died - substitute to NA)
      mutate(timestep_died = case_when(timestep_died == 0 ~ NA_real_,
                                       timestep_died != 0 ~ timestep_died)) %>%
      # Estimate their age at each time step (as days and as year rounding to lowest year)
      # Careful to use either timestep, or timestep_died after the person dies!!
      mutate(age_at_time = (pmin(timestep, coalesce(timestep_died, timestep)) - timestep_born)) %>%
      mutate(age_at_time_year = floor((pmin(timestep, coalesce(timestep_died, timestep)) - timestep_born)/365)) %>%
      # Maximum age reached (at death or end of sim)
      group_by(individual_index) %>% mutate(max_age = max(age_at_time_year)) %>% ungroup()
    
      # Finally, add some vars to signal whether the person ever received/removed net
      # or ever received treatment
    
    age_cohort <- age_cohort %>%
      group_by(individual_index) %>%
      mutate(recieved_treat = case_when("Tr" %in% state ~ 1,
                                        !("Tr" %in% state) ~ 0)) %>%
      mutate(recieved_net = case_when("recieved_net" %in% process ~ 1,
                                      !("recieved_net" %in% process) ~ 0)) %>%
      mutate(removed_net = case_when("removed_net" %in% process ~ 1,
                                     !("removed_net" %in% process) ~ 0)) %>%
      ungroup()
    
    return(age_cohort)

  } else {
    
    print("Oops! We cannot get age for everyone in the sim...")
    print("Revise age limit for age output is absurdly high (i.e. 100 years) and that snapshot timestep is 1!")
  }

}