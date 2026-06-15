
# INPUT: a df result of a verbose simulation with birth/death signalled.
# OUTPUT(S): a df with only a subset of indivudals to track.


# DESCRIPTION:

# Function to subset to a number of individuals from the verbose sim,
# to only "follow" them in a trial simulation.

# Mimic that transmission would occur in a wider community,
# but only a few would be enrolled in trial.

# Note, only going to subset among those alive at key timestep (i.e. trial start),
# thus need to run function after getting the timestep_died.

get_enrol_sample <- function(df, alive_by = min(df$timestep), trial_size,
                             age_min = 0, age_max = 100) {
  
  year <- 365
  
  require(dplyr)
  
  # First, filter to individuals alive at our timepoint of interest (e.g. trial start):
  # They must have not died but also they must have been born by then.
  # (By default, those alive at start of sim, smallest timestep)
  # Trials clear of infection at trial start, so need to filter
  # to individuals that are S on trial start.
  # Finally, if desired, sample only a particular age group at our timepoint
  # (by default, everyone 0-100 years)
  
  df <- df %>%
    filter(timestep_born < alive_by) %>%
    filter(is.na(timestep_died) | timestep_died > alive_by) %>%
    group_by(individual_index) %>%
    filter(any(timestep == alive_by & state == "S")) %>%
    ungroup() %>%
    # With the age group, we want people that at trial start have a particular age
    group_by(individual_index) %>%
    mutate(in_age_group = any(age_at_time[timestep == alive_by] >= (age_min * year) & age_at_time[timestep == alive_by] <= (age_max * year))) %>%
    ungroup() %>%
    filter(in_age_group) %>% select(-in_age_group)
  
  # Now randomly select x number of indv as our sample
  
  random_sample <- sample(unique(df$individual_index),
                          size = trial_size, # Our trial (or trial cluster) size
                          replace = FALSE)   # No duplicates
  
  df <- df %>%
    filter(individual_index %in% random_sample)
  
  return(df)
}