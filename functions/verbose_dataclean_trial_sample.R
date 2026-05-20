
# INPUT: a df result of a verbose simulation with birth/death signalled.
# OUTPUT(S): a df with only a subset of indivudals to track.


# DESCRIPTION:

# Function to subset to a number of individuals from the verbose sim,
# to only "follow" them in a trial simulation.

# Mimic that transmission would occur in a wider community,
# but only a few would be enrolled in trial.

# Note, only going to subset among those alive at key timestep (i.e. trial start),
# thus need to run function after getting the timestep_died.

get_enrol_sample <- function(df, alive_by = min(df$timestep), trial_size) {
  
  require(dplyr)
  
  # First, filter to individuals alive at our timepoint of interest (e.g. trial start):
  # They must have not died but also they must have been born by then.
  # (By default, those alive at start of sim, smallest timestep)
  
  df <- df %>%
    filter(timestep_born < alive_by) %>%
    filter(is.na(timestep_died) | timestep_died > alive_by)
  
  # Now randomly select x number of indv as our sample
  
  random_sample <- sample(unique(df$individual_index),
                          size = trial_size, # Our trial (or trial cluster) size
                          replace = FALSE)   # No duplicates
  
  df <- df %>%
    filter(individual_index %in% random_sample)
  
  return(df)
}