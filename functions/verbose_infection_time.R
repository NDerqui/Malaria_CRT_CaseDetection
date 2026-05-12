
# INPUT: a df with new infections or cases flagged (still one obs per time per person).
# OUTPUT: a df with time-to-event


# DESCRIPTION:

# Functions to get:

get_time_to_event <- function(df, time_inter) {
  
  require(dplyr)
  
  df <- df %>%
    # Censoring
    # Leaving this space to see which indv to censor
    ##
    # First, estimate time from intervention each timestep
    mutate(time_event = timestep - time_inter) %>%
    # Only interested in events occuring after intervention (no negative time_event)
    filter(timestep >= time_inter) %>%
    # Only interested on time to new_infections or new_case
    filter(new_infection_at_time | new_case_at_time) %>%
    # For each individual, get smallest time
    group_by(individual_index) %>%
    mutate(
      time_to_infection = if (any(new_infection_at_time)) {min(time_event[new_infection_at_time], na.rm = TRUE)} else {NA_real_},
      time_to_case = if (any(new_case_at_time)) {min(time_event[new_case_at_time], na.rm = TRUE)} else {NA_real_}) %>%
    filter(row_number() == 1) %>%
    ungroup() 
    
  return(df)
}

