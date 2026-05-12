
# INPUT: a df with new infections or cases flagged (still one obs per time per person).
# OUTPUT: a df with time-to-event


# DESCRIPTION:

# Functions to get:

get_time_to_event <- function(df, time_inter) {
  
  require(dplyr)
  
  df <- df %>%
    # Censoring
    # Get time from intervention until death
    mutate(time_to_death = timestep_died - time_inter) %>%
    # First, estimate time from intervention each timestep
    mutate(time_event = timestep - time_inter) %>%
    # Only interested in events occuring after intervention (no negative time_event)
    # which will automatically delete anybody dead at time_intervention
    filter(timestep >= time_inter) %>%
    # For each individual, get the smallest time to a new infection or new case
    group_by(individual_index) %>%
    mutate(
      time_to_infection = if (any(new_infection_at_time)) {min(time_event[new_infection_at_time], na.rm = TRUE)} else {NA_real_},
      time_to_case = if (any(new_case_at_time)) {min(time_event[new_case_at_time], na.rm = TRUE)} else {NA_real_}) %>%
    # Modify the ever_infected or ever_case to see individuals who have event AFTER time_inter
    mutate(
      ever_infected = if(any(new_infection_at_time)) {TRUE} else {FALSE},
      ever_case = if(any(new_case_at_time)) {TRUE} else {FALSE}) %>%
    filter(row_number() == 1) %>%
    ungroup() %>%
    # Clean
    select(
      # Basics of each individual
      individual_index, received_treat, received_net, removed_net,
      # Info on death (censoring?)
      ever_died, timestep_died, 
      # Info on time to event, and whether they overall had the event
      # (the new_* vars not relevant anymore as we already got timing to each new_*)
      ever_infected, ever_case, time_to_infection, time_to_case)
    
  return(df)
}

