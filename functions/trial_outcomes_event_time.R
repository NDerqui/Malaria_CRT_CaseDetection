
# INPUT: a df with new infections or cases flagged (still one obs per time per person).
# OUTPUT: a df with time-to-event


# DESCRIPTION:

# Functions to get time to event (first infection or first clinical case).
# Able to control from where to start counting time:
# e.g. start of trial or introduction of second intervention.

estimate_true_time_to_event <- function(df, time_inter) {
  
  require(dplyr)
  
  df <- df %>%
    # Get time from intervention until death
    mutate(time_to_death = timestep_died - time_inter) %>%
    # First, estimate time from intervention each timestep
    mutate(time_event = timestep - time_inter) %>%
    # Only interested in events occurring after intervention (no negative time_event)
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
    # Finally get the age at which each event (earliest infection or case) happens
    # (need to add the min because sometimes we get two lines on the same timestep)
    mutate(
      age_at_first_infection = min(age_at_time_year[time_event == time_to_infection]),
      age_at_first_case = min(age_at_time_year[time_event == time_to_case])) %>%
    filter(row_number() == 1) %>%
    ungroup() %>%
    # Clean
    select(
      # Basics of each individual
      individual_index, received_treat, received_net, removed_net,
      # Info on death (for censoring later)
      ever_died, timestep_died, 
      # Info on time to event, and whether they overall had the event
      # (the new_* vars not relevant anymore as we already got timing to each new_*)
      ever_infected, ever_case, time_to_infection, time_to_case,
      age_at_first_infection, age_at_first_case)
    
  return(df)
}

# A function to do as above but as if we were doing ACD visits:

estimate_acd_time_to_event <- function(df, time_inter,
                                   routine_visits_in_weeks,
                                   days_catchment) {
  
  require(dplyr)
  
  year <- 365
  
  # Subset to the timepoints of the Active Case Detection routine visits.
  # For incidence calcs we want a "catchment period" too:
  # include new infections in last x days to model "tested for malaria if fever <48 h"
  
  # IMP: routine visits passed as week measure (7 days)
  
  routine_visits_times <- (time_inter + routine_visits_in_weeks*7)
  # Get all events happening in interval (each_timepoint - days_catchment):each_timestep
  routine_visits_timesteps <- unlist( 
    lapply(routine_visits_times,
           function(routine_visits_times){
             return((routine_visits_times - days_catchment):routine_visits_times)}))
  # But then we want to aggregate all events happening in the interval to each routine visit time
  # Do so with an equivalence data frame (each visit timestep duplicated according to no of days in interval)
  time_eq <- data.frame(timestep = routine_visits_timesteps,
                        timestep_agg = unlist(lapply(routine_visits_times, rep, times = (1+days_catchment))))
  
  df <- df %>%
    filter(timestep %in% routine_visits_timesteps) %>%
    merge(time_eq) %>% select(-timestep) %>% rename(timestep = timestep_agg)
  
  df <- df %>%
    # Get time from intervention until death
    mutate(time_to_death = timestep_died - time_inter) %>%
    # First, estimate time from intervention each timestep
    mutate(time_event = timestep - time_inter) %>%
    # Only interested in events occurring after intervention (no negative time_event)
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
    # Finally get the age at which each event (earliest infection or case) happens
    # (need to add the min because sometimes we get two lines on the same timestep)
    mutate(
      age_at_first_infection = min(age_at_time_year[time_event == time_to_infection]),
      age_at_first_case = min(age_at_time_year[time_event == time_to_case])) %>%
    filter(row_number() == 1) %>%
    ungroup() %>%
    # Clean
    select(
      # Basics of each individual
      individual_index, received_treat, received_net, removed_net,
      # Info on death (for censoring later)
      ever_died, timestep_died, 
      # Info on time to event, and whether they overall had the event
      # (the new_* vars not relevant anymore as we already got timing to each new_*)
      ever_infected, ever_case, time_to_infection, time_to_case,
      age_at_first_infection, age_at_first_case)
  
  return(df)
}

# Prepare data for survival analysis:
# individuals without the event (infection or case) need a censoring time:
# either end of follow up period or death (censor if died).
# Follow up period is from start of time-to-event measure until sim end.

# (Keep this separate from above as it's counter-intuitive to have time-to-event
# if not having the event, only makes sense if performing a survival analysis)

prepare_time_to_event_survival <- function(df, time_inter, sim_length) {
  
  require(dplyr)
  
  df <- df %>%
    # Modify our time-to-event variables for infection
    mutate(time_to_infection = case_when(
      # keep our estimate if indv had the event
      ever_infected == TRUE ~ time_to_infection,
      # time until end of follow up if no infection or death
      ever_infected == FALSE & ever_died == 0 ~ (sim_length - time_inter), 
      # time until death if no infection but died
      ever_infected == FALSE & ever_died == 1 ~ (timestep_died - time_inter)
    )) %>%
    # Modify our time-to-event variables for cases
    mutate(time_to_case = case_when(
      # keep our estimate if indv had the event
      ever_case == TRUE ~ time_to_case,
      # time until end of follow up if no case or death
      ever_case == FALSE & ever_died == 0 ~ (sim_length - time_inter), 
      # time until death if no case but died
      ever_case == FALSE & ever_died == 1 ~ (timestep_died - time_inter)
    ))

  return(df)
}