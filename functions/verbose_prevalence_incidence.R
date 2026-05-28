
# INPUT: a df with new infections or cases flagged (still one obs per time per person).
# OUTPUT: a df with one observation per timestep, with information on prevalence and incidence.


# DESCRIPTION:

# Functions to get:
# infection prevalence at each timestep (infected at timestep/population),
# infection incidence at each timestep (newly infected at timestep/population at risk),
# case prevalence at each timestep (case at timestep/population),
# case incidence at each timestep (newly clinical case at timestep/population at risk).

# Option to get them for all timepoints (above),
# or preva/incidence as captured by a cross-sectional survey (specific timeoints).

get_prev_inc <- function(df) {
  
  require(dplyr)
  
  # Prepare to iter over all timesteps
  
  timesteps <- unique(df$timestep)
  result <- data.frame()
  
  # Loop over all timesteps to get prevalence and incidence at each timstep
  
  for (time in 1:length(timesteps)) {
    
    add <- df %>%
      # To ensure timings of transitions come okay...
      arrange(individual_index, timestep) %>%
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
      filter(is.na(timestep_died) | timestep_died > timesteps[time]) %>%
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

# Same as above adding age group

get_prev_inc_by_age <- function(df) {
  
  require(dplyr)
  
  # Prepare to iter over all timesteps
  
  timesteps <- unique(df$timestep)
  result <- data.frame()
  
  # Loop over all timesteps to get prevalence and incidence at each timstep
  
  for (time in 1:length(timesteps)) {
    
    add <- df %>%
      # To ensure timings of transitions come okay...
      arrange(individual_index, timestep) %>%
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
      filter(is.na(timestep_died) | timestep_died > timesteps[time]) %>%
      ## Age group
      group_by(age_at_time_year) %>%
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
      ungroup() %>%
      # Cleaning
      group_by(timestep, age_at_time_year) %>%
      filter(row_number() == 1) %>% ungroup() %>%
      select(timestep, age_at_time_year, n, at_risk,
             infections, cases, new_infections, new_cases,
             prevalence_infec, prevalence_case, incidence_infec, incidence_case)
    
    result <- rbind(result, add)
    rm(add)
  }
  
  return(result)
}


# ONLY FOR SELECTED TIMESTEPS (CROSS-SECTIONAL SURVEYS)

get_prev_survey <- function(df, cross_surveys, trial_start) {
  
  require(dplyr)
  year <- 365
  
  # Before getting any prevalence estimate,
  # subset to the timepoints of the cross-sectional survey
  
  # IMP: cross surveys passed as year measure, i.e. 6 month = 0.5 year
  
  cross_survey_times <- (trial_start + cross_surveys)*year
  
  df <- df %>%
    filter(timestep %in% cross_survey_times)
  
  # Prepare to iter over all timesteps
  
  timesteps <- unique(df$timestep)
  result <- data.frame()
  
  # Loop over all timesteps to get prevalence and incidence at each timstep
  
  for (time in 1:length(timesteps)) {
    
    add <- df %>%
      ## Ready for the each timestep calcs
      # Filter to each timestep and remove everyone dead by then (for denom)
      filter(timestep == timesteps[time]) %>%
      filter(is.na(timestep_died) | timestep_died > timesteps[time]) %>%
      # Get some basic counts
      mutate(n = n()) %>%
      mutate(infections = sum(infected_at_time)) %>%
      mutate(cases = sum(case_at_time)) %>%
      # Incidence and Prevalence calcs
      mutate(prevalence_infec = infections/n) %>%
      mutate(prevalence_case = cases/n) %>%
      # Cleaning
      group_by(timestep) %>% filter(row_number() == 1) %>% ungroup() %>%
      select(timestep, n,
             infections, cases, prevalence_infec, prevalence_case)
    
    result <- rbind(result, add)
    rm(add)
  }
  
  return(result)
}

get_inc_survey <- function(df, routine_visits, trial_start, days_catchment) {
  
  require(dplyr)
  year <- 365
  
  # Before getting any incidence estimate,
  # need to prepare data flagging for each indv and timestep their previous state.
  
  df <- df %>%
    # To ensure timings of transitions come okay...
    arrange(individual_index, timestep) %>%
    ## For incidence calcs:
    # get no. at risk at each timestep considering prev timstep.
    # No at risk of incident infection won't include infected in previous timestep.
    group_by(individual_index) %>%
    mutate(prev_state = lag(state)) %>%
    ungroup()
    
  # Subset to the timepoints of the cross-sectional survey.
  # For incidence calcs we want a "catchment period" too:
  # include new infections in last x days to model "tested for malaria if fever <48 h"
  
  # IMP: routine visits passed as week measure (7 days)
  
  routine_visits_times <- (trial_start*year + routine_visits*7)
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
  
  # Prepare to iter over all timesteps
  
  timesteps <- unique(df$timestep)
  result <- data.frame()
  
  # Loop over all timesteps to get prevalence and incidence at each timstep
  
  for (time in 1:length(timesteps)) {
    
    add <- df %>%
      ## Ready for the each timestep calcs
      # Filter to each timestep and remove everyone dead by then (for denom)
      filter(timestep == timesteps[time]) %>%
      filter(is.na(timestep_died) | timestep_died > timesteps[time]) %>%
      # Get some basic counts
      mutate(n = n()) %>%
      mutate(at_risk = sum(!(prev_state %in% c("U", "A", "D", "Tr")))) %>%
      mutate(new_infections = sum(new_infection_at_time)) %>%
      mutate(new_cases = sum(new_case_at_time)) %>%
      # Incidence and Prevalence calcs
      mutate(incidence_infec = new_infections/at_risk) %>%
      mutate(incidence_case = new_cases/at_risk) %>%
      # Cleaning
      group_by(timestep) %>% filter(row_number() == 1) %>% ungroup() %>%
      select(timestep, n, at_risk, new_infections, new_cases,
             incidence_infec, incidence_case)
    
    result <- rbind(result, add)
    rm(add)
  }
  
  return(result)
}


# AGGREGATED INCIDENCE OVER FOLLOW-UP PERIODS

# INPUT: a df with new infections or cases flagged (still one obs per time per person).
# OUTPUT: a df with one observation per follow-up period, with incidence rates.

# DESCRIPTION:

# Aggregates incident infections/cases over longer follow-up windows, such as
# 6 months or 1 year, instead of estimating incidence at each timestep.

get_incidence_period <- function(df,
                                 trial_start,
                                 period = 0.5,
                                 followup = NULL,
                                 by = NULL,
                                 infected_states = c("U", "A", "D", "Tr")) {

  require(dplyr)

  year <- 365
  start_time <- trial_start * year
  period_days <- period * year

  if (is.null(followup)) {
    end_time <- max(df$timestep, na.rm = TRUE)
  } else {
    end_time <- start_time + followup * year
  }

  grouping_vars <- c(by, "period")

  df_period <- df %>%
    arrange(individual_index, timestep) %>%
    group_by(individual_index) %>%
    mutate(prev_state = lag(state)) %>%
    ungroup() %>%
    filter(timestep >= start_time, timestep < end_time) %>%
    filter(is.na(timestep_died) | timestep_died > timestep) %>%
    mutate(
      period = floor((timestep - start_time) / period_days) + 1,
      period_start = start_time + (period - 1) * period_days,
      period_end = pmin(start_time + period * period_days, end_time),
      at_risk = !(prev_state %in% infected_states)
    )

  df_period %>%
    group_by(across(all_of(grouping_vars))) %>%
    summarise(
      period_start = min(period_start),
      period_end = max(period_end),
      person_days_alive = n(),
      person_days_at_risk = sum(at_risk, na.rm = TRUE),
      n_people = n_distinct(individual_index),
      new_infections = sum(new_infection_at_time, na.rm = TRUE),
      new_cases = sum(new_case_at_time, na.rm = TRUE),
      incidence_infection = new_infections / person_days_at_risk,
      incidence_case = new_cases / person_days_at_risk,
      incidence_infection_per_py = incidence_infection * year,
      incidence_case_per_py = incidence_case * year,
      .groups = "drop"
    )
}
