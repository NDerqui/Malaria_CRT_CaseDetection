
# INPUT: a df with new infections or cases flagged (still one obs per time per person).
# OUTPUT: a df with one observation per timestep, with information on prevalence and incidence.


# DESCRIPTION:

# Function to get:
# infection prevalence at each timestep (infected at timestep/population),
# infection incidence at each timestep (newly infected at timestep/population at risk),
# case prevalence at each timestep (case at timestep/population),
# case incidence at each timestep (newly clinical case at timestep/population at risk).

estimate_true_realtime_outcomes <- function(df) {
  
  require(dplyr)
  
  # Prepare data once before grouped calcs
  # then get prevalence and incidence each timestep
  
  result <- df %>%
    # To ensure timings of transitions come okay...
    arrange(individual_index, timestep) %>%
    ## Get no. at risk in this timestep considering prev timstep:
    # To later get no. at risk, get the state at our prior timestep.
    # No at risk of incident infection won't include infected in previous timestep.
    group_by(individual_index) %>%
    mutate(prev_state = lag(state)) %>%
    ungroup() %>%
    ## Ready for the each timestep calcs
    # Filter to each timestep and remove everyone dead by then (for denom)
    group_by(timestep) %>%
    filter(is.na(timestep_died) | timestep_died > timestep) %>%
    # Get some basic counts
    mutate(n = n()) %>%
    mutate(person_days_at_risk = sum(!(prev_state %in% c("U", "A", "D", "Tr")))) %>%
    mutate(infections = sum(infected_at_time)) %>%
    mutate(cases = sum(case_at_time)) %>%
    mutate(new_infections = sum(new_infection_at_time)) %>%
    mutate(new_cases = sum(new_case_at_time)) %>%
    # Incidence and Prevalence calcs
    mutate(prevalence_infection = infections/n) %>%
    mutate(prevalence_case = cases/n) %>%
    mutate(incidence_ppd_infection = new_infections/person_days_at_risk) %>%
    mutate(incidence_ppd_case = new_cases/person_days_at_risk) %>%
    mutate(incidence_ppy_infection = new_infections/(person_days_at_risk/365)) %>%
    mutate(incidence_ppy_case = new_cases/(person_days_at_risk/365)) %>%
    # Cleaning
    filter(row_number() == 1) %>% ungroup() %>%
    mutate(type_measure = "True Instantaneous") %>%
    select(timestep, type_measure, n, person_days_at_risk,
           infections, cases, new_infections, new_cases,
           prevalence_infection, prevalence_case,
           incidence_ppd_infection, incidence_ppd_case,
           incidence_ppy_infection, incidence_ppy_case)
  
  return(result)
}

# Same as above adding age group

estimate_true_realtime_outcomes_by_age <- function(df) {
  
  require(dplyr)
  
  # Prepare data once before grouped calcs
  # then get prevalence and incidence each timestep
  
  result <- df %>%
    # To ensure timings of transitions come okay...
    arrange(individual_index, timestep) %>%
    ## Get no. at risk in this timestep considering prev timstep:
    # To later get no. at risk, get the state at our prior timestep.
    # No at risk of incident infection won't include infected in previous timestep.
    group_by(individual_index) %>%
    mutate(prev_state = lag(state)) %>%
    ungroup() %>%
    ## Ready for the each timestep calcs
    # Filter to each timestep and remove everyone dead by then (for denom)
    group_by(timestep, age_at_time_year) %>%
    filter(is.na(timestep_died) | timestep_died > timestep) %>%
    # Get some basic counts
    mutate(n = n()) %>%
    mutate(person_days_at_risk = sum(!(prev_state %in% c("U", "A", "D", "Tr")))) %>%
    mutate(infections = sum(infected_at_time)) %>%
    mutate(cases = sum(case_at_time)) %>%
    mutate(new_infections = sum(new_infection_at_time)) %>%
    mutate(new_cases = sum(new_case_at_time)) %>%
    # Incidence and Prevalence calcs
    mutate(prevalence_infection = infections/n) %>%
    mutate(prevalence_case = cases/n) %>%
    mutate(incidence_ppd_infection = new_infections/person_days_at_risk) %>%
    mutate(incidence_ppd_case = new_cases/person_days_at_risk) %>%
    mutate(incidence_ppy_infection = new_infections/(person_days_at_risk/365)) %>%
    mutate(incidence_ppy_case = new_cases/(person_days_at_risk/365)) %>%
    # Cleaning
    filter(row_number() == 1) %>% ungroup() %>%
    mutate(type_measure = "True Instantaneous") %>%
    select(timestep, type_measure, age_at_time_year, n, person_days_at_risk,
           infections, cases, new_infections, new_cases,
           prevalence_infection, prevalence_case,
           incidence_ppd_infection, incidence_ppd_case,
           incidence_ppy_infection, incidence_ppy_case)
  
  return(result)
}

# AGGREGATED INCIDENCE OVER FOLLOW-UP PERIODS
# Main diff: instead each timestep, each row will represent aggregate new infections in a period

# Aggregates incident infections/cases over longer follow-up windows, such as
# 6 months or 1 year, instead of estimating incidence at each timestep.

estimate_true_aggregate_incidence <- function(df, trial_start,
                                       followup_period_year = 1,
                                       followup_end = max(df$timestep, na.rm = TRUE)) {
  
  require(dplyr)
  year <- 365
  
  # Define how many periods we will have considering length of followup_period:
  
  # IMP: followup period passed as year measure, i.e. 6 month = 0.5 year
  
  number_periods <- round((followup_end - trial_start*year) / (followup_period_year*year))
  
  periods <- data.frame(period = 1:number_periods,
                        period_label = paste0(
                          seq(0, by = followup_period_year*12, length.out = number_periods),
                          "-",
                          seq(followup_period_year*12, by = followup_period_year*12, length.out = number_periods),
                          " months"
                        ))
  
  # Identifying to which period does each timestep belong to
  
  followup_period_days <- round(followup_period_year * year)
  
  df <- df %>%
    # Only follow up from trial start
    filter(timestep >= trial_start*year) %>%
    # Getting period for each timestep
    mutate(
      period = pmin(
        floor((timestep - trial_start*year)/followup_period_days) + 1, # Add a one so that periods don't start at zero
        number_periods # pmin so that last timestep doesn't get assigned to other period out of a bad denominator
      )
    ) %>% merge(periods)
  
  # Prepare data once before grouped calcs
  # then get aggregated incidence for each period
  
  result <- df %>%
    # To ensure timings of transitions come okay...
    arrange(individual_index, timestep) %>%
    ## Get no. at risk in this timestep considering prev timstep:
    # To later get no. at risk, get the state at our prior timestep.
    # No at risk of incident infection won't include infected in previous timestep.
    group_by(individual_index) %>%
    mutate(prev_state = lag(state)) %>%
    ungroup() %>%
    ## Ready for the each period calcs
    # Filter to each timestep and remove everyone dead by then (for denom)
    group_by(period) %>%
    filter(is.na(timestep_died) | timestep_died > timestep) %>%
    # Get some basic counts
    mutate(n = n()) %>%
    mutate(person_days_at_risk = sum(!(prev_state %in% c("U", "A", "D", "Tr")))) %>%
    mutate(new_infections = sum(new_infection_at_time)) %>%
    mutate(new_cases = sum(new_case_at_time)) %>%
    # Incidence 
    mutate(incidence_ppd_infection = new_infections/person_days_at_risk) %>%
    mutate(incidence_ppd_case = new_cases/person_days_at_risk) %>%
    mutate(incidence_ppy_infection = new_infections/(person_days_at_risk/365)) %>%
    mutate(incidence_ppy_case = new_cases/(person_days_at_risk/365)) %>%
    # Get the last timestep of each period to plot
    mutate(timestep = max(timestep)) %>%
    # Cleaning
    filter(row_number() == 1) %>% ungroup() %>%
    mutate(type_measure = paste0("True, aggregate over ", followup_period_year*12, " mos.")) %>%
    select(timestep, type_measure, period, period_label,
           n, person_days_at_risk,
           new_infections, new_cases,
           incidence_ppd_infection, incidence_ppd_case,
           incidence_ppy_infection, incidence_ppy_case)
  
  return(result)
  
}


### ONLY FOR SELECTED TIMESTEPS (CROSS-SECTIONAL SURVEYS or ROUTINE VISITS)

estimate_survey_prevalence <- function(df, cross_surveys_in_years, trial_start) {
  
  require(dplyr)
  year <- 365
  
  # Before getting any prevalence estimate,
  # subset to the timepoints of the cross-sectional survey
  
  # IMP: cross surveys passed as year measure, i.e. 6 month = 0.5 year
  
  cross_survey_times <- round((trial_start + cross_surveys_in_years)*year)
  
  df <- df %>%
    filter(timestep %in% cross_survey_times)
  
  # Prepare data once before grouped calcs
  # then get prevalence each timestep
  
  result <- df %>%
    # To ensure timings of transitions come okay...
    arrange(individual_index, timestep) %>%
    ## Ready for the each timestep calcs
    # Filter to each timestep and remove everyone dead by then (for denom)
    group_by(timestep) %>%
    filter(is.na(timestep_died) | timestep_died > timestep) %>%
    # Get some basic counts
    mutate(n = n()) %>%
    mutate(infections = sum(infected_at_time)) %>%
    mutate(cases = sum(case_at_time)) %>%
    # Incidence and Prevalence calcs
    mutate(prevalence_infection = infections/n) %>%
    mutate(prevalence_case = cases/n) %>%
    # Cleaning
    filter(row_number() == 1) %>% ungroup() %>%
    mutate(type_measure = "Cross-sectional surveys") %>%
    select(timestep, type_measure, n,
           infections, cases, prevalence_infection, prevalence_case)
  
  return(result)
}

estimate_acd_incidence <- function(df, trial_start,
                             routine_visits_in_weeks, days_catchment,
                             followup_period_year = 1,
                             followup_end = max(df$timestep, na.rm = TRUE)) {
  
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
    
  # Subset to the timepoints of the Active Case Detection routine visits.
  # For incidence calcs we want a "catchment period" too:
  # include new infections in last x days to model "tested for malaria if fever <48 h"
  
  # IMP: routine visits passed as week measure (7 days)
  
  routine_visits_times <- (trial_start*year + routine_visits_in_weeks*7)
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
  
  # Because this is incidence calculation, we want to estimate incidence over period of time.
  
  # Define how many periods we will have considering length of followup_period.
  # IMP: followup period passed as year measure, i.e. 6 month = 0.5 year
  
  number_periods <- round((followup_end - trial_start*year) / (followup_period_year*year))
  
  periods <- data.frame(period = 1:number_periods,
                        period_label = paste0(
                          seq(0, by = followup_period_year*12, length.out = number_periods),
                          "-",
                          seq(followup_period_year*12, by = followup_period_year*12, length.out = number_periods),
                          " months"
                        ))
  
  # Identifying to which period does each timestep belong to
  
  followup_period_days <- followup_period_year * year
  
  df <- df %>%
    # Only follow up from trial start
    filter(timestep >= trial_start*year) %>%
    # Getting period for each timestep
    mutate(
      period = pmin(
        floor((timestep - trial_start*year)/followup_period_days) + 1, # Add a one so that periods don't start at zero
        number_periods # pmin so that last timestep doesn't get assigned to other period out of a bad denominator
      )
    ) %>% merge(periods)
  
  # Prepare data once before grouped calcs
  # then get aggregated incidence for each period
  
  result <- df %>%
    ## Ready for the each period calcs
    # Filter to each timestep and remove everyone dead by then (for denom)
    group_by(period) %>%
    filter(is.na(timestep_died) | timestep_died > timestep) %>%
    # Get some basic counts
    mutate(n = n()) %>%
    mutate(person_days_at_risk = sum(!(prev_state %in% c("U", "A", "D", "Tr")))) %>%
    mutate(new_infections = sum(new_infection_at_time)) %>%
    mutate(new_cases = sum(new_case_at_time)) %>%
    # Incidence 
    mutate(incidence_ppd_infection = new_infections/person_days_at_risk) %>%
    mutate(incidence_ppd_case = new_cases/person_days_at_risk) %>%
    mutate(incidence_ppy_infection = new_infections/(person_days_at_risk/365)) %>%
    mutate(incidence_ppy_case = new_cases/(person_days_at_risk/365)) %>%
    # Get the last timestep of each period to plot
    mutate(timestep = max(timestep)) %>%
    # Cleaning
    filter(row_number() == 1) %>% ungroup() %>%
    mutate(type_measure = paste(paste0("ACD every ", routine_visits_in_weeks[1], " weeks"),
                                paste0("with ", days_catchment, " days window,"),
                                paste("aggr. over ", followup_period_year*12, " mos."),
                                sep = "\n")) %>%
    select(timestep, type_measure, period, period_label,
           n, person_days_at_risk,
           new_infections, new_cases,
           incidence_ppd_infection, incidence_ppd_case,
           incidence_ppy_infection, incidence_ppy_case)
  
  return(result)
}


