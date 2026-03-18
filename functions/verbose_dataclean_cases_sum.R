
# INPUT: a df after the detect infections or detect cases functions (still one obs per time per person).
# OUTPUT: a summarised df with only one obs per person, or per person per age.

# DESCRIPTION:

# The functions estimate the overall infections/cases an individual has:
# (i) in the sim and (ii) at each age.

# Note the count_by_age funcs will only work if the detect infection/cases
# step has been done on the filtered age cohort.


## First, some functions to get one obs per individual or per indvidiual & age,
## with the final count of infections or clinical cases.

count_infections <- function(df) {
  
  require(dplyr)
  
  basic_info <- c("individual_index", "timestep_born", "timestep_died",
                  "max_age", "age_gr_final", "recieved_treat", "recieved_net", "removed_net")
  
  inc_count <- df %>%
    # Group by individual and only get one line
    group_by(individual_index) %>%
    # Create vars with total infections or clinical cases (in that individual over time)
    mutate(total_infection = max(infection_by_time)) %>%
    # Only one obs per individual
    filter(row_number() == 1) %>%
    ungroup() %>%
    select(any_of(basic_info), total_infection)
  
   return(inc_count)
  
}
count_clinical <- function(df) {
  
  require(dplyr)
  
  basic_info <- c("individual_index", "timestep_born", "timestep_died",
                  "max_age", "age_gr_final", "recieved_treat", "recieved_net", "removed_net")
  
  inc_count <- df %>%
    # Group by individual and only get one line
    group_by(individual_index) %>%
    # Create vars with total infections or clinical cases (in that individual over time)
    mutate(total_D = max(D_by_time),
           total_T = max(T_by_time),
           total_clinical = total_D + total_T) %>%
    # Only one obs per individual
    filter(row_number() == 1) %>%
    ungroup() %>%
    select(any_of(basic_info), total_D, total_T, total_clinical)
  
  return(inc_count)
  
}

count_infection_by_age <- function(df) {
  
  require(dplyr)
  
  basic_info <- c("individual_index", "timestep_born", "timestep_died",
                  "max_age", "age_gr_final", "recieved_treat", "recieved_net", "removed_net")
  
  inc_count <- df %>%
    # Group by individual and only get one line
    group_by(individual_index, age_year) %>%
    # Create vars with total infections or clinical cases (in that individual at that age)
    mutate(total_infection_year = case_when(new_infection ~ 1)) %>%
    mutate(total_infection_year = sum(total_infection_year, na.rm = TRUE)) %>%
    # Only one obs per individual
    filter(row_number() == 1) %>%
    ungroup() %>%
    select(any_of(basic_info), age_year, infection_by_time, total_infection_year)
    
  return(inc_count)
  
}

count_clinical_by_age <- function(df) {
  
  require(dplyr)
  
  basic_info <- c("individual_index", "timestep_born", "timestep_died",
                  "max_age", "age_gr_final", "recieved_treat", "recieved_net", "removed_net")
  
  inc_count <- df %>%
    # Group by individual and only get one line
    group_by(individual_index, age_year) %>%
    # Create vars with total infections or clinical cases (in that individual at that age)
    group_by(individual_index, age_year) %>%
    mutate(total_D_year = case_when(new_D ~ 1)) %>%
    mutate(total_D_year = sum(total_D_year, na.rm = TRUE)) %>%
    mutate(total_T_year = case_when(new_T ~ 1)) %>%
    mutate(total_T_year = sum(total_T_year, na.rm = TRUE)) %>%
    mutate(total_clin_year = case_when((new_D  | new_T) ~ 1)) %>%
    mutate(total_clin_year = sum(total_clin_year, na.rm = TRUE)) %>%
    # Only one obs per individual
    filter(row_number() == 1) %>%
    ungroup() %>%
    select(any_of(basic_info),
           age_year, D_by_time, T_by_time, total_D_year, total_T_year, total_clin_year)
  
  return(inc_count)
  
}