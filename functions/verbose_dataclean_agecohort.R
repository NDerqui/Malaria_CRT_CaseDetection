
# INPUT: a df result of a verbose simulation.
# OUTPUT: a df with verbose output only for individuals for whom we know age.

# DESCRIPTION:

# From the total output of a verbose sim, function subsets to the specific cohort
# for which we can know the age of each individual (age-cohort) .
# Age is estimated as timestep - timestep_born,
# thus this cohort is restricted to individuals born in the simulation.
# The function notes age at each timestep and age at final step.
# Note the final step can be death or end of sim.

get_age_cohort <- function(df) {
  
  require(dplyr)

  age_cohort <- df %>%
    # Flag and filter if individuals were born during the sim
    group_by(individual_index) %>%
    mutate(ever_born = case_when("born" %in% process ~ 1,
                                 !("born" %in% process) ~ 0)) %>%
    ungroup() %>% filter(ever_born == 1) %>% select(-ever_born) %>%
    # To get age, get the timestep at which the person is born
    mutate(timestep_born = case_when(process == "born" ~ timestep,
                                     process != "born" ~ 0)) %>%
    group_by(individual_index) %>% mutate(timestep_born = max(timestep_born)) %>% ungroup() %>%
    # Because we want to stop aging after death, get timestep of death
    mutate(timestep_died = case_when(process == "died" ~ timestep,
                                     process != "died" ~ 0)) %>%
    group_by(individual_index) %>% mutate(timestep_died = max(timestep_died)) %>% ungroup() %>%
    # If timestep_died = 0 (ie, individual not died - substitute to NA)
    mutate(timestep_died = case_when(timestep_died == 0 ~ NA_real_,
                                     timestep_died != 0 ~ timestep_died)) %>%
    # Estimate their age at each time step as a year rounding to lowest year
    # Careful to use either timestep, or timestep_died after the person dies!!
    mutate(age_year = floor((pmin(timestep, coalesce(timestep_died, timestep)) - timestep_born)/365)) %>%
    # Maximum age reached (at death or end of sim)
    group_by(individual_index) %>% mutate(max_age = max(age_year)) %>% ungroup() %>%
    mutate(age_gr_final = case_when(max_age > 15 ~ ">15 y",
                                    max_age <= 15 & max_age > 10 ~ "11-15 y",
                                    max_age <= 10 & max_age > 2 ~ "3-10 y",
                                    max_age <= 2 & max_age > 0 ~ "1-2 y",
                                    max_age == 0 ~ "<1 y")) %>%
    mutate(age_gr_final = factor(age_gr_final,
                                 levels = c("<1 y", "1-2 y", "3-10 y", "11-15 y", ">15 y"))) %>%
    # Finally, add some vars to signal whether the person ever received/removed net
    # or ever received treatment
    group_by(individual_index) %>%
    mutate(recieved_treat = case_when("Tr" %in% state ~ 1,
                                    !("Tr" %in% state) ~ 0)) %>%
    mutate(recieved_net = case_when("recieved_net" %in% process ~ 1,
                                    !("recieved_net" %in% process) ~ 0)) %>%
    mutate(removed_net = case_when("removed_net" %in% process ~ 1,
                                    !("removed_net" %in% process) ~ 0)) %>%
    ungroup()
  
  return(age_cohort)
  
}