
# INPUT: a df result of a verbose simulation.
# OUTPUT: a df with verbose output only for individuals for whom we can know age.

# DESCRIPTION:

# From the total output of a verbose sim,
# function that gets age-cohort (sub-group for which we can know their age).

# To estimate age, we use an age snapshot (gives the age of all individuals alive at that time),
# and the timestep individuals are born (for those who are born after the snapshot).

# Age is estimated as timestep - timestep_born.
# The function notes age at each timestep and age at final step.
# Note the final step can be death or end of sim.

get_age_cohort <- function(df, age_snapshot, snapshot_time) {
  
  require(dplyr)
  
  # First, signal if individuals are born during the sim
  
  df <- df %>%
    # Flag if individuals were born / died during the sim
    group_by(individual_index) %>%
    mutate(ever_born = case_when("born" %in% process ~ 1,
                                 !("born" %in% process) ~ 0)) %>%
    mutate(ever_died = case_when("died" %in% process ~ 1,
                                 !("died" %in% process) ~ 0)) %>%
    ungroup()
  
  # Second, check and filter to individuals who we can have age for:
  # Individuals alive at age snapshot + those born after.
  
  those_born <- unique(df$individual_index[df$ever_born == 1])
  those_died <- unique(df$individual_index[df$ever_died == 0])
  those_snapshot <- unique(age_snapshot$individual_index)
  
  age_cohort <- df %>%
    filter(individual_index %in% unique(c(those_born, those_snapshot)))
  
  print(paste0("We can have age for ", round(nrow(age_cohort)*100/nrow(df), 2), "% of sim pop."))
  
  # Third, add the age of our snapshot to the main df
  
  age_snapshot <- age_snapshot %>%
    filter(timestep == snapshot_time) %>% # Ensure we are on the snapshot time
    select(individual_index, age) %>%
    rename(age_at_snapshot = age)
  
  age_cohort <- merge(age_cohort, age_snapshot, all = TRUE)
  
  # Next move to calculate age:
  # To do this, get timestep at which individuals are born.
  
  # Straightforward for those born during the sim, and
  # calculate for those with snapshot considering how much they aged at time of snapshot.
  
  age_cohort <- age_cohort %>%
    # For those born, get the timestep for that process
    mutate(timestep_born = case_when(process == "born" ~ timestep,
                                     process != "born" ~ 0)) %>%
    group_by(individual_index) %>% mutate(timestep_born = max(timestep_born)) %>% ungroup() %>%
    # For those not born, timestep born is -(age_at_snapshot - snapshot_timestep)
    mutate(timestep_born = case_when(!is.na(age_at_snapshot) ~ -(age_at_snapshot - snapshot_time),
                                     is.na(age_at_snapshot) ~ timestep_born)) %>%
    # Clean
    select(-age_at_snapshot)
  
  # Now that we have timestep born, get age at each timestep and final
  
  age_cohort <- age_cohort %>%
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
    ungroup() %>%
    # Clean
    arrange(individual_index, timestep)
  
  return(age_cohort)

}