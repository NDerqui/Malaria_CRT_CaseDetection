
# INPUT: a df result of a verbose simulation, or from an age-cohort df.
# OUTPUT: a df with new infections or cases flagged (still one obs per time per person).

# DESCRIPTION:

# The functions detect new infections and new clinical infections (cases),
# but there is also a function to signal all infected states or clinical states
# (to be used with a cross sectional survey sim).


## For new infections, analyse possible processes changes first,
## and select new As, Ds and Ts.

detect_new_infection <- function(df) {
  
  require(dplyr)
  
  # New appearance of Ds and Ts are, by model spec, new infections.
  # But new infections that result in A are trickier, so ...
  # all movements to A that are NOT from D are new infections.
  
  # (Context: at each timestep, process shown with state from previous timestep process)
  
  new_A_infections <- c("Gone_to_A - S", "Gone_to_A - Tr",
                        "Gone_to_A - A", "Gone_to_A - U")

  agecohort_inc <- df %>%
    # To ensure timings of transitions come okay...
    arrange(individual_index, timestep) %>%
    # First get the transition pairing process to state,
    # because some processes will be new infection or not depending on previous state.
    mutate(transition = paste0(process, " - ", state)) %>%
    # Everything by indv_index as we want to follow
    group_by(individual_index) %>%
    # Detect appearances of a new infection as:
    # New D or T (by model definition, new Ds or Ts are new infections), or as
    # an accepted transition to A (see above)
    mutate(new_infection = (transition %in% new_A_infections) |
             (state == "D" & lag(state, default = first(state)) != "D") |
             (state == "Tr" & lag(state, default = first(state)) != "Tr")) %>%
    select(-transition) %>%
    # Count the number of infection in that individual (by each timestep)
    mutate(infection_by_time = cumsum(new_infection)) %>%
    ungroup()
  
  return(agecohort_inc)
  
}


## For new clinical cases, straightforward (select new T and new D).

detect_new_clinical <- function(df) {
  
  require(dplyr)
  
  agecohort_inc <- df %>%
    # To ensure timings of transitions come okay...
    arrange(individual_index, timestep) %>%
    # Everything by indv_index as we want to follow
    group_by(individual_index) %>%
    # Detect new appearances of D/T using lag function (TRUE / FALSE)
    mutate(new_D = state == "D" & lag(state, default = first(state)) != "D") %>%
    mutate(new_T = state == "Tr" & lag(state, default = first(state)) != "Tr") %>%
    # Count the number of each state in that individual (by each timestep)
    mutate(D_by_time = cumsum(new_D)) %>%
    mutate(T_by_time = cumsum(new_T)) %>%
    ungroup()
  
  return(agecohort_inc)
  
}


## For signalling all infections or cases, no matter the start,
## just flag Us, As, Ds and Ts.

detect_all_infection <- function(df) {
  
  require(dplyr)
  
  agecohort_inc <- df %>%
    # Flag all infection states
    mutate(infection = state %in% c("U", "A", "D", "Tr"))
  
  return(agecohort_inc)
  
}

detect_all_clinical <- function(df) {
  
  require(dplyr)
  
  agecohort_inc <- df %>%
    # Flag all infection states
    mutate(case = state %in% c("D", "Tr"))
  
  return(agecohort_inc)
  
}
