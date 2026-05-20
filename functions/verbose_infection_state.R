
# INPUT: a df result of a verbose simulation, or from an age-cohort df.
# OUTPUT: a df with new infections or cases flagged (still one obs per time per person).

# DESCRIPTION:

# First, functions to detect if an indv has ever been infected / had a clinical infection.

# Then, functions to detect if an indv has a new infection and new clinical case at time,
# or simply if the indv has an infection / case at each timetep.


## Signal an ever_infected or ever_case, as the individuals
## who overall ever had malaria or a clinical case.

detect_ever_malaria <- function(df) {
  
  require(dplyr)
  
  df <- df %>%
    # To ensure timings of transitions come okay...
    arrange(individual_index, timestep) %>%
    # Everything by indv_index as we want to see overall status
    group_by(individual_index) %>%
    mutate(ever_infected = any(state %in% c("U", "A", "D", "Tr"))) %>%
    mutate(ever_case = any(state %in% c("D", "Tr"))) %>%
    ungroup()
  
  return(df)
  
}


## For new clinical cases, straightforward (select new T and new D).
## For new infections, analyse possible processes changes first,
## and select new As, Ds and Ts.

## For signalling all infections or cases, no matter the start,
## at each timestep just flag Us, As, Ds and Ts.

detect_infection <- function(df) {
  
  require(dplyr)
  
  # New D and Ts are easy to track,
  # but new infections that result in A are trickier, so we operate on that
  # all movements to A that are NOT from D are new infections.
  
  # (Context: at each timestep, process shown with state from previous timestep process)
  
  new_A_infections <- c("Gone_to_A - S", "Gone_to_A - Tr",
                        "Gone_to_A - A", "Gone_to_A - U")
  
  df <- df %>%
    # To ensure timings of transitions come okay...
    arrange(individual_index, timestep) %>%
    # First get the transition pairing process to state,
    # because some processes will be new infection or not depending on previous state.
    mutate(transition = paste0(process, " - ", state)) %>%
    # Everything by indv_index as we do not want timstep lag from other indv_index
    # but no need of timestep (as it will go row_wise)
    group_by(individual_index) %>%
    # Detect appearances of a new infection as:
    # New D or T (by model definition, new Ds or Ts are new infections), or as
    # an accepted transition to A (see above)
    mutate(new_infection_at_time = (transition %in% new_A_infections) |
             (state == "D" & lag(state, default = first(state)) != "D") |
             (state == "Tr" & lag(state, default = first(state)) != "Tr")) %>%
    select(-transition) %>%
    ungroup()
  
  # New appearance of Ds and Ts are, by model spec, new infections.
  
  df <- df %>%
    # To ensure timings of transitions come okay...
    arrange(individual_index, timestep) %>%
    # Everything by indv_index as we do not want timstep lag from other indv_index
    # but no need of timestep (as it will go row_wise)
    group_by(individual_index) %>%
    # Detect new appearances of D/T using lag function (TRUE / FALSE)
    mutate(new_D = (state == "D" & lag(state, default = first(state)) != "D")) %>%
    mutate(new_T = (state == "Tr" & lag(state, default = first(state)) != "Tr")) %>%
    mutate(new_case_at_time = new_T | new_D) %>%
    select(-c(new_D, new_T)) %>%
    ungroup()
  
  # Signal infection or case status at each timepoint
  # (This will go rowwise, so no need to group)

  df <- df %>%
    # To ensure timings of transitions come okay...
    arrange(individual_index, timestep) %>%
    # Flag all infection states
    mutate(infected_at_time = state %in% c("U", "A", "D", "Tr")) %>%
    mutate(case_at_time = state %in% c("D", "Tr")) 
  
  return(df)
  
}