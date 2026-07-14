
# INPUT: just trial name
# OUTPUT: trial slug (short name) and (if not present already) folder structure


# DESCRIPTION:

# Functions to create consistent folder structures and file names

make_trial_slug <- function(trial_name) {
  
  trial_name %>%
    trimws() %>%
    tolower() %>%
    gsub("[^a-z0-9]+", "_", .) %>%
    gsub("^_|_$", "", .)

}

make_output_dirs <- function() {
  
  dirs <- c(
    "outputs/cohort_data",
    "outputs/metadata",
    
    "outputs/estimates/prevalence_incidence",
    "outputs/estimates/time_to_event",
    "outputs/estimates/effect_size",
    "outputs/estimates/protocol_tests",
    
    "outputs/plots/cohort",
    "outputs/plots/prevalence_incidence",
    "outputs/plots/time_to_event",
    "outputs/plots/effect_size",
    "outputs/plots/protocol_tests"
  )
  
  lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
}