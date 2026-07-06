
# INPUT: just trial name
# OUTPUT: trial slug (short name) and (if not present already) folder structure


# DESCRIPTION:

# Functions to create consistent folder structures and file names

make_trial_slug <- function(trial_name) {
  
  gsub(" ", "_", tolower(trial_name))

}

make_output_dirs <- function() {
  
  dirs <- c(
    "outputs/cohort_data",
    "outputs/estimates/prevalence_incidence",
    "outputs/estimates/time_to_event",
    "outputs/estimates/effect_size",
    "outputs/plots/prevalence_incidence",
    "outputs/plots/time_to_event",
    "outputs/plots/effect_size",
    "outputs/plots/protocol_tests"
  )
  
  lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
}