
# This function can be run at any time after a verbose sim:
# before or after getting the age cohort.

# The function subsets timesteps to those "observed" in a cross sectional survey.
# Idea being to then apply any of the other functions to detect infection or clinical cases.

# The period argument allows you to "observe" events in xx days prior to survey,
# e.g. "Did you have a fever today or in the last 7 days?"

cross_survey <- function(df, survey_time, period) {
  
  df <- df %>%
    # Filter to that timestep and the period included in survey
    filter(timestep %in% (survey_time-period):survey_time)
  
  return(df)
}