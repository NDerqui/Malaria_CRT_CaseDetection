#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#-#-#-#-#-#-#-#-# Modelling Malaria CRT - Case Detection #-#-#-#-#-#-#-#-#-#-#-
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#


## Modelling case detection in malaria CRTs using verbose simulations:
## statistical analyses from analyses cohort (read from previous step).



# SET-UP ------------------------------------------------------------------


rm(list = ls())

# Packages

library(tidyverse)
library(ggpubr)

# Functions

source("functions/trial_tidy_outputs.R")
source("functions/trial_metadata.R")
source("functions/trial_analysis.R")
source("functions/trial_visualisation.R")

# General options

year <- 365
month <- 30



# TRIAL CONDITIONS --------------------------------------------------------


# Get the general conditions from the metadata, which we load with trial name

trial_name <- "High transmission test for workflow"
metadata <- load_trial_metadata(trial_name)

trial_id <- metadata$trial_id
trial_slug <- metadata$trial_slug

sim_length <- metadata$simulation$sim_length
trial_start <- metadata$trial$trial_start
trial_second_intervention <- metadata$trial$trial_second_intervention
key_intervention_time <- c(trial_start, trial_start+trial_second_intervention)

human_population <- metadata$simulation$human_population
trial_size <- metadata$analysis_cohort$trial_size

# Ensure we have output dirs

make_output_dirs()



# TRIAL OUTCOMES ----------------------------------------------------------


# Define a PCD cross-sectional survey and an ACD routine visit protocol

survey_protocol <- list(
  cross_surveys_in_years = seq(0.5, 6, 0.5)
)

acd_protocol <- list(
  routine_visits_in_weeks = seq(4, 6 * 52, 4),
  days_catchment = 2
)

# Run function for all estimates just from the trial name and protocols defined above

trial_results <- analyse_two_arm_trial(
  trial_id = trial_id,
  trial_start = trial_start,
  trial_second_intervention = trial_second_intervention,
  sim_length = sim_length,
  survey_protocol = survey_protocol,
  acd_protocol = acd_protocol
)
gc()



# SAVE OUTPUTS ------------------------------------------------------------


# Save all estimates and effect sizes to corresponding folder/file.

save_two_arm_trial(
  trial_results = trial_results,
  trial_slug = trial_slug
)

# Save all plots to corresponding folder/file.

save_two_arm_trial_plots(
  trial_results = trial_results,
  trial_slug = trial_slug,
  key_intervention_time = key_intervention_time,
  sim_length = sim_length,
  trial_title = trial_name
)



# PROTOCOLS ---------------------------------------------------------------


#### ACD visits ####

# Testing different protocols for ACD routine visits.
# Only re-runs estimate_acd_incidence() per iteration — data read once outside the loop.

source("functions/trial_outcomes.R")
source("functions/trial_effects.R")

visits <- c(2, 4, 8)
visits_period <- c(2, 4, 7)

# Read cohort data once
infections_proto <- if (file.exists(paste0("outputs/cohort_data/", trial_id, ".csv"))) {
  read.csv(paste0("outputs/cohort_data/", trial_id, ".csv"))
} else {
  readRDS(paste0("outputs/cohort_data/", trial_id, ".rds"))
}
infections_proto_control      <- infections_proto %>% filter(run == "Control")
infections_proto_intervention  <- infections_proto %>% filter(run == "Intervention")
rm(infections_proto); gc()

results_protocol_test <- data.frame()

for (visit_window in seq_along(visits)) {

  for (visit_period_window in seq_along(visits_period)) {

    acd_protocol <- list(
      routine_visits_in_weeks = seq(visits[visit_window], 6 * 52, visits[visit_window]),
      days_catchment = visits_period[visit_period_window]
    )

    bind_rows(
      inject(infections_proto_control     %>% estimate_acd_incidence(trial_start = trial_start, !!!acd_protocol)) %>% mutate(run = "Control"),
      inject(infections_proto_intervention %>% estimate_acd_incidence(trial_start = trial_start, !!!acd_protocol)) %>% mutate(run = "Intervention")
    ) %>%
      tidy_outcomes_for_effect() %>%
      estimate_relative_effect() %>%
      filter(!is.na(mean), grepl("ACD", type_measure)) %>%
      mutate(
        type_measure = "ACD visits",
        protocol = paste0(visits[visit_window], " weeks - ", visits_period[visit_period_window], " days window")
      ) -> protocol_result

    results_protocol_test <- rbind(results_protocol_test, protocol_result)
    gc()

  }
}

rm(infections_proto_control, infections_proto_intervention); gc()
write.csv(results_protocol_test,
          paste0("outputs/estimates/protocol_tests/", trial_slug, "_acd_protocol.csv"),
          row.names = FALSE)

png(filename = paste0("outputs/plots/protocol_tests/", trial_slug, "_acd_protocol.png"),
    width = 18, height = 8, units = "in", res = 1200)
ggplot(data = results_protocol_test,
       aes(x = timestep, y = mean, color = protocol)) +
  geom_errorbar(aes(ymin = lower_95quant, ymax = upper_95quant), width = 0.2) +
  geom_jitter(size = 1.5) +
  geom_vline(xintercept = key_intervention_time*year, color = "firebrick", linetype = "dashed") +
  scale_color_manual(values = c(carto_pal(name = "Safe")[c(9, 10, 2, 3, 8, 4, 7, 1, 11, 5, 12)], "black")) +
  scale_x_continuous(breaks = seq(0, sim_length * year, by = year),
                     labels = (0:sim_length)) +
  scale_y_continuous(labels = scales::percent, limits = 0:1) +
  labs(x = "Year", y = NULL,
       title = paste0(trial_name)) +
  theme_bw() + theme(legend.position = "bottom", legend.title = element_blank()) +
  facet_nested(type_measure + measure ~ analysis_population, scales = "free")
dev.off()
