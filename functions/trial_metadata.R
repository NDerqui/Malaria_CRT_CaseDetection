
# INPUT: trial name and basic parameters
# OUTPUT: trial metadata


# DESCRIPTION:

# This script is designed to generate metadata for a given trial based on its name and basic parameters. The metadata includes information such as the trial's unique identifier, start and end dates, status, and any relevant notes or comments.

create_trial_metadata <- function(trial_name, simulation, trial,
                                  intervention, analysis_cohort, n_power) {
  
  # Basic systems
  
  created_at <- Sys.time()
  trial_slug <- make_trial_slug(trial_name)
  timestamp <- format(created_at, "%Y%m%dT%H%M%S", tz = "UTC")
  
  list(
    trial_name = trial_name,
    trial_slug = trial_slug,
    trial_id = paste(trial_slug, timestamp, sep = "_"),
    created_at = created_at,
    simulation = simulation,
    trial = trial,
    intervention = intervention,
    analysis_cohort = analysis_cohort,
    n_power = n_power
  )
}

save_trial_metadata <- function(metadata) {
  
  # Ensure folder system is there
  
  make_output_dirs()
  
  # To write a record for eahc trial, name the file and create is neeeded
  
  registry_file <- "outputs/metadata/trial_registry.csv"
  
  registry <- if (file.exists(registry_file)) {
    read.csv(registry_file)
  } else {
    data.frame()
  }
  
  # Check where there is a trial with the same name
  
  if (nrow(registry) > 0 &&
      metadata$trial_name %in% registry$trial_name) {
    stop("trial_name already exists: ", metadata$trial_name)
  }
  
  # Save a row for each trial simulation to keep a registry
  
  new_entry <- data.frame(
    trial_name = metadata$trial_name,
    trial_slug = metadata$trial_slug,
    trial_id = metadata$trial_id,
    created_at = as.character(metadata$created_at)
  )
  
  write.csv(
    dplyr::bind_rows(registry, new_entry),
    registry_file,
    row.names = FALSE
  )
  
  # Save the metada as Rds
  
  saveRDS(
    metadata,
    file.path("outputs/metadata", paste0(metadata$trial_id, ".rds"))
  )
}

load_trial_metadata <- function(trial_name, trial_id = NULL) {
  
  # Read our trials registry and check for trial with that name
  
  registry <- read.csv("outputs/metadata/trial_registry.csv")
  
  matches <- registry %>%
    filter(.data$trial_name == .env$trial_name)
  
  # Unless specific trial_id is given, get the latest trial_id with same name and slug
  
  if (!is.null(trial_id)) {
    matches <- matches %>%
      filter(.data$trial_id == .env$trial_id)
  } else {
    matches <- matches %>%
      arrange(desc(as.POSIXct(created_at))) %>%
      slice(1)
  }
  
  # Check for non-exiting trial
  
  if (nrow(matches) != 1) {
    stop("No matching trial metadata found.")
  }
  
  # Load the metada with our selected trial_id
  
  readRDS(file.path(
    "outputs/metadata",
    paste0(matches$trial_id, ".rds")
  ))
}