library(devtools)

load_all("../malariasimulation/")

library(individual)
library(malariasimulation)
library(ggplot2)

source("functions/verbose_set_parameters.R")
source("functions/verbose_runsim.R")

## Set up a couple of basic general options for simulations

year <- 365
month <- 30

start_time <- 3*year
init_EIR <- 15

init_EIRS <- exp(seq(log(1), log(100), length.out = 21))
init_EIRS <- c(20)
human_population <- 100000
trial_size <- 100


sim_length <- 23

trial_start <- 3
trial_second_intervention <- 3

key_intervention_time <- c(trial_start, trial_start+trial_second_intervention)

# Control when we get our age snapshot (best at start of trial, timestep = 1)

snapshot_time <- start_time

# Add a "trial name" to keep track of results

trial_name <- paste0("Seasonal_Init_EIR_", init_EIR)
# change_in_time_inf <- c()
# k <- 1
for(init_EIR in init_EIRS){
  
  baseline_parameters <- set_baseline_pars(
    sim_length = sim_length,
    init_EIR = init_EIR,
    human_population = human_population,
    seasonality = FALSE,
    treatment = TRUE,
    bednets = FALSE
    )
  
  ## Run sim
  
  # # Control arm
  out <- run_verbose_sim_pev_vaccines(simparams = baseline_parameters, sim_length = sim_length,
                         snapshot_time = snapshot_time, vaccine_bool = FALSE, run_note = "control")
  # 
  # # flop
  df_control <- read.csv("outputs_verbose_sims/verbose_dumping_control.csv")
  # 
  df_control$process <- out$process_vector[df_control$process_index]
  df_control$state <- out$state_list[df_control$state_index]
  # 
  # # Intervention arm
  
  
  out <- run_verbose_sim_pev_vaccines(
    simparams = baseline_parameters,
    sim_length = sim_length,
    snapshot_time = snapshot_time,
    vaccine_bool  = TRUE,
    run_note = "vaccine",
    profile = rtss_profile,
    vaccine_epi = TRUE,
    epi_timesteps = start_time,
    epi_coverages = 0.5,
    epi_min_wait = 6*month,
    epi_age = 6*month,
    epi_booster_spacing = year,
    epi_boost_coverage = matrix(0.5)
    )
  
  df_vaccine <- read.csv("outputs_verbose_sims/verbose_dumping_vaccine.csv")
  
  df_vaccine$process <- out$process_vector[df_vaccine$process_index]
  df_vaccine$state <- out$state_list[df_vaccine$state_index]
  
  df_vaccine_snapshot <- read.csv("outputs_verbose_sims/verbose_dumping_snapshot_vaccine.csv")
  df_control_snapshot <- read.csv("outputs_verbose_sims/verbose_dumping_snapshot_control.csv")
  
  df_vaccine_snapshot$state <- out$state_list[df_vaccine_snapshot$state_index]
  df_control_snapshot$state <- out$state_list[df_control_snapshot$state_index]
  
  vaccinated_people <- unique(df_vaccine[df_vaccine$process == "vaccinated_mass", "individual_index"])
  target_pop <- df_vaccine_snapshot$individual_index

  dft <- df_vaccine[df_vaccine$timestep >= start_time, ]
  time_till_first_infection <- c()
  for (i in seq(length(target_pop))){
    cat(i, "\r")
    person <- target_pop[i]
    states <- dft[dft$individual_index == person, "state"]
    for (j in seq(length(states))){
      if (states[j] != "S"){
        time_till_first_infection[i] <- j - 1
        break
      }
    }
    if (j == length(states)){
      time_till_first_infection[i] <- j
    }
  }
  
  num_first_born <- 10000
  df_vaccine_temp <- df_vaccine[df_vaccine$timestep >= start_time, ]
  first_born_vaccine <- df_vaccine_temp[df_vaccine_temp$process == "born", "individual_index"][1:num_first_born]
  df_vaccine_born <- df_vaccine[df_vaccine$individual_index %in% first_born_vaccine, ]
  df_vaccine_born$type <- "vaccinated"
  df_vaccine_born <- df_vaccine_born[df_vaccine_born$process == "state", ]
  
  df_control_temp <- df_control[df_control$timestep >= start_time, ]
  first_born_control <- df_control_temp[df_control_temp$process == "born", "individual_index"][1:num_first_born]
  df_control_born <- df_control[df_control$individual_index %in% first_born_control, ]
  df_control_born$type <- "control"
  df_control_born <- df_control_born[df_control_born$process == "state", ]
  
  df_plotting <- rbind(df_vaccine_born, df_control_born)
  
  ggplot(df_plotting) +
    geom_point(aes(x = individual_index, y = (timestep - start_time)/year, colour = state)) +
    facet_wrap(.~type) +
    theme_bw() +
    labs(x = "Unique Identifier", y = "Year since vaccine rollout")
  
  ggplot(df_plotting) +
    geom_point(aes(x = individual_index, y = (timestep - start_time)/year, colour = process)) +
    facet_wrap(.~type) +
    theme_bw() +
    labs(x = "Unique Identifier", y = "Year since vaccine rollout")
  
  time_infected_control <- length(df_control_born[df_control_born$state != "S", "individual_index"])
  time_infected_vaccine <- length(df_vaccine_born[df_vaccine_born$state != "S", "individual_index"])
  # change_in_time_inf[k] <- 100*(1 - time_infected_vaccine/time_infected_control)
  # k <- k + 1
  # times_infected_control <- c()
  # i <- 1
  # for(ind in first_born_control){
  #   dft <- df_control_born[df_control_born$individual_index == ind, ]
  #   time <- length(dft[dft$state != "S", "state"])
  #   times_infected_control[i] <- time
  #   i <- i + 1
  # }
  # 
  # times_infected_vaccine <- c()
  # i <- 1
  # for(ind in first_born_vaccine){
  #   dft <- df_vaccine_born[df_vaccine_born$individual_index == ind, ]
  #   time <- length(dft[dft$state != "S", "state"])
  #   times_infected_vaccine[i] <- time
  #   i <- i + 1
  # }
  # 
}

ggplot()+
  geom_point(aes(x = init_EIRS, y  = change_in_time_inf)) +
  geom_smooth(aes(x = init_EIRS, y  = change_in_time_inf)) +
  theme_bw() +
  labs(x = "Initial EIR value", y = "Percentage change in lifetime spent infected")

ages_inf_control_no_vacc <- c()
i <- 1
for (ind in first_born_control){
  dft <- df_control_born[df_control_born$individual_index == ind, ]
  born_time <- dft$timestep[1]
  dft <- dft[dft$timestep >= born_time + 6*month, ]
  first_inf_time <- dft[dft$state != "S", "timestep"][1]
  age_first_inf <- first_inf_time - born_time
  ages_inf_control_no_vacc[i] <- age_first_inf
  i <- i + 1
}

num_first_born <- 10000
df_vaccine_temp <- df_vaccine[df_vaccine$timestep >= start_time, ]
vaccinees <- unique(df_vaccine_temp[df_vaccine_temp$process == "vaccinated_epi", "individual_index"])
df_vaccine_temp <- df_vaccine_temp[df_vaccine_temp$individual_index %in% vaccinees, ]
first_born_vaccine_vaccinated <- df_vaccine_temp[df_vaccine_temp$process == "born", "individual_index"][1:num_first_born]
df_vaccine_vaccinated_born <- df_vaccine[df_vaccine$individual_index %in% first_born_vaccine_vaccinated, ]

df_vaccine_temp <- df_vaccine[df_vaccine$timestep >= start_time, ]
df_vaccine_temp <- df_vaccine_temp[!df_vaccine_temp$individual_index %in% vaccinees, ]
first_born_vaccine_not_vaccinated <- df_vaccine_temp[df_vaccine_temp$process == "born", "individual_index"][1:num_first_born]
df_vaccine_not_vaccinated_born <- df_vaccine[df_vaccine$individual_index %in% first_born_vaccine_not_vaccinated, ]

ages_inf_vaccinee_vacc <- c()
i <- 1
for (ind in first_born_vaccine_vaccinated){
  dft <- df_vaccine_vaccinated_born[df_vaccine_vaccinated_born$individual_index == ind, ]
  born_time <- dft$timestep[1]
  dft <- dft[dft$timestep >= born_time + 6*month, ]
  first_inf_time <- dft[dft$state != "S", "timestep"][1]
  age_first_inf <- first_inf_time - born_time
  ages_inf_vaccinee_vacc[i] <- age_first_inf
  i <- i + 1
}

ages_inf_not_vaccinee_vacc <- c()
i <- 1
for (ind in first_born_vaccine_not_vaccinated){
  dft <- df_vaccine_not_vaccinated_born[df_vaccine_not_vaccinated_born$individual_index == ind, ]
  born_time <- dft$timestep[1]
  dft <- dft[dft$timestep >= born_time + 6*month, ]
  first_inf_time <- dft[dft$state != "S", "timestep"][1]
  age_first_inf <- first_inf_time - born_time
  ages_inf_not_vaccinee_vacc[i] <- age_first_inf
  i <- i + 1
}

median(ages_inf_control_no_vacc/year, na.rm = TRUE)
median(ages_inf_not_vaccinee_vacc/year, na.rm = TRUE)
median(ages_inf_vaccinee_vacc/year, na.rm = TRUE)
sum(is.na(ages_inf_control_no_vacc))
sum(is.na(ages_inf_not_vaccinee_vacc))
sum(is.na(ages_inf_vaccinee_vacc))
