
# INPUT: a df result of a verbose simulation, or from an age-cohort df.
# OUTPUT: a plot tracking each individual's state over time.

# DESCRIPTION:

# Functions to do some basic plots of verbose sims, to understand
# what we can track (individual state over time).

plot_verbose <- function(df, sim_length, human_population, note) {
  
  month <- 30
  year <- 365
  
  dir.create("outputs_plots/", , showWarnings = FALSE)
  
  require(ggplot2)
  require(rcartocolor)
  
  # Colors 
  
  malaria_sim_colors <- carto_pal(name = "Safe")[c(4, 11, 1, 10, 2)]
  malaria_sim_states <- c("S", "A", "U", "D", "Tr")
  
  ## Heatmap-like plot
  
  plot <- ggplot(df, aes(x = timestep, y = individual_index, fill = state)) +
    geom_tile() +
    scale_x_continuous(breaks = seq(0, sim_length * year, by = year),
                       labels = (0:sim_length)) +
    scale_fill_manual(breaks = malaria_sim_states,
                      values = malaria_sim_colors) +
    labs(x = "Year", y = NULL,
         title = paste0(human_population, " ppl - ", sim_length, " years - ", note)) +
    theme_bw() + theme(legend.position = "bottom", legend.title = element_blank())
  
  return(plot)
  
}

plot_verbose_itn <- function(df, sim_length, human_population, bednetstimesteps, note) {
  
  month <- 30
  year <- 365
  
  require(ggplot2)
  require(rcartocolor)
  
  # Colors 
  
  malaria_sim_colors <- carto_pal(name = "Safe")[c(4, 11, 1, 10, 2)]
  malaria_sim_states <- c("S", "A", "U", "D", "Tr")
  
  ## Heatmap-like plot

  plot <- ggplot(df, aes(x = timestep, y = individual_index, fill = state)) +
    geom_vline(xintercept = bednetstimesteps, linetype = "dashed", color = "grey") +
    geom_tile() +
    scale_x_continuous(breaks = seq(0, sim_length * year, by = year),
                       labels = (0:sim_length)) +
    scale_fill_manual(breaks = malaria_sim_states,
                      values = malaria_sim_colors) +
    labs(x = "Year", y = NULL,
         title = paste0(human_population, " ppl - ", sim_length, " years - ", note)) +
    theme_bw() + theme(legend.position = "bottom", legend.title = element_blank())
  
  return(plot)
  
}