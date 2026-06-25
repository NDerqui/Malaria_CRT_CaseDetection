# functions/verbose_effect_plots.R

save_plot <- function(plot, filename, width = 12, height = 8) {
  png(filename = filename, width = width, height = height, units = "in", res = 1200)
  print(plot)
  dev.off()
}

plot_protective_effect <- function(protective_effect, key_intervention_time,
                                   sim_length, trial_title) {
  year <- 365
  
  ggplot(data = protective_effect,
         aes(x = timestep, y = effect)) +
    geom_point() + geom_line() +
    geom_vline(xintercept = key_intervention_time*year, color = "firebrick", linetype = "dashed") +
    scale_x_continuous(breaks = seq(0, sim_length * year, by = year),
                       labels = (0:sim_length)) +
    scale_y_continuous(labels = scales::percent, limits = 0:1) +
    labs(x = "Year", y = "Intervention Protective Effect",
         title = trial_title) +
    theme_bw() +
    theme(legend.position = "bottom", legend.title = element_blank()) +
    facet_nested(type_measure + measure ~ ., scales = "free")
}