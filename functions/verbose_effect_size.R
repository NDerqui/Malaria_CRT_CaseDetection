
# INPUT: trial outcome estimates (prevalence, incidence, etc.).
# OUTPUT: effect-size estimates comparing intervention with control.


# DESCRIPTION:

# Functions to estimate intervention effect sizes from trial outcomes:
# Prevalence, incidence and time-to-event.


# WRAP FOR PLOTTING AND EFFECT

# For ease later, use a function to wrap all the different outputs 
# (prevalence, incidence, whether true or observed by cross-sectional surveys, etc.),
# which is useful for plotting too.

wrap_for_plot_effect <- function(df) {
  
  require(dplyr)
  
  df <- df %>%
    select(-any_of(c("n", "person_days_at_risk", "infections", "cases",
                     "new_infections", "new_cases", "period", "period_label"))) %>%
    pivot_longer(-c(timestep, type_measure, run), names_to = "measure", values_to = "value") %>%
    mutate(measure = factor(measure, levels = measures, labels = measures_labels))
  
  return(df)
  
}

# RELATIVE EFFECT SIZE

# Getting the effect size of an intervention as: 1 - intervention/control.

get_relative_effect <- function(df,
                                value_col = "value", arm_col = "run",
                                control = "Control", intervention = "Intervention",
                                outcome_cols = c("type_measure", "measure", "timestep")) {

  require(dplyr)
  require(tidyr)
  
  # Make R recognise our columns for the numbers to compare and trial arms

  value_col <- rlang::sym(value_col)
  arm_col <- rlang::sym(arm_col)

  result <- df %>%
    # Select the cols we are interested in
    select(all_of(outcome_cols), arm = !!arm_col, value = !!value_col) %>%
    # Pivot to compare each
    tidyr::pivot_wider(names_from = arm, values_from = value) %>%
    mutate(
      effect = 1 - .data[[intervention]] / .data[[control]]
    ) %>%
    mutate(effect = as.numeric(effect))
  
  return(result)
}

# HAZARD RATIO EFFECT SIZE

# INPUT: a survival analysis df, usually after prepare_survival().
# OUTPUT: a one-row df with hazard ratio, confidence interval, and p-value.

get_hazard_ratio <- function(df,
                             time_col,
                             event_col,
                             arm_col = "run",
                             reference = "Control",
                             adjust_vars = NULL) {

  require(survival)

  df[[arm_col]] <- stats::relevel(as.factor(df[[arm_col]]), ref = reference)

  rhs <- paste(c(arm_col, adjust_vars), collapse = " + ")
  model_formula <- stats::as.formula(
    paste0("survival::Surv(", time_col, ", ", event_col, ") ~ ", rhs)
  )

  fit <- survival::coxph(model_formula, data = df)
  fit_summary <- summary(fit)

  arm_row <- grep(paste0("^", arm_col), rownames(fit_summary$coefficients))[1]

  return(data.frame(
    outcome = event_col,
    term = rownames(fit_summary$coefficients)[arm_row],
    hazard_ratio = fit_summary$coefficients[arm_row, "exp(coef)"],
    conf_low = fit_summary$conf.int[arm_row, "lower .95"],
    conf_high = fit_summary$conf.int[arm_row, "upper .95"],
    p_value = fit_summary$coefficients[arm_row, "Pr(>|z|)"],
    row.names = NULL
  ))
}
