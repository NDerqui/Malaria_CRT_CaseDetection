
# INPUT: trial outcome estimates (prevalence, incidence, etc.).
# OUTPUT: effect-size estimates comparing intervention with control.


# DESCRIPTION:

# Functions to estimate intervention effect sizes from trial outcomes:
# Prevalence, incidence and time-to-event.


# WRAP FOR PLOTTING AND EFFECT

# For ease later, use a function to wrap all the different outputs 
# (prevalence, incidence, whether true or observed by cross-sectional surveys, etc.),
# which is useful for plotting too.

tidy_outcomes_for_effect <- function(df,
                                 measures = c("prevalence_infection", "prevalence_case",
                                              "incidence_ppd_infection", "incidence_ppd_case",
                                              "incidence_ppy_infection", "incidence_ppy_case"),
                                 measures_labels = c("Infection Prevalence", "Case Prevalence",
                                                     "Infection Incidence p.p.day", "Case Incidence p.p.day",
                                                     "Infection Incidence p.p.year", "Case Incidence p.p.year"))
  {
  
  require(dplyr)
  require(tidyr)
  
  df <- df %>%
    select(-any_of(c("n", "person_days_at_risk", "infections", "cases",
                     "new_infections", "new_cases", "period", "period_label"))) %>%
    pivot_longer(-c(timestep, sim, type_measure, run),
                 names_to = "measure", values_to = "value") %>%
    mutate(measure = factor(measure, levels = measures, labels = measures_labels))
  
  return(df)
  
}

# RELATIVE EFFECT SIZE

# Getting the effect size of an intervention as: 1 - intervention/control.

estimate_relative_effect <- function(df,
                                value_col = "value", arm_col = "run",
                                control = "Control", intervention = "Intervention",
                                outcome_cols = c("type_measure", "measure", "timestep", "sim")) {

  require(dplyr)
  require(tidyr)
  require(rlang)
  
  # Make R recognise our columns for the numbers to compare and trial arms

  value_col <- rlang::sym(value_col)
  arm_col <- rlang::sym(arm_col)

  result <- df %>%
    # Select the cols we are interested in
    select(all_of(outcome_cols), arm = !!arm_col, value = !!value_col) %>%
    # Pivot to compare each
    tidyr::pivot_wider(names_from = arm, values_from = value) %>%
    mutate(
      effect = 1 - .data[[intervention]] / .data[[control]],
      effect = if_else(is.finite(effect), as.numeric(effect), NA_real_)
    ) %>%
    mutate(effect = as.numeric(effect))
  
  # Now that we have effect size by sim, get mean and quantiles across simulations
  
  result <- result %>%
    # Create mean and 95%CI across simulations
    group_by(across(all_of(outcome_cols[!(outcome_cols == "sim")]))) %>%
    mutate(mean = mean(effect, na.rm = TRUE),
           lower_95quant = quantile(effect, probs = 0.025, na.rm = TRUE),
           upper_95quant = quantile(effect, probs = 0.975, na.rm = TRUE)) %>%
    filter(row_number() == 1) %>% ungroup() %>%
    select(all_of(outcome_cols[!(outcome_cols == "sim")]), mean, lower_95quant, upper_95quant)
  
  return(result)
}

# HAZARD RATIO EFFECT SIZE

# INPUT: a survival analysis df, usually after estimate_true_time_to_event() and prepare_Survival().
# OUTPUT: a one-row df with hazard ratio, confidence interval, and p-value.

estimate_hazard_ratio <- function(df,
                             time_col, event_col,
                             arm_col = "run", control = "Control",
                             covariates = NULL) {

  require(survival)
  
  # Relevel our df so that our arm column has the control as reference
  df[[arm_col]] <- stats::relevel(as.factor(df[[arm_col]]), ref = control)
  
  # Make a vector with our arm_col + any other covariates we want to get HR from
  model_covars <- paste(c(arm_col, covariates), collapse = " + ")
  
  # Write our formula with our time to event col and covariates from before
  model_formula <- stats::as.formula(
    paste0("survival::Surv(", time_col, ", ", event_col, ") ~ ", model_covars))
  
  # Run our survival model
  fit <- survival::coxph(model_formula, data = df)
  fit_summary <- summary(fit)
  
  # Save an object

  arm_row <- grep(paste0("^", arm_col), rownames(fit_summary$coefficients))[1]

  return(data.frame(
    hazard_ratio = fit_summary$coefficients[arm_row, "exp(coef)"],
    conf_low = fit_summary$conf.int[arm_row, "lower .95"],
    conf_high = fit_summary$conf.int[arm_row, "upper .95"],
    p_value = fit_summary$coefficients[arm_row, "Pr(>|z|)"],
    row.names = NULL
  ) %>% # Add signaling for significance of the hazard ratio
    mutate(
      significant_any = conf_high < 1 | conf_low > 1,
      significant_benefit = conf_high < 1
    ))
}

# Above works when there are not multiple simulations,
# so use below wrapper when there are!

estimate_hazard_ratio_by_sim <- function(df,
                                         time_col, event_col,
                                         arm_col = "run", control = "Control",
                                         covariates = NULL) {
  require(dplyr)
  
  result <- df %>%
    group_by(sim) %>%
    group_modify(~ estimate_hazard_ratio(
      df = .x, time_col = time_col, event_col = event_col,
      arm_col = arm_col, control = control, covariates = covariates
    )) %>%
    ungroup()
  
  return(result)
}
