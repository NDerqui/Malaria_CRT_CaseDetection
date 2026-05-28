
# INPUT: trial outcome estimates with control and intervention arms.
# OUTPUT: effect-size estimates comparing intervention with control.

# DESCRIPTION:

# Functions to estimate intervention effect sizes from trial outcomes.
# Prevalence, incidence, and time-to-event functions estimate arm-specific
# outcomes first. These functions perform the final intervention comparison.


# RELATIVE EFFECT SIZE

# INPUT: a tidy outcome df with one value per arm and time/period.
# OUTPUT: a tidy df with intervention effect estimated as 1 - intervention/control.

get_relative_effect <- function(df,
                                value_col = "value",
                                arm_col = "run",
                                control = "Control",
                                intervention = "Intervention",
                                outcome_cols = c("measure"),
                                comparison_cols = NULL) {

  require(dplyr)
  require(tidyr)

  value_col <- rlang::sym(value_col)
  arm_col <- rlang::sym(arm_col)

  df %>%
    select(all_of(c(comparison_cols, outcome_cols)), arm = !!arm_col, value = !!value_col) %>%
    filter(arm %in% c(control, intervention)) %>%
    tidyr::pivot_wider(names_from = arm, values_from = value) %>%
    mutate(
      effect = 1 - .data[[intervention]] / .data[[control]]
    )
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

  data.frame(
    outcome = event_col,
    term = rownames(fit_summary$coefficients)[arm_row],
    hazard_ratio = fit_summary$coefficients[arm_row, "exp(coef)"],
    conf_low = fit_summary$conf.int[arm_row, "lower .95"],
    conf_high = fit_summary$conf.int[arm_row, "upper .95"],
    p_value = fit_summary$coefficients[arm_row, "Pr(>|z|)"],
    row.names = NULL
  )
}
