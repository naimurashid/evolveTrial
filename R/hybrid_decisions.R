#' hybrid_decisions.R
#' Decision rules for hybrid single-arm to between-arm trials
#'
#' Implements the decision logic for:
#' 1. Single-arm efficacy/futility decisions
#' 2. Conversion trigger evaluation
#' 3. Between-arm efficacy/futility decisions
#' 4. Overall trial conclusions

# ==============================================================================
# SINGLE-ARM DECISION RULES
# ==============================================================================

#' Evaluate single-arm efficacy for an arm
#'
#' @param p_single Posterior probability P(HR < c | data)
#' @param eff_threshold Efficacy threshold (default 0.90)
#' @param min_events Minimum events required (default 15)
#' @param current_events Current number of events
#'
#' @return List with decision and reason
#' @export
evaluate_sa_efficacy <- function(p_single, eff_threshold = 0.90,
                                  min_events = 15, current_events = 0) {

  if (is.na(p_single)) {
    return(list(
      efficacy = FALSE,
      reason = "Posterior probability not available"
    ))
  }

  if (current_events < min_events) {
    return(list(
      efficacy = FALSE,
      reason = sprintf("Insufficient events (%d < %d)", current_events, min_events)
    ))
  }

  if (p_single >= eff_threshold) {
    return(list(
      efficacy = TRUE,
      reason = sprintf("P(HR < c | data) = %.3f >= %.2f", p_single, eff_threshold)
    ))
  }

  list(
    efficacy = FALSE,
    reason = sprintf("P(HR < c | data) = %.3f < %.2f", p_single, eff_threshold)
  )
}

#' Evaluate single-arm futility for an arm
#'
#' @param p_single Posterior probability P(HR < c | data)
#' @param fut_threshold Futility threshold (default 0.10)
#' @param min_events Minimum events required (default 15)
#' @param current_events Current number of events
#'
#' @return List with decision and reason
#' @export
evaluate_sa_futility <- function(p_single, fut_threshold = 0.10,
                                  min_events = 15, current_events = 0) {

  if (is.na(p_single)) {
    return(list(
      futility = FALSE,
      reason = "Posterior probability not available"
    ))
  }

  if (current_events < min_events) {
    return(list(
      futility = FALSE,
      reason = sprintf("Insufficient events (%d < %d)", current_events, min_events)
    ))
  }

  if (p_single <= fut_threshold) {
    return(list(
      futility = TRUE,
      reason = sprintf("P(HR < c | data) = %.3f <= %.2f", p_single, fut_threshold)
    ))
  }

  list(
    futility = FALSE,
    reason = sprintf("P(HR < c | data) = %.3f > %.2f", p_single, fut_threshold)
  )
}

#' Apply futility action to trial state
#'
#' @param state Current trial state
#' @param arm Arm that hit futility
#' @param action Action to take: "stop_trial", "drop_arm", or "continue"
#'
#' @return Updated state and action taken
#' @export
apply_futility_action <- function(state, arm, action = "drop_arm") {

  action <- match.arg(action, c("stop_trial", "drop_arm", "continue"))

  result <- list(
    state = state,
    action_taken = action,
    arm = arm,
    trial_stopped = FALSE,
    arm_dropped = FALSE
  )

  switch(action,
    "stop_trial" = {
      state$current_state <- "STATE_STOP"
      state$trial_outcome <- "sa_futility_stop"
      state$stop_reason <- sprintf("Trial stopped due to futility in arm %s", arm)
      result$state <- state
      result$trial_stopped <- TRUE
    },
    "drop_arm" = {
      state$active_arms <- setdiff(state$active_arms, arm)
      state$dropped_arms <- c(state$dropped_arms, arm)
      result$state <- state
      result$arm_dropped <- TRUE
    },
    "continue" = {
      # No action needed, continue monitoring
      result$state <- state
    }
  )

  result
}

# ==============================================================================
# CONVERSION TRIGGER RULES
# ==============================================================================

#' Evaluate conversion trigger
#'
#' Checks if conditions are met to transition from SA to BA phase.
#'
#' @param sa_efficacy_reached Named logical vector of SA efficacy status per arm
#' @param active_arms Vector of active arm names
#' @param trigger Trigger type: "any_single_success", "all_single_success", or "k_of_K"
#' @param k_required Number required for k_of_K trigger
#'
#' @return List with triggered status and reason
#' @export
evaluate_conversion_trigger <- function(sa_efficacy_reached, active_arms,
                                         trigger = "any_single_success",
                                         k_required = 1) {

  trigger <- match.arg(trigger, c("any_single_success", "all_single_success", "k_of_K"))

  # Count efficacies among active arms
  efficacy_count <- sum(sa_efficacy_reached[active_arms], na.rm = TRUE)
  n_active <- length(active_arms)

  if (n_active == 0) {
    return(list(
      triggered = FALSE,
      reason = "No active arms remain"
    ))
  }

  result <- switch(trigger,
    "any_single_success" = {
      if (efficacy_count >= 1) {
        list(
          triggered = TRUE,
          reason = sprintf("At least 1 of %d active arms reached SA efficacy", n_active)
        )
      } else {
        list(
          triggered = FALSE,
          reason = sprintf("No active arms have reached SA efficacy (0 of %d)", n_active)
        )
      }
    },

    "all_single_success" = {
      if (efficacy_count == n_active) {
        list(
          triggered = TRUE,
          reason = sprintf("All %d active arms reached SA efficacy", n_active)
        )
      } else {
        list(
          triggered = FALSE,
          reason = sprintf("%d of %d active arms reached SA efficacy (need all)",
                           efficacy_count, n_active)
        )
      }
    },

    "k_of_K" = {
      if (efficacy_count >= k_required) {
        list(
          triggered = TRUE,
          reason = sprintf("%d of %d active arms reached SA efficacy (need >= %d)",
                           efficacy_count, n_active, k_required)
        )
      } else {
        list(
          triggered = FALSE,
          reason = sprintf("%d of %d active arms reached SA efficacy (need >= %d)",
                           efficacy_count, n_active, k_required)
        )
      }
    }
  )

  result$efficacy_count <- efficacy_count
  result$n_active <- n_active
  result$trigger_type <- trigger
  result
}

#' Select arms for between-arm comparison
#'
#' Determines which arms should proceed to BA phase.
#'
#' @param sa_efficacy_reached Named logical vector of SA efficacy status
#' @param active_arms Vector of active arm names
#' @param reference_arm Name of reference arm
#' @param selection_strategy Strategy: "all_successful", "best", or "reference_plus_best"
#'
#' @return Vector of arm names for BA phase
#' @export
select_arms_for_ba <- function(sa_efficacy_reached, active_arms, reference_arm,
                                selection_strategy = "all_successful") {

  selection_strategy <- match.arg(selection_strategy,
                                   c("all_successful", "best", "reference_plus_best"))

  # Ensure reference arm is included
  arms_to_include <- reference_arm

  switch(selection_strategy,
    "all_successful" = {
      # Include all arms that achieved SA efficacy
      successful <- names(sa_efficacy_reached[sa_efficacy_reached])
      arms_to_include <- union(arms_to_include, intersect(successful, active_arms))
    },

    "best" = {
      # Include only reference and the single best experimental arm
      # For now, take any successful experimental arm
      experimental <- setdiff(active_arms, reference_arm)
      successful_exp <- intersect(experimental, names(sa_efficacy_reached[sa_efficacy_reached]))
      if (length(successful_exp) > 0) {
        arms_to_include <- c(reference_arm, successful_exp[1])
      }
    },

    "reference_plus_best" = {
      # Same as "best" - reference plus best experimental
      experimental <- setdiff(active_arms, reference_arm)
      successful_exp <- intersect(experimental, names(sa_efficacy_reached[sa_efficacy_reached]))
      if (length(successful_exp) > 0) {
        arms_to_include <- c(reference_arm, successful_exp[1])
      }
    }
  )

  arms_to_include
}

# ==============================================================================
# BETWEEN-ARM DECISION RULES
# ==============================================================================

#' Evaluate between-arm efficacy
#'
#' @param p_between Posterior probability P(HR_AB < 1 | data)
#' @param eff_threshold Efficacy threshold (default 0.975)
#' @param min_events Minimum total events required (default 30)
#' @param current_events Current total events across arms
#'
#' @return List with decision and reason
#' @export
evaluate_ba_efficacy <- function(p_between, eff_threshold = 0.975,
                                  min_events = 30, current_events = 0) {

  if (is.na(p_between)) {
    return(list(
      efficacy = FALSE,
      reason = "Between-arm posterior probability not available"
    ))
  }

  if (current_events < min_events) {
    return(list(
      efficacy = FALSE,
      reason = sprintf("Insufficient total events (%d < %d)", current_events, min_events)
    ))
  }

  if (p_between >= eff_threshold) {
    return(list(
      efficacy = TRUE,
      reason = sprintf("P(HR_AB < 1 | data) = %.4f >= %.3f", p_between, eff_threshold)
    ))
  }

  list(
    efficacy = FALSE,
    reason = sprintf("P(HR_AB < 1 | data) = %.4f < %.3f", p_between, eff_threshold)
  )
}

#' Evaluate between-arm futility
#'
#' @param p_between Posterior probability P(HR_AB < 1 | data)
#' @param fut_threshold Futility threshold (default 0.05)
#' @param min_events Minimum total events required (default 30)
#' @param current_events Current total events across arms
#'
#' @return List with decision and reason
#' @export
evaluate_ba_futility <- function(p_between, fut_threshold = 0.05,
                                  min_events = 30, current_events = 0) {

  if (is.na(p_between)) {
    return(list(
      futility = FALSE,
      reason = "Between-arm posterior probability not available"
    ))
  }

  if (current_events < min_events) {
    return(list(
      futility = FALSE,
      reason = sprintf("Insufficient total events (%d < %d)", current_events, min_events)
    ))
  }

  if (p_between <= fut_threshold) {
    return(list(
      futility = TRUE,
      reason = sprintf("P(HR_AB < 1 | data) = %.4f <= %.3f", p_between, fut_threshold)
    ))
  }

  list(
    futility = FALSE,
    reason = sprintf("P(HR_AB < 1 | data) = %.4f > %.3f", p_between, fut_threshold)
  )
}

# ==============================================================================
# OVERALL TRIAL CONCLUSION
# ==============================================================================

#' Classify final trial outcome
#'
#' @param state Final trial state
#' @param theta Hybrid design parameters
#'
#' @return List with classification details
#' @export
classify_trial_outcome <- function(state, theta) {

  outcome <- state$trial_outcome

  classification <- list(
    outcome = outcome,
    phase_stopped = NA_character_,
    success = FALSE,
    futility = FALSE,
    max_n = FALSE,
    converted = state$conversion_decision == "GO",
    sa_success = any(state$sa_efficacy_reached),
    ba_success = state$ba_efficacy_reached
  )

  # Classify by outcome type
  if (grepl("^sa_", outcome)) {
    classification$phase_stopped <- "single_arm"
    if (outcome == "sa_futility_stop") {
      classification$futility <- TRUE
    }
  } else if (grepl("^ba_", outcome)) {
    classification$phase_stopped <- "between_arm"
    if (outcome == "ba_efficacy" || outcome == "ba_efficacy_at_max") {
      classification$success <- TRUE
    } else if (outcome == "ba_futility") {
      classification$futility <- TRUE
    }
  } else if (grepl("^conversion_", outcome)) {
    classification$phase_stopped <- "conversion"
    if (outcome == "conversion_nogo" || outcome == "conversion_ambiguous") {
      # SA was successful but no BA needed
      if (any(state$sa_efficacy_reached)) {
        classification$success <- TRUE  # Count SA success as overall success
      }
    }
  } else if (grepl("max_n", outcome)) {
    classification$max_n <- TRUE
    if (grepl("single", outcome)) {
      classification$phase_stopped <- "single_arm"
    } else {
      classification$phase_stopped <- "between_arm"
    }
  } else if (outcome == "all_arms_futile") {
    classification$phase_stopped <- "single_arm"
    classification$futility <- TRUE
  }

  classification
}

#' Determine if trial concluded with overall success
#'
#' Success is defined as:
#' - SA efficacy reached AND conversion not needed (NO_GO due to already effective)
#' - OR BA efficacy reached
#'
#' @param state Final trial state
#' @param theta Hybrid design parameters
#'
#' @return Logical
#' @export
is_trial_success <- function(state, theta) {

  # BA efficacy
  if (isTRUE(state$ba_efficacy_reached)) {
    return(TRUE)
  }

  # SA efficacy with conversion not proceeding
  # (implies SA conclusion was sufficient)
  if (any(state$sa_efficacy_reached) &&
      state$conversion_decision %in% c("NO_GO", "AMBIGUOUS")) {
    return(TRUE)
  }

  FALSE
}

#' Compute operating characteristics for a single trial
#'
#' @param state Final trial state
#' @param theta Hybrid design parameters
#' @param true_effect True effect scenario (for labeling)
#'
#' @return Data frame with operating characteristics
#' @export
compute_trial_oc <- function(state, theta, true_effect = "unknown") {

  classification <- classify_trial_outcome(state, theta)

  data.frame(
    true_effect = true_effect,
    outcome = state$trial_outcome,
    success = is_trial_success(state, theta),
    futility = classification$futility,
    max_n_reached = classification$max_n,
    converted = classification$converted,
    sa_efficacy = any(state$sa_efficacy_reached),
    sa_futility = any(state$sa_futility_reached),
    ba_efficacy = state$ba_efficacy_reached %||% FALSE,
    ba_futility = state$ba_futility_reached %||% FALSE,
    total_n = sum(state$n_enrolled),
    n_phase1 = sum(state$n_enrolled_phase1),
    n_phase2 = sum(state$n_enrolled) - sum(state$n_enrolled_phase1),
    duration = state$current_time,
    pp_at_conversion = state$pp_at_conversion %||% NA_real_,
    n_add_selected = state$n_add_selected %||% NA_integer_,
    stringsAsFactors = FALSE
  )
}

# ==============================================================================
# DECISION SUMMARY FUNCTIONS
# ==============================================================================

#' Generate human-readable decision summary
#'
#' @param state Final trial state
#' @param theta Hybrid design parameters
#'
#' @return Character string with trial summary
#' @export
generate_decision_summary <- function(state, theta) {

  lines <- character()

  # Header
  lines <- c(lines, "=== Hybrid Trial Decision Summary ===")
  lines <- c(lines, "")

  # Phase 1 summary
  lines <- c(lines, "PHASE 1 (Single-Arm):")
  for (arm in state$arm_names) {
    status <- if (state$sa_efficacy_reached[arm]) {
      "EFFICACY"
    } else if (state$sa_futility_reached[arm]) {
      "FUTILITY"
    } else if (arm %in% state$dropped_arms) {
      "DROPPED"
    } else {
      "INCONCLUSIVE"
    }
    p_val <- state$sa_posterior_prob[arm]
    p_str <- if (is.na(p_val)) "N/A" else sprintf("%.3f", p_val)
    lines <- c(lines, sprintf("  %s: %s (P = %s)", arm, status, p_str))
  }
  lines <- c(lines, sprintf("  Enrolled (Phase 1): %d", sum(state$n_enrolled_phase1)))
  lines <- c(lines, "")

  # Conversion summary
  lines <- c(lines, "CONVERSION:")
  lines <- c(lines, sprintf("  Decision: %s", state$conversion_decision %||% "N/A"))
  if (!is.na(state$pp_at_conversion)) {
    lines <- c(lines, sprintf("  PP at conversion: %.3f", state$pp_at_conversion))
  }
  if (!is.na(state$n_add_selected)) {
    lines <- c(lines, sprintf("  N_add selected: %d per arm", state$n_add_selected))
  }
  lines <- c(lines, "")

  # Phase 2 summary (if applicable)
  if (state$conversion_decision == "GO") {
    lines <- c(lines, "PHASE 2 (Between-Arm):")
    ba_status <- if (state$ba_efficacy_reached) {
      "EFFICACY"
    } else if (state$ba_futility_reached) {
      "FUTILITY"
    } else {
      "INCONCLUSIVE"
    }
    lines <- c(lines, sprintf("  Decision: %s", ba_status))
    if (!is.na(state$ba_posterior_prob)) {
      lines <- c(lines, sprintf("  P(HR < 1 | data): %.4f", state$ba_posterior_prob))
    }
    n_phase2 <- sum(state$n_enrolled) - sum(state$n_enrolled_phase1)
    lines <- c(lines, sprintf("  Enrolled (Phase 2): %d", n_phase2))
    lines <- c(lines, "")
  }

  # Final outcome
  lines <- c(lines, "FINAL:")
  lines <- c(lines, sprintf("  Outcome: %s", state$trial_outcome))
  lines <- c(lines, sprintf("  Reason: %s", state$stop_reason))
  lines <- c(lines, sprintf("  Total N: %d", sum(state$n_enrolled)))
  lines <- c(lines, sprintf("  Overall Success: %s",
                             if (is_trial_success(state, theta)) "YES" else "NO"))

  paste(lines, collapse = "\n")
}

#' Create structured decision report
#'
#' @param state Final trial state
#' @param theta Hybrid design parameters
#'
#' @return List with structured report components
#' @export
create_decision_report <- function(state, theta) {

  list(
    # Trial identification
    trial_id = state$trial_id %||% NA_character_,
    timestamp = Sys.time(),

    # Phase 1 results
    phase1 = list(
      arms = state$arm_names,
      active_arms = state$active_arms,
      dropped_arms = state$dropped_arms,
      sa_efficacy = state$sa_efficacy_reached,
      sa_futility = state$sa_futility_reached,
      sa_posterior_prob = state$sa_posterior_prob,
      n_enrolled = state$n_enrolled_phase1
    ),

    # Conversion results
    conversion = list(
      evaluated = state$conversion_evaluated,
      decision = state$conversion_decision,
      pp_at_decision = state$pp_at_conversion,
      n_add_selected = state$n_add_selected,
      pp_curve = state$pp_curve
    ),

    # Phase 2 results
    phase2 = list(
      active = state$conversion_decision == "GO",
      ba_efficacy = state$ba_efficacy_reached,
      ba_futility = state$ba_futility_reached,
      ba_posterior_prob = state$ba_posterior_prob,
      n_enrolled = sum(state$n_enrolled) - sum(state$n_enrolled_phase1)
    ),

    # Final summary
    final = list(
      outcome = state$trial_outcome,
      stop_reason = state$stop_reason,
      success = is_trial_success(state, theta),
      total_n = sum(state$n_enrolled),
      duration = state$current_time
    ),

    # Design parameters used
    design = theta
  )
}

# ==============================================================================
# DECISION COMPARISON UTILITIES
# ==============================================================================

#' Compare SA-only vs Hybrid decision
#'
#' For a given trial, determine what decision would have been made
#' under SA-only vs hybrid design.
#'
#' @param state Final hybrid trial state
#' @param theta Hybrid design parameters
#'
#' @return List with comparison
#' @export
compare_sa_vs_hybrid <- function(state, theta) {

  # SA-only decision: based purely on SA efficacy
  sa_decision <- if (any(state$sa_efficacy_reached)) {
    "SUCCESS"
  } else if (any(state$sa_futility_reached)) {
    "FUTILITY"
  } else {
    "INCONCLUSIVE"
  }

  # Hybrid decision
  hybrid_success <- is_trial_success(state, theta)
  hybrid_decision <- if (hybrid_success) {
    "SUCCESS"
  } else if (state$ba_futility_reached || all(state$sa_futility_reached)) {
    "FUTILITY"
  } else {
    "INCONCLUSIVE"
  }

  list(
    sa_only_decision = sa_decision,
    hybrid_decision = hybrid_decision,
    decisions_agree = sa_decision == hybrid_decision,
    hybrid_converted = state$conversion_decision == "GO",
    n_saved = if (!state$conversion_decision == "GO") {
      theta$nmax_ba * 2 - sum(state$n_enrolled)
    } else {
      0
    }
  )
}

# ==============================================================================
# NULL COALESCING OPERATOR
# ==============================================================================

if (!exists("%||%", mode = "function")) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
}
