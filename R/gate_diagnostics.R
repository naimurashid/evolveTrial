# --- INTERNAL HELPER: Resolve gate vectors with optional proportional scaling ---

#' Resolve gate parameter vector for specified arms with optional scaling
#'
#' Internal helper that converts gate parameters (which may be scalar, named vector,
#' or positional vector) into a named numeric vector for the specified arms.
#' Optionally applies proportional scaling based on randomization probabilities.
#'
#' @param raw Raw gate parameter value (scalar, named vector, or positional vector).
#' @param target_arms Character vector of arm names to resolve gates for.
#' @param all_arm_names Character vector of all arm names in the trial (for positional matching).
#' @param randomization_probs Optional named numeric vector of randomization probabilities.
#' @param default Default value if `raw` is NULL (default 0).
#' @param scale Logical; if TRUE and `raw` is scalar/NULL, apply proportional scaling
#'   by dividing each arm's randomization probability by the maximum probability (default FALSE).
#'
#' @return Named numeric vector with one element per arm in `target_arms`.
#' @keywords internal
resolve_gate_vec <- function(raw, target_arms, all_arm_names,
                              randomization_probs = NULL,
                              default = 0, scale = FALSE) {
  raw_len <- if (is.null(raw)) 0 else length(raw)

  # Step 1: Convert raw input to named vector for target arms
  if (is.null(raw)) {
    out <- rep(default, length(target_arms))
  } else if (!is.null(names(raw))) {
    # Named vector: extract by name
    out <- raw[target_arms]
    out[is.na(out)] <- default
  } else if (length(raw) == length(all_arm_names)) {
    # Positional vector matching all arms: name it then extract
    out <- raw
    names(out) <- all_arm_names
    out <- out[target_arms]
  } else if (length(raw) == length(target_arms)) {
    # Vector already matches target length
    out <- raw
  } else {
    # Scalar or mismatched length: recycle first element
    out <- rep(raw[1], length(target_arms))
  }
  names(out) <- target_arms

  # Step 2: Apply proportional scaling if requested
  should_scale <- scale && (raw_len <= 1 || is.null(raw))
  if (should_scale) {
    probs <- randomization_probs
    if (is.null(probs)) {
      probs <- rep(1 / length(all_arm_names), length(all_arm_names))
      names(probs) <- all_arm_names
    } else if (is.null(names(probs)) && length(probs) == length(all_arm_names)) {
      names(probs) <- all_arm_names
    }

    scale_vec <- rep(1, length(target_arms))
    names(scale_vec) <- target_arms
    if (!is.null(names(probs))) {
      max_prob <- max(probs, na.rm = TRUE)
      if (is.finite(max_prob) && max_prob > 0) {
        scale_vec <- probs[target_arms] / max_prob
        scale_vec[!is.finite(scale_vec)] <- 1
      }
    }
    out <- out * scale_vec
  }

  out
}


#' Estimate when the vs-reference interim gates can be satisfied
#'
#' This helper provides rough lower-bound heuristics for when the
#' between-arm (vs-reference) interim gating criteria in
#' [`interim_check()`] can all be satisfied under deterministic accrual.
#' It computes calendar-time approximations for the event, median
#' follow-up, and person-time requirements on a per-arm basis and reports
#' the largest of those times as a conservative bound on the earliest
#' informative interim look.
#'
#' The calculations assume constant accrual with rate
#' `overall_accrual_rate` split according to `randomization_probs`.
#' Person-time requirements are approximated using the relationship
#' PT ≈ rate × time² / 2 under steady accrual, while the median follow-up
#' requirement is approximated using the rule-of-thumb that the median
#' follow-up cannot exceed roughly half of the calendar time under uniform
#' accrual.
#'
#' @param args A named list following the structure passed into
#'   `run_simulation_pure()` / `run_scenarios()`.  The function uses the
#'   entries relevant to vs-reference gating, namely
#'   `arm_names`, `reference_arm_name`, `overall_accrual_rate`,
#'   `randomization_probs`, `max_total_patients_per_arm`,
#'   `max_follow_up_sim`, `min_events_per_arm`,
#'   `min_median_followup_per_arm`, and
#'   `min_person_time_frac_per_arm`.
#'
#' @return A list with two components:
#'   * `per_arm`: data frame containing the per-arm accrual rate and the
#'     heuristic lower-bound times (in months) for each gating component.
#'   * `joint_lower_bound`: the maximum of the per-arm lower bounds,
#'     representing the earliest calendar time at which all gates could
#'     plausibly be satisfied simultaneously under the heuristics.
#'
#' @examples
#' args <- list(
#'   arm_names = c("Doublet", "Triplet"),
#'  reference_arm_name = "Doublet",
#'   overall_accrual_rate = 3,
#'   randomization_probs = c(Doublet = 0.5, Triplet = 0.5),
#'   max_total_patients_per_arm = c(Doublet = 70, Triplet = 70),
#'   max_follow_up_sim = 24,
#'   min_events_per_arm = 8,
#'   min_median_followup_per_arm = 3,
#'   min_person_time_frac_per_arm = 0.15
#' )
#' estimate_vsref_gate_timing(args)
#'
#' @export
estimate_vsref_gate_timing <- function(args) {
  if (!is.list(args)) {
    stop("`args` must be a list of design parameters.")
  }

  arm_names <- args$arm_names
  if (is.null(arm_names) || length(arm_names) < 2) {
    stop("`args$arm_names` must contain at least two treatment arms.")
  }

  ref_arm <- args$reference_arm_name
  if (is.null(ref_arm) || !ref_arm %in% arm_names) {
    stop("`args$reference_arm_name` must name one of the arms.")
  }

  exp_arms <- setdiff(arm_names, ref_arm)
  if (length(exp_arms) == 0) {
    stop("At least one experimental arm is required for the vs-reference path.")
  }

  target_arms <- c(ref_arm, exp_arms)

  overall_rate <- args$overall_accrual_rate
  if (!is.numeric(overall_rate) || length(overall_rate) != 1 || !is.finite(overall_rate) || overall_rate <= 0) {
    stop("`overall_accrual_rate` must be a positive scalar numeric value.")
  }

  rand_probs <- args$randomization_probs
  if (is.null(rand_probs)) {
    rand_probs <- rep(1 / length(arm_names), length(arm_names))
    names(rand_probs) <- arm_names
  }
  if (is.null(names(rand_probs))) {
    if (length(rand_probs) != length(arm_names)) {
      stop("`randomization_probs` must either be named or match the arm ordering.")
    }
    names(rand_probs) <- arm_names
  }
  rand_probs <- rand_probs[target_arms]

  per_arm_rate <- overall_rate * rand_probs
  if (any(!is.finite(per_arm_rate)) || any(per_arm_rate <= 0)) {
    stop("Derived per-arm accrual rates must be positive and finite.")
  }

  resolve_per_arm <- function(x, default = 0) {
    if (is.null(x)) {
      return(rep(default, length(target_arms)))
    }
    if (length(x) == 1 && is.null(names(x))) {
      return(rep(x, length(target_arms)))
    }
    if (!is.null(names(x))) {
      vals <- x[target_arms]
      vals[is.na(vals)] <- default
      return(as.numeric(vals))
    }
    if (length(x) == length(target_arms)) {
      return(as.numeric(x))
    }
    stop("Unable to recycle values across arms for one of the gating parameters.")
  }

  min_events <- resolve_per_arm(args$min_events_per_arm, default = 0)
  min_mfu    <- resolve_per_arm(args$min_median_followup_per_arm, default = 0)
  min_pt_frac <- resolve_per_arm(args$min_person_time_frac_per_arm, default = 0)
  max_total <- resolve_per_arm(args$max_total_patients_per_arm)
  max_follow <- args$max_follow_up_sim
  if (is.null(max_follow) || !is.finite(max_follow) || max_follow < 0) {
    max_follow <- 0
  }

  pt_targets <- max_total * max_follow * min_pt_frac

  safe_div <- function(num, denom) {
    out <- rep(0, length(num))
    valid <- denom > 0
    out[valid] <- num[valid] / denom[valid]
    out
  }

  time_for_events <- safe_div(min_events, per_arm_rate)
  time_for_mfu    <- ifelse(min_mfu > 0, 2 * min_mfu, 0)
  time_for_pt     <- ifelse(pt_targets > 0, sqrt(2 * pt_targets / per_arm_rate), 0)

  per_arm_summary <- data.frame(
    Arm = target_arms,
    AccrualRate = per_arm_rate,
    MinEvents = min_events,
    TimeForEvents = time_for_events,
    MinMedianFollowup = min_mfu,
    TimeForMedianFollowup = time_for_mfu,
    MinPersonTimeMonths = pt_targets,
    TimeForPersonTime = time_for_pt,
    stringsAsFactors = FALSE
  )

  joint_lower_bound <- max(per_arm_summary$TimeForEvents,
                            per_arm_summary$TimeForMedianFollowup,
                            per_arm_summary$TimeForPersonTime,
                            na.rm = TRUE)

  list(
    per_arm = per_arm_summary,
    joint_lower_bound = joint_lower_bound
  )
}
