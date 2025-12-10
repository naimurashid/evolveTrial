
# --- UPDATED: Fast interval metrics with robust event indexing ---
#' Calculate interval-specific metrics from patient data
#'
#' Recalculates events and person-time using lightweight base-R operations.
#'
#' @param patient_data Data frame with columns `observed_time` and `event_status`.
#' @param interval_cutpoints Numeric vector of interval boundaries.
#'
#' @return A list with `events_per_interval` and `person_time_per_interval` vectors.
#' @keywords internal
calculate_interval_metrics_fast <- function(patient_data, interval_cutpoints) {
  num_intervals <- length(interval_cutpoints) - 1L
  if (nrow(patient_data) == 0) {
    return(list(
      events_per_interval = rep(0L, num_intervals),
      person_time_per_interval = rep(0, num_intervals)
    ))
  }

  observed_times <- patient_data$observed_time
  event_status   <- patient_data$event_status

  # 1) Events per interval using tabulate (avoids data.table/table overhead)
  if (any(event_status == 1L)) {
    event_times <- observed_times[event_status == 1L]
    idx <- findInterval(
      event_times,
      interval_cutpoints,
      left.open = FALSE,        # [t_i, t_{i+1})
      rightmost.closed = FALSE
    )
    idx <- idx[idx > 0L]
    if (length(idx) > 0) {
      idx <- pmin(idx, num_intervals)
      events_per_interval <- tabulate(idx, nbins = num_intervals)
    } else {
      events_per_interval <- integer(num_intervals)
    }
  } else {
    events_per_interval <- integer(num_intervals)
  }

  # 2) Person-time per interval (left-closed, right-open)
  lower <- interval_cutpoints[-length(interval_cutpoints)]
  upper <- interval_cutpoints[-1L]
  person_time_per_interval <- numeric(num_intervals)
  if (length(observed_times) > 0) {
    for (i in seq_len(num_intervals)) {
      lo <- lower[i]
      hi <- upper[i]
      at_risk <- observed_times > lo
      if (any(at_risk)) {
        contrib <- pmin(observed_times[at_risk], hi) - lo
        # guard against tiny negative values from floating point error
        contrib[contrib < 0] <- 0
        person_time_per_interval[i] <- sum(contrib)
      }
    }
  }

  list(
    events_per_interval = events_per_interval,
    person_time_per_interval = person_time_per_interval
  )
}


# 1) ---------- State container ----------
make_state <- function(arm_names, max_total_patients_per_arm) {
  # PERFORMANCE: Pre-allocate registries with capacity to avoid rbind copies
  # Use 1.2x max capacity as buffer for safety
  # FIX: as.integer() strips names, so preserve them explicitly
  capacity <- as.integer(ceiling(max_total_patients_per_arm * 1.2))
  names(capacity) <- names(max_total_patients_per_arm)

  # Create pre-allocated registry for one arm
  create_preallocated_registry <- function(cap) {
    data.frame(
      id = rep(NA_integer_, cap),
      enroll_time = rep(NA_real_, cap),
      true_event_time = rep(NA_real_, cap),
      random_censor_time = rep(NA_real_, cap)
    )
  }

  list(
    arm_status = setNames(rep("recruiting", length(arm_names)), arm_names),
    enrolled_counts = setNames(rep(0L, length(arm_names)), arm_names),
    registries = setNames(
      lapply(arm_names, function(a) create_preallocated_registry(capacity[a])),
      arm_names
    ),
    # Track next available row index for each arm's registry (for O(1) insertion)
    registry_row_idx = setNames(rep(1L, length(arm_names)), arm_names),
    # per-simulation outputs that interims mutate:
    stop_efficacy_per_sim_row = setNames(rep(0L, length(arm_names)), arm_names),
    stop_futility_per_sim_row = setNames(rep(0L, length(arm_names)), arm_names),
    sim_final_n_current_run   = setNames(rep(NA_integer_, length(arm_names)), arm_names),
    # Track calendar time when each arm stopped (for ET calculation)
    stop_time = setNames(rep(NA_real_, length(arm_names)), arm_names)
  )
}


# Helper: Get active (non-NA) rows from pre-allocated registry
# PERFORMANCE: Use this to trim pre-allocated registries before processing
get_active_registry <- function(registry_df) {
  if (nrow(registry_df) == 0) return(registry_df)
  # Filter to only rows with valid id (non-NA)
  registry_df[!is.na(registry_df$id), , drop = FALSE]
}

# --- slice_arm_data_at_time (UPDATED: tiny clarity tweak on nrow check) ---
slice_arm_data_at_time <- function(registry_df, calendar_time, max_follow_up, interval_cutpoints) {
  # PERFORMANCE: Handle pre-allocated registries by trimming NA rows first
  registry_df <- get_active_registry(registry_df)
  if (nrow(registry_df) == 0) {
    return(list(
      patient_data = data.frame(observed_time = numeric(0), event_status = integer(0)),
      metrics = list(
        events_per_interval = rep(0, length(interval_cutpoints) - 1),
        person_time_per_interval = rep(0, length(interval_cutpoints) - 1),
        events_total = 0L,
        person_time_total = 0,
        median_followup = 0
      )
    ))
  }
  time_since_enroll <- pmax(0, calendar_time - registry_df$enroll_time)
  time_available    <- pmin(time_since_enroll, max_follow_up)
  observed_time <- pmin(registry_df$true_event_time,
                        registry_df$random_censor_time,
                        time_available)
  event_status  <- as.integer(registry_df$true_event_time <= pmin(registry_df$random_censor_time, time_available))
  keep <- time_available > 0
  pd <- data.frame(
    observed_time = observed_time[keep],
    event_status  = event_status[keep]
  )
  metrics <- calculate_interval_metrics_fast(pd, interval_cutpoints)
  metrics$events_total      <- sum(metrics$events_per_interval)
  metrics$person_time_total <- sum(metrics$person_time_per_interval)
  metrics$median_followup   <- if (nrow(pd) > 0) stats::median(pd$observed_time) else 0
  list(patient_data = pd, metrics = metrics)
}


# --- PER-ARM GATES FOR vs-ref (UPDATED) --------------------------------------

gates_pass_for_both_arms <- function(slCtrl, slTrt, args, diagnostics = FALSE) {
  # Arm names: take control from args$reference_arm_name; treatment is "the other one"
  if (is.null(args$reference_arm_name)) {
    stop("args$reference_arm_name must be provided for vs-reference gates.")
  }
  ctrl_name <- args$reference_arm_name

  if (is.null(args$arm_names) || length(args$arm_names) == 0) {
    stop("args$arm_names must be provided for vs-reference gates.")
  }
  arm_names <- args$arm_names

  trt_candidates <- arm_names[arm_names != ctrl_name]
  if (length(trt_candidates) == 0) {
    stop("No experimental arms found for vs-reference gates (after excluding reference arm).")
  }
  # Assuming only one treatment arm is being evaluated against the reference at a time
  trt_name <- trt_candidates[1]

  # thresholds (per-arm gates) - use shared helper from gate_diagnostics.R
  arms <- c(ctrl_name, trt_name)
  min_ev_vec <- resolve_gate_vec(
    args$min_events_per_arm,
    target_arms = arms,
    all_arm_names = arm_names,
    randomization_probs = args$randomization_probs,
    default = 0,
    scale = TRUE
  )
  min_mfu_vec <- resolve_gate_vec(
    args$min_median_followup_per_arm,
    target_arms = arms,
    all_arm_names = arm_names,
    randomization_probs = args$randomization_probs,
    default = 0,
    scale = FALSE
  )
  min_pt_vec <- resolve_gate_vec(
    args$min_person_time_frac_per_arm,
    target_arms = arms,
    all_arm_names = arm_names,
    randomization_probs = args$randomization_probs,
    default = 0,
    scale = TRUE
  )
  min_ev_ctrl <- min_ev_vec[[ctrl_name]] %||% 0
  min_ev_trt  <- min_ev_vec[[trt_name]] %||% 0
  min_mfu_ctrl <- min_mfu_vec[[ctrl_name]] %||% 0
  min_mfu_trt  <- min_mfu_vec[[trt_name]] %||% 0
  min_pt_ctrl  <- min_pt_vec[[ctrl_name]] %||% 0
  min_pt_trt   <- min_pt_vec[[trt_name]] %||% 0

  # extract metrics
  evC  <- coalesce_num(slCtrl$metrics$events_total, 0)
  evT  <- coalesce_num(slTrt$metrics$events_total, 0)
  mfuC <- coalesce_num(slCtrl$metrics$median_followup, 0)
  mfuT <- coalesce_num(slTrt$metrics$median_followup, 0)
  ptC  <- coalesce_num(slCtrl$metrics$person_time_total, 0)
  ptT  <- coalesce_num(slTrt$metrics$person_time_total, 0)

  # denominators (person-time caps)
  # Use max_trial_time for denominator if available; otherwise fall back to max_follow_up_sim
  # This ensures person-time gates are achievable when fu_time is set very high (e.g., 120 months)
  # but max_trial_time limits actual trial duration (e.g., 72 months)
  mtn <- args$max_total_patients_per_arm
  max_follow <- coalesce_num(args$max_trial_time, coalesce_num(args$max_follow_up_sim, 0))

  if (is.null(mtn)) {
    stop("args$max_total_patients_per_arm must be provided for vs-reference gates.")
  }
  if (is.null(names(mtn)) || !all(c(ctrl_name, trt_name) %in% names(mtn))) {
    stop("args$max_total_patients_per_arm must be a named vector containing all arm names.")
  }

  maxPT_C <- coalesce_num(mtn[[ctrl_name]], 0) * max_follow
  maxPT_T <- coalesce_num(mtn[[trt_name]],  0) * max_follow

  fracC <- if (maxPT_C > 0) ptC / maxPT_C else 0
  fracT <- if (maxPT_T > 0) ptT / maxPT_T else 0

  passC <- (evC >= min_ev_ctrl) && (mfuC >= min_mfu_ctrl) && (fracC >= min_pt_ctrl)
  passT <- (evT >= min_ev_trt)  && (mfuT >= min_mfu_trt)  && (fracT >= min_pt_trt)

  pass <- passC && passT

  if (isTRUE(diagnostics)) {
    message(sprintf(
      "[vsREF gate] %s: ev=%d mFU=%.2f PT=%.1f frac=%.3f  |  %s: ev=%d mFU=%.2f PT=%.1f frac=%.3f  -> pass=%s",
      ctrl_name, evC, mfuC, ptC, fracC,
      trt_name,  evT, mfuT, ptT, fracT,
      as.character(pass)
    ))
  }
  pass
}


# --- BETWEEN-ARM POSTERIOR PROBABILITY: FUTILITY (NEW/UPDATED) ---------------
# Proper *futility* probability with the correct directionality:
#   P(Triplet <= Doublet - margin)
# DO NOT use 1 - P(Triplet > Doublet + margin); with a margin this is *not* equivalent.


# Coalesce for numerics
coalesce_num <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x)) y else x
}


# --- BETWEEN-ARM POSTERIOR PROBABILITY: EFFICACY (UPDATED) -------------------
# P(Triplet > Doublet + margin), where " > " means whatever clinical estimand you use.
# This version expects each slice to expose posterior *draws* for a scalar estimand
# on the same scale between arms (e.g., median PFS, or -log(hazard), etc.)
# If your code stores piecewise rates, call your existing aggregator that maps
# draws -> scalar estimand per arm before comparing.
