
# --- UPDATED: Fast interval metrics with robust event indexing ---
calculate_interval_metrics_fast <- function(patient_data, interval_cutpoints) {
  #' Recalculates events and person-time using a more efficient data.table approach.
  num_intervals <- length(interval_cutpoints) - 1
  if (nrow(patient_data) == 0) {
    return(list(
      events_per_interval = rep(0, num_intervals),
      person_time_per_interval = rep(0, num_intervals)
    ))
  }
  
  dt <- as.data.table(patient_data)
  results_template <- data.table(interval_num = 1:num_intervals)
  
  # 1) Events per interval with bullet-proof indexing
  if (nrow(dt[event_status == 1]) > 0) {
    ev_times <- dt[event_status == 1, observed_time]
    idx <- findInterval(
      ev_times,
      interval_cutpoints,
      left.open = FALSE,        # [t_i, t_{i+1})
      rightmost.closed = FALSE
    )
    # clamp into 1..num_intervals; drop zeros (events exactly at time 0)
    idx <- idx[idx > 0]
    if (length(idx) > 0) {
      idx <- pmin(pmax(idx, 1L), num_intervals)
      event_counts <- as.data.table(table(factor(idx, levels = 1:num_intervals)))
      setnames(event_counts, c("interval_num", "events"))
      event_counts[, interval_num := as.integer(as.character(interval_num))]
    } else {
      event_counts <- data.table(interval_num = 1:num_intervals, events = 0L)
    }
  } else {
    event_counts <- data.table(interval_num = 1:num_intervals, events = 0L)
  }
  
  # 2) Person-time per interval (left-closed, right-open)
  pt_list <- lapply(1:num_intervals, function(i) {
    lower_bound <- interval_cutpoints[i]
    upper_bound <- interval_cutpoints[i + 1]
    at_risk_dt <- dt[observed_time >= lower_bound]
    if (nrow(at_risk_dt) == 0) {
      return(data.table(interval_num = i, person_time = 0.0))
    }
    time_spent <- pmin(at_risk_dt$observed_time, upper_bound) - lower_bound
    time_spent[time_spent < 0] <- 0
    data.table(interval_num = i, person_time = sum(time_spent))
  })
  pt_summary <- rbindlist(pt_list)
  
  # 3) Merge to full vectors
  merged <- merge(results_template, event_counts, by = "interval_num", all.x = TRUE)
  merged <- merge(merged, pt_summary,   by = "interval_num", all.x = TRUE)
  merged[is.na(events), events := 0L]
  merged[is.na(person_time), person_time := 0]
  
  list(
    events_per_interval = merged$events,
    person_time_per_interval = merged$person_time
  )
}


# 1) ---------- State container ----------
make_state <- function(arm_names, max_total_patients_per_arm) {
  list(
    arm_status = setNames(rep("recruiting", length(arm_names)), arm_names),
    enrolled_counts = setNames(rep(0L, length(arm_names)), arm_names),
    registries = setNames(
      replicate(length(arm_names),
                data.frame(id = integer(0),
                           enroll_time = numeric(0),
                           true_event_time = numeric(0),
                           random_censor_time = numeric(0)),
                simplify = FALSE),
      arm_names
    ),
    # per-simulation outputs that interims mutate:
    stop_efficacy_per_sim_row = setNames(rep(0L, length(arm_names)), arm_names),
    stop_futility_per_sim_row = setNames(rep(0L, length(arm_names)), arm_names),
    sim_final_n_current_run   = setNames(rep(NA_integer_, length(arm_names)), arm_names)
  )
}


# --- slice_arm_data_at_time (UPDATED: tiny clarity tweak on nrow check) ---
slice_arm_data_at_time <- function(registry_df, calendar_time, max_follow_up, interval_cutpoints) {
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
  event_status  <- as.integer(registry_df$true_event_time <= pmin(registry_df$random_censor_time,
                                                                  time_available))
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

gates_pass_for_both_arms <- function(slA, slB, args, armA_name, armB_name) {
  stopifnot(!missing(armA_name), !missing(armB_name))
  min_ev   <- coalesce_num(args$min_events_per_arm, 0)
  min_mfu  <- coalesce_num(args$min_median_followup_per_arm, 0)
  min_pt_f <- coalesce_num(args$min_person_time_frac_per_arm, 0)
  
  evA  <- coalesce_num(slA$metrics$events_total, 0)
  evB  <- coalesce_num(slB$metrics$events_total, 0)
  mfuA <- coalesce_num(slA$metrics$median_followup, 0)
  mfuB <- coalesce_num(slB$metrics$median_followup, 0)
  
  ptA  <- coalesce_num(slA$metrics$person_time_total, 0)
  ptB  <- coalesce_num(slB$metrics$person_time_total, 0)
  
  lookup_max_pt <- function(arm_name) {
    max_pt <- NULL
    if (!is.null(args$max_PT_per_arm)) {
      max_pt <- args$max_PT_per_arm[[arm_name]]
    }
    if (is.null(max_pt)) {
      max_total <- args$max_total_patients_per_arm[[arm_name]]
      if (is.null(max_total)) return(0)
      max_pt <- max_total * coalesce_num(args$max_follow_up_sim, 0)
    }
    coalesce_num(max_pt, 0)
  }

  maxPT_A <- lookup_max_pt(armA_name)
  maxPT_B <- lookup_max_pt(armB_name)
  
  fracA <- if (maxPT_A > 0) ptA / maxPT_A else 0
  fracB <- if (maxPT_B > 0) ptB / maxPT_B else 0
  
  (evA >= min_ev  && evB >= min_ev) &&
    (mfuA >= min_mfu && mfuB >= min_mfu) &&
    (fracA >= min_pt_f && fracB >= min_pt_f)
}


# --- BETWEEN-ARM POSTERIOR PROBABILITY: FUTILITY (NEW/UPDATED) ---------------
# Proper *futility* probability with the correct directionality:
#   P(Triplet <= Doublet - margin)
# DO NOT use 1 - P(Triplet > Doublet + margin); with a margin this is *not* equivalent.


# Coalesce for numerics
coalesce_num <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x)) y else x
}


# Adapter: slice a single arm at a given calendar time
# Replace body with your project's existing function if names differ.
slice_arm_at_time <- function(arm_name, current_time, arm_data, args) {
  # Your internal slicer likely already exists; keep this as a thin wrapper.
  # Example placeholder:
  calculate_arm_slice(arm_name = arm_name,
                      current_time = current_time,
                      data = arm_data,
                      args = args)
}


# Adapter: map posterior of the piecewise model to a scalar estimand per arm
# Replace the body with your *existing* code that returns posterior draws
# for the median PFS (or any monotone transform you used in calibration).
posterior_scalar_draws <- function(arm_slice, args) {
  # Example placeholder: you likely have something like:
  #   out <- sample_posterior_piecewise_exp(arm_slice, args, n_draws = args$num_posterior_draws)
  #   out$median_pfs_draws
  #
  # For now, fail loud if not implemented:
  if (!is.null(arm_slice$posterior_scalar_draws)) {
    return(arm_slice$posterior_scalar_draws)
  }
  stop("posterior_scalar_draws(): please connect to your posterior draw helper (median PFS etc.)")
}


# Single-arm interim (unchanged). Keep your original implementation.
run_single_arm_interim <- function(current_time, data_by_arm, args) {
  # ... your existing HC interim logic ...
  list(decision = "continue", path = "hc")
}




# --- BETWEEN-ARM POSTERIOR PROBABILITY: EFFICACY (UPDATED) -------------------
# P(Triplet > Doublet + margin), where " > " means whatever clinical estimand you use.
# This version expects each slice to expose posterior *draws* for a scalar estimand
# on the same scale between arms (e.g., median PFS, or -log(hazard), etc.)
# If your code stores piecewise rates, call your existing aggregator that maps
# draws -> scalar estimand per arm before comparing.

