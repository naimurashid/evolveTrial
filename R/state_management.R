
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
  ctrl_name <- args$reference_arm_name %||% "Doublet"
  arm_names <- args$arm_names
  if (is.null(arm_names) || length(arm_names) == 0) {
    arm_names <- c("Doublet", "Triplet")
  }
  trt_candidates <- arm_names[arm_names != ctrl_name]
  trt_name <- if (length(trt_candidates) > 0) trt_candidates[1] else "Triplet"

  # thresholds (per-arm gates)
  resolve_gate_vec <- function(raw, default = 0, scale = FALSE) {
    arms <- c(ctrl_name, trt_name)
    raw_len <- if (is.null(raw)) 0 else length(raw)
    if (is.null(raw)) {
      out <- rep(default, length(arms))
    } else if (!is.null(names(raw))) {
      out <- raw[arms]
      out[is.na(out)] <- default
    } else if (length(raw) == length(arm_names)) {
      out <- raw
      names(out) <- arm_names
      out <- out[arms]
    } else if (length(raw) == length(arms)) {
      out <- raw
    } else {
      out <- rep(raw[1], length(arms))
    }
    names(out) <- arms
    should_scale <- scale && (raw_len <= 1 || is.null(raw))
    if (should_scale) {
      probs <- args$randomization_probs
      if (is.null(probs)) {
        probs <- rep(1 / length(arm_names), length(arm_names))
        names(probs) <- arm_names
      } else if (is.null(names(probs)) && length(probs) == length(arm_names)) {
        names(probs) <- arm_names
      }
      scale_vec <- rep(1, length(arms))
      names(scale_vec) <- arms
      if (!is.null(names(probs))) {
        max_prob <- max(probs, na.rm = TRUE)
        if (is.finite(max_prob) && max_prob > 0) {
          scale_vec <- probs[arms] / max_prob
          scale_vec[!is.finite(scale_vec)] <- 1
        }
      }
      out <- out * scale_vec
    }
    out
  }

  min_ev_vec   <- resolve_gate_vec(args$min_events_per_arm, 0, scale = TRUE)
  min_mfu_vec  <- resolve_gate_vec(args$min_median_followup_per_arm, 0, scale = FALSE)
  min_pt_vec   <- resolve_gate_vec(args$min_person_time_frac_per_arm, 0, scale = TRUE)
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

  # denominators (person-time caps) — handle unnamed vectors safely
  mtn <- args$max_total_patients_per_arm
  max_follow <- coalesce_num(args$max_follow_up_sim, 0)
  if (!is.null(mtn) && (is.null(names(mtn)) || any(!nzchar(names(mtn))))) {
    if (!is.null(arm_names) && length(mtn) == length(arm_names)) {
      names(mtn) <- arm_names
    }
  }
  if (is.null(mtn)) {
    maxPT_C <- 0; maxPT_T <- 0
  } else if (!is.null(names(mtn)) &&
             all(c(ctrl_name, trt_name) %in% names(mtn))) {
    maxPT_C <- coalesce_num(mtn[[ctrl_name]], 0) * max_follow
    maxPT_T <- coalesce_num(mtn[[trt_name]],  0) * max_follow
  } else {
    idx_ctrl <- match(ctrl_name, arm_names)
    idx_trt  <- match(trt_name,  arm_names)
    if (is.na(idx_ctrl) || idx_ctrl > length(mtn)) idx_ctrl <- 1L
    if (is.na(idx_trt)  || idx_trt  > length(mtn)) {
      idx_trt <- if (length(mtn) >= 2L) 2L else idx_ctrl
    }
    maxPT_C <- coalesce_num(as.numeric(mtn[idx_ctrl]), 0) * max_follow
    maxPT_T <- coalesce_num(as.numeric(mtn[idx_trt]),  0) * max_follow
  }

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


# Adapter: slice a single arm at a given calendar time
# Replace body with your project's existing function if names differ.
slice_arm_at_time <- function(arm_name, current_time, arm_data, args) {
  interval_cutpoints <- args$interval_cutpoints_sim
  if (is.null(interval_cutpoints)) {
    stop("slice_arm_at_time(): args$interval_cutpoints_sim must be supplied.")
  }
  max_follow <- coalesce_num(args$max_follow_up_sim, 0)
  if (max_follow <= 0) {
    stop("slice_arm_at_time(): args$max_follow_up_sim must be positive.")
  }
  if (is.null(arm_data)) {
    stop(sprintf("slice_arm_at_time(): data for arm '%s' is NULL.", arm_name))
  }
  slice_arm_data_at_time(
    registry_df = arm_data,
    calendar_time = current_time,
    max_follow_up = max_follow,
    interval_cutpoints = interval_cutpoints
  )
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
  if (is.null(names(data_by_arm))) {
    stop("run_single_arm_interim(): data_by_arm must be a named list of registries.")
  }
  decisions <- setNames(rep("continue", length(data_by_arm)), names(data_by_arm))
  pr_eff <- pr_fut <- setNames(rep(NA_real_, length(data_by_arm)), names(data_by_arm))

  for (arm in names(data_by_arm)) {
    slice <- slice_arm_at_time(
      arm_name = arm,
      current_time = current_time,
      arm_data = data_by_arm[[arm]],
      args = args
    )

    n_pat <- nrow(slice$patient_data)
    events_total <- sum(slice$patient_data$event_status)
    median_fu <- if (n_pat > 0) stats::median(slice$patient_data$observed_time) else 0
    pt_total <- sum(slice$metrics$person_time_per_interval)

    min_pat <- coalesce_num(args$min_patients_for_analysis, 0)
    min_events <- coalesce_num(args$min_events_for_analysis, 0)
    min_median_fu <- coalesce_num(args$min_median_followup, 0)

    min_pt_frac <- 0
    if (!is.null(args$min_person_time_frac_per_arm)) {
      gate_vec <- args$min_person_time_frac_per_arm
      if (!is.null(names(gate_vec)) && arm %in% names(gate_vec)) {
        min_pt_frac <- coalesce_num(gate_vec[[arm]], 0)
      } else if (length(gate_vec) == length(args$arm_names)) {
        min_pt_frac <- coalesce_num(gate_vec[match(arm, args$arm_names)], 0)
      } else if (length(gate_vec) >= 1) {
        min_pt_frac <- coalesce_num(gate_vec[1], 0)
      }
    }

    maxPT_arm <- 0
    mtpa <- args$max_total_patients_per_arm
    if (!is.null(mtpa)) {
      if (!is.null(names(mtpa)) && arm %in% names(mtpa)) {
        maxPT_arm <- coalesce_num(mtpa[[arm]], 0) * coalesce_num(args$max_follow_up_sim, 0)
      } else if (length(mtpa) == length(args$arm_names)) {
        maxPT_arm <- coalesce_num(mtpa[match(arm, args$arm_names)], 0) *
          coalesce_num(args$max_follow_up_sim, 0)
      } else if (length(mtpa) >= 1) {
        maxPT_arm <- coalesce_num(mtpa[1], 0) * coalesce_num(args$max_follow_up_sim, 0)
      }
    }
    pt_frac <- if (maxPT_arm > 0) pt_total / maxPT_arm else 0

    if (n_pat < min_pat ||
        events_total < min_events ||
        median_fu < min_median_fu ||
        pt_frac < min_pt_frac) {
      next
    }

    probs <- calculate_current_probs_hc(slice, args, arm)
    pr_eff[arm] <- probs$pr_eff
    pr_fut[arm] <- probs$pr_fut

    if (!is.null(args$efficacy_threshold_current_prob_hc) &&
        is.finite(args$efficacy_threshold_current_prob_hc) &&
        probs$pr_eff >= args$efficacy_threshold_current_prob_hc) {
      decisions[arm] <- "stop_efficacy"
      next
    }

    if (!is.null(args$posterior_futility_threshold_hc) &&
        is.finite(args$posterior_futility_threshold_hc) &&
        probs$pr_fut >= args$posterior_futility_threshold_hc) {
      decisions[arm] <- "stop_futility"
    }
  }

  list(
    decisions = decisions,
    pr_eff = pr_eff,
    pr_fut = pr_fut,
    path = "hc"
  )
}




# --- BETWEEN-ARM POSTERIOR PROBABILITY: EFFICACY (UPDATED) -------------------
# P(Triplet > Doublet + margin), where " > " means whatever clinical estimand you use.
# This version expects each slice to expose posterior *draws* for a scalar estimand
# on the same scale between arms (e.g., median PFS, or -log(hazard), etc.)
# If your code stores piecewise rates, call your existing aggregator that maps
# draws -> scalar estimand per arm before comparing.
