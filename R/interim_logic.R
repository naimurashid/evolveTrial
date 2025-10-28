calculate_current_probs_vs_ref <- function(slCtrl, slTrt, args) {
  draws <- sample_vs_ref_medians(
    slCtrl = slCtrl,
    slTrt = slTrt,
    args = args,
    num_samples = args$num_posterior_draws
  )
  if (isTRUE(args$use_ph_model_vs_ref) && !is.null(draws$logHR)) {
    log_margin <- 0
    if (!is.null(args$compare_arms_hr_margin)) {
      log_margin <- log1p(max(0, args$compare_arms_hr_margin))
    }
    list(
      pr_eff = mean(draws$logHR < 0),
      pr_fut = mean(draws$logHR >= log_margin)
    )
  } else {
    diff_med <- draws$medTrt - draws$medCtrl
    margin <- coalesce_num(args$compare_arms_futility_margin, 0)
    list(
      pr_eff = mean(diff_med > 0),
      pr_fut = mean(diff_med <= -margin)
    )
  }
}

calculate_current_probs_hc <- function(slArm, args, arm_name) {
  interval_lengths <- diff(args$interval_cutpoints_sim)
  lam <- draw_posterior_hazard_samples(
    num_intervals = length(interval_lengths),
    events_per_interval = slArm$metrics$events_per_interval,
    person_time_per_interval = slArm$metrics$person_time_per_interval,
    prior_alpha_params = args$prior_alpha_params_model,
    prior_beta_params  = args$prior_beta_params_model,
    num_samples = args$num_posterior_draws
  )
  med_draws <- apply(lam, 1, calculate_median_survival_piecewise,
                     interval_lengths = interval_lengths)
  success_thr_vec <- args$median_pfs_success_threshold_arms %||% args$null_median_arms
  futility_thr_vec <- args$median_pfs_futility_threshold_arms %||%
    args$futility_median_arms %||% success_thr_vec
  success_thr <- success_thr_vec[[arm_name]] %||% success_thr_vec[1]
  futility_thr <- futility_thr_vec[[arm_name]] %||% futility_thr_vec[1]
  list(
    pr_eff = mean(med_draws > success_thr),
    pr_fut = mean(med_draws < futility_thr)
  )
}

# 3) ---------- Interim checker: pure function (reads state, returns updated state)
# --- INTERIM DECISION ENGINE (UPDATED) ----------------------------------------
# This function assumes you already have utilities that can slice data for a given
# calendar time per arm (what you were calling `arm_slice`), and that those slices
# contain:
#   $metrics$events, $metrics$median_followup, $metrics$person_time_per_interval
# It *short-circuits* to the vs-ref branch whenever compare_arms_option=TRUE.
# Set args$diagnostics=TRUE to print the active path and probabilities.

interim_check <- function(state, current_time, args, diagnostics = FALSE) {
  if (isTRUE(args$compare_arms_option)) {
    reference_arm <- args$reference_arm_name
    if (is.null(reference_arm) || !reference_arm %in% names(state$registries)) {
      stop("interim_check(): args$reference_arm_name must reference an existing arm.")
    }

    experimental_arms <- args$arm_names[args$arm_names != reference_arm]
    if (length(experimental_arms) == 0) {
      if (diagnostics) message(sprintf("[t=%.2f] vsREF skipped: no experimental arms provided", current_time))
      return(state)
    }

    active_experimental <- experimental_arms[
      experimental_arms %in% names(state$arm_status) &
        state$arm_status[experimental_arms] == "recruiting"
    ]
    if (length(active_experimental) == 0) {
      if (diagnostics) message(sprintf("[t=%.2f] vsREF skipped: no recruiting experimental arms", current_time))
      return(state)
    }

    target_override <- args$compare_arms_target_arm
    if (!is.null(target_override)) {
      target_override <- target_override[target_override %in% active_experimental]
    }
    eval_arms <- if (length(target_override) > 0) target_override else active_experimental
    # slice the control arm once per look
    slC <- slice_arm_data_at_time(state$registries[[reference_arm]], current_time,
                                  args$max_follow_up_sim, args$interval_cutpoints_sim)

    for (trt_name in eval_arms) {
      if (state$arm_status[trt_name] != "recruiting") next

      slT <- slice_arm_data_at_time(state$registries[[trt_name]], current_time,
                                    args$max_follow_up_sim, args$interval_cutpoints_sim)

      args_gate <- args
      args_gate$arm_names <- c(reference_arm, trt_name)

      if (!gates_pass_for_both_arms(slC, slT, args_gate, diagnostics = diagnostics)) {
        if (diagnostics) {
          message(sprintf("[t=%.2f] vsREF gated out for %s vs %s", current_time, trt_name, reference_arm))
        }
        next
      }

      probs_vs_ref <- calculate_current_probs_vs_ref(slC, slT, args)
      pr_eff <- probs_vs_ref$pr_eff
      pr_fut <- probs_vs_ref$pr_fut

      if (diagnostics) {
        message(sprintf("[t=%.2f] vsREF %s>%s PrEff>=%.2f: %.3f | PrFut>=%.2f (m=%.2f): %.3f",
                        current_time, trt_name, reference_arm,
                        coalesce_num(args$efficacy_threshold_vs_ref_prob, NA_real_), pr_eff,
                        coalesce_num(args$futility_threshold_vs_ref_prob, NA_real_), coalesce_num(args$compare_arms_futility_margin, 0), pr_fut))
      }

      # Success?
      if (!is.null(args$efficacy_threshold_vs_ref_prob) && is.finite(args$efficacy_threshold_vs_ref_prob) &&
          pr_eff >= args$efficacy_threshold_vs_ref_prob) {
        state$arm_status[trt_name] <- "stopped_efficacy"
        state$stop_efficacy_per_sim_row[trt_name] <- 1L
        state$sim_final_n_current_run[trt_name]   <- state$enrolled_counts[trt_name]
        next
      }

      # Futility?
      if (!is.null(args$futility_threshold_vs_ref_prob) && is.finite(args$futility_threshold_vs_ref_prob) &&
          pr_fut >= args$futility_threshold_vs_ref_prob) {
        state$arm_status[trt_name] <- "stopped_futility"
        state$stop_futility_per_sim_row[trt_name] <- 1L
        state$sim_final_n_current_run[trt_name]   <- state$enrolled_counts[trt_name]
        next
      }
    }

    return(state)
  }

  # HC path (single-arm)
  for (arm in args$arm_names) {
    if (state$arm_status[arm] != "recruiting") next

    slA <- slice_arm_data_at_time(state$registries[[arm]], current_time,
                                  args$max_follow_up_sim, args$interval_cutpoints_sim)
    n_pat <- nrow(slA$patient_data)
    events_total <- sum(slA$patient_data$event_status)
    median_fu <- if (n_pat > 0) stats::median(slA$patient_data$observed_time) else 0
    pt_total <- sum(slA$metrics$person_time_per_interval)

    min_pat <- coalesce_num(args$min_patients_for_analysis, 0)
    min_events <- coalesce_num(args$min_events_for_analysis, 0)
    min_median_fu <- coalesce_num(args$min_median_followup, 0)

    min_pt_frac <- 0
    if (!is.null(args$min_person_time_frac_per_arm)) {
      gate_vec <- args$min_person_time_frac_per_arm
      if (!is.null(names(gate_vec))) {
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

    probs_hc <- calculate_current_probs_hc(slA, args, arm)
    pr_eff <- probs_hc$pr_eff
    pr_fut <- probs_hc$pr_fut

    if (diagnostics) {
      message(sprintf("[t=%.2f] HC %s PrEff>=%.3f: %.3f | PrFut>=%.3f: %.3f",
                      current_time, arm,
                      coalesce_num(args$efficacy_threshold_current_prob_hc, NA_real_),
                      pr_eff,
                      coalesce_num(args$posterior_futility_threshold_hc, NA_real_),
                      pr_fut))
    }

    if (!is.null(args$efficacy_threshold_current_prob_hc) &&
        is.finite(args$efficacy_threshold_current_prob_hc) &&
        pr_eff >= args$efficacy_threshold_current_prob_hc) {
      state$arm_status[arm] <- "stopped_efficacy"
      state$stop_efficacy_per_sim_row[arm] <- 1L
      state$sim_final_n_current_run[arm] <- state$enrolled_counts[arm]
      next
    }

    if (!is.null(args$posterior_futility_threshold_hc) &&
        is.finite(args$posterior_futility_threshold_hc) &&
        pr_fut >= args$posterior_futility_threshold_hc) {
      state$arm_status[arm] <- "stopped_futility"
      state$stop_futility_per_sim_row[arm] <- 1L
      state$sim_final_n_current_run[arm] <- state$enrolled_counts[arm]
      next
    }
  }

  state
}
