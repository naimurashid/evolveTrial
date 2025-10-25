calculate_current_probs_vs_ref <- function(slCtrl, slTrt, args) {
  medians <- sample_vs_ref_medians(
    slCtrl = slCtrl,
    slTrt = slTrt,
    args = args,
    num_samples = args$num_posterior_draws
  )
  diff_med <- medians$medTrt - medians$medCtrl
  margin <- coalesce_num(args$compare_arms_futility_margin, 0)
  list(
    pr_eff = mean(diff_med > 0),
    pr_fut = mean(diff_med <= -margin)
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
    trt_name <- if (length(target_override) > 0) {
      target_override[1]
    } else {
      active_experimental[1]
    }

    slC <- slice_arm_data_at_time(state$registries[[reference_arm]], current_time,
                                  args$max_follow_up_sim, args$interval_cutpoints_sim)
    slT <- slice_arm_data_at_time(state$registries[[trt_name]], current_time,
                                  args$max_follow_up_sim, args$interval_cutpoints_sim)

    args_gate <- args
    args_gate$arm_names <- c(reference_arm, trt_name)

    if (!gates_pass_for_both_arms(slC, slT, args_gate, diagnostics = diagnostics)) {
      if (diagnostics) message(sprintf("[t=%.2f] vsREF gated out for %s vs %s", current_time, trt_name, reference_arm))
      return(state)
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
      return(state)
    }

    # Futility?
    if (!is.null(args$futility_threshold_vs_ref_prob) && is.finite(args$futility_threshold_vs_ref_prob) &&
        pr_fut >= args$futility_threshold_vs_ref_prob) {
      state$arm_status[trt_name] <- "stopped_futility"
      state$stop_futility_per_sim_row[trt_name] <- 1L
      state$sim_final_n_current_run[trt_name]   <- state$enrolled_counts[trt_name]
      return(state)
    }

    return(state)
  }

  # HC path (single-arm) — leave as is / your existing logic
  state
}
