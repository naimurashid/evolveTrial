calculate_current_prob_vs_ref_futility <- function(slCtrl, slTrt, args) {
  K <- length(args$interval_cutpoints_sim) - 1
  L  <- diff(args$interval_cutpoints_sim)
  lamC <- draw_posterior_hazard_samples(K,
                                        slCtrl$metrics$events_per_interval,
                                        slCtrl$metrics$person_time_per_interval,
                                        args$prior_alpha_params_model,
                                        args$prior_beta_params_model,
                                        num_samples = args$num_posterior_draws)
  lamT <- draw_posterior_hazard_samples(K,
                                        slTrt$metrics$events_per_interval,
                                        slTrt$metrics$person_time_per_interval,
                                        args$prior_alpha_params_model,
                                        args$prior_beta_params_model,
                                        num_samples = args$num_posterior_draws)
  medC <- apply(lamC, 1, calculate_median_survival_piecewise, interval_lengths = L)
  medT <- apply(lamT, 1, calculate_median_survival_piecewise, interval_lengths = L)
  
  margin <- coalesce_num(args$compare_arms_futility_margin, 0)
  mean((medT - medC) <= -margin)
}

calculate_current_prob_vs_ref <- function(slCtrl, slTrt, args) {
  K <- length(args$interval_cutpoints_sim) - 1
  L  <- diff(args$interval_cutpoints_sim)
  # posterior hazards
  lamC <- draw_posterior_hazard_samples(K,
                                        slCtrl$metrics$events_per_interval,
                                        slCtrl$metrics$person_time_per_interval,
                                        args$prior_alpha_params_model,
                                        args$prior_beta_params_model,
                                        num_samples = args$num_posterior_draws)
  lamT <- draw_posterior_hazard_samples(K,
                                        slTrt$metrics$events_per_interval,
                                        slTrt$metrics$person_time_per_interval,
                                        args$prior_alpha_params_model,
                                        args$prior_beta_params_model,
                                        num_samples = args$num_posterior_draws)
  medC <- apply(lamC, 1, calculate_median_survival_piecewise, interval_lengths = L)
  medT <- apply(lamT, 1, calculate_median_survival_piecewise, interval_lengths = L)
  
  margin <- coalesce_num(args$compare_arms_futility_margin, 0) # see item 4 below
  mean((medT - medC) > margin)
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
    experimental_arms <- setdiff(args$arm_names, reference_arm)
    target_arm <- args$compare_arms_target_arm %||% experimental_arms[1]
    if (length(experimental_arms) == 0 || is.null(target_arm) || is.na(target_arm)) {
      stop("compare_arms_option requires at least one non-reference arm to evaluate.")
    }
    if (!target_arm %in% args$arm_names) {
      stop(sprintf("Arm '%s' is not present in args$arm_names.", target_arm))
    }

    slC <- slice_arm_data_at_time(state$registries[[reference_arm]], current_time,
                                  args$max_follow_up_sim, args$interval_cutpoints_sim)
    slT <- slice_arm_data_at_time(state$registries[[target_arm]],  current_time,
                                  args$max_follow_up_sim, args$interval_cutpoints_sim)
    
    if (!gates_pass_for_both_arms(slC, slT, args, reference_arm, target_arm)) {
      if (diagnostics) message(sprintf("[t=%.2f] vsREF gated out for %s vs %s", current_time, target_arm, reference_arm))
      return(state)
    }
    
    pr_eff <- calculate_current_prob_vs_ref(slC, slT, args)
    pr_fut <- calculate_current_prob_vs_ref_futility(slC, slT, args)
    
    if (diagnostics) {
      message(sprintf("[t=%.2f] vsREF %s>%s PrEff>=%.2f: %.3f | PrFut>=%.2f (m=%.2f): %.3f",
                      current_time, target_arm, reference_arm,
                      coalesce_num(args$efficacy_threshold_vs_ref_prob, NA_real_), pr_eff,
                      coalesce_num(args$futility_threshold_vs_ref_prob,  NA_real_), coalesce_num(args$compare_arms_futility_margin, 0), pr_fut))
    }
    
    # Success?
    if (!is.null(args$efficacy_threshold_vs_ref_prob) && is.finite(args$efficacy_threshold_vs_ref_prob) &&
        pr_eff >= args$efficacy_threshold_vs_ref_prob) {
      state$arm_status[target_arm] <- "stopped_efficacy"
      state$stop_efficacy_per_sim_row[target_arm] <- 1L
      state$sim_final_n_current_run[target_arm]   <- state$enrolled_counts[target_arm]
      return(state)
    }
    
    # Futility?
    if (!is.null(args$futility_threshold_vs_ref_prob) && is.finite(args$futility_threshold_vs_ref_prob) &&
        pr_fut >= args$futility_threshold_vs_ref_prob) {
      state$arm_status[target_arm] <- "stopped_futility"
      state$stop_futility_per_sim_row[target_arm] <- 1L
      state$sim_final_n_current_run[target_arm]   <- state$enrolled_counts[target_arm]
      return(state)
    }
    
    return(state)
  }
  
  # HC path (not compare_arms)
  # ... your existing single-arm logic that mutates state ...
  state
}

