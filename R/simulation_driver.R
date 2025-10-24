
# 4) ---------- Main simulation using pure-state threading ----------
run_simulation_pure <- function(
    num_simulations,
    arm_names,
    reference_arm_name,
    compare_arms_option,
    weibull_shape_true_arms,
    weibull_median_true_arms,
    null_median_arms,
    futility_median_arms,
    interval_cutpoints_sim,
    max_follow_up_sim,
    censor_max_time_sim,
    prior_alpha_params_model,
    prior_beta_params_model,
    num_posterior_draws,
    cohort_size_per_arm,
    max_total_patients_per_arm,
    min_patients_for_analysis,
    # stopping rule params
    efficacy_stopping_rule_hc,
    efficacy_threshold_current_prob_hc,
    posterior_futility_threshold_hc,
    futility_stopping_rule_hc,
    efficacy_stopping_rule_vs_ref,
    futility_stopping_rule_vs_ref,
    efficacy_threshold_vs_ref_prob,
    futility_threshold_vs_ref_prob,
    compare_arms_futility_margin,
    # final analysis params
    median_pfs_success_threshold_arms,
    final_success_posterior_prob_threshold,
    median_pfs_futility_threshold_arms,
    final_futility_posterior_prob_threshold,
    # accrual & randomization
    overall_accrual_rate,
    randomization_probs,
    min_follow_up_at_final = 0,
    # legacy info gates + calendar-beat interims
    min_events_for_analysis = 0,
    min_median_followup   = 0,
    interim_calendar_beat = 2,
    diagnostics = FALSE,
    pred_success_pp_threshold_hc,
    pred_futility_pp_threshold_hc,
    num_posterior_draws_pred,
    predictive_fast = FALSE,
    # NEW: per-arm info gates
    min_events_per_arm = NULL,
    min_median_followup_per_arm = NULL,
    min_person_time_frac_per_arm = 0,
    min_events_ratio_arm_vs_ref = 0.0,
    # NEW: person-time milestones
    person_time_milestones = NULL,
    latest_calendar_look = Inf
) {
  weibull_scale_true_arms <- sapply(arm_names, function(arm) {
    calculate_weibull_scale(weibull_median_true_arms[arm], weibull_shape_true_arms[arm])
  })
  names(weibull_scale_true_arms) <- arm_names
  
  results_data <- data.frame(
    Arm_Name = arm_names,
    True_Median = round(weibull_median_true_arms, 2),
    Type_I_Error_or_Power = 0,
    PET_Efficacy = 0,
    PET_Futility = 0,
    Pr_Reach_Max_N = 0,
    Pr_Final_Efficacy = 0,
    Pr_Final_Futility = 0,
    Pr_Final_Inconclusive = 0,
    Exp_N = 0,
    stringsAsFactors = FALSE
  )
  
  final_n_per_sim <- matrix(0, nrow = num_simulations, ncol = length(arm_names),
                            dimnames = list(NULL, arm_names))
  stop_efficacy_per_sim <- matrix(0, nrow = num_simulations, ncol = length(arm_names),
                                  dimnames = list(NULL, arm_names))
  stop_futility_per_sim <- matrix(0, nrow = num_simulations, ncol = length(arm_names),
                                  dimnames = list(NULL, arm_names))
  final_efficacy_per_sim <- matrix(0, nrow = num_simulations, ncol = length(arm_names),
                                   dimnames = list(NULL, arm_names))
  final_futility_per_sim <- matrix(0, nrow = num_simulations, ncol = length(arm_names),
                                   dimnames = list(NULL, arm_names))
  final_inconclusive_per_sim <- matrix(0, nrow = num_simulations, ncol = length(arm_names),
                                       dimnames = list(NULL, arm_names))
  
  pb <- progress::progress_bar$new(
    format = "  Sims [:bar] :percent in :elapsed",
    total = num_simulations, clear = FALSE, width = 60
  )
  
  # total max person-time across arms (for milestone schedule)
  max_PT_per_arm <- setNames(as.numeric(max_total_patients_per_arm) * max_follow_up_sim,
                             names(max_total_patients_per_arm))
  total_max_PT <- sum(max_PT_per_arm)
  
  # Build absolute PT milestones (if any)
  pt_targets_abs <- NULL
  if (!is.null(person_time_milestones)) {
    stopifnot(all(person_time_milestones > 0 & person_time_milestones <= 1))
    pt_targets_abs <- sort(unique(person_time_milestones)) * total_max_PT
  }
  
  # pack args for interim_check
  args <- list(
    arm_names = arm_names,
    reference_arm_name = reference_arm_name,
    compare_arms_option = compare_arms_option,
    interval_cutpoints_sim = interval_cutpoints_sim,
    max_follow_up_sim = max_follow_up_sim,
    prior_alpha_params_model = prior_alpha_params_model,
    prior_beta_params_model  = prior_beta_params_model,
    num_posterior_draws = num_posterior_draws,
    min_patients_for_analysis = min_patients_for_analysis,
    min_events_for_analysis = min_events_for_analysis,
    min_median_followup = min_median_followup,
    null_median_arms = null_median_arms,
    futility_median_arms = futility_median_arms,
    efficacy_threshold_current_prob_hc = efficacy_threshold_current_prob_hc,
    posterior_futility_threshold_hc = posterior_futility_threshold_hc,
    efficacy_threshold_vs_ref_prob = efficacy_threshold_vs_ref_prob,
    futility_threshold_vs_ref_prob = futility_threshold_vs_ref_prob,
    compare_arms_futility_margin = compare_arms_futility_margin,
    max_total_patients_per_arm = max_total_patients_per_arm,
    median_pfs_success_threshold_arms = median_pfs_success_threshold_arms,
    final_success_posterior_prob_threshold = final_success_posterior_prob_threshold,
    final_futility_posterior_prob_threshold = final_futility_posterior_prob_threshold,
    censor_max_time_sim = censor_max_time_sim,
    predictive_fast = predictive_fast,
    num_posterior_draws_pred = num_posterior_draws_pred,
    predictive_boundary_band  = 0.12,
    run_pred_every_k_beats    = 2,
    pred_success_pp_threshold_hc  = pred_success_pp_threshold_hc,
    pred_futility_pp_threshold_hc = pred_futility_pp_threshold_hc,
    # per-arm gates + PT metadata
    min_events_per_arm = min_events_per_arm,
    min_median_followup_per_arm = min_median_followup_per_arm,
    min_person_time_frac_per_arm = min_person_time_frac_per_arm,
    min_events_ratio_arm_vs_ref = min_events_ratio_arm_vs_ref,
    max_PT_per_arm = max_PT_per_arm
  )
  
  for (s in 1:num_simulations) {
    state <- make_state(arm_names, max_total_patients_per_arm)
    
    current_time <- 0.0
    next_calendar_look <- interim_calendar_beat
    next_pt_idx <- 1L
    patient_id <- 0L
    
    is_eligible <- function(st) {
      which((st$arm_status == "recruiting") & (st$enrolled_counts < max_total_patients_per_arm))
    }
    
    # accrual loop
    while (length(is_eligible(state)) > 0) {
      interarrival <- rexp(1, rate = overall_accrual_rate)
      current_time <- current_time + interarrival
      
      # ---- LOOK SCHEDULING ----
      if (!is.null(pt_targets_abs)) {
        # PT-based schedule (with calendar backstop)
        total_PT_now <- cum_person_time_all_arms(state, current_time, max_follow_up_sim,
                                                 interval_cutpoints_sim, arm_names)
        if (is.null(state$.__backstop_fired__)) state$.__backstop_fired__ <- FALSE
        
        while (!is.null(pt_targets_abs) &&
               next_pt_idx <= length(pt_targets_abs) &&
               (total_PT_now >= pt_targets_abs[next_pt_idx] ||
                (!state$.__backstop_fired__ && is.finite(latest_calendar_look) && current_time >= latest_calendar_look))) {
          
          state <- interim_check(state, current_time, args, diagnostics = diagnostics)
          
          if (total_PT_now >= pt_targets_abs[next_pt_idx]) {
            next_pt_idx <- next_pt_idx + 1L
          }
          if (!state$.__backstop_fired__ && is.finite(latest_calendar_look) && current_time >= latest_calendar_look) {
            state$.__backstop_fired__ <- TRUE
          }
          
          if (all(state$arm_status != "recruiting")) break
          total_PT_now <- cum_person_time_all_arms(state, current_time, max_follow_up_sim,
                                                   interval_cutpoints_sim, arm_names)
        }
      } else {
        # calendar-beat schedule
        if (current_time >= next_calendar_look) {
          state <- interim_check(state, current_time, args, diagnostics = diagnostics)
          next_calendar_look <- next_calendar_look + interim_calendar_beat
        }
      }
      
      # if all arms stopped at a look, quit
      if (all(state$arm_status != "recruiting")) break
      
      # randomize next patient
      elig_idx <- is_eligible(state)
      if (length(elig_idx) == 0) break
      elig_arms <- arm_names[elig_idx]
      
      probs <- randomization_probs[elig_arms]
      probs <- probs / sum(probs)
      chosen_arm <- sample(elig_arms, size = 1, prob = probs)
      
      patient_id <- patient_id + 1L
      t_event_true <- rweibull(1, shape = weibull_shape_true_arms[chosen_arm],
                               scale = weibull_scale_true_arms[chosen_arm])
      t_random_censor <- runif(1, min = 0, max = censor_max_time_sim)
      
      state$registries[[chosen_arm]] <- rbind(state$registries[[chosen_arm]],
                                              data.frame(
                                                id = patient_id,
                                                enroll_time = current_time,
                                                true_event_time = t_event_true,
                                                random_censor_time = t_random_censor
                                              ))
      state$enrolled_counts[chosen_arm] <- state$enrolled_counts[chosen_arm] + 1L
    }
    
    # last interim if we ended between looks
    state <- interim_check(state, current_time, args, diagnostics = diagnostics)
    
    # final analysis
    last_enroll_time <- max(c(0, unlist(lapply(state$registries, function(df) df$enroll_time))), na.rm = TRUE)
    final_time <- last_enroll_time + min_follow_up_at_final
    interval_lengths <- diff(interval_cutpoints_sim)
    num_intervals <- length(interval_lengths)
    
    for (arm in arm_names) {
      if (state$arm_status[arm] != "recruiting") next
      
      arm_slice <- slice_arm_data_at_time(state$registries[[arm]], final_time,
                                          max_follow_up_sim, interval_cutpoints_sim)
      post_arm <- draw_posterior_hazard_samples(
        num_intervals = num_intervals,
        events_per_interval = arm_slice$metrics$events_per_interval,
        person_time_per_interval = arm_slice$metrics$person_time_per_interval,
        prior_alpha_params = prior_alpha_params_model,
        prior_beta_params  = prior_beta_params_model,
        num_samples = num_posterior_draws
      )
      med_arm <- apply(post_arm, 1, function(h) {
        calculate_median_survival_piecewise(h, interval_lengths)
      })
      
      if (!compare_arms_option) {
        # Single-arm final vs absolute thresholds
        p_eff <- mean(med_arm > median_pfs_success_threshold_arms[arm])
        if (p_eff >= final_success_posterior_prob_threshold) {
          final_efficacy_per_sim[s, arm] <- 1L
        } else {
          p_fut <- mean(med_arm < median_pfs_futility_threshold_arms[arm])
          if (p_fut >= final_futility_posterior_prob_threshold) {
            final_futility_per_sim[s, arm] <- 1L
          } else {
            final_inconclusive_per_sim[s, arm] <- 1L
          }
        }
      } else {
        # Between-arm final (Triplet vs Doublet) with ABSOLUTE margin
        if (arm == reference_arm_name) {
          final_inconclusive_per_sim[s, arm] <- 1L
        } else {
          ref_slice <- slice_arm_data_at_time(state$registries[[reference_arm_name]], final_time,
                                              max_follow_up_sim, interval_cutpoints_sim)
          post_ref <- draw_posterior_hazard_samples(
            num_intervals = num_intervals,
            events_per_interval = ref_slice$metrics$events_per_interval,
            person_time_per_interval = ref_slice$metrics$person_time_per_interval,
            prior_alpha_params = prior_alpha_params_model,
            prior_beta_params  = prior_beta_params_model,
            num_samples = num_posterior_draws
          )
          med_ref <- apply(post_ref, 1, function(h) {
            calculate_median_survival_piecewise(h, interval_lengths)
          })
          
          margin_abs <- coalesce_num(compare_arms_futility_margin, 0)
          pr <- final_vsref_probs_abs(med_arm, med_ref, margin_abs)
          p_eff_ref <- pr$p_eff_ref
          p_fut_ref <- pr$p_fut_ref
          
          if (p_eff_ref >= efficacy_threshold_vs_ref_prob) {
            final_efficacy_per_sim[s, arm] <- 1L
          } else if (p_fut_ref >= futility_threshold_vs_ref_prob) {
            final_futility_per_sim[s, arm] <- 1L
          } else {
            final_inconclusive_per_sim[s, arm] <- 1L
          }
        }
      }
    }
    
    for (arm in arm_names) {
      final_n_per_sim[s, arm] <- if (is.na(state$sim_final_n_current_run[arm])) {
        state$enrolled_counts[arm]
      } else {
        state$sim_final_n_current_run[arm]
      }
      stop_efficacy_per_sim[s, arm] <- state$stop_efficacy_per_sim_row[arm]
      stop_futility_per_sim[s, arm] <- state$stop_futility_per_sim_row[arm]
    }
    
    pb$tick()
  } # sims
  
  for (j in seq_along(arm_names)) {
    arm <- arm_names[j]
    early_eff <- mean(stop_efficacy_per_sim[, arm])
    early_fut <- mean(stop_futility_per_sim[, arm])
    final_eff <- mean(final_efficacy_per_sim[, arm])
    final_fut <- mean(final_futility_per_sim[, arm])
    final_inc <- mean(final_inconclusive_per_sim[, arm])
    results_data$Type_I_Error_or_Power[j]   <- early_eff + final_eff
    results_data$PET_Efficacy[j]            <- early_eff
    results_data$PET_Futility[j]            <- early_fut
    results_data$Pr_Reach_Max_N[j]          <- 1 - early_eff - early_fut
    results_data$Pr_Final_Efficacy[j]       <- final_eff
    results_data$Pr_Final_Futility[j]       <- final_fut
    results_data$Pr_Final_Inconclusive[j]   <- final_inc
    results_data$Exp_N[j]                   <- mean(final_n_per_sim[, arm])
  }
  
  results_data
}


# ---------- 1) Build scenarios from a grid of options ----------
# Pass a named list where each element is either:
#   - a vector of scalar options (e.g., c(60, 70)), OR
#   - a list of per-arm options (e.g., list(c(Arm1=6,Arm2=6,Arm3=9), c(...)))
# Returns: a list of "override" lists, one per scenario.
scenarios_from_grid <- function(choices) {
  stopifnot(is.list(choices), length(choices) > 0)
  keys <- names(choices)
  # For list-valued choices (e.g., per-arm vectors), expand over their indices
  grid_input <- lapply(choices, function(v) if (is.list(v)) seq_along(v) else v)
  raw_grid <- do.call(expand.grid, c(grid_input, stringsAsFactors = FALSE))
  # Turn each grid row into an override list, replacing list indices by actual list elements
  scen_list <- lapply(seq_len(nrow(raw_grid)), function(i) {
    row <- raw_grid[i, , drop = FALSE]
    over <- setNames(vector("list", length(keys)), keys)
    for (k in keys) {
      v <- choices[[k]]
      if (is.list(v)) {
        over[[k]] <- v[[ as.integer(row[[k]]) ]]
      } else {
        over[[k]] <- row[[k]]
      }
    }
    over
  })
  # Attach a compact label for readability (optional)
  attr(scen_list, "grid") <- raw_grid
  scen_list
}

run_scenarios <- function(base_args, scens, parallel = FALSE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # keep only args that run_simulation_pure actually accepts
  rs_formals <- names(formals(run_simulation_pure))
  
  run_one <- function(i) {
    over   <- scens[[i]]
    args_i <- utils::modifyList(base_args, over, keep.null = TRUE)
    args_i <- args_i[intersect(names(args_i), rs_formals)]
    res <- do.call(run_simulation_pure, args_i)
    res$scenario <- i
    res
  }
  
  if (isTRUE(parallel)) {
    cores <- max(1L, parallel::detectCores() - 1L)
    out <- parallel::mclapply(seq_along(scens), run_one, mc.cores = cores)
  } else {
    out <- lapply(seq_along(scens), run_one)
  }
  
  # validate all returned items are tabular
  ok_types <- vapply(out, function(x) is.data.frame(x) || is.list(x) || data.table::is.data.table(x), logical(1))
  if (!all(ok_types)) stop("run_scenarios: one or more scenarios did not return a tabular result.")
  
  data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
}

