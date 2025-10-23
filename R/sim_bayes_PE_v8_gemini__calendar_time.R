
`%||%` <- function(a, b) if (!is.null(a)) a else b

cum_person_time_arm <- function(state, arm, current_time, max_follow_up, interval_cutpoints) {
  sl <- slice_arm_data_at_time(state$registries[[arm]], current_time, max_follow_up, interval_cutpoints)
  sum(sl$patient_data$observed_time)
}

cum_events_arm <- function(state, arm, current_time, max_follow_up, interval_cutpoints) {
  sl <- slice_arm_data_at_time(state$registries[[arm]], current_time, max_follow_up, interval_cutpoints)
  sum(sl$patient_data$event_status)
}

median_followup_arm <- function(state, arm, current_time, max_follow_up, interval_cutpoints) {
  sl <- slice_arm_data_at_time(state$registries[[arm]], current_time, max_follow_up, interval_cutpoints)
  if (nrow(sl$patient_data) == 0) 0 else stats::median(sl$patient_data$observed_time)
}

cum_person_time_all_arms <- function(state, current_time, max_follow_up, interval_cutpoints, arm_names) {
  sum(vapply(arm_names, function(a)
    cum_person_time_arm(state, a, current_time, max_follow_up, interval_cutpoints), numeric(1)))
}

# --- HELPER: Weibull scale from median (R parameterization) ---
# R's rweibull(n, shape = k, scale = λ) uses S(t) = exp(-(t/λ)^k).
# The population median t~ satisfies 0.5 = exp(-(t~/λ)^k)  ⇒  t~ = λ * (ln 2)^(1/k).
# Solve for λ:  λ = t~ / (ln 2)^(1/k).
calculate_weibull_scale <- function(desired_median, weibull_shape) {
  stopifnot(is.numeric(desired_median), is.numeric(weibull_shape),
            length(desired_median) == 1, length(weibull_shape) == 1,
            desired_median > 0, weibull_shape > 0)
  desired_median / (log(2)^(1 / weibull_shape))
}

# --- HELPER FUNCTION: Simulate PFS data from a Piecewise Exponential Distribution ---
simulate_piecewise_exponential_data <- function(
    num_patients,
    hazard_rates,
    interval_cutpoints,
    max_follow_up,
    censor_min_time = 0,
    censor_max_time = NULL,
    start_id = 1
) {
  # Simulates PFS times from a piecewise exponential distribution and applies censoring.
  # The LAST interval is treated as open-ended (extends to infinity).

  num_intervals <- length(interval_cutpoints) - 1
  interval_lengths <- diff(interval_cutpoints)

  if (length(hazard_rates) != num_intervals) {
    stop("hazard_rates must have length num_intervals (length(interval_cutpoints) - 1).")
  }

  if (is.null(censor_max_time)) censor_max_time <- max_follow_up
  if (censor_min_time < 0 || censor_max_time < 0 || censor_min_time > censor_max_time) {
    stop("Require 0 <= censor_min_time <= censor_max_time.")
  }

  true_pfs_times <- numeric(num_patients)
  cumulative_hazards_at_cutpoints <- c(0, cumsum(hazard_rates * interval_lengths))

  for (p in 1:num_patients) {
    U <- runif(1)
    target_H <- -log(U)

    # Find first interval where cumulative hazard exceeds target
    event_interval_idx <- which(target_H <= cumulative_hazards_at_cutpoints[-1])[1]

    if (is.na(event_interval_idx)) {
      # Spill into an open-ended last interval
      event_interval_idx <- num_intervals
      lambda_j <- hazard_rates[event_interval_idx]
      interval_start_time <- interval_cutpoints[event_interval_idx]
      H_start <- cumulative_hazards_at_cutpoints[event_interval_idx]

      if (lambda_j == 0) {
        true_pfs_times[p] <- Inf
      } else {
        time_in_interval <- (target_H - H_start) / lambda_j
        true_pfs_times[p] <- interval_start_time + time_in_interval
      }
    } else {
      # Event occurs within a finite interval
      lambda_j <- hazard_rates[event_interval_idx]
      interval_start_time <- interval_cutpoints[event_interval_idx]
      H_start <- cumulative_hazards_at_cutpoints[event_interval_idx]

      if (lambda_j == 0) {
        if (H_start >= target_H) {
          true_pfs_times[p] <- interval_start_time
        } else {
          true_pfs_times[p] <- Inf
        }
      } else {
        time_in_interval <- (target_H - H_start) / lambda_j
        true_pfs_times[p] <- interval_start_time + time_in_interval
      }
    }
  }

  # Apply censoring
  random_censor_times <- runif(num_patients, min = censor_min_time, max = censor_max_time)
  admin_censor_times  <- rep(max_follow_up, num_patients)
  effective_censor_times <- pmin(random_censor_times, admin_censor_times)

  observed_time <- pmin(true_pfs_times, effective_censor_times)
  event_status  <- as.numeric(true_pfs_times <= effective_censor_times)

  data.frame(
    id = start_id:(start_id + num_patients - 1),
    true_pfs_time = true_pfs_times,
    effective_censor_time = effective_censor_times,
    observed_time = observed_time,
    event_status = event_status
  )
}

# --- HELPER FUNCTION: Calculate Median Survival for a Piecewise Exponential Model ---
# --- UPDATED: Calculate Median Survival for a Piecewise Exponential Model ---
calculate_median_survival_piecewise <- function(hazard_rates, interval_lengths) {
  #' Computes the median survival time for a piecewise exponential model.
  #' If the 0.5 survival is not reached by the end of the last interval,
  #' we *continue past the last cutpoint* with the last interval's hazard (open-ended tail).
  #' Only return Inf if the last hazard is exactly zero.

  if (length(hazard_rates) != length(interval_lengths)) {
    stop("hazard_rates and interval_lengths must have the same length.")
  }

  cumulative_hazard <- 0.0
  current_time <- 0.0

  for (i in seq_along(hazard_rates)) {
    lambda_j <- hazard_rates[i]
    delta_t_j <- interval_lengths[i]
    cumulative_hazard_end_interval <- cumulative_hazard + (lambda_j * delta_t_j)
    # Survival at end of interval:
    # S(t_end) = exp(- cumulative_hazard_end_interval)
    if (cumulative_hazard_end_interval >= log(2)) {
      # Median occurs *within* this interval
      if (lambda_j == 0) {
        # If hazard is zero but we've already crossed log(2) (numerically), clamp to start
        # (this path is practically avoided by the check above unless due to rounding)
        return(current_time)
      }
      time_in_interval_to_median <- (log(2) - cumulative_hazard) / lambda_j
      return(current_time + time_in_interval_to_median)
    }
    # Move to next interval
    cumulative_hazard <- cumulative_hazard_end_interval
    current_time <- current_time + delta_t_j
  }

  # If we get here, the end-of-grid survival is still > 0.5.
  # Extend beyond last cutpoint with last hazard.
  lambda_last <- tail(hazard_rates, 1)
  if (lambda_last > 0) {
    time_beyond_grid <- (log(2) - cumulative_hazard) / lambda_last
    return(current_time + time_beyond_grid)
  } else {
    # Truly zero tail hazard => median is infinite
    return(Inf)
  }
}



# --- HELPER FUNCTION: Draw Posterior Hazard Samples ---
draw_posterior_hazard_samples <- function(
    num_intervals,
    events_per_interval,
    person_time_per_interval,
    prior_alpha_params,
    prior_beta_params,
    num_samples = 1000
) {
  #' Draws samples from the posterior distribution of hazard rates for each interval
  #' in a Bayesian piecewise exponential model with Gamma priors.

  if (!all(sapply(list(events_per_interval, person_time_per_interval, prior_alpha_params, prior_beta_params),
                  length) == num_intervals)) {
    stop("All input vectors must have length equal to num_intervals.")
  }

  posterior_hazard_samples <- matrix(NA, nrow = num_samples, ncol = num_intervals)

  for (j in 1:num_intervals) {
    posterior_alpha_j <- prior_alpha_params[j] + events_per_interval[j]
    posterior_beta_j <- prior_beta_params[j] + person_time_per_interval[j]

    if (posterior_beta_j <= 0) {
      warning(paste("Posterior beta parameter for interval", j, "is non-positive. Setting to a small value (1e-6)."))
      posterior_beta_j <- 1e-6
    }

    posterior_hazard_samples[, j] <- rgamma(num_samples, shape = posterior_alpha_j, rate = posterior_beta_j)
  }

  return(posterior_hazard_samples)
}

final_vsref_probs_abs <- function(med_arm, med_ref, margin_abs) {
  list(
    p_eff_ref = mean(med_arm >  med_ref + margin_abs),
    p_fut_ref = mean(med_arm <= med_ref - margin_abs)
  )
}

# --- HELPER FUNCTION: Calculate Predicted Probability of Success (vs HC) ---
calculate_predicted_success_prob_vs_hc <- function(
    current_patient_data,
    max_total_patients,
    interval_cutpoints,
    prior_alpha_params,
    prior_beta_params,
    num_posterior_draws = 1000,
    median_pfs_success_threshold,
    final_success_posterior_prob_threshold,
    max_follow_up_sim,
    censor_max_time_sim,
    inner_final_draws = 250,
    target_success_threshold = NULL
) {
  current_num_patients <- nrow(current_patient_data)
  num_intervals <- length(interval_cutpoints) - 1
  interval_lengths <- diff(interval_cutpoints)

  cur <- calculate_interval_metrics_fast(current_patient_data, interval_cutpoints)
  cur_post <- draw_posterior_hazard_samples(
    num_intervals,
    cur$events_per_interval, cur$person_time_per_interval,
    prior_alpha_params, prior_beta_params,
    num_samples = num_posterior_draws
  )

  s <- 0L
  N <- num_posterior_draws
  k <- 0L

  for (kk in 1:num_posterior_draws) {
    k <- kk
    lambda_true_k <- cur_post[kk, ]
    n_add <- max_total_patients - current_num_patients

    if (n_add > 0) {
      fut <- simulate_piecewise_exponential_data(
        num_patients = n_add,
        hazard_rates = lambda_true_k,
        interval_cutpoints = interval_cutpoints,
        max_follow_up = max_follow_up_sim,
        censor_max_time = censor_max_time_sim,
        start_id = current_num_patients + 1
      )
      dat <- rbind(current_patient_data, fut)
    } else {
      dat <- current_patient_data
    }

    fin <- calculate_interval_metrics_fast(dat, interval_cutpoints)
    fin_post <- draw_posterior_hazard_samples(
      num_intervals,
      fin$events_per_interval, fin$person_time_per_interval,
      prior_alpha_params, prior_beta_params,
      num_samples = min(inner_final_draws, num_posterior_draws)
    )

    fin_meds <- apply(fin_post, 1, function(h) {
      calculate_median_survival_piecewise(h, interval_lengths)
    })

    succ_k <- mean(fin_meds > median_pfs_success_threshold) >= final_success_posterior_prob_threshold
    s <- s + as.integer(succ_k)

    if (!is.null(target_success_threshold)) {
      n_rem <- N - k
      ub <- (s + n_rem) / N
      lb <- s / N
      if (ub < target_success_threshold || lb > target_success_threshold) break
    }
  }

  # unbiased MC estimate given possible early stop
  s / k
}



# --- NEW HELPER FUNCTION: Calculate Predicted Probability of Efficacy (vs a Reference Arm) ---
# --- Predicted probability vs-ref (UPDATED: fix censoring arg + empty rbind safety) ---
calculate_predicted_prob_vs_ref <- function(
    current_patient_data_arm,
    current_patient_data_ref,
    max_total_patients_arm,
    max_total_patients_ref,
    interval_cutpoints,
    prior_alpha_params,
    prior_beta_params,
    num_posterior_draws = 1000,
    final_efficacy_posterior_prob_threshold,
    final_futility_posterior_prob_threshold,
    max_follow_up_sim,
    censor_max_time_sim,
    compare_arms_futility_margin,
    inner_final_draws = 250,
    target_pred_efficacy_threshold = NULL,
    target_pred_futility_threshold = NULL
) {
  num_intervals <- length(interval_cutpoints) - 1
  interval_lengths <- diff(interval_cutpoints)
  stopifnot(length(prior_alpha_params) == num_intervals,
            length(prior_beta_params)  == num_intervals)

  curA <- calculate_interval_metrics_fast(current_patient_data_arm, interval_cutpoints)
  curR <- calculate_interval_metrics_fast(current_patient_data_ref, interval_cutpoints)

  postA <- draw_posterior_hazard_samples(
    num_intervals, curA$events_per_interval, curA$person_time_per_interval,
    prior_alpha_params, prior_beta_params, num_samples = num_posterior_draws
  )
  postR <- draw_posterior_hazard_samples(
    num_intervals, curR$events_per_interval, curR$person_time_per_interval,
    prior_alpha_params, prior_beta_params, num_samples = num_posterior_draws
  )

  s_eff <- 0L; s_fut <- 0L; N <- num_posterior_draws
  margin_abs <- coalesce_num(compare_arms_futility_margin, 0)

  for (k in 1:num_posterior_draws) {
    lamA <- postA[k, ]; lamR <- postR[k, ]
    nA <- max_total_patients_arm - nrow(current_patient_data_arm)
    nR <- max_total_patients_ref - nrow(current_patient_data_ref)

    # Use *named* args so censor_max_time is correct; keep censor_min_time at 0
    if (nA > 0) {
      futA <- simulate_piecewise_exponential_data(
        num_patients = nA,
        hazard_rates = lamA,
        interval_cutpoints = interval_cutpoints,
        max_follow_up = max_follow_up_sim,
        censor_min_time = 0,
        censor_max_time = censor_max_time_sim,
        start_id = nrow(current_patient_data_arm) + 1
      )
    } else {
      futA <- current_patient_data_arm[FALSE, ]  # empty with same columns
    }

    if (nR > 0) {
      futR <- simulate_piecewise_exponential_data(
        num_patients = nR,
        hazard_rates = lamR,
        interval_cutpoints = interval_cutpoints,
        max_follow_up = max_follow_up_sim,
        censor_min_time = 0,
        censor_max_time = censor_max_time_sim,
        start_id = nrow(current_patient_data_ref) + 1
      )
    } else {
      futR <- current_patient_data_ref[FALSE, ]
    }

    datA <- rbind(current_patient_data_arm, futA)
    datR <- rbind(current_patient_data_ref, futR)

    finA <- calculate_interval_metrics_fast(datA, interval_cutpoints)
    finR <- calculate_interval_metrics_fast(datR, interval_cutpoints)

    finPostA <- draw_posterior_hazard_samples(
      num_intervals, finA$events_per_interval, finA$person_time_per_interval,
      prior_alpha_params, prior_beta_params, num_samples = min(inner_final_draws, num_posterior_draws)
    )
    finPostR <- draw_posterior_hazard_samples(
      num_intervals, finR$events_per_interval, finR$person_time_per_interval,
      prior_alpha_params, prior_beta_params, num_samples = min(inner_final_draws, num_posterior_draws)
    )

    medA <- apply(finPostA, 1, function(h) calculate_median_survival_piecewise(h, interval_lengths))
    medR <- apply(finPostR, 1, function(h) calculate_median_survival_piecewise(h, interval_lengths))

    p_eff_final <- mean(medA >  medR + margin_abs)
    p_fut_final <- mean(medA <= medR - margin_abs)

    s_eff <- s_eff + as.integer(p_eff_final >= final_efficacy_posterior_prob_threshold)
    s_fut <- s_fut + as.integer(p_fut_final >= final_futility_posterior_prob_threshold)

    n_rem <- N - k
    if (!is.null(target_pred_efficacy_threshold)) {
      ub_eff <- (s_eff + n_rem) / N; lb_eff <- s_eff / N
      if (ub_eff < target_pred_efficacy_threshold || lb_eff > target_pred_efficacy_threshold) {
        if (is.null(target_pred_futility_threshold)) break
      }
    }
    if (!is.null(target_pred_futility_threshold)) {
      ub_fut <- (s_fut + n_rem) / N; lb_fut <- s_fut / N
      if (ub_fut < target_pred_futility_threshold || lb_fut > target_pred_futility_threshold) {
        if (is.null(target_pred_efficacy_threshold)) break
      }
    }
    if (!is.null(target_pred_efficacy_threshold) && !is.null(target_pred_futility_threshold)) {
      decided_eff <- ((s_eff + n_rem)/N < target_pred_efficacy_threshold) || (s_eff/N > target_pred_efficacy_threshold)
      decided_fut <- ((s_fut + n_rem)/N < target_pred_futility_threshold) || (s_fut/N > target_pred_futility_threshold)
      if (decided_eff && decided_fut) break
    }
  }

  list(
    predicted_prob_efficacy = s_eff / N,
    predicted_prob_futility = s_fut / N
  )
}



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

gates_pass_for_both_arms <- function(slA, slB, args) {
  min_ev   <- coalesce_num(args$min_events_per_arm, 0)
  min_mfu  <- coalesce_num(args$min_median_followup_per_arm, 0)
  min_pt_f <- coalesce_num(args$min_person_time_frac_per_arm, 0)

  evA  <- coalesce_num(slA$metrics$events_total, 0)
  evB  <- coalesce_num(slB$metrics$events_total, 0)
  mfuA <- coalesce_num(slA$metrics$median_followup, 0)
  mfuB <- coalesce_num(slB$metrics$median_followup, 0)

  ptA  <- coalesce_num(slA$metrics$person_time_total, 0)
  ptB  <- coalesce_num(slB$metrics$person_time_total, 0)

  maxPT_A <- coalesce_num(args$max_total_patients_per_arm[["Doublet"]], 0) * coalesce_num(args$max_follow_up_sim, 0)
  maxPT_B <- coalesce_num(args$max_total_patients_per_arm[["Triplet"]],  0) * coalesce_num(args$max_follow_up_sim, 0)

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
    slC <- slice_arm_data_at_time(state$registries[["Doublet"]], current_time,
                                  args$max_follow_up_sim, args$interval_cutpoints_sim)
    slT <- slice_arm_data_at_time(state$registries[["Triplet"]],  current_time,
                                  args$max_follow_up_sim, args$interval_cutpoints_sim)

    if (!gates_pass_for_both_arms(slC, slT, args)) {
      if (diagnostics) message(sprintf("[t=%.2f] vsREF gated out", current_time))
      return(state)
    }

    pr_eff <- calculate_current_prob_vs_ref(slC, slT, args)
    pr_fut <- calculate_current_prob_vs_ref_futility(slC, slT, args)

    if (diagnostics) {
      message(sprintf("[t=%.2f] vsREF PrEff>=%.2f: %.3f | PrFut>=%.2f (m=%.2f): %.3f",
                      current_time,
                      coalesce_num(args$efficacy_threshold_vs_ref_prob, NA_real_), pr_eff,
                      coalesce_num(args$futility_threshold_vs_ref_prob,  NA_real_), coalesce_num(args$compare_arms_futility_margin, 0), pr_fut))
    }

    # Success?
    if (!is.null(args$efficacy_threshold_vs_ref_prob) && is.finite(args$efficacy_threshold_vs_ref_prob) &&
        pr_eff >= args$efficacy_threshold_vs_ref_prob) {
      state$arm_status["Triplet"] <- "stopped_efficacy"
      state$stop_efficacy_per_sim_row["Triplet"] <- 1L
      state$sim_final_n_current_run["Triplet"]   <- state$enrolled_counts["Triplet"]
      return(state)
    }

    # Futility?
    if (!is.null(args$futility_threshold_vs_ref_prob) && is.finite(args$futility_threshold_vs_ref_prob) &&
        pr_fut >= args$futility_threshold_vs_ref_prob) {
      state$arm_status["Triplet"] <- "stopped_futility"
      state$stop_futility_per_sim_row["Triplet"] <- 1L
      state$sim_final_n_current_run["Triplet"]   <- state$enrolled_counts["Triplet"]
      return(state)
    }

    return(state)
  }

  # HC path (not compare_arms)
  # ... your existing single-arm logic that mutates state ...
  state
}




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

# Approximate future person-time per interval (no patient sim)
approx_future_exposure <- function(current_time,
                                   enroll_times,
                                   max_total_patients,
                                   overall_accrual_rate,
                                   interval_cutpoints,
                                   max_follow_up) {
  # Expected future exposure for remaining patients, on the patient-time clock,
  # capped at max_follow_up and binned into [cuts[j], cuts[j+1)).
  cuts <- interval_cutpoints
  K <- length(cuts) - 1
  n_enrolled <- length(enroll_times)
  remaining_expected <- max(0, max_total_patients - n_enrolled)
  if (remaining_expected == 0) return(numeric(K))

  # Each future patient contributes the same shape of follow-up up to max_follow_up
  interval_follow_up <- numeric(K)
  for (j in 1:K) {
    lo <- cuts[j]; hi <- cuts[j+1]
    interval_follow_up[j] <- max(0, min(hi, max_follow_up) - lo)
  }

  remaining_expected * interval_follow_up
}



calculate_predicted_success_prob_vs_hc_fast <- function(
    current_patient_data,
    current_enroll_times,
    max_total_patients,
    interval_cutpoints,
    prior_alpha_params,
    prior_beta_params,
    num_posterior_draws = 400,
    median_pfs_success_threshold,
    final_success_posterior_prob_threshold,
    max_follow_up_sim,
    overall_accrual_rate,
    current_time
) {
  K <- length(interval_cutpoints) - 1
  stopifnot(length(prior_alpha_params) == K, length(prior_beta_params) == K)

  m <- calculate_interval_metrics_fast(current_patient_data, interval_cutpoints)
  alpha0 <- prior_alpha_params + m$events_per_interval
  beta0  <- prior_beta_params  + m$person_time_per_interval

  T_future <- approx_future_exposure(
    current_time,
    enroll_times = current_enroll_times,
    max_total_patients = max_total_patients,
    overall_accrual_rate = overall_accrual_rate,
    interval_cutpoints = interval_cutpoints,
    max_follow_up = max_follow_up_sim
  )

  interval_lengths <- diff(interval_cutpoints)
  hit <- 0L

  for (k in seq_len(num_posterior_draws)) {
    lambda <- rgamma(K, shape = alpha0, rate = beta0)
    ev_future <- ifelse(T_future > 0, rpois(K, lambda * T_future), 0L)

    alpha_fin <- alpha0 + ev_future
    beta_fin  <- beta0  + T_future

    lambda_fin <- rgamma(K, shape = alpha_fin, rate = beta_fin)
    med_fin <- calculate_median_survival_piecewise(lambda_fin, interval_lengths)

    hit <- hit + as.integer(med_fin > median_pfs_success_threshold)
  }

  hit / num_posterior_draws
}



near_boundary <- function(prob_now, target, band = 0.10) {
  abs(prob_now - target) <= band
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




pretty_scenario_matrix <- function(results_df) {
  # Validate
  if (!"scenario" %in% names(results_df)) stop("results_df must include 'scenario'")
  if (!"Arm_Name" %in% names(results_df)) stop("results_df must include 'Arm_Name'")

  # Aggregate summary by scenario + arm
  tbl <- data.table::as.data.table(results_df)[, .(
    True_Median         = mean(True_Median, na.rm = TRUE),
    Max_Planned_N       = max(Exp_N, na.rm = TRUE),
    Exp_N               = mean(Exp_N, na.rm = TRUE),
    Pr_Reach_Max_N      = mean(Pr_Reach_Max_N, na.rm = TRUE),
    Type_I_Error_or_Power = mean(Type_I_Error_or_Power, na.rm = TRUE),
    PET_Efficacy        = mean(PET_Efficacy, na.rm = TRUE),
    PET_Futility        = mean(PET_Futility, na.rm = TRUE),
    Pr_Final_Efficacy   = mean(Pr_Final_Efficacy, na.rm = TRUE),
    Pr_Final_Futility   = mean(Pr_Final_Futility, na.rm = TRUE)
  ), by = .(scenario, Arm_Name)]

  # Pivot wide: one row per scenario, arms side-by-side
  wide <- data.table::dcast(
    tbl,
    scenario ~ Arm_Name,
    value.var = c("True_Median", "Exp_N", "Pr_Reach_Max_N",
                  "Type_I_Error_or_Power", "PET_Efficacy", "PET_Futility",
                  "Pr_Final_Efficacy", "Pr_Final_Futility"),
    fun.aggregate = mean
  )

  # Optional rounding for cleaner display
  num_cols <- names(wide)[sapply(wide, is.numeric)]
  wide[, (num_cols) := lapply(.SD, function(x) round(x, 2)), .SDcols = num_cols]

  # Return formatted data.frame (for printing)
  as.data.frame(wide)
}

export_scenario_table_to_excel <- function(pretty_tbl, file_path = "scenario_summary.xlsx") {
  library(openxlsx)
  wb <- createWorkbook()
  addWorksheet(wb, "Scenario Summary")
  writeDataTable(wb, sheet = 1, x = pretty_tbl, tableStyle = "TableStyleMedium9")
  setColWidths(wb, sheet = 1, cols = 1:ncol(pretty_tbl), widths = "auto")
  saveWorkbook(wb, file_path, overwrite = TRUE)
  message("✅ Exported formatted table to: ", normalizePath(file_path))
}

export_scenario_table_to_png <- function(results_df,
                                         file_path = "scenario_summary.png",
                                         title = "Bayesian Adaptive Design Summary",
                                         subtitle = NULL,
                                         highlight_arm = "Triplet",
                                         snapshot_engine = c("auto", "webshot2", "webshot"),
                                         vwidth = 1400,
                                         vheight = 600,
                                         zoom = 1) {
  if (!requireNamespace("gt", quietly = TRUE)) stop("Please install.packages('gt')")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please install.packages('dplyr')")

  snapshot_engine <- match.arg(snapshot_engine)

  # engine auto-detect
  if (snapshot_engine == "auto") {
    if (requireNamespace("webshot2", quietly = TRUE)) {
      snapshot_engine <- "webshot2"
    } else if (requireNamespace("webshot", quietly = TRUE)) {
      snapshot_engine <- "webshot"
    } else {
      stop("Install either 'webshot2' (Chrome/Chromium) or 'webshot' (PhantomJS).")
    }
  }

  # build table data (expects you already defined pretty_scenario_matrix)
  tbl <- pretty_scenario_matrix(results_df)

  # minimal formatting
  library(gt)
  library(dplyr)

  tbl <- tbl %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 2)))

  gt_tbl <- gt::gt(tbl) %>%
    gt::tab_header(
      title = gt::md(paste0("**", title, "**")),
      subtitle = subtitle
    ) %>%
    gt::fmt_number(columns = where(is.numeric), decimals = 2) %>%
    gt::opt_align_table_header("left") %>%
    gt::tab_options(
      table.font.size = 12,
      data_row.padding = gt::px(4),
      heading.background.color = "#f0f0f0",
      column_labels.background.color = "#f7f7f7",
      table.border.top.width = gt::px(1),
      table.border.bottom.width = gt::px(1),
      table.border.bottom.color = "gray"
    )

  if (!is.null(highlight_arm)) {
    highlight_cols <- grep(highlight_arm, names(tbl), value = TRUE)
    if (length(highlight_cols) > 0) {
      gt_tbl <- gt_tbl %>%
        gt::tab_style(
          style = gt::cell_fill(color = "#e8f4f8"),
          locations = gt::cells_body(columns = dplyr::all_of(highlight_cols))
        )
    }
  }

  # branch by engine
  if (snapshot_engine == "webshot2") {
    # Directly save PNG (gt uses webshot2/chromote internally)
    gt::gtsave(
      data = gt_tbl,
      filename = file_path,
      vwidth = vwidth,
      vheight = vheight,
      zoom = zoom
    )
  } else {
    # Save HTML then rasterize with webshot (PhantomJS)
    if (!requireNamespace("webshot", quietly = TRUE)) {
      stop("snapshot_engine='webshot' requires the 'webshot' package.")
    }
    html_tmp <- sub("\\.png$", ".html", file_path, ignore.case = TRUE)
    gt::gtsave(
      data = gt_tbl,
      filename = html_tmp,
      vwidth = vwidth,
      vheight = vheight,
      zoom = zoom
    )
    webshot::webshot(
      url = html_tmp,
      file = file_path,
      vwidth = vwidth,
      vheight = vheight,
      zoom = zoom
    )
  }

  message("✅ Image saved: ", normalizePath(file_path))
}


calibrate_alpha <- function(base_args, scens_null, thr_grid_interim = c(0.9, 0.95, 0.975),
                            thr_grid_final = c(0.95, 0.975, 0.99), sims = 300) {
  best <- NULL
  base_args$num_simulations <- sims
  for (ti in thr_grid_interim) for (tf in thr_grid_final) {
    args <- base_args
    args$efficacy_threshold_current_prob_hc <- ti
    args$final_success_posterior_prob_threshold <- tf
    res <- run_scenarios(args, scens_null, parallel = TRUE, seed = 123)
    # Type I = early+final efficacy for the experimental arm(s) under null scenario
    typeI <- res[res$Arm_Name != args$reference_arm_name & res$scenario == 1, "Type_I_Error_or_Power"]
    alpha_hat <- mean(typeI)
    cand <- list(ti = ti, tf = tf, alpha = alpha_hat)
    if (is.null(best) || (alpha_hat < best$alpha && alpha_hat <= 0.10)) best <- cand
    cat(sprintf("Interim=%.3f Final=%.3f => alpha=%.3f\n", ti, tf, alpha_hat))
  }
  best
}

# === Grid search for Type I / Power vs thresholds & margin ===
# - Sweeps interim success threshold (vs HC), final success threshold, and Triplet success margin
# - Evaluates Type I under all-null (Doublet=6, Triplet=6) and Power under alt (Doublet=6, Triplet=9)
# - Ranks designs that satisfy alpha <= target_alpha by highest power, then lowest Exp_N under alt

library(data.table)

grid_calibrate <- function(base_args,
                           null_med = 6,
                           alt_med  = 9,
                           margins_abs = c(1, 2, 3),          # success margin for Triplet in months
                           interim_thr_grid = c(0.90, 0.95),  # Pr(success|data) at interim vs HC
                           final_thr_grid   = c(0.95, 0.975, 0.99), # Pr(success|final posterior)
                           sims = 400,                         # increase for stability
                           target_alpha = 0.10,
                           seed = 123,
                           parallel = TRUE) {

  # Build two scenarios: (1) all-null, (2) Triplet=alt
  scens <- scenarios_from_grid(list(
    weibull_median_true_arms = list(
      c(Doublet = null_med, Triplet = null_med),
      c(Doublet = null_med, Triplet = alt_med)
    )
  ))

  combos <- CJ(margin = margins_abs,
               interim_thr = interim_thr_grid,
               final_thr   = final_thr_grid)

  results <- vector("list", nrow(combos))

  for (i in seq_len(nrow(combos))) {
    m  <- combos$margin[i]
    ti <- combos$interim_thr[i]
    tf <- combos$final_thr[i]

    args_i <- base_args

    # --- thresholds being calibrated ---
    args_i$efficacy_threshold_current_prob_hc <- ti
    args_i$final_success_posterior_prob_threshold <- tf

    # --- require superiority margin for Triplet only (Doublet stays at null) ---
    args_i$median_pfs_success_threshold_arms <- c(Doublet = null_med, Triplet = null_med + m)

    # --- keep predictive disabled during calibration (safer for Type I) ---
    args_i$predictive_fast <- FALSE

    # --- set sims & seed ---
    args_i$num_simulations <- sims

    # run
    res <- run_scenarios(args_i, scens, parallel = parallel, seed = seed)

    # Pull alpha (scenario 1, Triplet arm) & power (scenario 2, Triplet arm)
    r_null <- res[res$scenario == 1 & res$Arm_Name == "Triplet", ]
    r_alt  <- res[res$scenario == 2 & res$Arm_Name == "Triplet", ]

    alpha_hat <- mean(r_null$Type_I_Error_or_Power)
    power_hat <- mean(r_alt$Type_I_Error_or_Power)

    out <- data.table(
      margin_abs = m,
      interim_thr = ti,
      final_thr = tf,
      alpha = alpha_hat,
      power = power_hat,
      ExpN_null = mean(r_null$Exp_N),
      ExpN_alt  = mean(r_alt$Exp_N),
      PET_Eff_null = mean(r_null$PET_Efficacy),
      PET_Eff_alt  = mean(r_alt$PET_Efficacy),
      PET_Fut_null = mean(r_null$PET_Futility),
      PET_Fut_alt  = mean(r_alt$PET_Futility)
    )

    results[[i]] <- out
  }

  grid <- rbindlist(results)

  # Designs meeting the alpha target
  ok <- grid[alpha <= target_alpha]
  setorder(ok, -power, ExpN_alt, margin_abs, interim_thr, final_thr)

  list(all = grid[order(margin_abs, interim_thr, final_thr)],
       feasible = ok,
       top = head(ok, 10))
}

# =========================
# Power vs Type I plots
# =========================
# Requires: data.table, ggplot2
# Optional: ggrepel (for nicer labels), scales (percent axes)

library(data.table)
library(ggplot2)

plot_calibration <- function(cal,
                             target_alpha = 0.10,
                             label_top_n = 3) {
  stopifnot(is.list(cal), !is.null(cal$all))
  df <- as.data.table(cal$all)
  if (nrow(df) == 0L) stop("cal$all is empty")

  # Facet label helpers
  df[, margin_lab := sprintf("Δ = %s", format(margin_abs, trim = TRUE))]
  df[, final_lab  := paste0("Final thr = ", final_thr)]
  df[, interim_lab := paste0("Interim thr=", interim_thr)]

  # Pareto frontier within each facet (margin x final_thr)
  setorder(df, margin_abs, final_lab, alpha)
  df[, frontier := power == cummax(power), by = .(margin_abs, final_lab)]

  # best feasible (alpha <= target) to annotate
  best_tbl <- if (!is.null(cal$feasible) && nrow(cal$feasible) > 0L) {
    bt <- as.data.table(cal$feasible)
    setorder(bt, -power, ExpN_alt)
    bt[, margin_lab := sprintf("Δ = %s", format(margin_abs, trim = TRUE))]
    bt[, final_lab  := paste0("Final thr = ", final_thr)]
    bt[1:min(label_top_n, .N)]
  } else data.table()

  p <- ggplot(df, aes(x = alpha, y = power,
                      color = margin_lab,
                      shape = factor(interim_thr))) +
    geom_hline(yintercept = 0, linewidth = 0.2, color = "grey70") +
    geom_vline(xintercept = target_alpha, linetype = 2, linewidth = 0.4) +
    geom_point(size = 2, alpha = 0.9) +
    # frontier lines
    geom_line(data = df[frontier == TRUE],
              aes(group = interaction(margin_lab, final_lab)),
              linewidth = 0.7, alpha = 0.8) +
    facet_wrap(~ final_lab) +
    labs(
      x = "Type I error (alpha)",
      y = "Power",
      color = "Success margin",
      shape = "Interim threshold",
      title = "Design calibration: Power vs Type I",
      subtitle = sprintf("Vertical dashed line = target alpha (%g)", target_alpha)
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")

  # Annotate top feasible designs (if any)
  if (nrow(best_tbl) > 0L) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_label_repel(
        data = best_tbl,
        aes(x = alpha, y = power,
            label = paste0("Δ=", margin_abs,
                           "; Int=", interim_thr,
                           "; Fin=", final_thr,
                           "; N~", round(ExpN_alt, 1))),
        size = 3, label.size = 0.15, alpha = 0.9,
        fill = "white", label.padding = unit(0.15, "lines"),
        max.overlaps = Inf
      )
    } else {
      p <- p + geom_text(
        data = best_tbl,
        aes(x = alpha, y = power,
            label = paste0("Δ=", margin_abs,
                           ", Int=", interim_thr,
                           ", Fin=", final_thr)),
        size = 3, vjust = -0.8
      )
    }
  }

  # Optional percent formatting if 'scales' is available
  if (requireNamespace("scales", quietly = TRUE)) {
    p <- p +
      scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                         limits = c(0, 1))
  } else {
    p <- p + coord_cartesian(ylim = c(0, 1))
  }

  p
}
# --- Pick best calibration row and bake into args ---
adopt_calibration <- function(cal, base_args, null_med, alt_med, which = 1L) {
  stopifnot(is.list(cal), !is.null(cal$top), nrow(cal$top) >= which)
  best <- data.table::as.data.table(cal$top)[which]

  args_star <- base_args
  # success thresholds from calibration
  args_star$efficacy_threshold_current_prob_hc     <- best$interim_thr
  args_star$final_success_posterior_prob_threshold <- best$final_thr
  # set Triplet superiority margin (absolute months)
  args_star$median_pfs_success_threshold_arms <- c(
    Doublet = null_med,
    Triplet = null_med + best$margin_abs
  )
  # keep predictive off during calibration/exploration (safer for Type I)
  args_star$predictive_fast <- FALSE

  # 2-scenario grid: all-null vs Triplet=alt
  scens2 <- scenarios_from_grid(list(
    weibull_median_true_arms = list(
      c(Doublet = null_med, Triplet = null_med),
      c(Doublet = null_med, Triplet = alt_med)
    )
  ))

  list(args_star = args_star, pick = best, scens2 = scens2)
}
# --- Explore early-stopping & information gates around a calibrated design ---
# Sweeps: futility threshold, min events, min median FU, and calendar look beat
# Explore early stopping while varying the futility comparator for *all* experimental arms
# Comparator options:
#   base = "null+delta": uses P(median < null + delta | data) >= fut_thr  => stop
#   base = "alt"       : uses P(median < alt         | data) >= fut_thr  => stop
#
# Notes:
# - Uses adopt_calibration() to import the calibrated success/margin from `cal` (your prior step).
# - Keeps predictive_fast = FALSE while exploring early-stopping knobs (clean Type I accounting).
#
# --- Explore early-stopping & information gates around a calibrated design ---
# Adds sweeps over:
#   * per-arm gates (min_events_per_arm, min_median_followup_per_arm, min_person_time_frac_per_arm)
#   * schedule type: "calendar" (beats) OR "persontime" (milestones with a calendar backstop)
#
# For person-time schedules, pass a LIST of milestone vectors via `pt_milestones_choices`,
# e.g. list(c(0.30,0.45,0.60,0.80,1.00), c(0.30,0.60,0.90)).
#
explore_early_stopping_from_cal <- function(
    cal,
    base_args,
    null_med,
    alt_med,
    base                  = c("null+delta","alt"),
    futility_delta_grid   = c(0, 1, 2, 3),    # only used when base = "null+delta"
    fut_thr_grid          = c(0.6, 0.7, 0.8, 0.9),

    # legacy/global gates
    min_events_grid       = c(12, 18),
    min_medFU_grid        = c(3, 4.5),

    # schedule sweep
    schedule_modes        = c("calendar","persontime"),
    beat_grid             = c(3, 6),          # used when schedule="calendar"

    # person-time schedule knobs
    pt_milestones_choices = list(c(0.30,0.45,0.60,0.80,1.00)),
    latest_calendar_look_grid = c(Inf),       # e.g., c(12, 18) months as a backstop

    # per-arm gates (NEW)
    min_events_per_arm_grid          = c(8, 12),
    min_median_followup_per_arm_grid = c(0, 4.5),
    min_person_time_frac_per_arm_grid= c(0.00, 0.25),

    sims                 = 400,
    seed                 = 123,
    parallel             = (.Platform$OS.type == "unix")
) {
  base <- match.arg(base)
  stopifnot(is.list(cal), !is.null(cal$top) || !is.null(cal$feasible))

  # 1) import calibrated success thresholds/margin & the two scenarios (null, alt)
  adopted  <- adopt_calibration(cal, base_args, null_med = null_med, alt_med = alt_med, which = 1L)
  args0    <- adopted$args_star
  scens2   <- adopted$scens2
  exp_arms <- exp_arms_from_args(args0)

  # helper to stringify milestone vectors for the output table
  fmt_frac_vec <- function(x) if (is.null(x)) "NULL" else paste0(sprintf("%.2f", x), collapse = ",")

  # build the list of combinations manually (because milestones is a list-column)
  combos <- list()
  for (ft in fut_thr_grid) {
    # choose the futility comparator set for exp arms
    fut_deltas <- if (base == "null+delta") futility_delta_grid else 0
    for (fd in fut_deltas) {
      for (ev in min_events_grid) {
        for (mu in min_medFU_grid) {
          for (me in min_events_per_arm_grid) {
            for (mfu in min_median_followup_per_arm_grid) {
              for (mpt in min_person_time_frac_per_arm_grid) {

                # schedule: calendar beats
                if ("calendar" %in% schedule_modes) {
                  for (bt in beat_grid) {
                    combos[[length(combos) + 1L]] <- data.table::data.table(
                      schedule = "calendar",
                      fut_base = base,
                      fut_delta = fd,
                      fut_thr = ft,
                      min_events = ev,
                      min_medFU = mu,
                      beat = bt,
                      pt_milestones = NA_character_,
                      latest_calendar_look = NA_real_,
                      min_events_per_arm = me,
                      min_median_followup_per_arm = mfu,
                      min_person_time_frac_per_arm = mpt
                    )
                  }
                }

                # schedule: person-time milestones
                if ("persontime" %in% schedule_modes) {
                  for (ml in seq_along(pt_milestones_choices)) {
                    mlv <- pt_milestones_choices[[ml]]
                    stopifnot(is.numeric(mlv), all(mlv > 0 & mlv <= 1))
                    for (lc in latest_calendar_look_grid) {
                      combos[[length(combos) + 1L]] <- data.table::data.table(
                        schedule = "persontime",
                        fut_base = base,
                        fut_delta = fd,
                        fut_thr = ft,
                        min_events = ev,
                        min_medFU = mu,
                        beat = NA_real_,
                        pt_milestones = fmt_frac_vec(mlv),
                        latest_calendar_look = lc,
                        min_events_per_arm = me,
                        min_median_followup_per_arm = mfu,
                        min_person_time_frac_per_arm = mpt
                      )
                    }
                  }
                }

              } # mpt
            } # mfu
          } # me
        } # mu
      } # ev
    } # fd
  } # ft

  combos_dt <- data.table::rbindlist(combos, use.names = TRUE)
  out <- vector("list", nrow(combos_dt))

  # 2) iterate & simulate
  for (i in seq_len(nrow(combos_dt))) {
    row <- combos_dt[i]

    args_i <- args0

    # set the *interim futility* comparator for experimental arms
    args_i <- set_futility_medians(
      args    = args_i,
      null_med = null_med,
      alt_med  = alt_med,
      base     = row$fut_base,
      delta    = row$fut_delta
    )
    args_i$posterior_futility_threshold_hc <- row$fut_thr

    # global gates
    args_i$min_events_for_analysis <- row$min_events
    args_i$min_median_followup     <- row$min_medFU

    # per-arm gates
    args_i$min_events_per_arm             <- row$min_events_per_arm
    args_i$min_median_followup_per_arm    <- row$min_median_followup_per_arm
    args_i$min_person_time_frac_per_arm   <- row$min_person_time_frac_per_arm

    # schedule
    if (row$schedule == "calendar") {
      args_i$interim_calendar_beat   <- row$beat
      args_i$person_time_milestones  <- NULL
      args_i$latest_calendar_look    <- Inf
    } else {
      # person-time milestones
      mlv <- as.numeric(strsplit(row$pt_milestones, ",")[[1]])
      args_i$person_time_milestones <- mlv
      args_i$latest_calendar_look   <- row$latest_calendar_look
      # keep a (small) calendar beat as a guard if you like, but we’ll let the backstop rule dominate
      args_i$interim_calendar_beat  <- args0$interim_calendar_beat %||% 2
    }

    # clean exploration (predictive off)
    args_i$predictive_fast <- FALSE
    args_i$num_simulations <- sims

    # run null vs alt
    res <- run_scenarios(args_i, adopted$scens2, parallel = parallel, seed = seed)

    r_null <- res[res$scenario == 1 & res$Arm_Name %in% exp_arms, ]
    r_alt  <- res[res$scenario == 2 & res$Arm_Name %in% exp_arms, ]

    out[[i]] <- data.table::data.table(
      schedule     = row$schedule,
      fut_base     = row$fut_base,
      fut_delta    = row$fut_delta,
      fut_thr      = row$fut_thr,
      min_events   = row$min_events,
      min_medFU    = row$min_medFU,
      beat         = ifelse(row$schedule == "calendar", row$beat, NA_real_),
      pt_milestones = ifelse(row$schedule == "persontime", row$pt_milestones, NA_character_),
      latest_calendar_look = ifelse(row$schedule == "persontime", row$latest_calendar_look, NA_real_),

      min_events_per_arm           = row$min_events_per_arm,
      min_median_followup_per_arm  = row$min_median_followup_per_arm,
      min_person_time_frac_per_arm = row$min_person_time_frac_per_arm,

      alpha        = mean(r_null$Type_I_Error_or_Power),
      power        = mean(r_alt$Type_I_Error_or_Power),
      ExpN_null    = mean(r_null$Exp_N),
      ExpN_alt     = mean(r_alt$Exp_N),
      PET_Eff_null = mean(r_null$PET_Efficacy),
      PET_Eff_alt  = mean(r_alt$PET_Efficacy),
      PET_Fut_null = mean(r_null$PET_Futility),
      PET_Fut_alt  = mean(r_alt$PET_Futility)
    )
  }

  early <- data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
  data.table::setorder(
    early,
    schedule, fut_base, fut_delta, fut_thr,
    min_events, min_medFU,
    min_events_per_arm, min_median_followup_per_arm, min_person_time_frac_per_arm,
    beat, pt_milestones, latest_calendar_look
  )
  early[]
}

# --- Quick plot (optional) ---
plot_early_tradeoff <- function(early_df,
                                target_alpha = 0.10,
                                fix_min_ev = NULL,
                                fix_mfu    = NULL,
                                fix_beat   = NULL) {
  df <- data.table::as.data.table(early_df)
  if (!is.null(fix_min_ev)) df <- df[min_events == fix_min_ev]
  if (!is.null(fix_mfu))    df <- df[min_medFU  == fix_mfu]
  if (!is.null(fix_beat))   df <- df[beat       == fix_beat]
  if (nrow(df) == 0L) stop("No rows after filtering—relax the fixed settings.")

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    ggplot2::ggplot(df, ggplot2::aes(alpha, power, color = factor(fut_thr))) +
      ggplot2::geom_vline(xintercept = target_alpha, linetype = 2, linewidth = 0.4) +
      ggplot2::geom_point(size = 2) +
      ggplot2::geom_path(ggplot2::aes(group = fut_thr), linewidth = 0.6, alpha = 0.7) +
      ggplot2::facet_grid(min_events ~ min_medFU, labeller = "label_both") +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::labs(x = "Type I error", y = "Power", color = "Futility thr")
  } else {
    message("Install ggplot2 to plot. Returning data.")
    df
  }
}

# --- Filter + rank early-stopping designs by constraints ---
filter_early_grid <- function(early_df,
                              alpha_cap   = 0.10,
                              power_floor = 0.70) {
  df <- data.table::as.data.table(early_df)
  ok <- df[alpha <= alpha_cap & power >= power_floor]
  if (nrow(ok) == 0L) return(ok)

  # Rank: smallest N under alt, then more futility under null, then less efficacy under null
  data.table::setorder(ok, ExpN_alt, -PET_Fut_null, PET_Eff_null, fut_thr, min_events, min_medFU, beat)
  ok[]
}

# --- Pick the single recommended design under constraints ---
recommend_design_from_early <- function(df,
                                        alpha_cap = 0.10,
                                        power_floor = 0.80,
                                        pet_fut_cap = NULL) {
  df_filt <- df[alpha <= alpha_cap & power >= power_floor]
  if (!is.null(pet_fut_cap)) {
    df_filt <- df_filt[PET_Fut_alt <= pet_fut_cap]
  }
  if (nrow(df_filt) == 0)
    stop("No designs meet specified criteria.")

  df_filt[order(-power, PET_Fut_alt, ExpN_alt)][1]
}

# --- Bake the recommended early-stopping knobs back into args ---
apply_recommended_to_args <- function(args_star, rec_row) {
  stopifnot(nrow(rec_row) == 1L)
  a <- args_star
  a$posterior_futility_threshold_hc <- rec_row$fut_thr
  a$min_events_for_analysis         <- rec_row$min_events
  a$min_median_followup             <- rec_row$min_medFU
  a$interim_calendar_beat           <- rec_row$beat

  # NEW: carry per-arm gates if present
  if ("min_events_per_arm" %in% names(rec_row)) {
    a$min_events_per_arm <- rec_row$min_events_per_arm
  }
  if ("min_median_followup_per_arm" %in% names(rec_row)) {
    a$min_median_followup_per_arm <- rec_row$min_median_followup_per_arm
  }
  if ("min_person_time_frac_per_arm" %in% names(rec_row)) {
    a$min_person_time_frac_per_arm <- rec_row$min_person_time_frac_per_arm
  }
  a
}

# All non-reference arms are considered "experimental"
exp_arms_from_args <- function(args) {
  setdiff(args$arm_names, args$reference_arm_name)
}

# Modify args so interim *futility* compares to either:
#   - null + delta  (base = "null+delta"; delta is a single number for all exp arms)
#   - alt           (base = "alt")
# The HC/reference arm keeps its "null" futility threshold.
set_futility_medians <- function(args, null_med, alt_med, base = c("null+delta","alt"), delta = 0) {
  base <- match.arg(base)
  arms_exp <- exp_arms_from_args(args)
  fut <- args$futility_median_arms
  if (base == "null+delta") {
    fut[arms_exp] <- null_med + delta
  } else if (base == "alt") {
    fut[arms_exp] <- alt_med
  }
  args$futility_median_arms <- fut
  args
}
# --- Predicted probability vs-ref (FAST) (UPDATED: add argument checks) ---
calculate_predicted_prob_vs_ref_fast <- function(
    current_patient_data_arm,
    current_patient_data_ref,
    enroll_times_arm,
    enroll_times_ref,
    max_total_patients_arm,
    max_total_patients_ref,
    interval_cutpoints,
    prior_alpha_params,
    prior_beta_params,
    num_posterior_draws = 400,
    final_efficacy_posterior_prob_threshold,
    final_futility_posterior_prob_threshold,
    compare_arms_futility_margin,
    max_follow_up_sim,
    overall_accrual_rate,
    current_time
) {
  K <- length(interval_cutpoints) - 1
  L <- diff(interval_cutpoints)
  stopifnot(length(prior_alpha_params) == K,
            length(prior_beta_params)  == K)

  ma <- calculate_interval_metrics_fast(current_patient_data_arm, interval_cutpoints)
  mr <- calculate_interval_metrics_fast(current_patient_data_ref, interval_cutpoints)

  alpha_a0 <- prior_alpha_params + ma$events_per_interval
  beta_a0  <- prior_beta_params  + ma$person_time_per_interval
  alpha_r0 <- prior_alpha_params + mr$events_per_interval
  beta_r0  <- prior_beta_params  + mr$person_time_per_interval

  T_future_a <- approx_future_exposure(current_time, enroll_times_arm,
                                       max_total_patients_arm, overall_accrual_rate,
                                       interval_cutpoints, max_follow_up_sim)
  T_future_r <- approx_future_exposure(current_time, enroll_times_ref,
                                       max_total_patients_ref, overall_accrual_rate,
                                       interval_cutpoints, max_follow_up_sim)

  eff_hits <- 0L; fut_hits <- 0L
  margin_abs <- coalesce_num(compare_arms_futility_margin, 0)

  for (k in seq_len(num_posterior_draws)) {
    lam_a <- rgamma(K, shape = alpha_a0, rate = beta_a0)
    lam_r <- rgamma(K, shape = alpha_r0, rate = beta_r0)

    ev_a <- ifelse(T_future_a > 0, rpois(K, lam_a * T_future_a), 0L)
    ev_r <- ifelse(T_future_r > 0, rpois(K, lam_r * T_future_r), 0L)

    alpha_af <- alpha_a0 + ev_a; beta_af <- beta_a0 + T_future_a
    alpha_rf <- alpha_r0 + ev_r; beta_rf <- beta_r0 + T_future_r

    lam_af <- rgamma(K, shape = alpha_af, rate = beta_af)
    lam_rf <- rgamma(K, shape = alpha_rf, rate = beta_rf)

    med_a <- calculate_median_survival_piecewise(lam_af, L)
    med_r <- calculate_median_survival_piecewise(lam_rf, L)

    p_eff <- as.integer( mean(med_a >  med_r + margin_abs) >= final_efficacy_posterior_prob_threshold )
    p_fut <- as.integer( mean(med_a <= med_r - margin_abs) >= final_futility_posterior_prob_threshold )

    eff_hits <- eff_hits + p_eff
    fut_hits <- fut_hits + p_fut
  }

  list(
    predicted_prob_efficacy = eff_hits / num_posterior_draws,
    predicted_prob_futility = fut_hits / num_posterior_draws
  )
}
