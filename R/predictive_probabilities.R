
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
  
  s_eff <- 0L; s_fut <- 0L; N <- num_posterior_draws; k_eff <- 0L
  margin_abs <- coalesce_num(compare_arms_futility_margin, 0)
  
  for (k in 1:num_posterior_draws) {
    k_eff <- k
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
  
  if (k_eff == 0L) {
    k_eff <- 1L
  }

  list(
    predicted_prob_efficacy = s_eff / k_eff,
    predicted_prob_futility = s_fut / k_eff
  )
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

