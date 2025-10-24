library(survminer)
library(survival)
library(data.table)
library(progress)

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

