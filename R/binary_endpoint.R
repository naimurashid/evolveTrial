# binary_endpoint.R
# Binary endpoint support for evolveTrial
# Enables two-stage Simon-type designs with binary (response/no-response) outcomes

# =============================================================================
# Data Generation
# =============================================================================

#' Simulate binary response data for a cohort
#'
#' Generates binary response outcomes for n patients with true response probability p.
#'
#' @param n Number of patients
#' @param p True response probability (0 to 1)
#' @param start_id Starting patient ID (default 1)
#'
#' @return Data frame with columns: id, response (0/1)
#' @keywords internal
simulate_binary_response_data <- function(n, p, start_id = 1L) {

  if (n <= 0) {
    return(data.frame(id = integer(0), response = integer(0)))
  }

  responses <- rbinom(n, size = 1, prob = p)

  data.frame(
    id = seq(from = start_id, length.out = n),
    response = as.integer(responses)
  )
}

# =============================================================================
# Posterior Sampling (Beta-Binomial Conjugacy)
# =============================================================================

#' Draw posterior samples for binary response rate
#'
#' Uses Beta-Binomial conjugacy to sample from the posterior distribution
#' of the response rate.
#'
#' Prior: Beta(alpha_prior, beta_prior)
#' Likelihood: Binomial(n, p)
#' Posterior: Beta(alpha_prior + successes, beta_prior + failures)
#'
#' @param n_responses Number of responders
#' @param n_total Total number of patients
#' @param alpha_prior Beta prior shape1 parameter (default 1 for uniform)
#' @param beta_prior Beta prior shape2 parameter (default 1 for uniform)
#' @param num_samples Number of posterior samples to draw
#'
#' @return Numeric vector of posterior samples for response rate
#' @export
draw_posterior_response_rate <- function(
    n_responses,
    n_total,
    alpha_prior = 1,
    beta_prior = 1,
    num_samples = 1000
) {
  # Posterior parameters (Beta-Binomial conjugacy)
  alpha_post <- alpha_prior + n_responses
  beta_post <- beta_prior + (n_total - n_responses)

  rbeta(num_samples, shape1 = alpha_post, shape2 = beta_post)
}

#' Calculate posterior probability that response rate exceeds threshold
#'
#' P(p > threshold | data) using Beta posterior
#'
#' @param n_responses Number of responders
#' @param n_total Total number of patients
#' @param threshold Response rate threshold
#' @param alpha_prior Beta prior shape1 parameter (default 1)
#' @param beta_prior Beta prior shape2 parameter (default 1)
#'
#' @return Posterior probability P(p > threshold)
#' @export
prob_response_exceeds <- function(
    n_responses,
    n_total,
    threshold,
    alpha_prior = 1,
    beta_prior = 1
) {
  alpha_post <- alpha_prior + n_responses
  beta_post <- beta_prior + (n_total - n_responses)

  # P(p > threshold) = 1 - P(p <= threshold) = 1 - pbeta(threshold)
  1 - pbeta(threshold, shape1 = alpha_post, shape2 = beta_post)
}

#' Calculate posterior probability that response rate is below threshold
#'
#' P(p < threshold | data) using Beta posterior
#'
#' @param n_responses Number of responders
#' @param n_total Total number of patients
#' @param threshold Response rate threshold
#' @param alpha_prior Beta prior shape1 parameter (default 1)
#' @param beta_prior Beta prior shape2 parameter (default 1)
#'
#' @return Posterior probability P(p < threshold)
#' @export
prob_response_below <- function(
    n_responses,
    n_total,
    threshold,
    alpha_prior = 1,
    beta_prior = 1
) {
  alpha_post <- alpha_prior + n_responses
  beta_post <- beta_prior + (n_total - n_responses)

  pbeta(threshold, shape1 = alpha_post, shape2 = beta_post)
}

# =============================================================================
# Binary Interim Logic
# =============================================================================

#' Calculate binary endpoint interim probabilities
#'
#' Computes posterior probabilities for efficacy and futility decisions
#' in binary endpoint trials.
#'
#' @param n_responses Number of responders observed
#' @param n_total Total patients enrolled
#' @param p0 Null hypothesis response rate (for futility)
#' @param p1 Alternative hypothesis response rate (for efficacy target)
#' @param alpha_prior Beta prior shape1 (default 1)
#' @param beta_prior Beta prior shape2 (default 1)
#'
#' @return List with pr_eff (P(p > p0)) and pr_fut (P(p < p1))
#' @keywords internal
calculate_binary_probs <- function(
    n_responses,
    n_total,
    p0,
    p1,
    alpha_prior = 1,
    beta_prior = 1
) {
  list(
    pr_eff = prob_response_exceeds(n_responses, n_total, p0, alpha_prior, beta_prior),
    pr_fut = prob_response_below(n_responses, n_total, p1, alpha_prior, beta_prior)
  )
}

#' Binary interim decision check
#'
#' Evaluates whether to stop for efficacy or futility at an interim look
#' in a binary endpoint trial.
#'
#' @param n_responses Number of responders
#' @param n_total Total enrolled
#' @param args Trial arguments containing thresholds
#' @param diagnostics Print diagnostic messages
#'
#' @return List with decision ("continue", "stop_efficacy", "stop_futility")
#'   and probabilities
#' @keywords internal
binary_interim_decision <- function(
    n_responses,
    n_total,
    args,
    diagnostics = FALSE
) {
  # Get thresholds
  p0 <- args$binary_p0  # null response rate
  p1 <- args$binary_p1 %||% p0  # alternative (futility reference)

  # Posterior probability thresholds for decisions
  eff_prob_threshold <- args$efficacy_threshold_binary_prob %||% 0.95
  fut_prob_threshold <- args$futility_threshold_binary_prob %||% 0.95

  # Prior parameters
  alpha_prior <- args$binary_alpha_prior %||% 1
  beta_prior <- args$binary_beta_prior %||% 1

  # Calculate probabilities
  probs <- calculate_binary_probs(
    n_responses = n_responses,
    n_total = n_total,
    p0 = p0,
    p1 = p1,
    alpha_prior = alpha_prior,
    beta_prior = beta_prior
  )

  if (diagnostics) {
    message(sprintf(
      "[Binary] n=%d, r=%d, P(p>%.2f)=%.3f, P(p<%.2f)=%.3f",
      n_total, n_responses, p0, probs$pr_eff, p1, probs$pr_fut
    ))
  }

  decision <- "continue"

  # Efficacy: P(p > p0) >= threshold
  if (probs$pr_eff >= eff_prob_threshold) {
    decision <- "stop_efficacy"
  }
  # Futility: P(p < p1) >= threshold (only if not already stopping for efficacy)
  else if (probs$pr_fut >= fut_prob_threshold) {
    decision <- "stop_futility"
  }

  list(
    decision = decision,
    pr_eff = probs$pr_eff,
    pr_fut = probs$pr_fut
  )
}

# =============================================================================
# Simon Two-Stage Design Helpers
# =============================================================================

#' Check Simon-style stage 1 futility rule
#'
#' In classical Simon design, stop at stage 1 if responses <= r1
#'
#' @param n_responses Number of responders at stage 1
#' @param r1 Stage 1 futility boundary (stop if responses <= r1)
#'
#' @return TRUE if should stop for futility
#' @keywords internal
simon_stage1_futility <- function(n_responses, r1) {

  n_responses <= r1
}

#' Check Simon-style final efficacy rule
#'
#' In classical Simon design, declare efficacy if total responses > r
#'
#' @param n_responses Total number of responders
#' @param r Total response threshold (success if responses > r)
#'
#' @return TRUE if efficacious
#' @keywords internal
simon_final_efficacy <- function(n_responses, r) {
  n_responses > r
}

#' Calculate Simon design operating characteristics analytically
#'
#' Computes exact operating characteristics for a Simon two-stage design.
#'
#' @param n1 Stage 1 sample size
#' @param r1 Stage 1 futility boundary (stop if X1 <= r1)
#' @param n Total sample size (n1 + n2)
#' @param r Total response threshold for efficacy (success if X > r)
#' @param p True response probability
#'
#' @return List with:
#'   - reject_prob: Probability of rejecting null (power or type I error)
#'   - pet: Probability of early termination at stage 1
#'   - en: Expected sample size
#' @export
simon_oc_exact <- function(n1, r1, n, r, p) {
  n2 <- n - n1


  # Probability of early termination (PET) = P(X1 <= r1)
  pet <- pbinom(r1, size = n1, prob = p)


  # Probability of rejection = P(X1 > r1 AND X1 + X2 > r)
  # = sum over x1 from (r1+1) to n1 of:
  #     P(X1 = x1) * P(X2 > r - x1)
  reject_prob <- 0
  for (x1 in (r1 + 1):n1) {
    p_x1 <- dbinom(x1, size = n1, prob = p)
    # Need X2 > r - x1, i.e., X2 >= r - x1 + 1
    # P(X2 > r - x1) = 1 - P(X2 <= r - x1)
    threshold_x2 <- r - x1
    if (threshold_x2 < 0) {
      # Always reject if x1 alone exceeds r
      p_reject_given_x1 <- 1
    } else if (threshold_x2 >= n2) {
      # Can never get enough responses in stage 2
      p_reject_given_x1 <- 0
    } else {
      p_reject_given_x1 <- 1 - pbinom(threshold_x2, size = n2, prob = p)
    }
    reject_prob <- reject_prob + p_x1 * p_reject_given_x1
  }

  # Expected sample size
  # E[N] = n1 + (1 - PET) * n2 = n1 + P(X1 > r1) * n2
  en <- n1 + (1 - pet) * n2

  list(
    reject_prob = reject_prob,
    pet = pet,
    en = en,
    n1 = n1,
    r1 = r1,
    n = n,
    r = r,
    p = p
  )
}

#' Find optimal Simon design
#'
#' Searches for the Simon optimal or minimax design meeting constraints.
#' Uses optimized search with early termination.
#'
#' @param p0 Null response rate
#' @param p1 Alternative response rate
#' @param alpha Maximum type I error
#' @param beta Maximum type II error (1 - power)
#' @param n_max Maximum total sample size to search
#' @param criterion "optimal" (minimize E[N] under null) or "minimax" (minimize max N)
#'
#' @return Data frame with design parameters and operating characteristics
#' @export
find_simon_design <- function(
    p0,
    p1,
    alpha = 0.10,
    beta = 0.20,
    n_max = 100,
    criterion = c("optimal", "minimax")
) {
  criterion <- match.arg(criterion)

  best_design <- NULL
  best_metric <- Inf

  # For minimax, search from small n upward and stop when we find a feasible design
  # For optimal, need to search more broadly but can use pruning

  for (n in 2:n_max) {
    # For minimax: if we already have a design with smaller n, skip
    if (criterion == "minimax" && !is.null(best_design)) break

    # For optimal: if best E[N] < n1_min, no point searching (E[N] >= n1)
    # Skip if no chance of improvement
    if (criterion == "optimal" && best_metric < 1) break

    for (n1 in 1:(n - 1)) {
      n2 <- n - n1

      # Pruning: if n1 alone > best_metric, can't improve (for optimal)
      if (criterion == "optimal" && n1 >= best_metric) next

      for (r1 in 0:(min(n1 - 1, floor(p0 * n1) + 3))) {
        # Pruning: check if stage 1 PET is reasonable
        pet1_null <- pbinom(r1, n1, p0)

        # If PET under null is too low and we're looking for optimal, E[N] will be high
        # Skip some obviously bad r1 values

        for (r in max(r1, floor(p0 * n) - 2):(min(n - 1, ceiling(p1 * n) + 5))) {
          # Calculate OC under null (fast)
          oc_null <- simon_oc_exact(n1, r1, n, r, p0)

          # Check type I error constraint
          if (oc_null$reject_prob > alpha) next

          # Early termination for optimal: if EN_null >= best_metric, skip
          if (criterion == "optimal" && oc_null$en >= best_metric) next

          # Calculate OC under alternative
          oc_alt <- simon_oc_exact(n1, r1, n, r, p1)

          # Check power constraint
          if (oc_alt$reject_prob < 1 - beta) next

          # This design is feasible!
          metric <- if (criterion == "optimal") oc_null$en else n

          if (metric < best_metric) {
            best_metric <- metric
            best_design <- data.frame(
              n1 = n1,
              r1 = r1,
              n = n,
              r = r,
              n2 = n2,
              EN_null = oc_null$en,
              EN_alt = oc_alt$en,
              PET_null = oc_null$pet,
              PET_alt = oc_alt$pet,
              type1_error = oc_null$reject_prob,
              power = oc_alt$reject_prob,
              criterion = criterion
            )

            # For minimax, first feasible at this n is the best at this n
            if (criterion == "minimax") break
          }
        }
        if (criterion == "minimax" && !is.null(best_design) && best_design$n == n) break
      }
      if (criterion == "minimax" && !is.null(best_design) && best_design$n == n) break
    }
  }

  if (is.null(best_design)) {
    warning("No feasible design found within n_max = ", n_max)
    return(NULL)
  }

  best_design
}

# =============================================================================
# State Management for Binary Endpoints
# =============================================================================

#' Create state container for binary endpoint trial
#'
#' @param arm_names Character vector of arm names
#' @param max_total_patients_per_arm Named integer vector of max N per arm
#'
#' @return State list for binary trial simulation
#' @keywords internal
make_state_binary <- function(arm_names, max_total_patients_per_arm) {
  # Pre-allocate registries for binary outcomes
  capacity <- as.integer(ceiling(max_total_patients_per_arm * 1.2))
  names(capacity) <- names(max_total_patients_per_arm)

  create_binary_registry <- function(cap) {
    data.frame(
      id = rep(NA_integer_, cap),
      stage = rep(NA_integer_, cap),  # 1 or 2
      response = rep(NA_integer_, cap)  # 0 or 1
    )
  }

  list(
    arm_status = setNames(rep("recruiting", length(arm_names)), arm_names),
    enrolled_counts = setNames(rep(0L, length(arm_names)), arm_names),
    stage_counts = setNames(
      lapply(arm_names, function(a) c(stage1 = 0L, stage2 = 0L)),
      arm_names
    ),
    current_stage = setNames(rep(1L, length(arm_names)), arm_names),
    registries = setNames(
      lapply(arm_names, function(a) create_binary_registry(capacity[a])),
      arm_names
    ),
    registry_row_idx = setNames(rep(1L, length(arm_names)), arm_names),
    response_counts = setNames(rep(0L, length(arm_names)), arm_names),
    # Per-simulation outputs
    stop_efficacy_per_sim_row = setNames(rep(0L, length(arm_names)), arm_names),
    stop_futility_per_sim_row = setNames(rep(0L, length(arm_names)), arm_names),
    sim_final_n_current_run = setNames(rep(NA_integer_, length(arm_names)), arm_names)
  )
}

#' Get binary trial metrics from registry
#'
#' @param registry_df Binary registry data frame
#' @param stage Optional: filter to specific stage (1, 2, or NULL for all)
#'
#' @return List with n_enrolled, n_responses
#' @keywords internal
get_binary_metrics <- function(registry_df, stage = NULL) {
  # Filter to active rows
  active <- registry_df[!is.na(registry_df$id), , drop = FALSE]

  if (!is.null(stage)) {
    active <- active[active$stage == stage, , drop = FALSE]
  }

  list(
    n_enrolled = nrow(active),
    n_responses = sum(active$response, na.rm = TRUE)
  )
}

# =============================================================================
# Binary Simulation Driver
# =============================================================================

#' Run binary endpoint trial simulation
#'
#' Simulates a two-stage binary endpoint trial with Bayesian decision rules.
#' Supports Simon-style designs with cohort-based enrollment.
#'
#' @param num_simulations Number of Monte Carlo replicates
#' @param arm_names Character vector of arm names (typically single arm)
#' @param true_response_prob Named numeric vector of true response probabilities
#' @param n1_per_arm Named integer vector of stage 1 sample sizes
#' @param n_total_per_arm Named integer vector of total sample sizes
#' @param p0 Null hypothesis response rate
#' @param p1 Alternative hypothesis response rate
#' @param efficacy_threshold_binary_prob Posterior probability threshold for efficacy
#' @param futility_threshold_binary_prob Posterior probability threshold for futility
#' @param r1_per_arm Optional: Simon-style stage 1 boundaries (stop if X1 <= r1)
#' @param r_per_arm Optional: Simon-style total boundaries (success if X > r)
#' @param use_simon_rules Use exact Simon counting rules instead of Bayesian
#' @param alpha_prior Beta prior shape1 (default 1)
#' @param beta_prior Beta prior shape2 (default 1)
#' @param disable_interim_eff_stop If TRUE, do not stop for efficacy at interim
#'   (only futility stopping at interim, efficacy only at final). This makes

#'   the Bayesian design comparable to Simon's two-stage design which only
#'   stops for futility at the interim.
#' @param diagnostics Print diagnostic messages
#'
#' @return Data frame with operating characteristics per arm
#' @export
run_simulation_binary <- function(
    num_simulations,
    arm_names,
    true_response_prob,
    n1_per_arm,
    n_total_per_arm,
    p0,
    p1 = NULL,
    efficacy_threshold_binary_prob = 0.95,
    futility_threshold_binary_prob = 0.95,
    r1_per_arm = NULL,
    r_per_arm = NULL,
    use_simon_rules = FALSE,
    alpha_prior = 1,
    beta_prior = 1,
    disable_interim_eff_stop = FALSE,
    diagnostics = FALSE,
    progress = interactive()
) {
  if (is.null(p1)) p1 <- p0

  # Ensure named vectors
  if (is.null(names(true_response_prob))) {
    names(true_response_prob) <- arm_names
  }
  if (is.null(names(n1_per_arm))) {
    names(n1_per_arm) <- arm_names
  }
  if (is.null(names(n_total_per_arm))) {
    names(n_total_per_arm) <- arm_names
  }

  # Initialize results
  results_data <- data.frame(
    Arm_Name = arm_names,
    True_Response_Prob = true_response_prob[arm_names],
    Type_I_Error_or_Power = 0,
    PET_Efficacy = 0,
    PET_Futility = 0,
    Pr_Reach_Stage2 = 0,
    Pr_Final_Efficacy = 0,
    Pr_Final_Futility = 0,
    Exp_N = 0,
    Exp_Responses = 0,
    stringsAsFactors = FALSE
  )

  # Progress bar
  show_progress <- isTRUE(progress)
  if (show_progress) {
    pb <- progress::progress_bar$new(
      format = "  Binary sims [:bar] :percent in :elapsed",
      total = num_simulations, clear = FALSE, width = 60
    )
  }

  # Accumulators
  sum_final_n <- setNames(numeric(length(arm_names)), arm_names)
  sum_final_responses <- setNames(numeric(length(arm_names)), arm_names)
  sum_stop_efficacy <- setNames(numeric(length(arm_names)), arm_names)
  sum_stop_futility <- setNames(numeric(length(arm_names)), arm_names)
  sum_final_efficacy <- setNames(numeric(length(arm_names)), arm_names)
  sum_final_futility <- setNames(numeric(length(arm_names)), arm_names)
  sum_reach_stage2 <- setNames(numeric(length(arm_names)), arm_names)

  # Run simulations
  for (sim_idx in seq_len(num_simulations)) {
    if (show_progress) pb$tick()

    for (arm in arm_names) {
      p_true <- true_response_prob[arm]
      n1 <- n1_per_arm[arm]
      n_total <- n_total_per_arm[arm]
      n2 <- n_total - n1

      # Stage 1: Enroll n1 patients
      stage1_data <- simulate_binary_response_data(n1, p_true, start_id = 1)
      x1 <- sum(stage1_data$response)

      # Stage 1 decision
      if (use_simon_rules && !is.null(r1_per_arm)) {
        # Simon counting rule: stop if X1 <= r1
        r1 <- r1_per_arm[arm]
        stop_stage1 <- (x1 <= r1)
        decision1 <- if (stop_stage1) "stop_futility" else "continue"
      } else {
        # Bayesian rule
        args <- list(
          binary_p0 = p0,
          binary_p1 = p1,
          efficacy_threshold_binary_prob = efficacy_threshold_binary_prob,
          futility_threshold_binary_prob = futility_threshold_binary_prob,
          binary_alpha_prior = alpha_prior,
          binary_beta_prior = beta_prior
        )
        result1 <- binary_interim_decision(x1, n1, args, diagnostics = diagnostics)
        decision1 <- result1$decision

        # If disable_interim_eff_stop is TRUE, convert interim efficacy to continue
        # (efficacy will still be evaluated at final analysis)
        if (disable_interim_eff_stop && decision1 == "stop_efficacy") {
          decision1 <- "continue"
        }
      }

      # Process stage 1 outcome
      if (decision1 == "stop_efficacy") {
        sum_stop_efficacy[arm] <- sum_stop_efficacy[arm] + 1
        sum_final_n[arm] <- sum_final_n[arm] + n1
        sum_final_responses[arm] <- sum_final_responses[arm] + x1
        next
      }

      if (decision1 == "stop_futility") {
        sum_stop_futility[arm] <- sum_stop_futility[arm] + 1
        sum_final_n[arm] <- sum_final_n[arm] + n1
        sum_final_responses[arm] <- sum_final_responses[arm] + x1
        next
      }

      # Continue to stage 2
      sum_reach_stage2[arm] <- sum_reach_stage2[arm] + 1

      # Stage 2: Enroll remaining n2 patients
      stage2_data <- simulate_binary_response_data(n2, p_true, start_id = n1 + 1)
      x2 <- sum(stage2_data$response)
      x_total <- x1 + x2

      # Final decision
      if (use_simon_rules && !is.null(r_per_arm)) {
        # Simon counting rule: success if X > r
        r <- r_per_arm[arm]
        final_efficacy <- (x_total > r)
      } else {
        # Bayesian rule at final
        args <- list(
          binary_p0 = p0,
          binary_p1 = p1,
          efficacy_threshold_binary_prob = efficacy_threshold_binary_prob,
          futility_threshold_binary_prob = futility_threshold_binary_prob,
          binary_alpha_prior = alpha_prior,
          binary_beta_prior = beta_prior
        )
        result_final <- binary_interim_decision(x_total, n_total, args, diagnostics = diagnostics)
        final_efficacy <- (result_final$decision == "stop_efficacy" ||
                            result_final$pr_eff >= efficacy_threshold_binary_prob)
      }

      if (final_efficacy) {
        sum_final_efficacy[arm] <- sum_final_efficacy[arm] + 1
      } else {
        sum_final_futility[arm] <- sum_final_futility[arm] + 1
      }

      sum_final_n[arm] <- sum_final_n[arm] + n_total
      sum_final_responses[arm] <- sum_final_responses[arm] + x_total
    }
  }

  # Compute summary statistics
  inv_num_sims <- 1 / num_simulations

  for (j in seq_along(arm_names)) {
    arm <- arm_names[j]
    early_eff <- sum_stop_efficacy[arm] * inv_num_sims
    early_fut <- sum_stop_futility[arm] * inv_num_sims
    final_eff <- sum_final_efficacy[arm] * inv_num_sims
    final_fut <- sum_final_futility[arm] * inv_num_sims

    results_data$Type_I_Error_or_Power[j] <- early_eff + final_eff
    results_data$PET_Efficacy[j] <- early_eff
    results_data$PET_Futility[j] <- early_fut
    results_data$Pr_Reach_Stage2[j] <- sum_reach_stage2[arm] * inv_num_sims
    results_data$Pr_Final_Efficacy[j] <- final_eff
    results_data$Pr_Final_Futility[j] <- final_fut
    results_data$Exp_N[j] <- sum_final_n[arm] * inv_num_sims
    results_data$Exp_Responses[j] <- sum_final_responses[arm] * inv_num_sims
  }

  results_data
}

# =============================================================================
# Utility: Compare BO calibration to Simon enumeration
# =============================================================================

#' Validate BO calibration against Simon enumeration
#'
#' Runs BO calibration on binary simulator and compares to exact Simon design.
#'
#' @param p0 Null response rate
#' @param p1 Alternative response rate
#' @param alpha Type I error constraint
#' @param beta Type II error constraint (1 - power)
#' @param n_max_search Maximum N for Simon search
#' @param bo_fit Optional: pre-computed BO fit object
#' @param num_sims Number of simulations for Monte Carlo validation
#'
#' @return Data frame comparing Simon (exact) vs BO (calibrated) designs
#' @export
compare_simon_to_bo <- function(
    p0,
    p1,
    alpha = 0.10,
    beta = 0.20,
    n_max_search = 100,
    bo_fit = NULL,
    num_sims = 10000
) {
  # Get Simon optimal design
  simon_opt <- find_simon_design(p0, p1, alpha, beta, n_max_search, "optimal")
  simon_mm <- find_simon_design(p0, p1, alpha, beta, n_max_search, "minimax")

  # Validate Simon designs with simulation
  if (!is.null(simon_opt)) {
    # Under null
    oc_null_opt <- run_simulation_binary(
      num_simulations = num_sims,
      arm_names = "Arm",
      true_response_prob = c(Arm = p0),
      n1_per_arm = c(Arm = simon_opt$n1),
      n_total_per_arm = c(Arm = simon_opt$n),
      p0 = p0,
      p1 = p1,
      r1_per_arm = c(Arm = simon_opt$r1),
      r_per_arm = c(Arm = simon_opt$r),
      use_simon_rules = TRUE,
      progress = FALSE
    )
    # Under alternative
    oc_alt_opt <- run_simulation_binary(
      num_simulations = num_sims,
      arm_names = "Arm",
      true_response_prob = c(Arm = p1),
      n1_per_arm = c(Arm = simon_opt$n1),
      n_total_per_arm = c(Arm = simon_opt$n),
      p0 = p0,
      p1 = p1,
      r1_per_arm = c(Arm = simon_opt$r1),
      r_per_arm = c(Arm = simon_opt$r),
      use_simon_rules = TRUE,
      progress = FALSE
    )

    simon_opt$sim_type1 <- oc_null_opt$Type_I_Error_or_Power
    simon_opt$sim_power <- oc_alt_opt$Type_I_Error_or_Power
    simon_opt$sim_EN_null <- oc_null_opt$Exp_N
    simon_opt$sim_EN_alt <- oc_alt_opt$Exp_N
    simon_opt$sim_PET_null <- oc_null_opt$PET_Futility
  }

  list(
    simon_optimal = simon_opt,
    simon_minimax = simon_mm,
    p0 = p0,
    p1 = p1,
    alpha = alpha,
    beta = beta
  )
}
