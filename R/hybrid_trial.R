#' hybrid_trial.R
#' Core hybrid single-arm to between-arm Bayesian adaptive trial simulator
#'
#' Implements a 4-state machine:
#' STATE_SINGLE -> STATE_CONSIDER_CONVERSION -> STATE_BETWEEN -> STATE_STOP
#'
#' @description
#' This module provides the main simulation engine for hybrid trials that:
#' 1. Begin with single-arm monitoring vs historical control

#' 2. Evaluate conversion to between-arm comparison based on predictive probability
#' 3. Continue with between-arm monitoring using seamless data

# ==============================================================================
# CONSTANTS
# ==============================================================================

#' Trial states for hybrid design
HYBRID_STATES <- list(
  SINGLE = "STATE_SINGLE",
  CONSIDER_CONVERSION = "STATE_CONSIDER_CONVERSION",
  BETWEEN = "STATE_BETWEEN",
 STOP = "STATE_STOP"
)

#' Conversion triggers
CONVERSION_TRIGGERS <- c("any_single_success", "all_single_success", "k_of_K")

#' Futility actions
FUTILITY_ACTIONS <- c("stop_trial", "drop_arm", "continue")

# ==============================================================================
# STATE INITIALIZATION
# ==============================================================================

#' Create initial hybrid trial state
#'
#' @param arm_names Character vector of arm names (e.g., c("Control", "Experimental"))
#' @param reference_arm Name of reference arm for between-arm comparison
#' @param theta Hybrid design parameters (see create_hybrid_theta)
#' @param base_args evolveTrial base configuration
#'
#' @return List containing trial state
#' @export
create_hybrid_state <- function(arm_names, reference_arm, theta, base_args) {
  K <- length(arm_names)

  state <- list(
    # Current state
    current_state = HYBRID_STATES$SINGLE,

    # Arm tracking
    arm_names = arm_names,
    reference_arm = reference_arm,
    active_arms = arm_names,
    dropped_arms = character(),

    # Per-arm single-arm status
    sa_efficacy_reached = setNames(rep(FALSE, K), arm_names),
    sa_futility_reached = setNames(rep(FALSE, K), arm_names),
    sa_decision_time = setNames(rep(NA_real_, K), arm_names),
    sa_posterior_prob = setNames(rep(NA_real_, K), arm_names),

    # Conversion tracking
    conversion_evaluated = FALSE,
    conversion_decision = NA_character_,
    conversion_time = NA_real_,
    pp_at_conversion = NA_real_,
    n_add_selected = NA_integer_,
    pp_curve = NULL,

    # Between-arm tracking
    ba_efficacy_reached = FALSE,
    ba_futility_reached = FALSE,
    ba_decision_time = NA_real_,
    ba_posterior_prob = NA_real_,

    # Enrollment tracking
    n_enrolled = setNames(rep(0L, K), arm_names),
    n_enrolled_phase1 = setNames(rep(0L, K), arm_names),
    n_enrolled_phase2 = setNames(rep(0L, K), arm_names),

    # Posterior parameters (gamma: shape=a, rate=b)
    # For PWE: list of vectors per arm
    posterior_a = setNames(vector("list", K), arm_names),
    posterior_b = setNames(vector("list", K), arm_names),

    # Patient registries (accumulated data)
    registries = setNames(vector("list", K), arm_names),

    # Timeline
    current_time = 0,
    interim_count = 0,

    # Final outcome
    trial_outcome = NA_character_,
    stop_reason = NA_character_
  )

  # Initialize posterior parameters with priors
  n_intervals <- base_args$n_intervals %||% 8
  prior_a <- rep(theta$prior_strength %||% 0.5, n_intervals)
  prior_b <- rep(theta$prior_strength %||% 0.5, n_intervals)

  for (arm in arm_names) {
    state$posterior_a[[arm]] <- prior_a
    state$posterior_b[[arm]] <- prior_b
    state$registries[[arm]] <- data.frame(
      patient_id = integer(),
      enrollment_time = numeric(),
      event_time = numeric(),
      observed_time = numeric(),
      event = logical(),
      stringsAsFactors = FALSE
    )
  }

  class(state) <- c("hybrid_trial_state", "list")
  state
}

# ==============================================================================
# HYBRID DESIGN PARAMETERS
# ==============================================================================

#' Create hybrid design parameter structure
#'
#' @param eff_sa SA efficacy threshold (default 0.90)
#' @param fut_sa SA futility threshold (default 0.10)
#' @param hr_threshold_sa Target HR vs historical (default 0.80)
#' @param ev_sa Minimum events for SA interim (default 15)
#' @param nmax_sa Maximum N in SA phase per arm (default 40)
#' @param conversion_trigger Trigger type: "any_single_success", "all_single_success", "k_of_K"
#' @param k_required For k_of_K trigger, number required (default 1)
#' @param pp_go PP threshold to proceed to BA (default 0.70)
#' @param pp_nogo PP threshold to stop (default 0.20)
#' @param ss_method SSR method: "predictive" or "posterior" (default "predictive")
#' @param max_additional_n Maximum additional patients for BA (default 60)
#' @param n_add_candidates Candidate N values for PP curve (default seq(10, 100, 10))
#' @param eff_ba BA efficacy threshold (default 0.975)
#' @param fut_ba BA futility threshold (default 0.05)
#' @param ev_ba Minimum events for BA interim (default 15)
#' @param nmax_ba Maximum N per arm in BA phase (default 80)
#' @param futility_action Action on SA futility (default "drop_arm")
#' @param prior_strength Gamma prior concentration (default 0.5)
#' @param n_outer MC outer samples for PP (default 1000)
#' @param n_inner MC inner samples for PP (default 250)
#'
#' @return Named list of hybrid design parameters
#' @export
create_hybrid_theta <- function(
  # Phase 1 (Single-Arm)
  eff_sa = 0.90,
  fut_sa = 0.10,
  hr_threshold_sa = 0.80,
  ev_sa = 15,
  nmax_sa = 40,

  # Conversion
  conversion_trigger = c("any_single_success", "all_single_success", "k_of_K"),
  k_required = 1,
  pp_go = 0.70,
  pp_nogo = 0.20,
  ss_method = c("predictive", "posterior"),
  max_additional_n = 60,
  n_add_candidates = seq(10, 100, by = 10),

  # Phase 2 (Between-Arm)
  eff_ba = 0.975,
  fut_ba = 0.05,
  ev_ba = 15,
  nmax_ba = 80,

  # Futility handling
  futility_action = c("drop_arm", "stop_trial", "continue"),

  # Structural
  prior_strength = 0.5,
  n_outer = 1000,
  n_inner = 250
) {
  conversion_trigger <- match.arg(conversion_trigger)
  ss_method <- match.arg(ss_method)
  futility_action <- match.arg(futility_action)

  theta <- list(
    # Phase 1
    eff_sa = eff_sa,
    fut_sa = fut_sa,
    hr_threshold_sa = hr_threshold_sa,
    ev_sa = as.integer(ev_sa),
    nmax_sa = as.integer(nmax_sa),

    # Conversion
    conversion_trigger = conversion_trigger,
    k_required = as.integer(k_required),
    pp_go = pp_go,
    pp_nogo = pp_nogo,
    ss_method = ss_method,
    max_additional_n = as.integer(max_additional_n),
    n_add_candidates = as.integer(n_add_candidates),

    # Phase 2
    eff_ba = eff_ba,
    fut_ba = fut_ba,
    ev_ba = as.integer(ev_ba),
    nmax_ba = as.integer(nmax_ba),

    # Futility
    futility_action = futility_action,

    # Structural
    prior_strength = prior_strength,
    n_outer = as.integer(n_outer),
    n_inner = as.integer(n_inner)
  )

  class(theta) <- c("hybrid_theta", "list")
  theta
}

# ==============================================================================
# MAIN STATE MACHINE
# ==============================================================================

#' Update hybrid trial state based on current conditions
#'
#' This is the main state machine driver. It evaluates the current state
#' and applies appropriate transitions.
#'
#' @param state Current trial state
#' @param theta Hybrid design parameters
#' @param base_args evolveTrial base configuration
#' @param scenario_params Scenario parameters (true hazards, etc.)
#'
#' @return Updated trial state
#' @export
update_hybrid_state <- function(state, theta, base_args, scenario_params) {

  switch(state$current_state,

    "STATE_SINGLE" = {
      state <- handle_state_single(state, theta, base_args, scenario_params)
    },

    "STATE_CONSIDER_CONVERSION" = {
      state <- handle_state_consider_conversion(state, theta, base_args, scenario_params)
    },

    "STATE_BETWEEN" = {
      state <- handle_state_between(state, theta, base_args, scenario_params)
    },

    "STATE_STOP" = {
      # Terminal state - no transitions
    }
  )

  state
}

# ==============================================================================
# STATE HANDLERS
# ==============================================================================

#' Handle STATE_SINGLE phase
#'
#' @param state Current trial state
#' @param theta Hybrid design parameters
#' @param base_args evolveTrial base configuration
#' @param scenario_params Scenario parameters
#'
#' @return Updated state
handle_state_single <- function(state, theta, base_args, scenario_params) {

  # Update posteriors for each arm
  state <- update_posteriors(state, base_args)

  # Compute single-arm posterior probabilities vs historical
  for (arm in state$active_arms) {
    if (state$sa_efficacy_reached[arm] || state$sa_futility_reached[arm]) next

    # Get historical hazard for this arm
    hist_hazard <- get_historical_hazard(arm, base_args, scenario_params)

    # Compute P(HR < c | data) using current posterior
    p_single <- compute_p_single_arm(
      state$posterior_a[[arm]],
      state$posterior_b[[arm]],
      hist_hazard,
      theta$hr_threshold_sa,
      base_args
    )

    state$sa_posterior_prob[arm] <- p_single

    # Count events for this arm
    total_events <- count_arm_events(state, arm, state$current_time)

    # Check information gate
    if (total_events < theta$ev_sa) next

    # Check efficacy
    if (p_single >= theta$eff_sa) {
      state$sa_efficacy_reached[arm] <- TRUE
      state$sa_decision_time[arm] <- state$current_time
    }

    # Check futility
    if (p_single <= theta$fut_sa) {
      state$sa_futility_reached[arm] <- TRUE
      state$sa_decision_time[arm] <- state$current_time
      state <- handle_sa_futility(state, arm, theta)
    }
  }

  # Check if all arms dropped
  if (length(state$active_arms) == 0) {
    state$current_state <- HYBRID_STATES$STOP
    state$trial_outcome <- "all_arms_futile"
    state$stop_reason <- "All arms dropped for futility in SA phase"
    return(state)
  }

  # Check transition trigger
  if (check_conversion_trigger(state, theta)) {
    state$current_state <- HYBRID_STATES$CONSIDER_CONVERSION
    return(state)
  }

  # Check max N for SA phase
  total_enrolled <- sum(state$n_enrolled[state$active_arms])
  max_sa_n <- theta$nmax_sa * length(state$arm_names)

  if (total_enrolled >= max_sa_n) {
    state$current_state <- HYBRID_STATES$STOP
    state$trial_outcome <- "max_n_single_phase"
    state$stop_reason <- "Reached max N in SA phase without trigger"
    return(state)
  }

  state
}

#' Handle STATE_CONSIDER_CONVERSION phase
#'
#' @param state Current trial state
#' @param theta Hybrid design parameters
#' @param base_args evolveTrial base configuration
#' @param scenario_params Scenario parameters
#'
#' @return Updated state
handle_state_consider_conversion <- function(state, theta, base_args, scenario_params) {

  if (state$conversion_evaluated) return(state)

  # Compute PP curve for candidate N values
  pp_results <- compute_pp_curve(state, theta, base_args, scenario_params)

  state$pp_curve <- pp_results
  state$conversion_evaluated <- TRUE
  state$conversion_time <- state$current_time

  # Find viable N achieving pp_go
  viable_idx <- which(pp_results$pp >= theta$pp_go)

  if (length(viable_idx) > 0) {
    # GO: Proceed to between-arm phase
    state$n_add_selected <- pp_results$n_add[min(viable_idx)]
    state$pp_at_conversion <- pp_results$pp[min(viable_idx)]
    state$conversion_decision <- "GO"
    state$current_state <- HYBRID_STATES$BETWEEN

    # Record phase 1 enrollment
    state$n_enrolled_phase1 <- state$n_enrolled

  } else if (max(pp_results$pp) < theta$pp_nogo) {
    # NO-GO: Not worth continuing
    state$pp_at_conversion <- max(pp_results$pp)
    state$conversion_decision <- "NO_GO"
    state$current_state <- HYBRID_STATES$STOP
    state$trial_outcome <- "conversion_nogo"
    state$stop_reason <- sprintf("PP (%.2f) below no-go threshold (%.2f)",
                                  max(pp_results$pp), theta$pp_nogo)

  } else {
    # AMBIGUOUS: Default to no-go
    state$pp_at_conversion <- max(pp_results$pp)
    state$conversion_decision <- "AMBIGUOUS_NOGO"
    state$current_state <- HYBRID_STATES$STOP
    state$trial_outcome <- "conversion_ambiguous"
    state$stop_reason <- sprintf("PP (%.2f) in ambiguous region [%.2f, %.2f)",
                                  max(pp_results$pp), theta$pp_nogo, theta$pp_go)
  }

  state
}

#' Handle STATE_BETWEEN phase
#'
#' @param state Current trial state
#' @param theta Hybrid design parameters
#' @param base_args evolveTrial base configuration
#' @param scenario_params Scenario parameters
#'
#' @return Updated state
handle_state_between <- function(state, theta, base_args, scenario_params) {

  # Update posteriors (using all accumulated data - seamless)
  state <- update_posteriors(state, base_args)

  # Compute between-arm posterior probability
  p_between <- compute_p_between_arm(state, theta, base_args)
  state$ba_posterior_prob <- p_between

  # Count events across arms
  total_events <- sum(sapply(state$active_arms, function(arm) {
    count_arm_events(state, arm, state$current_time)
  }))

  # Check information gate
  if (total_events >= theta$ev_ba * 2) {  # Events across both arms

    # Check efficacy
    if (p_between >= theta$eff_ba) {
      state$ba_efficacy_reached <- TRUE
      state$ba_decision_time <- state$current_time
      state$current_state <- HYBRID_STATES$STOP
      state$trial_outcome <- "ba_efficacy"
      state$stop_reason <- sprintf("BA efficacy: P(HR<1|data) = %.3f >= %.3f",
                                    p_between, theta$eff_ba)
      return(state)
    }

    # Check futility
    if (p_between <= theta$fut_ba) {
      state$ba_futility_reached <- TRUE
      state$ba_decision_time <- state$current_time
      state$current_state <- HYBRID_STATES$STOP
      state$trial_outcome <- "ba_futility"
      state$stop_reason <- sprintf("BA futility: P(HR<1|data) = %.3f <= %.3f",
                                    p_between, theta$fut_ba)
      return(state)
    }
  }

  # Check max N for BA phase
  total_enrolled <- sum(state$n_enrolled[state$active_arms])
  max_ba_n <- theta$nmax_ba * length(state$arm_names)

  if (total_enrolled >= max_ba_n) {
    state$current_state <- HYBRID_STATES$STOP
    state$trial_outcome <- "max_n_between_phase"
    state$stop_reason <- "Reached max N in BA phase"

    # Final decision based on posterior
    if (p_between > 0.5) {
      state$trial_outcome <- "ba_efficacy_at_max"
    } else {
      state$trial_outcome <- "ba_no_efficacy_at_max"
    }
  }

  state
}

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Handle single-arm futility for an arm
#'
#' @param state Current trial state
#' @param arm Arm name that hit futility
#' @param theta Hybrid design parameters
#'
#' @return Updated state
handle_sa_futility <- function(state, arm, theta) {

  switch(theta$futility_action,
    "stop_trial" = {
      state$current_state <- HYBRID_STATES$STOP
      state$trial_outcome <- "sa_futility_stop"
      state$stop_reason <- sprintf("Arm %s hit SA futility, trial stopped", arm)
    },
    "drop_arm" = {
      state$active_arms <- setdiff(state$active_arms, arm)
      state$dropped_arms <- c(state$dropped_arms, arm)
    },
    "continue" = {
      # Do nothing, continue monitoring
    }
  )

  state
}

#' Check if conversion trigger is met
#'
#' @param state Current trial state
#' @param theta Hybrid design parameters
#'
#' @return Logical
check_conversion_trigger <- function(state, theta) {

  efficacy_count <- sum(state$sa_efficacy_reached[state$active_arms])
  n_active <- length(state$active_arms)

  switch(theta$conversion_trigger,
    "any_single_success" = efficacy_count >= 1,
    "all_single_success" = efficacy_count == n_active && n_active > 0,
    "k_of_K" = efficacy_count >= theta$k_required,
    FALSE
  )
}

#' Update posterior parameters from trial data
#'
#' @param state Current trial state
#' @param base_args evolveTrial base configuration
#'
#' @return Updated state with posterior parameters
update_posteriors <- function(state, base_args) {

  n_intervals <- base_args$n_intervals %||% 8
  interval_cutpoints <- base_args$interval_cutpoints_sim %||%
    seq(0, 24, length.out = n_intervals + 1)

  for (arm in state$arm_names) {
    registry <- state$registries[[arm]]

    if (nrow(registry) == 0) next

    # Compute interval-specific sufficient statistics
    metrics <- compute_interval_metrics(
      registry,
      state$current_time,
      interval_cutpoints
    )

    # Update posterior: a = a_prior + events, b = b_prior + exposure
    prior_a <- rep(base_args$prior_alpha_params_model[1] %||% 0.5, n_intervals)
    prior_b <- rep(base_args$prior_beta_params_model[1] %||% 0.5, n_intervals)

    state$posterior_a[[arm]] <- prior_a + metrics$events_per_interval
    state$posterior_b[[arm]] <- prior_b + metrics$exposure_per_interval
  }

  state
}

#' Compute interval-specific metrics from registry
#'
#' @param registry Patient registry data frame
#' @param analysis_time Current analysis time
#' @param interval_cutpoints PWE interval boundaries
#'
#' @return List with events_per_interval and exposure_per_interval
compute_interval_metrics <- function(registry, analysis_time, interval_cutpoints) {

  n_intervals <- length(interval_cutpoints) - 1
  events <- numeric(n_intervals)
  exposure <- numeric(n_intervals)

  for (i in seq_len(nrow(registry))) {
    patient <- registry[i, ]

    # Time since enrollment
    entry_time <- 0
    exit_time <- min(
      patient$event_time - patient$enrollment_time,
      analysis_time - patient$enrollment_time
    )
    exit_time <- max(exit_time, 0)

    had_event <- patient$event &&
      (patient$event_time <= analysis_time)

    event_time_since_entry <- if (had_event) {
      patient$event_time - patient$enrollment_time
    } else {
      Inf
    }

    # Allocate to intervals
    for (j in seq_len(n_intervals)) {
      int_start <- interval_cutpoints[j]
      int_end <- interval_cutpoints[j + 1]

      # Exposure in this interval
      exp_start <- max(entry_time, int_start)
      exp_end <- min(exit_time, int_end)

      if (exp_end > exp_start) {
        exposure[j] <- exposure[j] + (exp_end - exp_start)
      }

      # Event in this interval
      if (had_event && event_time_since_entry >= int_start &&
          event_time_since_entry < int_end) {
        events[j] <- events[j] + 1
      }
    }
  }

  list(
    events_per_interval = events,
    exposure_per_interval = exposure
  )
}

#' Get historical hazard for an arm
#'
#' @param arm Arm name
#' @param base_args evolveTrial base configuration
#' @param scenario_params Scenario parameters
#'
#' @return Numeric hazard rate (or vector for PWE)
get_historical_hazard <- function(arm, base_args, scenario_params) {

  # Use historical median to compute hazard
  hist_median <- scenario_params$historical_median %||%
    base_args$historical_median %||% 6

  # For exponential: lambda = log(2) / median
  log(2) / hist_median
}

#' Count events for an arm up to analysis time
#'
#' @param state Current trial state
#' @param arm Arm name
#' @param analysis_time Analysis time
#'
#' @return Integer event count
count_arm_events <- function(state, arm, analysis_time) {

  registry <- state$registries[[arm]]
  if (nrow(registry) == 0) return(0L)

  sum(registry$event & registry$event_time <= analysis_time)
}

#' Compute single-arm posterior probability
#'
#' P(HR < c | data) where HR = lambda_arm / lambda_hist
#'
#' @param post_a Posterior shape parameters
#' @param post_b Posterior rate parameters
#' @param hist_hazard Historical hazard
#' @param hr_threshold Target HR threshold
#' @param base_args evolveTrial base configuration
#'
#' @return Posterior probability
compute_p_single_arm <- function(post_a, post_b, hist_hazard, hr_threshold, base_args) {

  # For exponential (single interval or use median)
  if (length(post_a) == 1) {
    # P(lambda < c * lambda_hist) = Gamma CDF
    threshold <- hr_threshold * hist_hazard
    return(pgamma(threshold, shape = post_a, rate = post_b))
  }

  # For PWE: Monte Carlo sampling
  n_samples <- 1000
  samples <- matrix(0, n_samples, length(post_a))

  for (j in seq_along(post_a)) {
    samples[, j] <- rgamma(n_samples, shape = post_a[j], rate = post_b[j])
  }

  # Compute median survival for each sample
  interval_cutpoints <- base_args$interval_cutpoints_sim %||%
    seq(0, 24, length.out = length(post_a) + 1)

  medians <- apply(samples, 1, function(lambda) {
    compute_pwe_median(lambda, interval_cutpoints)
  })

  # Historical median
  hist_median <- log(2) / hist_hazard

  # HR < c means median_arm > hist_median / c (longer survival)
  # This is equivalent to better performance
  mean(medians > hist_median / hr_threshold)
}

#' Compute between-arm posterior probability
#'
#' P(HR_AB < 1 | data) where A is experimental, B is reference
#'
#' @param state Current trial state
#' @param theta Hybrid design parameters
#' @param base_args evolveTrial base configuration
#'
#' @return Posterior probability
compute_p_between_arm <- function(state, theta, base_args) {

  # Identify experimental and reference arms
  exp_arm <- setdiff(state$active_arms, state$reference_arm)[1]
  ref_arm <- state$reference_arm

  if (is.na(exp_arm) || !exp_arm %in% state$active_arms) {
    return(NA_real_)
  }

  a_exp <- state$posterior_a[[exp_arm]]
  b_exp <- state$posterior_b[[exp_arm]]
  a_ref <- state$posterior_a[[ref_arm]]
  b_ref <- state$posterior_b[[ref_arm]]

  # For exponential (single interval)
  if (length(a_exp) == 1 && length(a_ref) == 1) {
    # Closed-form F-distribution
    # P(lambda_exp / lambda_ref < 1)
    return(pf(
      q = (b_exp / b_ref) * (a_ref / a_exp),
      df1 = 2 * a_exp,
      df2 = 2 * a_ref
    ))
  }

  # For PWE: Monte Carlo sampling
  n_samples <- 2000
  interval_cutpoints <- base_args$interval_cutpoints_sim %||%
    seq(0, 24, length.out = length(a_exp) + 1)

  # Sample hazards for each arm
  samples_exp <- matrix(0, n_samples, length(a_exp))
  samples_ref <- matrix(0, n_samples, length(a_ref))

  for (j in seq_along(a_exp)) {
    samples_exp[, j] <- rgamma(n_samples, shape = a_exp[j], rate = b_exp[j])
    samples_ref[, j] <- rgamma(n_samples, shape = a_ref[j], rate = b_ref[j])
  }

  # Compute median survival for each sample
  medians_exp <- apply(samples_exp, 1, function(lambda) {
    compute_pwe_median(lambda, interval_cutpoints)
  })

  medians_ref <- apply(samples_ref, 1, function(lambda) {
    compute_pwe_median(lambda, interval_cutpoints)
  })

  # P(median_exp > median_ref) = P(HR < 1)
  mean(medians_exp > medians_ref)
}

#' Compute median survival from PWE hazards
#'
#' @param lambda Vector of hazard rates per interval
#' @param interval_cutpoints Interval boundaries
#'
#' @return Median survival time
compute_pwe_median <- function(lambda, interval_cutpoints) {

  # Find time where S(t) = 0.5
  # S(t) = exp(-cumulative hazard)

  cum_haz <- 0
  for (j in seq_along(lambda)) {
    int_start <- interval_cutpoints[j]
    int_end <- interval_cutpoints[j + 1]
    int_width <- int_end - int_start

    # Cumulative hazard at end of this interval
    cum_haz_end <- cum_haz + lambda[j] * int_width

    # Check if median is in this interval
    # S(t) = 0.5 means cum_haz = log(2)
    if (cum_haz_end >= log(2)) {
      # Median is in this interval
      # cum_haz + lambda[j] * (t - int_start) = log(2)
      remaining <- log(2) - cum_haz
      return(int_start + remaining / lambda[j])
    }

    cum_haz <- cum_haz_end
  }

  # Median beyond last interval - extrapolate
  remaining <- log(2) - cum_haz
  last_idx <- length(lambda)
  interval_cutpoints[last_idx + 1] + remaining / lambda[last_idx]
}

# ==============================================================================
# PP CURVE COMPUTATION (placeholder - detailed in hybrid_ssr.R)
# ==============================================================================

#' Compute predictive probability curve
#'
#' @param state Current trial state
#' @param theta Hybrid design parameters
#' @param base_args evolveTrial base configuration
#' @param scenario_params Scenario parameters
#'
#' @return Data frame with n_add and pp columns
compute_pp_curve <- function(state, theta, base_args, scenario_params) {

  n_candidates <- theta$n_add_candidates

  if (theta$ss_method == "posterior") {
    # Fast posterior method
    pp_values <- sapply(n_candidates, function(n) {
      compute_pp_posterior(state, n, theta, base_args)
    })
  } else {
    # Predictive probability method (Monte Carlo)
    pp_values <- sapply(n_candidates, function(n) {
      compute_pp_predictive(state, n, theta, base_args, scenario_params)
    })
  }

  data.frame(
    n_add = n_candidates,
    pp = pp_values
  )
}

#' Compute PP using posterior method (fast approximation)
#'
#' @param state Current trial state
#' @param n_add Additional patients per arm
#' @param theta Hybrid design parameters
#' @param base_args evolveTrial base configuration
#'
#' @return Predictive probability
compute_pp_posterior <- function(state, n_add, theta, base_args) {

  # Use current posterior point estimate to calculate conditional power
  # This is a simplified approximation

  p_current <- compute_p_between_arm(state, theta, base_args)

  if (is.na(p_current)) return(0)

  # Rough approximation: PP increases with sample size
  # Using logistic scaling
  current_n <- sum(state$n_enrolled[state$active_arms])
  future_n <- current_n + n_add * 2  # Both arms

  # Information ratio
  info_ratio <- sqrt(future_n / max(current_n, 1))

  # Transform current probability
  # If p_current > 0.5, PP should increase toward 1
  # If p_current < 0.5, PP should stay low

  if (p_current > 0.5) {
    # Positive effect - PP increases with sample size
    effect_size <- qnorm(p_current)
    future_z <- effect_size * info_ratio
    pp <- pnorm(future_z - qnorm(theta$eff_ba))
  } else {
    # Negative or null effect - PP stays low
    pp <- p_current^info_ratio
  }

  min(max(pp, 0), 1)
}

#' Compute PP using predictive method (Monte Carlo)
#'
#' @param state Current trial state
#' @param n_add Additional patients per arm
#' @param theta Hybrid design parameters
#' @param base_args evolveTrial base configuration
#' @param scenario_params Scenario parameters
#'
#' @return Predictive probability
compute_pp_predictive <- function(state, n_add, theta, base_args, scenario_params) {

  # This will be fully implemented in hybrid_ssr.R
  # For now, use a simplified Monte Carlo approach

  n_outer <- min(theta$n_outer, 500)  # Reduce for speed during development
  success_count <- 0

  exp_arm <- setdiff(state$active_arms, state$reference_arm)[1]
  ref_arm <- state$reference_arm

  if (is.na(exp_arm)) return(0)

  a_exp <- state$posterior_a[[exp_arm]]
  b_exp <- state$posterior_b[[exp_arm]]
  a_ref <- state$posterior_a[[ref_arm]]
  b_ref <- state$posterior_b[[ref_arm]]

  interval_cutpoints <- base_args$interval_cutpoints_sim %||%
    seq(0, 24, length.out = length(a_exp) + 1)
  n_intervals <- length(a_exp)

  accrual_rate <- base_args$overall_accrual_rate %||% 2.0
  followup <- base_args$max_follow_up_sim %||% 24

  for (i in seq_len(n_outer)) {
    # Step 1: Draw "true" hazards from current posterior
    lambda_exp_true <- rgamma(n_intervals, shape = a_exp, rate = b_exp)
    lambda_ref_true <- rgamma(n_intervals, shape = a_ref, rate = b_ref)

    # Step 2: Simulate future events/exposure
    future_exp <- simulate_future_arm_pwe(
      n_add, lambda_exp_true, interval_cutpoints, accrual_rate, followup
    )
    future_ref <- simulate_future_arm_pwe(
      n_add, lambda_ref_true, interval_cutpoints, accrual_rate, followup
    )

    # Step 3: Update posteriors
    a_exp_final <- a_exp + future_exp$events
    b_exp_final <- b_exp + future_exp$exposure
    a_ref_final <- a_ref + future_ref$events
    b_ref_final <- b_ref + future_ref$exposure

    # Step 4: Compute P(HR < 1 | final data)
    # For simplicity, use exponential approximation (sum over intervals)
    a_exp_sum <- sum(a_exp_final)
    b_exp_sum <- sum(b_exp_final)
    a_ref_sum <- sum(a_ref_final)
    b_ref_sum <- sum(b_ref_final)

    p_between <- pf(
      q = (b_exp_sum / b_ref_sum) * (a_ref_sum / a_exp_sum),
      df1 = 2 * a_exp_sum,
      df2 = 2 * a_ref_sum
    )

    # Step 5: Check success criterion
    if (p_between >= theta$eff_ba) {
      success_count <- success_count + 1
    }
  }

  success_count / n_outer
}

#' Simulate future arm data under PWE model
#'
#' @param n_patients Number of patients to simulate
#' @param lambda True hazard rates per interval
#' @param interval_cutpoints Interval boundaries
#' @param accrual_rate Patients per month
#' @param followup Follow-up time in months
#'
#' @return List with events and exposure per interval
simulate_future_arm_pwe <- function(n_patients, lambda, interval_cutpoints,
                                     accrual_rate, followup) {

  n_intervals <- length(lambda)
  events <- numeric(n_intervals)
  exposure <- numeric(n_intervals)

  # Enrollment times
  enrollment_duration <- n_patients / accrual_rate
  enrollment_times <- sort(runif(n_patients, 0, enrollment_duration))

  # Analysis time
  analysis_time <- enrollment_duration + followup

  for (p in seq_len(n_patients)) {
    enroll_time <- enrollment_times[p]

    # Simulate survival time under PWE
    surv_time <- simulate_pwe_survival(lambda, interval_cutpoints)

    # Observed time
    observed_time <- min(surv_time, analysis_time - enroll_time)
    had_event <- surv_time <= (analysis_time - enroll_time)

    # Allocate to intervals
    for (j in seq_len(n_intervals)) {
      int_start <- interval_cutpoints[j]
      int_end <- interval_cutpoints[j + 1]

      # Exposure
      exp_start <- max(0, int_start)
      exp_end <- min(observed_time, int_end)

      if (exp_end > exp_start) {
        exposure[j] <- exposure[j] + (exp_end - exp_start)
      }

      # Event
      if (had_event && surv_time >= int_start && surv_time < int_end) {
        events[j] <- events[j] + 1
      }
    }
  }

  list(events = events, exposure = exposure)
}

#' Simulate survival time from PWE model
#'
#' @param lambda Hazard rates per interval
#' @param interval_cutpoints Interval boundaries
#'
#' @return Survival time
simulate_pwe_survival <- function(lambda, interval_cutpoints) {

  u <- runif(1)
  cum_haz <- 0

  for (j in seq_along(lambda)) {
    int_start <- interval_cutpoints[j]
    int_end <- interval_cutpoints[j + 1]
    int_width <- int_end - int_start

    # Survival probability to end of interval
    # S(t) = exp(-cum_haz)
    cum_haz_end <- cum_haz + lambda[j] * int_width

    # Check if event occurs in this interval
    # u < S(t) at start but u >= S(t) at end means event in interval
    s_start <- exp(-cum_haz)
    s_end <- exp(-cum_haz_end)

    if (u >= s_end && u < s_start) {
      # Event in this interval
      # Solve: u = exp(-cum_haz - lambda[j] * (t - int_start))
      # t = int_start - (log(u) + cum_haz) / lambda[j]
      return(int_start - (log(u) + cum_haz) / lambda[j])
    }

    cum_haz <- cum_haz_end
  }

  # Event beyond last interval - extrapolate
  last_idx <- length(lambda)
  int_end <- interval_cutpoints[last_idx + 1]
  return(int_end - (log(u) + cum_haz) / lambda[last_idx])
}

# ==============================================================================
# TRIAL COMPILATION
# ==============================================================================

#' Compile hybrid trial results
#'
#' @param state Final trial state
#' @param theta Hybrid design parameters
#'
#' @return List of trial metrics
#' @export
compile_hybrid_results <- function(state, theta) {

  # Determine final decision
  final_decision <- state$trial_outcome

  # Total N
  total_n <- sum(state$n_enrolled)

  # Phase-specific N
  n_phase1 <- sum(state$n_enrolled_phase1)
  n_phase2 <- sum(state$n_enrolled) - n_phase1

  # Was conversion reached?
  converted <- state$current_state == HYBRID_STATES$BETWEEN ||
    state$conversion_decision == "GO"

  # SA efficacy/futility rates
  sa_efficacy_any <- any(state$sa_efficacy_reached)
  sa_futility_any <- any(state$sa_futility_reached)

  # BA efficacy/futility
  ba_efficacy <- state$ba_efficacy_reached
  ba_futility <- state$ba_futility_reached

  list(
    # Final outcome
    decision = final_decision,
    stop_reason = state$stop_reason,

    # Sample sizes
    total_n = total_n,
    n_phase1 = n_phase1,
    n_phase2 = n_phase2,
    n_per_arm = state$n_enrolled,

    # Conversion
    converted = converted,
    conversion_decision = state$conversion_decision,
    pp_at_conversion = state$pp_at_conversion,
    n_add_selected = state$n_add_selected,

    # SA phase
    sa_efficacy = state$sa_efficacy_reached,
    sa_futility = state$sa_futility_reached,
    sa_efficacy_any = sa_efficacy_any,
    sa_futility_any = sa_futility_any,
    sa_posterior_prob = state$sa_posterior_prob,

    # BA phase
    ba_efficacy = ba_efficacy,
    ba_futility = ba_futility,
    ba_posterior_prob = state$ba_posterior_prob,

    # Timeline
    duration = state$current_time,
    conversion_time = state$conversion_time,
    interim_count = state$interim_count
  )
}

# ==============================================================================
# NULL COALESCING OPERATOR
# ==============================================================================

if (!exists("%||%", mode = "function")) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
}
