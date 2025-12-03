#!/usr/bin/env Rscript
# Profile evolveTrial::run_simulation_pure() sequentially (no PSOCK workers)
# This allows Rprof to see the actual function-level bottlenecks

suppressPackageStartupMessages({
  library(evolveTrial)
})

cat("\n")
cat("===============================================================================\n")
cat("  EVOLVETRIAL SEQUENTIAL PROFILING\n")
cat("  Running without parallel workers for accurate Rprof\n")
cat("===============================================================================\n\n")

# Create output directory
dir.create("profiling_results", showWarnings = FALSE, recursive = TRUE)

# Test parameters (typical single-arm design)
arm_names <- c("Treatment")
reference_arm_name <- NULL
compare_arms_option <- FALSE

weibull_shape_true_arms <- c(Treatment = 1.0)
weibull_median_true_arms <- c(Treatment = 8.0)  # Alternative (H1)
null_median_arms <- c(Treatment = 5.0)

# 8-interval piecewise exponential
num_intervals <- 8
max_follow_up <- 24
interval_cutpoints <- seq(0, max_follow_up, length.out = num_intervals + 1)

# Priors
prior_alpha <- rep(0.5, num_intervals)
prior_beta <- rep(0.1, num_intervals)

# Design parameters
max_n <- 60
efficacy_threshold <- 0.95
futility_threshold <- 0.10

cat("Test configuration:\n")
cat(sprintf("  Arms: %s\n", paste(arm_names, collapse = ", ")))
cat(sprintf("  True median: %.1f months\n", weibull_median_true_arms["Treatment"]))
cat(sprintf("  Null median: %.1f months\n", null_median_arms["Treatment"]))
cat(sprintf("  Max N: %d\n", max_n))
cat(sprintf("  Intervals: %d\n", num_intervals))
cat("\n")

#-------------------------------------------------------------------------------
# SECTION 1: Timing single replicate (no profiling)
#-------------------------------------------------------------------------------
cat("===============================================================================\n")
cat("  SECTION 1: BASELINE TIMING (10 reps, no profiling)\n")
cat("===============================================================================\n\n")

t0 <- Sys.time()
result_baseline <- evolveTrial::run_simulation_pure(
  num_simulations = 10,
  arm_names = arm_names,
  reference_arm_name = reference_arm_name,
  compare_arms_option = compare_arms_option,
  weibull_shape_true_arms = weibull_shape_true_arms,
  weibull_median_true_arms = weibull_median_true_arms,
  null_median_arms = null_median_arms,
  interval_cutpoints_sim = interval_cutpoints,
  max_follow_up_sim = max_follow_up,
  censor_max_time_sim = max_follow_up * 2,
  prior_alpha_params_model = prior_alpha,
  prior_beta_params_model = prior_beta,
  num_posterior_draws = 1000,
  cohort_size_per_arm = 1,
  max_total_patients_per_arm = c(Treatment = max_n),
  efficacy_stopping_rule_hc = TRUE,
  efficacy_threshold_hc_prob = efficacy_threshold,
  futility_stopping_rule_hc = TRUE,
  futility_threshold_hc_prob = futility_threshold,
  median_pfs_success_threshold_arms = null_median_arms,
  final_success_posterior_prob_threshold = efficacy_threshold,
  overall_accrual_rate = 3,
  randomization_probs = c(Treatment = 1.0),
  min_follow_up_at_final = 6,
  interim_calendar_beat = 2,
  parallel_replicates = FALSE,  # CRITICAL: Sequential execution
  diagnostics = FALSE,
  progress = FALSE
)
t1 <- Sys.time()
baseline_time <- as.numeric(t1 - t0)
cat(sprintf("  10 reps completed in %.2f seconds (%.1f ms/rep)\n", baseline_time, baseline_time * 100))
cat(sprintf("  Power estimate: %.3f\n", result_baseline$Type_I_Error_or_Power[1]))
cat("\n")

#-------------------------------------------------------------------------------
# SECTION 2: Profile with 50 reps (short interval for high resolution)
#-------------------------------------------------------------------------------
cat("===============================================================================\n")
cat("  SECTION 2: DETAILED PROFILING (50 reps, 5ms interval)\n")
cat("===============================================================================\n\n")

cat("Starting Rprof...\n")
Rprof("profiling_results/sequential_profile.out", interval = 0.005, memory.profiling = FALSE)

t0 <- Sys.time()
result_profiled <- evolveTrial::run_simulation_pure(
  num_simulations = 50,
  arm_names = arm_names,
  reference_arm_name = reference_arm_name,
  compare_arms_option = compare_arms_option,
  weibull_shape_true_arms = weibull_shape_true_arms,
  weibull_median_true_arms = weibull_median_true_arms,
  null_median_arms = null_median_arms,
  interval_cutpoints_sim = interval_cutpoints,
  max_follow_up_sim = max_follow_up,
  censor_max_time_sim = max_follow_up * 2,
  prior_alpha_params_model = prior_alpha,
  prior_beta_params_model = prior_beta,
  num_posterior_draws = 1000,
  cohort_size_per_arm = 1,
  max_total_patients_per_arm = c(Treatment = max_n),
  efficacy_stopping_rule_hc = TRUE,
  efficacy_threshold_hc_prob = efficacy_threshold,
  futility_stopping_rule_hc = TRUE,
  futility_threshold_hc_prob = futility_threshold,
  median_pfs_success_threshold_arms = null_median_arms,
  final_success_posterior_prob_threshold = efficacy_threshold,
  overall_accrual_rate = 3,
  randomization_probs = c(Treatment = 1.0),
  min_follow_up_at_final = 6,
  interim_calendar_beat = 2,
  parallel_replicates = FALSE,
  diagnostics = FALSE,
  progress = FALSE
)
t1 <- Sys.time()

Rprof(NULL)

profiled_time <- as.numeric(t1 - t0)
cat(sprintf("  50 reps completed in %.2f seconds (%.1f ms/rep)\n", profiled_time, profiled_time * 20))

# Parse profile
prof <- summaryRprof("profiling_results/sequential_profile.out")

cat("\n--- TOP 30 FUNCTIONS BY SELF TIME ---\n\n")
self_time <- prof$by.self[, c("self.time", "self.pct", "total.time", "total.pct")]
print(head(self_time, 30))

cat("\n--- TOP 30 FUNCTIONS BY TOTAL TIME ---\n\n")
total_time <- prof$by.total[, c("total.time", "total.pct", "self.time", "self.pct")]
print(head(total_time, 30))

#-------------------------------------------------------------------------------
# SECTION 3: Categorize bottlenecks
#-------------------------------------------------------------------------------
cat("\n===============================================================================\n")
cat("  SECTION 3: BOTTLENECK CATEGORIZATION\n")
cat("===============================================================================\n\n")

all_funcs <- rownames(prof$by.self)

# Find evolveTrial-specific functions
evolve_funcs <- all_funcs[grepl("evolveTrial|run_simulation|interim|posterior|hazard|median|slice|calculate",
                                 all_funcs, ignore.case = TRUE)]

if (length(evolve_funcs) > 0) {
  cat("evolveTrial-related functions:\n\n")
  evolve_prof <- prof$by.self[evolve_funcs, , drop = FALSE]
  evolve_prof <- evolve_prof[order(-evolve_prof$self.time), ]
  print(evolve_prof)
} else {
  cat("No evolveTrial-prefixed functions in profile (may be called via .Call)\n")
}

# Find RNG functions
rng_funcs <- all_funcs[grepl("rbinom|rnorm|runif|rpois|rexp|rgamma|rweibull|sample",
                              all_funcs)]
if (length(rng_funcs) > 0) {
  cat("\n\nRandom number generation functions:\n\n")
  rng_prof <- prof$by.self[rng_funcs, c("self.time", "self.pct"), drop = FALSE]
  rng_prof <- rng_prof[order(-rng_prof$self.time), ]
  print(rng_prof)

  rng_total <- sum(rng_prof$self.time)
  total_time_profiled <- sum(prof$by.self$self.time)
  cat(sprintf("\nRNG overhead: %.2f sec (%.1f%% of total)\n", rng_total, 100 * rng_total / total_time_profiled))
}

# Find data manipulation functions
data_funcs <- all_funcs[grepl("data\\.frame|rbind|cbind|\\$<-|\\[<-|\\[\\[<-|merge|nrow|ncol",
                               all_funcs)]
if (length(data_funcs) > 0) {
  cat("\n\nData manipulation functions:\n\n")
  data_prof <- prof$by.self[data_funcs, c("self.time", "self.pct"), drop = FALSE]
  data_prof <- data_prof[order(-data_prof$self.time), ]
  print(data_prof)

  data_total <- sum(data_prof$self.time)
  total_time_profiled <- sum(prof$by.self$self.time)
  cat(sprintf("\nData manipulation overhead: %.2f sec (%.1f%% of total)\n",
              data_total, 100 * data_total / total_time_profiled))
}

# Find apply functions
apply_funcs <- all_funcs[grepl("^apply$|^sapply$|^lapply$|^mapply$|^vapply$|^FUN$",
                                all_funcs)]
if (length(apply_funcs) > 0) {
  cat("\n\nApply family functions:\n\n")
  apply_prof <- prof$by.self[apply_funcs, c("self.time", "self.pct"), drop = FALSE]
  apply_prof <- apply_prof[order(-apply_prof$self.time), ]
  print(apply_prof)
}

# Find C++ calls
cpp_funcs <- all_funcs[grepl("\\.Call|_cpp", all_funcs)]
if (length(cpp_funcs) > 0) {
  cat("\n\nC++/Rcpp functions:\n\n")
  cpp_prof <- prof$by.self[cpp_funcs, c("self.time", "self.pct"), drop = FALSE]
  cpp_prof <- cpp_prof[order(-cpp_prof$self.time), ]
  print(cpp_prof)
}

#-------------------------------------------------------------------------------
# SUMMARY
#-------------------------------------------------------------------------------
cat("\n===============================================================================\n")
cat("  SUMMARY\n")
cat("===============================================================================\n\n")

total_time_profiled <- sum(prof$by.self$self.time)
cat(sprintf("Total profiled time: %.2f seconds\n", total_time_profiled))
cat(sprintf("Time per rep: %.1f ms\n", total_time_profiled * 1000 / 50))

# Top 10 hotspots
cat("\nTop 10 hotspots:\n")
for (i in 1:min(10, nrow(prof$by.self))) {
  func_name <- rownames(prof$by.self)[i]
  self_time <- prof$by.self$self.time[i]
  self_pct <- prof$by.self$self.pct[i]
  cat(sprintf("  %2d. %-45s %6.2f sec (%5.1f%%)\n",
              i, substr(func_name, 1, 45), self_time, self_pct))
}

cat("\nProfile saved to: profiling_results/sequential_profile.out\n")
cat("Done!\n")
