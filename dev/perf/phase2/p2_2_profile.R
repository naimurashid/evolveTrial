# p2_2_profile.R -- step P2-2 profile of the phase 2 engine at one worker (decides conditional step P2-4).
# Needs a scratch build of the package with timers (scope timers around slice_arm, the exact single-arm probability, the
# PH Laplace step, the per-arrival rebalancing scan and the whole trial; built from the worktree by the instructions in
# P2-2_REPORT.md Section 7, exported .exact_profile). Usage: Rscript dev/perf/phase2/p2_2_profile.R <scratch lib> [R]
a <- commandArgs(TRUE); lib <- a[1]; R <- as.integer(if (length(a) >= 2) a[2] else 2000)
library(evolveTrial, lib.loc = lib)
ns <- asNamespace("evolveTrial")
options(evolveTrial.looks_stop_at_final = TRUE)
frame <- function(d, s) {
  rsp <- get("run_simulation_pure", ns)
  args_i <- utils::modifyList(d$base_args, d$scens[[s]], keep.null = TRUE)
  args_i$return_percentiles <- FALSE; args_i$percentile_probs <- c(0, 0.25, 0.5, 0.75, 0.9, 1.0); args_i$return_variance <- TRUE
  args_i <- args_i[intersect(names(args_i), names(formals(rsp)))]; args_i$num_simulations <- R
  old <- options(evolveTrial.chunk_probe = function(env) env); on.exit(options(old)); do.call(rsp, args_i)
}
lab <- c("look slicing (slice_arm)", "exact single-arm probability", "PH Laplace (ph_beta_mode_var_cpp)", "rebalancing scan", "trial total")
for (D in c("R01", "R23")) {
  d <- readRDS(sprintf("dev/perf/evidence/args/args_%s.rds", D))
  ns$.exact_profile(TRUE)
  wall <- 0
  for (s in seq_along(d$scens)) {
    env <- frame(d, s)
    t0 <- proc.time()[["elapsed"]]; invisible(ns$.simulate_chunk_exact(seq_len(R), env$cfg, c(1L, 2L))); wall <- wall + proc.time()[["elapsed"]] - t0
  }
  p <- ns$.exact_profile(FALSE)
  tot <- p$seconds[5]
  cat(sprintf("%s: %d trials (both scenarios), wall %.3f s, %.3f s per 1,000 trials\n", D, R * length(d$scens), wall, wall / (R * length(d$scens) / 1000)))
  for (i in 1:5) cat(sprintf("  %-36s %8.3f s  %5.1f%% of trial time  %9d calls\n", lab[i], p$seconds[i], 100 * p$seconds[i] / tot, p$calls[i]))
}
