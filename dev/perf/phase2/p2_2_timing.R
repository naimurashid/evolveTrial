# p2_2_timing.R -- step P2-2 rough timing (not G6): seconds per 1,000 trials at one worker (parallel_replicates = FALSE,
# in-process), run_scenarios on both scenarios of a capture, median of 3 repetitions after one warm-up.
# Usage: Rscript dev/perf/phase2/p2_2_timing.R <lib> <cpp|exact> <R> <designs>
a <- commandArgs(TRUE); lib <- a[1]; eng <- a[2]; R <- as.integer(a[3]); dn <- strsplit(a[4], ",")[[1]]
library(evolveTrial, lib.loc = lib)
options(evolveTrial.looks_stop_at_final = TRUE, evolveTrial.engine = eng)
for (D in dn) {
  d <- readRDS(sprintf("dev/perf/evidence/args/args_%s.rds", D))
  ba <- d$base_args; ba$num_simulations <- R; ba$parallel_replicates <- FALSE
  go <- function() { t0 <- proc.time()[["elapsed"]]; run_scenarios(ba, d$scens, seed = 900001); proc.time()[["elapsed"]] - t0 }
  go()
  s <- replicate(3, go())
  cat(sprintf("%s %s %s R=%d per_1000_trials_s=%.4f (reps %s) load=%s\n", basename(lib), eng, D, R,
              median(s) / (R * length(d$scens) / 1000), paste(round(s, 3), collapse = ","), readLines("/proc/loadavg")))
}
