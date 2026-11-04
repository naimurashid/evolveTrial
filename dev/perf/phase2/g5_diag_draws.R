# g5_diag_draws.R <lib> <engine cpp|exact> <draws> <R> <workers> <seed> <out.rds>
# G5 failure diagnostic (R01 EN_null): does the phase 1 engine's null expected sample size move toward the exact engine's
# value as the number of posterior draws grows? Same production path as p2_t7.R; only the draw counts and R change.
a <- commandArgs(TRUE); LIB <- normalizePath(a[1]); ENG <- a[2]; D <- as.integer(a[3]); R <- as.integer(a[4])
W <- as.integer(a[5]); SEED <- as.integer(a[6])
suppressPackageStartupMessages(library(evolveTrial, lib.loc = LIB))
stopifnot(startsWith(find.package("evolveTrial"), LIB))
options(evolveTrial.looks_stop_at_final = TRUE, evolveTrial.engine = ENG)
x <- readRDS("/home/naimrashid/research/evolveTrial-perf/dev/perf/evidence/args/args_R01.rds")
ba <- x$base_args; ba$num_simulations <- R; ba$parallel_replicates <- TRUE; ba$num_workers <- W
ba$num_posterior_draws <- D; ba$num_posterior_draws_interim <- D
pt <- system.time(res <- run_scenarios(ba, x$scens, parallel = FALSE, seed = SEED, return_variance = TRUE))
out <- list(engine = ENG, draws = D, R = R, seed = SEED, elapsed = pt[["elapsed"]], res = res)
print(as.data.frame(res)[, 1:min(14, ncol(as.data.frame(res)))]); saveRDS(out, a[7])
