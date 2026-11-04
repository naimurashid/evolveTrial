# g6_supp_scaling.R <lib> <design> <R per scenario> <workers> -- supplementary to G6: the step 7 protocol (1,000 trials
# per scenario) is dominated by fixed start-up cost for the exact engine, so 16-worker efficiency is also measured at a
# realistic R. Same production path as p2_t7.R, engine "exact", seed 900001. Prints elapsed and CPU seconds.
a <- commandArgs(TRUE); LIB <- normalizePath(a[1]); D <- a[2]; R <- as.integer(a[3]); W <- as.integer(a[4])
suppressPackageStartupMessages(library(evolveTrial, lib.loc = LIB)); stopifnot(startsWith(find.package("evolveTrial"), LIB))
options(evolveTrial.looks_stop_at_final = TRUE, evolveTrial.engine = "exact")
x <- readRDS(sprintf("/home/naimrashid/research/evolveTrial-perf/dev/perf/evidence/args/args_%s.rds", D))
ba <- x$base_args; ba$num_simulations <- R; ba$parallel_replicates <- TRUE; ba$num_workers <- W
pt <- system.time(run_scenarios(ba, x$scens, parallel = FALSE, seed = 900001, return_variance = TRUE))
cat(sprintf("SUPP %s R %d workers %d elapsed %.3f cpu %.3f load %s\n", D, R, W, pt[["elapsed"]],
            sum(pt[c("user.self", "sys.self", "user.child", "sys.child")], na.rm = TRUE), readLines("/proc/loadavg")))
