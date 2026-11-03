# profile.R -- single-worker (in-process, parallel_replicates = FALSE) timing and Rprof of the captured call.
#   Rscript profile.R <design> <R> [noprof]
a <- commandArgs(TRUE); DES <- a[1]; R <- as.integer(a[2]); PROF <- !identical(a[3], "noprof")
S <- "/tmp/claude-1000/-home-naimrashid-Downloads-adaptive-trial-bo-paper/537dc1ca-2d2f-4a33-9259-7d770244bb0a/scratchpad/perf"
suppressPackageStartupMessages(library(evolveTrial))
stopifnot(startsWith(find.package("evolveTrial"), file.path(S, "Rlib")))
options(evolveTrial.looks_stop_at_final = TRUE)
x <- readRDS(file.path(S, sprintf("args_%s.rds", DES)))
ba <- x$base_args; ba$num_simulations <- R; ba$parallel_replicates <- FALSE; ba$num_workers <- 1L
stopifnot(identical(ba$min_follow_up_at_final, 6), identical(ba$random_dropout, FALSE))
pf <- file.path(S, sprintf("rprof_%s_R%d.out", DES, R))
if (PROF) Rprof(pf, interval = 0.01, line.profiling = TRUE)
pt0 <- proc.time()
res <- run_scenarios(ba, x$scens, parallel = FALSE, seed = 900001, return_variance = TRUE)
pt <- proc.time() - pt0
if (PROF) Rprof(NULL)
res <- as.data.frame(res)
n_trials <- R * length(x$scens)
cat(sprintf("%s R=%d per scenario x %d scenarios = %d trials: elapsed %.2f s, user %.2f, sys %.2f; %.2f s per 1000 trials\n",
            DES, R, length(x$scens), n_trials, pt[["elapsed"]], pt[["user.self"]], pt[["sys.self"]], 1000 * pt[["elapsed"]] / n_trials))
print(res[, intersect(c("scenario", "Arm_Name", "Type_I_Error_or_Power", "PET_Futility", "Expected_N", "Expected_Sample_Size"), names(res))])
saveRDS(list(pt = pt, res = res, n_trials = n_trials), file.path(S, sprintf("timing_seq_%s_R%d.rds", DES, R)))
if (PROF) {
  sf <- summaryRprof(pf)
  sl <- summaryRprof(pf, lines = "show")
  saveRDS(list(fun = sf, line = sl), file.path(S, sprintf("summaryRprof_%s_R%d.rds", DES, R)))
  sink(file.path(S, sprintf("summaryRprof_%s_R%d.txt", DES, R)))
  cat("sampling.time", sf$sampling.time, "\n\n== by.self (function, top 40)\n"); print(head(sf$by.self, 40))
  cat("\n== by.total (function, top 40)\n"); print(head(sf$by.total, 40))
  cat("\n== by.self (line, top 40)\n"); print(head(sl$by.self, 40))
  cat("\n== by.line (top 60 by self)\n"); bl <- sl$by.line; print(head(bl[order(-bl$self.time), ], 60))
  sink()
}
