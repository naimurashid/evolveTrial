# P3-2 gate S follow-up: does the R body's type I error on D01/D03 (null) move toward the exact engine as the posterior
# draw count grows? Phase 1 engine (bit-identical to the R body, gate B), null scenario only, fresh seeds 3240xxx.
# Usage: Rscript p32_draws_check.R <lib> <frame> <draw multiplier> <seed> <n> <workers> <out.rds>
a <- commandArgs(TRUE)
setwd(normalizePath("~/research/evolveTrial-merge"))
suppressMessages(library(evolveTrial, lib.loc = normalizePath(a[1])))
source("dev/perf/phase3/p32_frames.R")
frame <- a[2]; mult <- as.numeric(a[3]); seed <- as.integer(a[4]); n <- as.integer(a[5]); W <- as.integer(a[6])
inp <- p31_frame_inputs(p32_frame_specs()[[frame]])
options(evolveTrial.engine = "cpp", evolveTrial.use_cpp_engine = TRUE)
ba <- inp$base_args
for (k in c("num_posterior_draws", "num_posterior_draws_interim")) if (!is.null(ba[[k]])) ba[[k]] <- round(ba[[k]] * mult)
ba$num_simulations <- n; ba$parallel_replicates <- W > 1L; ba$num_workers <- W; ba$progress <- FALSE
t0 <- Sys.time()
tab <- as.data.frame(run_scenarios(ba, inp$scens[1], seed = seed, return_variance = TRUE))
wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
saveRDS(list(frame = frame, mult = mult, draws = ba$num_posterior_draws, draws_int = ba$num_posterior_draws_interim,
             seed = seed, n = n, tab = tab, wall = wall), a[7])
cat(sprintf("%s x%g draws %s/%s n %d: type I %.4f EN %.3f, %.0f s\n", frame, mult, ba$num_posterior_draws,
            format(ba$num_posterior_draws_interim), n, tab$Type_I_Error_or_Power[2], tab$Exp_N[2], wall))
