# P3-1 gate S, one run: a P3-1 frame (p31_frames.R) through run_scenarios() with return_variance = TRUE, either the R body
# (options(evolveTrial.use_cpp_engine = FALSE), FORK workers) or the exact engine (the default selection, checked to be the
# exact engine). Fresh seeds only (P3-1 range 3_120_000 + k, not used elsewhere in dev/perf).
# Usage: Rscript p31_gateS_run.R <lib> <frame> <r|exact> <seed> <n> <workers> <out.rds>
a <- commandArgs(TRUE)
setwd(normalizePath("~/research/evolveTrial-merge"))
library(evolveTrial, lib.loc = normalizePath(a[1]))
source("dev/perf/phase3/p31_frames.R")
frame <- a[2]; eng <- a[3]; seed <- as.integer(a[4]); n <- as.integer(a[5]); W <- as.integer(a[6]); out <- a[7]
spec <- p31_frame_specs()[[frame]]
stopifnot(!is.null(spec))
inp <- p31_frame_inputs(spec)
options(evolveTrial.looks_stop_at_final = NULL, evolveTrial.engine = NULL, evolveTrial.exact_block_size = NULL,
        evolveTrial.use_cpp_engine = !identical(eng, "r"))
# the exact engine must be the one selected: the frame is admitted for every scenario
for (s in seq_along(inp$scens)) {
  env <- p31_chunk_frame(spec, s, 10L)
  if (!isTRUE(env$cfg$supported)) stop(frame, " scenario ", s, " not admitted")
}
ba <- inp$base_args
ba$num_simulations <- n
ba$parallel_replicates <- W > 1L
ba$num_workers <- W
ba$progress <- FALSE
t0 <- Sys.time()
tab <- as.data.frame(run_scenarios(ba, inp$scens, seed = seed, return_variance = TRUE))
wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
saveRDS(list(frame = frame, engine = eng, seed = seed, n = n, workers = W, tab = tab, wall = wall,
             lib = normalizePath(a[1]), use_cpp_engine = getOption("evolveTrial.use_cpp_engine")), out)
cat(sprintf("%s %s seed %d n %d: %.0f s\n", frame, eng, seed, n, wall))
