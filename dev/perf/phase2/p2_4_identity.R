# p2_4_identity.R -- step P2-4 identity run: the phase 2 chunk engine of one library on all 10 acceptance designs.
# For each design and scenario: the raw per-trial chunk result of .simulate_chunk_exact(1:R, cfg, c(seed, scenario))
# (frame built as in tests/testthat/test-simulate-chunk-exact.R), plus run_scenarios() at one worker
# (parallel_replicates = FALSE) with the capture's seed, and the .Random.seed it leaves. Saved to <out.rds>; compare two
# libraries with p2_4_compare.R.
# Usage: Rscript dev/perf/phase2/p2_4_identity.R <lib> <R> <out.rds>
a <- commandArgs(TRUE); lib <- normalizePath(a[1]); R <- as.integer(a[2]); out <- a[3]
library(evolveTrial, lib.loc = lib)
stopifnot(normalizePath(dirname(getNamespaceInfo("evolveTrial", "path"))) == lib)
options(evolveTrial.looks_stop_at_final = TRUE, evolveTrial.engine = "exact")
ns <- function(name) getFromNamespace(name, "evolveTrial")
build_frame_x <- function(d, scenario, R) {
  rsp <- ns("run_simulation_pure")
  args_i <- utils::modifyList(d$base_args, d$scens[[scenario]], keep.null = TRUE)
  args_i$return_percentiles <- FALSE
  args_i$percentile_probs <- c(0, 0.25, 0.5, 0.75, 0.9, 1.0)
  args_i$return_variance <- TRUE
  args_i <- args_i[intersect(names(args_i), names(formals(rsp)))]
  args_i$num_simulations <- as.integer(R)
  old <- options(evolveTrial.chunk_probe = function(env) env)
  on.exit(options(old), add = TRUE)
  env <- do.call(rsp, args_i)
  stopifnot(is.environment(env), isTRUE(env$looks_stop_at_final))
  env
}
dn <- sub("^args_(R[0-9]+)\\.rds$", "\\1", list.files("dev/perf/evidence/args", pattern = "^args_R[0-9]+\\.rds$"))
res <- list()
for (D in dn) {
  d <- readRDS(sprintf("dev/perf/evidence/args/args_%s.rds", D))
  seed <- as.integer(d$cert_seed %||% d$dots$seed %||% 900001L)
  chunks <- lapply(seq_along(d$scens), function(s) {
    env <- build_frame_x(d, s, R)
    ns(".simulate_chunk_exact")(seq_len(R), env$cfg, c(seed, as.integer(s)))
  })
  ba <- d$base_args; ba$num_simulations <- R; ba$parallel_replicates <- FALSE
  set.seed(seed)
  tab <- run_scenarios(ba, d$scens, seed = seed)
  res[[D]] <- list(seed = seed, chunks = chunks, tab = tab, rng = get(".Random.seed", envir = globalenv()))
  cat(D, "done\n")
}
saveRDS(res, out)
