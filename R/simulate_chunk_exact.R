# simulate_chunk_exact.R -- R side of the phase 2 chunk engine (engine speed-up phase 2, step P2-2; plan
# dev/plans/2026-11-04-phase2-plan.md Section 1 items 2 and 3, report dev/perf/phase2/P2-2_REPORT.md).
#
# Selection: options(evolveTrial.engine = "exact") (phase 2, src/simulate_chunk_exact.cpp; the default since 0.0.0.9002,
# PI ruling 2026-11-04) or "cpp" (phase 1, bit-identical to the R body). evolveTrial.use_cpp_engine = FALSE still selects
# the R body whatever the engine option.
# The option unset (NULL) means the default. Under the default, a frame outside the phase 1 supported set falls back
# silently to the phase 1 dispatcher (C++ body or R body), exactly as "cpp" does, and the L2 notice stays silent. Only an
# explicit options(evolveTrial.engine = "exact") makes an unsupported frame stop and makes the L2 notice fire.
# In "exact" mode (explicit, or the default on a supported frame) run_simulation_pure():
#   1. when the option is explicit, stops unless the frame is in the phase 1 supported set (.chunk_cpp_assert_supported;
#      never falls back);
#   2. draws the scenario seed with one call sample.int(.Machine$integer.max, 2L, replace = TRUE), so set.seed() before
#      the call reproduces the result; nothing else is drawn from R's generator;
#   3. runs trials 1..num_simulations through .simulate_chunk_exact in blocks of getOption("evolveTrial.exact_block_size")
#      trials (default .exact_default_block_size()), dispatched dynamically over the cluster (parallel::clusterApplyLB),
#      or in one call when the run is sequential;
#   4. combines the blocks with .exact_reduce (trial-index order), which returns the phase 1 chunk list, so
#      aggregate_results(), run_scenarios() and BATON see the same names and meanings.
# Each trial's result depends only on the scenario seed and its index, so the returned list is identical() for any
# worker count and block size.

# TRUE when the phase 2 engine is selected. Reads options only (draws nothing).
.exact_engine_selected <- function() {
  eng <- getOption("evolveTrial.engine", "exact")
  if (!(is.character(eng) && length(eng) == 1L && !is.na(eng) && eng %in% c("cpp", "exact"))) {
    stop("options(evolveTrial.engine) must be \"cpp\" or \"exact\"", call. = FALSE)
  }
  identical(eng, "exact") && isTRUE(getOption("evolveTrial.use_cpp_engine", TRUE))
}

# TRUE when the caller set options(evolveTrial.engine) to any value; FALSE when it is unset (the default applies).
.exact_engine_explicit <- function() !is.null(getOption("evolveTrial.engine"))

# One-time notice (review L2) for a caller who set options(evolveTrial.engine = "exact") on a path that does not use the
# exact engine: run_simulation_pure with use_cpp_engine = FALSE, and the other simulators (hybrid, binary). Reads options
# and prints one message per R session; it draws nothing, so R's random stream is untouched. Silent when the option is
# unset (the default "exact" is not a request), in "cpp" mode, and whenever the exact engine is in use.
.exact_note_env <- new.env(parent = emptyenv())
.exact_requested <- function() {
  eng <- getOption("evolveTrial.engine")   # no default on purpose: unset is not an explicit request
  is.character(eng) && length(eng) == 1L && !is.na(eng) && identical(eng, "exact")
}
.exact_note_unused <- function(reason) {
  if (!.exact_requested() || isTRUE(.exact_note_env$noted)) return(invisible(FALSE))
  .exact_note_env$noted <- TRUE
  message("options(evolveTrial.engine = \"exact\") is set, but ", reason, ", so the phase 1 engine or R body runs here. ",
          "This notice appears once per session.")
  invisible(TRUE)
}

# The one draw from R's generator per scenario dispatch: two integers in 1..(2^31 - 1), 62 bits of seed.
.exact_draw_scenario_seed <- function() sample.int(.Machine$integer.max, 2L, replace = TRUE)

# Default block size: about eight blocks per worker, between 1 and 250 trials.
.exact_default_block_size <- function(n, workers) {
  as.integer(max(1L, min(250L, ceiling(n / (8 * max(1L, workers))))))
}

.exact_block_size <- function(n, workers) {
  bs <- getOption("evolveTrial.exact_block_size", NULL)
  if (is.null(bs)) return(.exact_default_block_size(n, workers))
  if (!(is.numeric(bs) && length(bs) == 1L && is.finite(bs) && bs >= 1 && bs == floor(bs))) {
    stop("options(evolveTrial.exact_block_size) must be one positive integer", call. = FALSE)
  }
  as.integer(bs)
}

# Contiguous blocks of trial indices, in index order.
.exact_blocks <- function(n, bs) {
  starts <- seq.int(1L, n, by = bs)
  lapply(starts, function(s) seq.int(s, min(n, s + bs - 1L)))
}

# The function the workers run. A namespace function, so the cluster ships a reference to the namespace rather than a
# closure environment, and .simulate_chunk_exact is looked up in the worker's namespace at call time.
.exact_block_worker <- function(ix, cfg, scenario_seed) .simulate_chunk_exact(ix, cfg, scenario_seed)

# Sequential run (one block) -> phase 1 chunk list.
.exact_run_sequential <- function(n, cfg, scenario_seed) {
  .exact_reduce(list(.exact_block_worker(seq_len(n), cfg, scenario_seed)), cfg, as.integer(n))
}

# Dynamic dispatch over a cluster -> phase 1 chunk list.
.exact_run_cluster <- function(cl, n, cfg, scenario_seed, workers) {
  blocks <- .exact_blocks(as.integer(n), .exact_block_size(n, workers))
  parts <- tryCatch(
    parallel::clusterApplyLB(cl, blocks, .exact_block_worker, cfg = cfg, scenario_seed = scenario_seed),
    error = function(e) stop("Parallel simulation failed: ", conditionMessage(e), call. = FALSE)
  )
  .exact_reduce(parts, cfg, as.integer(n))
}
