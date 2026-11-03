# Harness self-tests H1-H7 (plan 2b.5). Run against dev/perf/lib/dev:
#   Rscript -e 'testthat::test_file("dev/perf/harness/tests/test-harness.R")'
# Environment: EVT_LIB overrides the library (default <root>/dev/perf/lib/dev).
library(testthat)
ROOT <- normalizePath(testthat::test_path("..", "..", "..", ".."))
Sys.setenv(EVT_ROOT = ROOT)
LIB  <- Sys.getenv("EVT_LIB", file.path(ROOT, "dev/perf/lib/dev"))
source(file.path(ROOT, "dev/perf/harness/chunk_harness.R"))
EV <- file.path(ROOT, "dev/perf/evidence")
options(evolveTrial.looks_stop_at_final = TRUE)   # I7: for the whole session, before any frame is built
d01 <- load_design(file.path(EV, "args", "args_R01.rds"))
d23 <- load_design(file.path(EV, "args", "args_R23.rds"))

# A candidate that runs the R body one trial at a time and recombines (the composition of plan H2).
# `hook(i, res)` may modify the k-th trial's result; `pre(i)` runs just before trial i.
composed_cand <- function(env, hook = function(i, res) res, pre = function(i) NULL) {
  function(sim_indices, seed = NULL, tick = function() {}) {
    if (!is.null(seed)) set.seed(seed)
    parts <- vector("list", length(sim_indices))
    for (j in seq_along(sim_indices)) {
      pre(j)
      parts[[j]] <- hook(j, env$simulate_chunk_r(sim_indices[j], seed = NULL))
    }
    combine_trial_results(parts)
  }
}

test_that("H6: chunk_frame errors if evolveTrial.looks_stop_at_final is unset", {
  old <- options(evolveTrial.looks_stop_at_final = NULL)
  on.exit(options(old))
  expect_error(chunk_frame(d01, 1, 20, LIB), "looks_stop_at_final")
})

test_that("chunk_frame builds a frame with the P6 branch on, without touching the random stream", {
  set.seed(1); s0 <- .Random.seed
  env <- chunk_frame(d01, 1, 20, LIB)
  expect_identical(.Random.seed, s0)
  expect_true(isTRUE(env$looks_stop_at_final))
  expect_true(is.function(env$simulate_chunk_r))
  expect_identical(env$num_simulations, 20L)
})

test_that("H1: reference against itself is identical on R01 and R23 (20 trials)", {
  for (d in list(d01, d23)) {
    env <- chunk_frame(d, 1, 20, LIB)
    cmp <- compare_chunk(env, seq_len(20), 4242L, ref = "simulate_chunk_r", cand = env$simulate_chunk_r)
    expect_true(cmp$same_result)
    expect_true(cmp$same_rng)
  }
})

test_that("H2: trial-by-trial composition of the R body equals the R body (30 trials, R01 and R23)", {
  for (d in list(d01, d23)) {
    env <- chunk_frame(d, 2, 30, LIB)
    cand <- composed_cand(env)
    cmp <- compare_chunk(env, seq_len(30), 31337L, ref = "simulate_chunk_r", cand = cand)
    expect_true(cmp$same_result)
    expect_true(cmp$same_rng)
  }
})

test_that("H3: chunk plan plus R body reproduces the run_scenarios summary columns (R01, R = 200, 4 workers)", {
  r <- plan_emulation_check(d01, LIB, R = 200L, workers = 4L, master_seed = 900001)
  expect_gt(length(r$columns), 10L)
  bad <- names(r$columns)[!unlist(r$columns)]
  expect_identical(bad, character(0))
  expect_true(r$same_rng_after_run)
})

test_that("H4: bisection finds a stream fault at trial 7 of 20", {
  env <- chunk_frame(d01, 1, 20, LIB)
  cand <- composed_cand(env, pre = function(j) if (j == 7L) runif(1))
  b <- bisect_chunk(env, seq_len(20), 777L, ref = "simulate_chunk_r", cand = cand)
  expect_identical(b$k, 7L)
  expect_identical(b$kind, "rng")
  expect_identical(b$sim_index, 7L)
})

test_that("H5: bisection finds a value fault at trial 12 of 20", {
  env <- chunk_frame(d01, 1, 20, LIB)
  cand <- composed_cand(env, hook = function(j, res) {
    if (j == 12L) { res$sum_final_n <- res$sum_final_n + 1; res$final_n_raw[[1]] <- res$final_n_raw[[1]] + 1 }
    res
  })
  b <- bisect_chunk(env, seq_len(20), 777L, ref = "simulate_chunk_r", cand = cand)
  expect_identical(b$k, 12L)
  expect_identical(b$kind, "result")
  expect_equal(unname(b$cand_trial$sum_delta$sum_final_n - b$ref_trial$sum_delta$sum_final_n), 1)
})

test_that("reference cache returns the stored reference and keys on the frame", {
  env <- chunk_frame(d01, 1, 20, LIB)
  tmp <- tempfile("cache"); dir.create(tmp)
  a <- cached_reference(env, seq_len(10), 5L, cache_dir = tmp)
  expect_false(a$from_cache)
  b <- cached_reference(env, seq_len(10), 5L, cache_dir = tmp)
  expect_true(b$from_cache)
  expect_identical(a$value, b$value)
  c2 <- cached_reference(env, seq_len(10), 6L, cache_dir = tmp)
  expect_false(c2$from_cache)
})
