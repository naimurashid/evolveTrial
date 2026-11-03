# chunk_harness.R -- layer 1 harness for the C++ engine port (plan 2b.4).
#
# Compares a candidate chunk body with the R chunk body (`simulate_chunk_r`, the renamed closure inside
# run_simulation_pure) chunk by chunk, on the returned list and on .Random.seed, and bisects a mismatch to the first
# diverging trial. Source this file; every function takes the library to test from an explicit argument or from a
# frame that was built from it.
#
# A candidate is any function callable as fn(sim_indices, seed = <int or NULL>, tick = function() {}) that returns the
# same list as simulate_chunk_r. `cand` may be the name of an object in the frame (env) or in the evolveTrial
# namespace, or the function itself.

`%||%` <- function(x, y) if (is.null(x)) y else x

# ---- paths -----------------------------------------------------------------------------------------------------------
harness_root <- function() {
  r <- Sys.getenv("EVT_ROOT", "")
  if (nzchar(r)) return(normalizePath(r))
  p <- normalizePath(getwd())
  repeat {
    if (file.exists(file.path(p, "DESCRIPTION")) &&
        any(grepl("^Package: evolveTrial", readLines(file.path(p, "DESCRIPTION"), n = 3L)))) return(p)
    if (identical(dirname(p), p)) stop("cannot find the evolveTrial worktree from ", getwd(), "; set EVT_ROOT")
    p <- dirname(p)
  }
}

# ---- loading ---------------------------------------------------------------------------------------------------------
# Every entry point asserts the RNG kind (plan Section 0) and that evolveTrial comes from `lib`.
harness_load <- function(lib) {
  lib <- normalizePath(lib)
  stopifnot(identical(RNGkind(), c("Mersenne-Twister", "Inversion", "Rejection")))
  if (!"evolveTrial" %in% loadedNamespaces()) {
    suppressPackageStartupMessages(library(evolveTrial, lib.loc = lib))
  }
  got <- dirname(getNamespaceInfo("evolveTrial", "path"))
  if (!identical(got, lib)) stop("evolveTrial is loaded from ", got, ", not from ", lib)
  invisible(lib)
}

load_design <- function(path) {
  d <- readRDS(path)
  stopifnot(is.list(d), !is.null(d$base_args), !is.null(d$scens), !is.null(d$dots))
  d
}

# ---- the frame -------------------------------------------------------------------------------------------------------
# Builds args_i exactly as run_scenarios' run_one does (simulation_driver.R:1280-1288), then calls run_simulation_pure
# with a probe that returns its own environment before any chunk runs (the seam at simulation_driver.R, after
# aggregate_results). The returned environment holds the closures and the argument list the chunk body sees.
chunk_frame <- function(design, scenario, R, lib) {
  harness_load(lib)
  # I7: looks_stop_at_final is a lazy default read from this option. Require it before building the frame.
  if (!isTRUE(getOption("evolveTrial.looks_stop_at_final"))) {
    stop("options(evolveTrial.looks_stop_at_final = TRUE) must be set for the whole session before chunk_frame() ",
         "(plan I7); the frozen designs run the looks_stop_at_final branch.")
  }
  rsp <- get("run_simulation_pure", envir = asNamespace("evolveTrial"))
  over <- design$scens[[scenario]]
  args_i <- utils::modifyList(design$base_args, over, keep.null = TRUE)
  args_i$return_percentiles <- FALSE
  args_i$percentile_probs   <- c(0, 0.25, 0.5, 0.75, 0.9, 1.0)
  args_i$return_variance    <- design$dots$return_variance
  args_i <- args_i[intersect(names(args_i), names(formals(rsp)))]
  args_i$num_simulations <- as.integer(R)
  old <- options(evolveTrial.chunk_probe = function(env) env)
  on.exit(options(old), add = TRUE)
  env <- do.call(rsp, args_i)
  if (!is.environment(env)) {
    stop("run_simulation_pure did not return its frame: the chunk_probe seam is missing from the installed package")
  }
  if (!is.function(env$simulate_chunk_r)) stop("frame has no simulate_chunk_r closure (seam missing?)")
  if (!isTRUE(env$looks_stop_at_final)) stop("frame has looks_stop_at_final != TRUE (plan I7)")
  env
}

# ---- the chunk plan --------------------------------------------------------------------------------------------------
# Mirrors simulation_driver.R:976-1043 (splitIndices, empty chunks dropped, one sample.int seed per chunk) and
# run_scenarios:1272 (one set.seed, then the scenarios in order), for R at or above the sequential threshold.
chunk_plan <- function(master_seed, R, workers, n_scen = 2L) {
  set.seed(master_seed)
  plan <- vector("list", n_scen)
  for (s in seq_len(n_scen)) {
    w <- max(1L, min(as.integer(workers), as.integer(R)))
    chunks <- parallel::splitIndices(R, w)
    chunks <- chunks[lengths(chunks) > 0]
    seeds <- replicate(length(chunks), sample.int(.Machine$integer.max, 1L), simplify = TRUE)
    plan[[s]] <- list(chunks = chunks, seeds = seeds)
  }
  attr(plan, "rng_after") <- .Random.seed
  plan
}

# ---- running and comparing -------------------------------------------------------------------------------------------
resolve_fn <- function(fn, env) {
  if (is.function(fn)) return(fn)
  stopifnot(is.character(fn), length(fn) == 1L)
  if (exists(fn, envir = env, inherits = FALSE)) return(get(fn, envir = env, inherits = FALSE))
  ns <- asNamespace("evolveTrial")
  if (exists(fn, envir = ns, inherits = FALSE)) return(get(fn, envir = ns, inherits = FALSE))
  stop("candidate '", fn, "' not found in the frame or the evolveTrial namespace")
}

run_chunk <- function(fn, indices, seed) {
  res <- fn(indices, seed = seed)
  list(result = res, rng = if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL)
}

# ---- reference cache -------------------------------------------------------------------------------------------------
# Keyed by md5 of deparse(simulate_chunk_r), the serialized argument list, the indices and the seed, plus the state the
# reference body reads from outside the closure: looks_stop_at_final, EVOLVETRIAL_USE_CPP, the R version, and the md5 of
# every package source file the R body reaches (all R/*.R except simulation_driver.R and RcppExports.R, all src/*.cpp
# except RcppExports.cpp and the files that belong to the port). Delete the cache directory to force a rerun.
.cache_state <- new.env()
cache_code_md5 <- function(root = harness_root()) {
  if (!is.null(.cache_state[[root]])) return(.cache_state[[root]])
  rf <- list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE)
  rf <- rf[!basename(rf) %in% c("simulation_driver.R", "RcppExports.R")]
  cf <- list.files(file.path(root, "src"), pattern = "\\.(cpp|h)$", full.names = TRUE)
  cf <- cf[!grepl("^(RcppExports|r_semantics|simulate_chunk)", basename(cf))]
  v <- unname(tools::md5sum(c(rf, cf)))
  names(v) <- basename(c(rf, cf))
  .cache_state[[root]] <- v
  v
}

reference_key <- function(env, indices, seed) {
  digest::digest(list(
    body = deparse(env$simulate_chunk_r),
    args = env$args,
    indices = indices, seed = seed,
    lsaf = isTRUE(env$looks_stop_at_final),
    use_cpp = Sys.getenv("EVOLVETRIAL_USE_CPP", ""),
    rver = R.version.string,
    code = cache_code_md5()
  ), algo = "md5")
}

cached_reference <- function(env, indices, seed, cache_dir) {
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  key <- reference_key(env, indices, seed)
  f <- file.path(cache_dir, paste0(key, ".rds"))
  if (file.exists(f)) return(list(value = readRDS(f), from_cache = TRUE, key = key))
  v <- run_chunk(env$simulate_chunk_r, indices, seed)
  saveRDS(v, f)
  list(value = v, from_cache = FALSE, key = key)
}

default_cache_dir <- function() file.path(harness_root(), "dev/perf/harness/cache")

# Runs reference and candidate from the same seed. Returns the two runs so that a caller can look inside.
compare_chunk <- function(env, indices, seed, ref = "simulate_chunk_r", cand = "simulate_chunk_cpp", cache_dir = NULL) {
  r <- if (!is.null(cache_dir) && identical(ref, "simulate_chunk_r")) {
    cached_reference(env, indices, seed, cache_dir)$value
  } else run_chunk(resolve_fn(ref, env), indices, seed)
  c1 <- run_chunk(resolve_fn(cand, env), indices, seed)
  list(same_result = identical(r$result, c1$result), same_rng = identical(r$rng, c1$rng), ref = r, cand = c1)
}

# ---- combining per-trial results (the composition used by the self-tests) ---------------------------------------------
# Recombines the chunk lists of single-trial runs the way one chunk accumulates them: plain double sums from zero in trial
# order, raw vectors concatenated, list order and names as in simulate_chunk_r (simulation_driver.R:885-898).
SUM_FIELDS <- c("sum_final_n", "sum_final_events", "sum_stop_eff", "sum_stop_fut", "sum_final_eff", "sum_final_fut",
                "sum_final_inc")
combine_trial_results <- function(parts) {
  acc <- lapply(SUM_FIELDS, function(f) { v <- parts[[1]][[f]]; setNames(numeric(length(v)), names(v)) })
  names(acc) <- SUM_FIELDS
  tm <- 0
  for (p in parts) {
    for (f in SUM_FIELDS) acc[[f]] <- acc[[f]] + p[[f]]
    tm <- tm + p$sum_final_time
  }
  cat_arm <- function(field) {
    if (is.null(parts[[1]][[field]])) return(NULL)
    arms <- names(parts[[1]][[field]])
    out <- lapply(arms, function(a) unlist(lapply(parts, function(p) p[[field]][[a]]), use.names = FALSE))
    names(out) <- arms
    out
  }
  list(
    sum_final_n = acc$sum_final_n, sum_final_events = acc$sum_final_events,
    sum_stop_eff = acc$sum_stop_eff, sum_stop_fut = acc$sum_stop_fut,
    sum_final_eff = acc$sum_final_eff, sum_final_fut = acc$sum_final_fut, sum_final_inc = acc$sum_final_inc,
    sum_final_time = tm,
    n_sims = sum(vapply(parts, function(p) p$n_sims, integer(1))),
    final_n_raw = cat_arm("final_n_raw"),
    final_events_raw = cat_arm("final_events_raw"),
    final_time_raw = if (is.null(parts[[1]]$final_time_raw)) NULL else
      unlist(lapply(parts, function(p) p$final_time_raw), use.names = FALSE)
  )
}

# ---- bisection -------------------------------------------------------------------------------------------------------
# Smallest k such that running indices[1:k] from `seed` differs between reference and candidate in .Random.seed or in the
# result. Each probe is a fresh run from the seed. Invariant: the prefix lo = 0 is the same, the prefix hi = n differs.
# `kind` is "rng" when the stream state after the first diverging trial differs (a value difference in that trial is then a
# consequence of the shifted stream and is reported in `result_differs`), and "result" when the stream agrees and only the
# returned values differ.
.trial_values <- function(res, k) {
  list(
    final_events_raw = if (is.null(res$final_events_raw)) NULL else vapply(res$final_events_raw, function(v) v[k], numeric(1)),
    final_n_raw = if (is.null(res$final_n_raw)) NULL else vapply(res$final_n_raw, function(v) v[k], numeric(1)),
    final_time_raw = if (is.null(res$final_time_raw)) NULL else res$final_time_raw[k]
  )
}
.sum_delta <- function(now, before) {
  out <- lapply(c(SUM_FIELDS, "sum_final_time"), function(f) now[[f]] - (if (is.null(before)) 0 else before[[f]]))
  names(out) <- c(SUM_FIELDS, "sum_final_time")
  out
}

bisect_chunk <- function(env, indices, seed, ref = "simulate_chunk_r", cand = "simulate_chunk_cpp", cache_dir = NULL) {
  n <- length(indices)
  probe <- function(k) compare_chunk(env, indices[seq_len(k)], seed, ref, cand, cache_dir)
  full <- probe(n)
  if (full$same_result && full$same_rng) return(NULL)
  lo <- 0L; hi <- n
  hi_cmp <- full
  while (hi - lo > 1L) {
    mid <- (lo + hi) %/% 2L
    m <- probe(mid)
    if (m$same_result && m$same_rng) lo <- mid else { hi <- mid; hi_cmp <- m }
  }
  prev <- if (lo >= 1L) probe(lo) else NULL
  k <- hi
  rng_d <- !hi_cmp$same_rng
  res_d <- !hi_cmp$same_result
  list(
    k = as.integer(k), sim_index = indices[k],
    kind = if (rng_d) "rng" else "result",
    rng_differs = rng_d, result_differs = res_d,
    ref_trial  = c(list(values = .trial_values(hi_cmp$ref$result, k)),
                   list(sum_delta = .sum_delta(hi_cmp$ref$result, prev$ref$result))),
    cand_trial = c(list(values = .trial_values(hi_cmp$cand$result, k)),
                   list(sum_delta = .sum_delta(hi_cmp$cand$result, prev$cand$result)))
  )
}

# ---- layer 1 ---------------------------------------------------------------------------------------------------------
# For every design, scenario, master seed and worker count, every chunk of the plan; then one extra chunk seq_len(long_chunk)
# per design and scenario at the first master seed. Returns a data.frame with one row per chunk compared.
run_layer1 <- function(designs, lib, master_seeds = c(900001, 900002, 900003), workers = c(16L, 19L), R = 50L,
                       long_chunk = 200L, cand = "simulate_chunk_cpp", cache_dir = default_cache_dir(), verbose = TRUE) {
  rows <- list()
  add <- function(design, scen, seed, w, chunk, n, cmp, bis) {
    rows[[length(rows) + 1L]] <<- data.frame(
      design = design, scenario = scen, master_seed = seed, workers = w, chunk = chunk, n = n,
      same_result = cmp$same_result, same_rng = cmp$same_rng,
      bisect_k = if (is.null(bis)) NA_integer_ else bis$k,
      bisect_kind = if (is.null(bis)) NA_character_ else bis$kind, stringsAsFactors = FALSE)
  }
  for (dpath in designs) {
    d <- load_design(dpath)
    nsc <- length(d$scens)
    envs <- lapply(seq_len(nsc), function(s) chunk_frame(d, s, R, lib))
    envs_long <- lapply(seq_len(nsc), function(s) chunk_frame(d, s, long_chunk, lib))
    for (ms in master_seeds) for (w in workers) {
      plan <- chunk_plan(ms, R, w, nsc)
      for (s in seq_len(nsc)) for (j in seq_along(plan[[s]]$chunks)) {
        idx <- plan[[s]]$chunks[[j]]; sd <- plan[[s]]$seeds[[j]]
        cmp <- compare_chunk(envs[[s]], idx, sd, "simulate_chunk_r", cand, cache_dir)
        bis <- if (cmp$same_result && cmp$same_rng) NULL else bisect_chunk(envs[[s]], idx, sd, "simulate_chunk_r", cand, cache_dir)
        add(d$design, s, ms, w, j, length(idx), cmp, bis)
      }
      if (verbose) message(sprintf("  %s master %d workers %d done (%d chunks so far)", d$design, ms, w, length(rows)))
    }
    plan <- chunk_plan(master_seeds[1], R, workers[1], nsc)
    for (s in seq_len(nsc)) {
      sd <- plan[[s]]$seeds[[1]]
      idx <- seq_len(long_chunk)
      cmp <- compare_chunk(envs_long[[s]], idx, sd, "simulate_chunk_r", cand, cache_dir)
      bis <- if (cmp$same_result && cmp$same_rng) NULL else bisect_chunk(envs_long[[s]], idx, sd, "simulate_chunk_r", cand, cache_dir)
      add(d$design, s, master_seeds[1], workers[1], 0L, length(idx), cmp, bis)   # chunk 0 = the long chunk
    }
  }
  do.call(rbind, rows)
}

# ---- H3: the chunk plan plus the R body reproduces run_scenarios ------------------------------------------------------
# Aggregates chunk results in chunk order with the arithmetic of simulation_driver.R:856-863, 901-944 and 1058-1135, and
# compares each summary column with the table run_scenarios returns for the same master seed and worker count.
plan_emulation_check <- function(design, lib, R = 200L, workers = 4L, master_seed = 900001) {
  harness_load(lib)
  base_args <- design$base_args
  base_args$num_simulations <- as.integer(R)
  base_args$num_workers <- as.integer(workers)
  tab <- evolveTrial::run_scenarios(base_args, design$scens, parallel = FALSE, seed = master_seed,
                                    return_variance = isTRUE(design$dots$return_variance))
  rng_run <- get(".Random.seed", envir = globalenv())
  nsc <- length(design$scens)
  plan <- chunk_plan(master_seed, R, workers, nsc)
  same_rng_after_run <- identical(rng_run, attr(plan, "rng_after"))
  columns <- list()
  inv <- 1 / R
  for (s in seq_len(nsc)) {
    env <- chunk_frame(design, s, R, lib)
    arms <- env$arm_names
    parts <- Map(function(idx, sd) env$simulate_chunk_r(idx, seed = sd), plan[[s]]$chunks, plan[[s]]$seeds)
    zero <- function(f) setNames(numeric(length(arms)), arms)
    S <- list(n = zero(), ev = zero(), se = zero(), sf = zero(), fe = zero(), ff = zero(), fi = zero())
    tm <- 0
    for (p in parts) {   # aggregate_results, in chunk order
      S$n <- S$n + p$sum_final_n; S$ev <- S$ev + p$sum_final_events
      S$se <- S$se + p$sum_stop_eff; S$sf <- S$sf + p$sum_stop_fut
      S$fe <- S$fe + p$sum_final_eff; S$ff <- S$ff + p$sum_final_fut; S$fi <- S$fi + p$sum_final_inc
      tm <- tm + p$sum_final_time
    }
    rawcat <- function(field, a) unlist(lapply(parts, function(p) p[[field]][[a]]), use.names = FALSE)
    obs <- tab[tab$scenario == s, ]
    emu <- list(
      PET_Efficacy = S$se * inv, PET_Futility = S$sf * inv,
      Pr_Final_Efficacy = S$fe * inv, Pr_Final_Futility = S$ff * inv, Pr_Final_Inconclusive = S$fi * inv,
      Exp_N = S$n * inv, Exp_Events = S$ev * inv)
    emu$Type_I_Error_or_Power <- emu$PET_Efficacy + emu$Pr_Final_Efficacy
    emu$Pr_Reach_Max_N <- 1 - emu$PET_Efficacy - emu$PET_Futility
    emu$Exp_Time <- rep(tm * inv, length(arms))
    if (isTRUE(design$dots$return_variance)) {
      emu$Var_Exp_N <- vapply(arms, function(a) stats::var(rawcat("final_n_raw", a)) / R, numeric(1))
      emu$Var_Exp_Events <- vapply(arms, function(a) stats::var(rawcat("final_events_raw", a)) / R, numeric(1))
      emu$Var_Exp_Time <- rep(stats::var(unlist(lapply(parts, function(p) p$final_time_raw), use.names = FALSE)) / R, length(arms))
      pv <- function(p) p * (1 - p) / R
      emu$Var_Type_I_Error_or_Power <- pv(emu$Type_I_Error_or_Power)
      emu$Var_PET_Efficacy <- pv(emu$PET_Efficacy); emu$Var_PET_Futility <- pv(emu$PET_Futility)
      emu$Var_Pr_Reach_Max_N <- pv(emu$Pr_Reach_Max_N)
      emu$Var_Pr_Final_Efficacy <- pv(emu$Pr_Final_Efficacy); emu$Var_Pr_Final_Futility <- pv(emu$Pr_Final_Futility)
      emu$Var_Pr_Final_Inconclusive <- pv(emu$Pr_Final_Inconclusive)
    }
    for (nm in names(emu)) {
      columns[[sprintf("scen%d.%s", s, nm)]] <- identical(unname(as.numeric(obs[[nm]])), unname(as.numeric(emu[[nm]])))
    }
  }
  list(columns = columns, same_rng_after_run = same_rng_after_run, table = tab)
}
