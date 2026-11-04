# P3-2 gate B (copy of p31_gateB.R on the P3-2 frames of p32_frames.R, counters of p32_instrument.R): the phase 1 engine
# (simulate_chunk_cpp) against the R body (simulate_chunk_r) on the P3-2 acceptance frames, identical() chunk result and .Random.seed after every chunk, chunk plans as in layer 1
# (master seed x worker count, one sample.int seed per chunk, plus one long chunk per scenario). On every chunk also:
#   - the instrumented R body (p31_instrument.R) gives the same result and stream as the R body, and its branch counts
#     are accumulated per frame;
#   - the exact file's replay instantiation (.simulate_chunk_exact_replay, R stream x Monte Carlo posteriors) gives the
#     same result and stream as phase 1, which checks the exact engine's copy of the trial loop.
# Frames marked expect_error (min_follow_up_at_final < 0): every engine must stop with the R body's message, and the
# stream at the stop must be the same.
# Frames with same_as: the R body on the frame must equal the R body on the frame with same_as applied (result and
# stream), chunk by chunk (rule TRUE with a NULL or Inf probability threshold equals rule FALSE).
# Per frame it also reports interim efficacy and futility stops (trials, from the R body's chunk sums).
# Usage: Rscript dev/perf/phase3/p32_gateB.R <lib> <out.rds> [frame regex]
a <- commandArgs(TRUE)
setwd(normalizePath("~/research/evolveTrial-merge"))
library(evolveTrial, lib.loc = normalizePath(a[1]))
options(evolveTrial.looks_stop_at_final = NULL, evolveTrial.engine = NULL)
source("dev/perf/phase3/p32_frames.R")
source("dev/perf/phase3/p32_instrument.R")
source("dev/perf/harness/chunk_harness.R")   # chunk_plan()
replay <- getFromNamespace(".simulate_chunk_exact_replay", "evolveTrial")
pat <- if (length(a) >= 3) a[3] else "."
masters <- c(3210001, 3210002)
workers <- c(16L, 19L)
rng_now <- function() if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
run1 <- function(f, idx, seed) {
  r <- tryCatch(list(result = f(idx, seed = seed), error = NULL), error = function(e) list(result = NULL, error = conditionMessage(e)))
  r$rng <- rng_now()
  r
}
specs <- p32_frame_specs()
specs <- specs[grepl(pat, names(specs))]
rows <- list(); counts <- list(); supported <- list(); stops <- list()
t0 <- Sys.time()
for (fn in names(specs)) {
  spec <- specs[[fn]]
  R <- if (spec$D == "R23") 24L else 40L
  long <- if (spec$D == "R23") 60L else 120L
  nsc <- length(p31_frame_inputs(spec)$scens)
  p31_counter_reset()
  stops[[fn]] <- c(eff = 0, fut = 0)
  for (s in seq_len(nsc)) {
    for (sz in c(R, long)) {
      env <- p31_chunk_frame(spec, s, sz)
      supported[[paste(fn, s)]] <- list(cfg = isTRUE(env$cfg$supported),
                                        violations = getFromNamespace(".chunk_cpp_violations", "evolveTrial")(env))
      if (!isTRUE(env$cfg$supported)) stop(fn, " scenario ", s, " is not admitted: ",
                                           paste(supported[[paste(fn, s)]]$violations, collapse = "; "))
      inst <- p31_instrument(env)
      env2 <- if (!is.null(spec$same_as)) {
        sp2 <- spec; sp2$mod <- utils::modifyList(spec$mod, spec$same_as, keep.null = TRUE); sp2$same_as <- NULL
        p31_chunk_frame(sp2, s, sz)
      } else NULL
      plans <- list()
      if (sz == R) {
        for (ms in masters) for (w in workers) {
          p <- chunk_plan(ms, R, w, nsc)[[s]]
          for (j in seq_along(p$chunks)) plans[[length(plans) + 1L]] <- list(ms = ms, w = w, j = j, idx = p$chunks[[j]], seed = p$seeds[[j]])
        }
      } else {
        p <- chunk_plan(masters[1], R, workers[1], nsc)[[s]]
        plans[[1]] <- list(ms = masters[1], w = workers[1], j = 0L, idx = seq_len(long), seed = p$seeds[[1]])
      }
      for (pl in plans) {
        ref <- run1(env$simulate_chunk_r, pl$idx, pl$seed)
        ins <- p32_with_trace(run1(inst, pl$idx, pl$seed))
        if (!is.null(ref$result)) stops[[fn]] <- stops[[fn]] + c(sum(ref$result$sum_stop_eff), sum(ref$result$sum_stop_fut))
        same_as_ok <- if (is.null(env2)) NA else {
          r2 <- run1(env2$simulate_chunk_r, pl$idx, pl$seed)
          identical(r2$result, ref$result) && identical(r2$rng, ref$rng)
        }
        cpp <- run1(env$simulate_chunk_cpp, pl$idx, pl$seed)
        rep <- run1(function(idx, seed) { set.seed(seed); replay(idx, env$cfg) }, pl$idx, pl$seed)
        same_err <- function(x, y) identical(is.null(x$error), is.null(y$error)) &&
          (is.null(x$error) || identical(sub("^Error[^:]*: ", "", x$error), sub("^Error[^:]*: ", "", y$error)))
        rows[[length(rows) + 1L]] <- data.frame(
          frame = fn, design = spec$D, scenario = s, master = pl$ms, workers = pl$w, chunk = pl$j, n = length(pl$idx),
          ref_error = !is.null(ref$error), ref_msg = ref$error %||% "",
          cpp_same_result = identical(cpp$result, ref$result), cpp_same_rng = identical(cpp$rng, ref$rng),
          cpp_same_error = same_err(cpp, ref), cpp_msg = cpp$error %||% "",
          inst_same = identical(ins$result, ref$result) && identical(ins$rng, ref$rng) && same_err(ins, ref),
          replay_same_result = identical(rep$result, cpp$result), replay_same_rng = identical(rep$rng, cpp$rng),
          replay_same_error = same_err(rep, ref),
          nan_free = is.null(cpp$result) || !anyNA(unlist(cpp$result)),
          same_as_ok = same_as_ok,
          stringsAsFactors = FALSE)
      }
    }
  }
  counts[[fn]] <- p31_counts()
  r <- do.call(rbind, rows[vapply(rows, function(x) x$frame == fn, logical(1))])
  ok <- with(r, cpp_same_result & cpp_same_rng & cpp_same_error & inst_same & replay_same_result & replay_same_rng &
               replay_same_error & nan_free & (is.na(same_as_ok) | same_as_ok))
  ev <- counts[[fn]]$events
  message(sprintf("%-26s chunks %4d trials %5d errors %3d mismatches %d same_as %s stops eff %d fut %d  [%s] {looks: %s} %.0f s",
                  fn, nrow(r), sum(r$n), sum(r$ref_error), sum(!ok),
                  if (all(is.na(r$same_as_ok))) "-" else sprintf("%d/%d", sum(r$same_as_ok), nrow(r)),
                  as.integer(stops[[fn]][1]), as.integer(stops[[fn]][2]),
                  paste(sprintf("%s=%d", names(counts[[fn]]$trial_hits), counts[[fn]]$trial_hits), collapse = " "),
                  paste(sprintf("%s=%d", names(ev)[grepl("^(hc|vr)_", names(ev))], ev[grepl("^(hc|vr)_", names(ev))]), collapse = " "),
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
res <- do.call(rbind, rows)
saveRDS(list(rows = res, counts = counts, supported = supported, stops = stops), a[2])
ok <- with(res, cpp_same_result & cpp_same_rng & cpp_same_error & inst_same & replay_same_result & replay_same_rng &
             replay_same_error & nan_free & (is.na(same_as_ok) | same_as_ok))
cat(sprintf("gate B: %d frames, %d chunks, %d trials, %d chunks ending in the R body's error, %d mismatches, %.0f s\n",
            length(specs), nrow(res), sum(res$n), sum(res$ref_error), sum(!ok),
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
if (any(!ok)) print(res[!ok, ])
