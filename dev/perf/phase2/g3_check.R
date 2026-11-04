# g3_check.R -- gate G3 (determinism) of phase 2 (dev/plans/2026-11-04-phase2-plan.md Section 3; step P2-2).
#
# For each acceptance design (dev/perf/evidence/args/args_R??.rds, the step 2a captures, loaded with the phase 1 harness
# loader), run_scenarios() with options(evolveTrial.engine = "exact") at R = 2,000 and the capture's cert seed, at 1, 4,
# 16 and 19 workers (parallel_replicates = TRUE) and two block sizes, plus one sequential run (parallel_replicates =
# FALSE). Every returned table and the .Random.seed left behind must be identical() to the first; a run with seed + 1 must
# give a different table. Writes dev/perf/phase2/g3_results.csv (one row per run) and exits 1 on any failure.
#
# Usage: Rscript dev/perf/phase2/g3_check.R [lib = dev/perf/lib/dev_p2] [R = 2000] [designs = all] [block sizes = 13,128]
root <- local({ p <- normalizePath(getwd()); while (!file.exists(file.path(p, "DESCRIPTION"))) p <- dirname(p); p })
Sys.setenv(EVT_ROOT = root)
source(file.path(root, "dev/perf/harness/chunk_harness.R"))
a <- commandArgs(TRUE)
lib <- normalizePath(file.path(root, if (length(a) >= 1) a[1] else "dev/perf/lib/dev_p2"))
R <- as.integer(if (length(a) >= 2) a[2] else 2000)
evdir <- file.path(root, "dev/perf/evidence/args")
dn <- if (length(a) >= 3 && a[3] != "all") strsplit(a[3], ",")[[1]] else
  sub("^args_(R[0-9]+)\\.rds$", "\\1", list.files(evdir, pattern = "^args_R[0-9]+\\.rds$"))
bss <- as.integer(strsplit(if (length(a) >= 4) a[4] else "13,128", ",")[[1]])
workers <- c(1L, 4L, 16L, 19L)

options(evolveTrial.looks_stop_at_final = TRUE, evolveTrial.engine = "exact")
harness_load(lib)
stopifnot(exists(".simulate_chunk_exact", envir = asNamespace("evolveTrial"), inherits = FALSE))

# Count the trials through the exact entry (per process, files, as p2_cert.R does), so a silent phase 1 fallback shows.
pdir <- tempfile("g3probe"); dir.create(pdir)
ns <- asNamespace("evolveTrial"); orig_ex <- get(".simulate_chunk_exact", envir = ns); orig_cpp <- get(".simulate_chunk_cpp", envir = ns)
utils::assignInNamespace(".simulate_chunk_exact", function(sim_indices, cfg, scenario_seed) {
  cat(length(sim_indices), "\n", file = file.path(pdir, paste0("ex_", Sys.getpid())), append = TRUE)
  orig_ex(sim_indices, cfg, scenario_seed)
}, ns = "evolveTrial")
utils::assignInNamespace(".simulate_chunk_cpp", function(sim_indices, cfg) stop("phase 1 entry reached in exact mode"), ns = "evolveTrial")
trials_seen <- function() {
  f <- list.files(pdir, full.names = TRUE)
  n <- sum(unlist(lapply(f, function(x) scan(x, quiet = TRUE))))
  unlink(f)
  c(trials = n, pids = length(f))
}

run_one <- function(d, w, bs, seq, seed) {
  ba <- d$base_args
  ba$num_simulations <- R
  ba$num_workers <- w
  ba$parallel_replicates <- !seq
  options(evolveTrial.exact_block_size = bs)
  t0 <- Sys.time()
  tab <- evolveTrial::run_scenarios(ba, d$scens, parallel = FALSE, seed = seed,
                                    return_variance = isTRUE(d$dots$return_variance))
  list(tab = tab, rng = get(".Random.seed", envir = globalenv()), tr = trials_seen(),
       wall = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}

rows <- list()
for (D in dn) {
  d <- load_design(file.path(evdir, paste0("args_", D, ".rds")))
  seed <- as.integer(d$cert_seed %||% d$dots$seed)
  nsc <- length(d$scens)
  runs <- list()
  for (w in workers) for (bs in bss) runs[[sprintf("w%d_bs%d", w, bs)]] <- c(w = w, bs = bs, seq = 0)
  runs[["seq"]] <- c(w = 1, bs = bss[1], seq = 1)
  base <- NULL
  for (nm in names(runs)) {
    r <- runs[[nm]]
    o <- run_one(d, r[["w"]], r[["bs"]], r[["seq"]] == 1, seed)
    if (is.null(base)) base <- o
    rows[[length(rows) + 1L]] <- data.frame(
      design = D, run = nm, workers = r[["w"]], block_size = r[["bs"]], sequential = r[["seq"]] == 1, seed = seed,
      trials_exact = o$tr[["trials"]], pids = o$tr[["pids"]], trials_expected = nsc * R,
      identical_table = identical(o$tab, base$tab), identical_rng = identical(o$rng, base$rng),
      power_or_t1_s1 = o$tab$Type_I_Error_or_Power[1], wall_s = round(o$wall, 2), check = "same", stringsAsFactors = FALSE)
  }
  o <- run_one(d, 19L, bss[1], FALSE, seed + 1L)
  rows[[length(rows) + 1L]] <- data.frame(
    design = D, run = "w19_seed+1", workers = 19, block_size = bss[1], sequential = FALSE, seed = seed + 1L,
    trials_exact = o$tr[["trials"]], pids = o$tr[["pids"]], trials_expected = nsc * R,
    identical_table = identical(o$tab, base$tab), identical_rng = identical(o$rng, base$rng),
    power_or_t1_s1 = o$tab$Type_I_Error_or_Power[1], wall_s = round(o$wall, 2), check = "differs", stringsAsFactors = FALSE)
  message(sprintf("%s done", D))
}
res <- do.call(rbind, rows)
res$pass <- with(res, trials_exact == trials_expected &
                   ifelse(check == "same", identical_table & identical_rng, !identical_table))
out <- file.path(root, "dev/perf/phase2/g3_results.csv")
write.csv(res, out, row.names = FALSE)
cat(sprintf("G3: %d runs, %d designs, %d failed -> %s\n", nrow(res), length(dn), sum(!res$pass), if (all(res$pass)) "PASS" else "FAIL"))
if (!all(res$pass)) print(res[!res$pass, ])
quit(status = if (all(res$pass)) 0L else 1L)
