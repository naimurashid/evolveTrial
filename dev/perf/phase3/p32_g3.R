# p32_g3.R -- P3-2 (copy of p31_g3.R): the phase 2 G3 determinism check (dev/perf/phase2/g3_check.R, same runs and pass rule) on two P3-2
# frames from dev/perf/phase3/p32_frames.R. Seeds are fresh (3240001, 3240002).
# Usage: Rscript dev/perf/phase3/p32_g3.R <lib> <R> <frames, comma separated> <out.csv>
root <- normalizePath("~/research/evolveTrial-merge")
Sys.setenv(EVT_ROOT = root)
setwd(root)
source(file.path(root, "dev/perf/harness/chunk_harness.R"))
source(file.path(root, "dev/perf/phase3/p32_frames.R"))
a <- commandArgs(TRUE)
lib <- normalizePath(a[1])
R <- as.integer(a[2])
dn <- strsplit(a[3], ",")[[1]]
bss <- c(13L, 128L)
workers <- c(1L, 4L, 16L, 19L)
seeds <- c(3240001L, 3240002L)

options(evolveTrial.looks_stop_at_final = NULL, evolveTrial.engine = "exact")
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
  ba$progress <- FALSE
  tab <- evolveTrial::run_scenarios(ba, d$scens, parallel = FALSE, seed = seed, return_variance = TRUE)
  list(tab = tab, rng = get(".Random.seed", envir = globalenv()), tr = trials_seen(),
       wall = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}

rows <- list()
for (D in dn) {
  d <- p31_frame_inputs(p32_frame_specs()[[D]])
  seed <- seeds[match(D, dn)]
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
out <- a[4]
write.csv(res, out, row.names = FALSE)
cat(sprintf("G3: %d runs, %d designs, %d failed -> %s\n", nrow(res), length(dn), sum(!res$pass), if (all(res$pass)) "PASS" else "FAIL"))
if (!all(res$pass)) print(res[!res$pass, ])
quit(status = if (all(res$pass)) 0L else 1L)
