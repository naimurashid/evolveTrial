# Layer 1 runner (plan 2b.4):  Rscript dev/perf/harness/run_layer1.R [--designs R01,R23] [--cand simulate_chunk_cpp]
#   [--lib dev/perf/lib/dev] [--R 50] [--long 200] [--workers 16,19] [--masters 900001,900002,900003] [--no-cache]
# Designs default to every args_R??.rds in dev/perf/evidence/args/ (step 2a output), the only inputs accepted.
# Writes dev/perf/harness/out/layer1_<timestamp>.csv, prints the number of chunks compared and of mismatches,
# and exits 1 on any mismatch.
root <- local({ p <- normalizePath(getwd()); while (!file.exists(file.path(p, "DESCRIPTION"))) p <- dirname(p); p })
Sys.setenv(EVT_ROOT = root)
source(file.path(root, "dev/perf/harness/chunk_harness.R"))
a <- commandArgs(TRUE)
opt <- function(name, default) { i <- match(paste0("--", name), a); if (is.na(i)) default else a[i + 1L] }
lib     <- normalizePath(file.path(root, opt("lib", "dev/perf/lib/dev")))
cand    <- opt("cand", "simulate_chunk_cpp")
R       <- as.integer(opt("R", "50"))
long    <- as.integer(opt("long", "200"))
workers <- as.integer(strsplit(opt("workers", "16,19"), ",")[[1]])
masters <- as.numeric(strsplit(opt("masters", "900001,900002,900003"), ",")[[1]])
cache   <- if ("--no-cache" %in% a) NULL else default_cache_dir()
evdir   <- file.path(root, "dev/perf/evidence")
dn <- opt("designs", NA)
if (is.na(dn)) {
  f <- list.files(file.path(evdir, "args"), pattern = "^args_R[0-9]+\\.rds$")
  dn <- paste(sub("^args_(R[0-9]+)\\.rds$", "\\1", f), collapse = ",")
}
dn <- strsplit(dn, ",")[[1]]
paths <- vapply(dn, function(d) {
  p1 <- file.path(evdir, "args", paste0("args_", d, ".rds"))
  if (file.exists(p1)) p1 else stop("no args file for ", d, " in dev/perf/evidence/args (step 2a captures only)")
}, character(1))
options(evolveTrial.looks_stop_at_final = TRUE)
harness_load(lib)
cat(sprintf("layer 1: designs %s, candidate %s, lib %s, R %d, workers %s, masters %s, long chunk %d\n",
            paste(dn, collapse = ","), cand, lib, R, paste(workers, collapse = ","), paste(masters, collapse = ","), long))
t0 <- Sys.time()
res <- run_layer1(unname(paths), lib, masters, workers, R, long, cand, cache)
dir.create(file.path(root, "dev/perf/harness/out"), showWarnings = FALSE)
f <- file.path(root, "dev/perf/harness/out", paste0("layer1_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"))
write.csv(res, f, row.names = FALSE)
bad <- !(res$same_result & res$same_rng)
cat(sprintf("layer 1: %d chunks compared, %d mismatches (%d with result differing, %d with stream differing), %.0f s, %s\n",
            nrow(res), sum(bad), sum(!res$same_result), sum(!res$same_rng),
            as.numeric(difftime(Sys.time(), t0, units = "secs")), f))
if (any(bad)) print(res[bad, ])
quit(status = if (any(bad)) 1L else 0L)
