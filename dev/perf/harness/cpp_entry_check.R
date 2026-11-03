# Gate 3 of step 3: show that the C++ chunk body actually ran, by counting entries into the compiled entry point
# evolveTrial:::.simulate_chunk_cpp (traced in the namespace; FORK workers inherit the trace and append to one file).
# The tracer draws nothing and prints nothing, so the traced run must reproduce the untraced gate output exactly.
#
#   Rscript dev/perf/harness/cpp_entry_check.R layer1 <R01,R02,...>
#       reruns layer 1 (cached reference) with the trace on; expects entries == chunks compared (no bisection probes).
#   Rscript dev/perf/harness/cpp_entry_check.R layer0 <lib> <args.rds> <use_cpp_engine> <gate_out.rds>
#       reruns the layer 0 call (run_scenarios, R = 200, 4 workers, seed 900001) with the trace on and the option set,
#       checks identical() to the saved gate output (table and .Random.seed), and reports the entry count.
root <- local({ p <- normalizePath(getwd()); while (!file.exists(file.path(p, "DESCRIPTION"))) p <- dirname(p); p })
a <- commandArgs(TRUE)
trace_file <- tempfile("cpp_entries_")
invisible(file.create(trace_file))
Sys.setenv(EVT_CPP_TRACE = trace_file)
install_trace <- function() {
  suppressMessages(trace(".simulate_chunk_cpp", where = asNamespace("evolveTrial"), print = FALSE,
                         tracer = quote(cat("1\n", file = Sys.getenv("EVT_CPP_TRACE"), append = TRUE))))
}
entries <- function() length(readLines(trace_file))

if (identical(a[1], "layer1")) {
  Sys.setenv(EVT_ROOT = root)
  source(file.path(root, "dev/perf/harness/chunk_harness.R"))
  lib <- normalizePath(file.path(root, "dev/perf/lib/dev"))
  options(evolveTrial.looks_stop_at_final = TRUE)
  harness_load(lib)
  install_trace()
  dn <- strsplit(a[2], ",")[[1]]
  paths <- file.path(root, "dev/perf/evidence/args", paste0("args_", dn, ".rds"))
  res <- run_layer1(paths, lib, verbose = FALSE)
  bad <- sum(!(res$same_result & res$same_rng))
  cat(sprintf("cpp_entry_check layer1 %s: %d chunks compared, %d mismatches, %d entries into .simulate_chunk_cpp\n",
              a[2], nrow(res), bad, entries()))
  quit(status = if (bad == 0L && entries() == nrow(res)) 0L else 1L)
}

if (identical(a[1], "layer0")) {
  lib <- normalizePath(a[2]); args_rds <- a[3]; use_cpp <- as.logical(a[4]); gate_out <- a[5]
  stopifnot(!is.na(use_cpp), identical(RNGkind(), c("Mersenne-Twister", "Inversion", "Rejection")))
  options(evolveTrial.use_cpp_engine = use_cpp)
  suppressPackageStartupMessages(library(evolveTrial, lib.loc = lib))
  stopifnot(dirname(getNamespaceInfo("evolveTrial", "path")) == lib)
  options(evolveTrial.looks_stop_at_final = TRUE)
  install_trace()
  d <- readRDS(args_rds)
  base_args <- d$base_args
  base_args$num_simulations <- 200L
  base_args$num_workers <- 4L
  res <- evolveTrial::run_scenarios(base_args, d$scens, parallel = FALSE, seed = 900001, return_variance = TRUE)
  rng <- .Random.seed
  # round trip through saveRDS/readRDS as the gate files were, so that data.table's .internal.selfref external pointer
  # (live here, nil after readRDS) does not make identical() fail on an otherwise identical table
  tf <- tempfile(fileext = ".rds"); saveRDS(list(res = res), tf); res <- readRDS(tf)$res
  g <- readRDS(gate_out)
  same <- identical(g$res, res) && identical(g$rng, rng)
  cat(sprintf("cpp_entry_check layer0 %s use_cpp_engine=%s: identical to gate output=%s, %d entries into .simulate_chunk_cpp\n",
              d$design, use_cpp, same, entries()))
  quit(status = if (same) 0L else 1L)
}
stop("usage: cpp_entry_check.R layer1 <designs> | layer0 <lib> <args.rds> <use_cpp_engine> <gate_out.rds>")
