# p2_t7.R <lib> <design> <workers> <engine cpp|exact> <out.rds> -- adapted from step7_timing/scripts/t7.R (P2: engine argument).
# Both scenarios, 1,000 trials, seed 900001, production path (internal), as in t7.R. P2: sets options(evolveTrial.engine)
# after loading and records it, the elapsed seconds per 1,000 trials and the user+system CPU seconds (self and children).
a <- commandArgs(TRUE); LIB <- normalizePath(a[1]); DES <- a[2]; W <- as.integer(a[3]); ENG <- a[4]
stopifnot(ENG %in% c("cpp", "exact"))
if (!exists("%||%", mode = "function")) `%||%` <- function(x, y) if (is.null(x)) y else x
suppressPackageStartupMessages(library(evolveTrial))
stopifnot(startsWith(find.package("evolveTrial"), LIB))
options(evolveTrial.looks_stop_at_final = TRUE, evolveTrial.engine = ENG)                    # P2
if (ENG == "exact" && !exists(".simulate_chunk_exact", envir = asNamespace("evolveTrial"), inherits = FALSE)) stop("no exact engine in ", LIB)
x <- readRDS(file.path("/home/naimrashid/research/evolveTrial-perf/dev/perf/evidence/args", sprintf("args_%s.rds", DES)))
ba <- x$base_args; ba$num_simulations <- 1000L; ba$parallel_replicates <- TRUE; ba$num_workers <- W
pt <- system.time(res <- run_scenarios(ba, x$scens, parallel = FALSE, seed = 900001, return_variance = TRUE))
row <- data.frame(design = DES, lib = basename(LIB), engine = ENG, workers = W, elapsed = pt[["elapsed"]],
                  cpu = sum(pt[c("user.self", "sys.self", "user.child", "sys.child")], na.rm = TRUE),
                  digest = digest::digest(list(as.data.frame(res), .Random.seed)),
                  sha = unname(packageDescription("evolveTrial")$EngineCppSHA %||% NA), load1 = as.numeric(strsplit(readLines("/proc/loadavg"), " ")[[1]][1]))
print(row); saveRDS(row, a[5])
