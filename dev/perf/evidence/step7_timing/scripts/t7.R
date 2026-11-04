# t7.R <lib> <design> <workers> <out.rds>: both scenarios, 1,000 trials, seed 900001, production path (internal)
a <- commandArgs(TRUE); LIB <- normalizePath(a[1]); DES <- a[2]; W <- as.integer(a[3])
suppressPackageStartupMessages(library(evolveTrial))
stopifnot(startsWith(find.package("evolveTrial"), LIB))
options(evolveTrial.looks_stop_at_final = TRUE)
x <- readRDS(file.path("/home/naimrashid/research/evolveTrial-perf/dev/perf/evidence/args", sprintf("args_%s.rds", DES)))
ba <- x$base_args; ba$num_simulations <- 1000L; ba$parallel_replicates <- TRUE; ba$num_workers <- W
pt <- system.time(res <- run_scenarios(ba, x$scens, parallel = FALSE, seed = 900001, return_variance = TRUE))
row <- data.frame(design = DES, lib = basename(LIB), workers = W, elapsed = pt[["elapsed"]],
                  digest = digest::digest(list(as.data.frame(res), .Random.seed)),
                  sha = unname(packageDescription("evolveTrial")$EngineCppSHA %||% NA))
print(row); saveRDS(row, a[4])
