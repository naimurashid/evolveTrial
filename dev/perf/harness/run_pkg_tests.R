# Run the package test suite against a private library (absolute path needed: testthat changes directory).
# Usage: Rscript dev/perf/harness/run_pkg_tests.R <lib> <out.rds> [filter]
a <- commandArgs(TRUE)
lib <- normalizePath(a[1])
library(evolveTrial, lib.loc = lib)
stopifnot(dirname(getNamespaceInfo("evolveTrial", "path")) == lib)
filt <- if (length(a) >= 3) a[3] else NULL
t0 <- Sys.time()
res <- testthat::test_dir("tests/testthat", package = "evolveTrial", load_package = "installed",
                          filter = filt, reporter = "summary", stop_on_failure = FALSE)
df <- as.data.frame(res)
cat("PASSED", sum(df$passed), "FAILED", sum(df$failed), "SKIPPED", sum(df$skipped),
    "ERRORS", sum(df$error), "WARN", sum(df$warning), "\n")
print(aggregate(cbind(passed, failed, skipped, error) ~ file, df, sum))
cat("wall seconds", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "\n")
saveRDS(df, a[2])
