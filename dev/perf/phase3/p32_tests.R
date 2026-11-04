# P3-2 gate T: the package test suite (tests/testthat) against a built scratch library, as in P3-1.
# Usage: Rscript dev/perf/phase3/p32_tests.R <lib> <out.rds> [filter]   (EVOLVETRIAL_LONG_TESTS=1 for the long tier)
a <- commandArgs(TRUE)
setwd(normalizePath("~/research/evolveTrial-merge"))
library(evolveTrial, lib.loc = normalizePath(a[1]))
library(testthat)
cat("evolveTrial", as.character(packageVersion("evolveTrial", lib.loc = normalizePath(a[1]))), "from", a[1], "\n")
t0 <- Sys.time()
res <- testthat::test_dir("tests/testthat", filter = if (length(a) >= 3) a[3] else NULL, reporter = "summary",
                          stop_on_failure = FALSE, load_package = "none",
                          env = new.env(parent = asNamespace("evolveTrial")))   # internal functions visible, as R CMD check
df <- as.data.frame(res)
saveRDS(df, a[2])
agg <- aggregate(cbind(nb, failed, skipped = as.integer(skipped), error = as.integer(error), warning, passed) ~ file, df, sum)
print(agg)
cat(sprintf("TOTAL passed %d failed %d skipped %d errors %d warnings %d\n", sum(df$passed), sum(df$failed),
            sum(df$skipped), sum(df$error), sum(df$warning)))
cat("wall seconds", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "\n")
