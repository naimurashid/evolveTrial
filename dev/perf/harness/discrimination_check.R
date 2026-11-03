# How often would a naive implementation differ from R on the unit-test inputs? (shows the tests can fail)
# Run: Rscript dev/perf/harness/discrimination_check.R   (uses the test generator from tests/testthat/test-r-semantics.R)
src <- readLines("tests/testthat/test-r-semantics.R")
i0 <- grep("^gen_inputs <- function", src); i1 <- grep("^test_that", src)[1] - 1L
eval(parse(text = src[i0:i1]))
cs <- gen_inputs(2000L, 11L); cm <- gen_inputs(2000L, 12L)
naive_sum  <- vapply(cs, function(x) Reduce(`+`, x, 0), 0)
naive_mean <- vapply(cm, function(x) if (length(x)) Reduce(`+`, x, 0) / length(x) else NaN, 0)
cat("sum:  naive double accumulation differs from sum() on", sum(!mapply(identical, as.list(naive_sum), lapply(cs, sum))), "of 2000 inputs\n")
cat("mean: naive sum/n differs from mean() on", sum(!mapply(identical, as.list(naive_mean), lapply(cm, mean))), "of 2000 inputs\n")
K <- 1:40
cat("seq:  i/K differs from seq(0,1,length.out=K+1) for", sum(!vapply(K, function(k) identical((0:k)/k, seq(0,1,length.out=k+1)), TRUE)), "of 40 K\n")
mlg <- 0; for (n in 1:2000) { cn <- 0:n; v <- vapply(cn, function(c) mean(c(rep(TRUE, c), rep(FALSE, n - c))), 0); mlg <- mlg + sum(v != cn / n) }
cat("mean_logical: plain count/n differs from mean() in", mlg, "of 2003000-2000 cases (n = 1..2000)\n")
