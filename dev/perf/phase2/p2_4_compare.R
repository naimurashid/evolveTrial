# p2_4_compare.R -- compare two p2_4_identity.R outputs: identical() per design (chunks, table, .Random.seed), and for
# any chunk that differs, the number of trials whose per-trial outcomes differ and the largest relative double gap.
# Usage: Rscript dev/perf/phase2/p2_4_compare.R <old.rds> <new.rds>
a <- commandArgs(TRUE); o <- readRDS(a[1]); n <- readRDS(a[2])
stopifnot(identical(names(o), names(n)))
ok_all <- TRUE
for (D in names(o)) {
  ch <- mapply(identical, o[[D]]$chunks, n[[D]]$chunks)
  tb <- identical(o[[D]]$tab, n[[D]]$tab); rg <- identical(o[[D]]$rng, n[[D]]$rng)
  ntr <- sapply(o[[D]]$chunks, function(x) NROW(x$final_n))
  cat(sprintf("%s seed=%d trials_per_scenario=%s chunks_identical=%s table_identical=%s rng_identical=%s\n", D,
              o[[D]]$seed, paste(ntr, collapse = "/"), paste(ch, collapse = "/"), tb, rg))
  if (!all(ch) || !tb || !rg) ok_all <- FALSE
  for (s in which(!ch)) {
    x <- o[[D]]$chunks[[s]]; y <- n[[D]]$chunks[[s]]
    for (f in names(x)) if (!identical(x[[f]], y[[f]])) {
      dx <- as.matrix(x[[f]]); dy <- as.matrix(y[[f]])
      rows <- if (nrow(dx) == ntr[s]) sum(rowSums(dx != dy | xor(is.na(dx), is.na(dy)), na.rm = TRUE) > 0) else NA
      rel <- suppressWarnings(max(abs(dx - dy) / pmax(abs(dx), 1e-300), na.rm = TRUE))
      cat(sprintf("   scenario %d field %s: trials differing %s, max rel gap %.3g\n", s, f, rows, rel))
    }
  }
}
cat(if (ok_all) "ALL IDENTICAL\n" else "DIFFERENCES FOUND\n")
