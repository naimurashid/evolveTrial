# Compare two layer 0 outputs. Usage: Rscript layer0_compare.R <ref.rds> <dev.rds>
a <- commandArgs(TRUE)
x <- readRDS(a[1]); y <- readRDS(a[2])
t_ok <- identical(x$res, y$res); r_ok <- identical(x$rng, y$rng)
cat(sprintf("layer0 %s: table identical=%s  .Random.seed identical=%s  (%d rows x %d cols)\n",
            x$design, t_ok, r_ok, nrow(x$res), ncol(x$res)))
quit(status = if (t_ok && r_ok) 0L else 1L)
