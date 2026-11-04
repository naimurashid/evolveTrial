# P3-2 gate S (copy of p31_gateS_compare.R, unchanged test): exact engine against fresh R-body runs, per frame, with the phase 2 test (dev/perf/phase2/p2_compare.R,
# G5): two one-sided tests at level 0.05 each, margins 0.01 on power, 0.005 on type I error and 1% of the pooled R-body
# mean on expected sample size (null and alternative), plus the classification checks (power on the same side of 0.80,
# type I error on the same side of 0.10). Inputs: the p32_gateS_run.R outputs under <dir>. The experimental arm is the
# arm the table reports for single-arm designs and the non-reference arm for between-arm designs. Scenario 1 is the null
# (type I error, EN_null), scenario 2 the alternative (power, EN_alt), as in the captures.
# Also reports the power of each test: the probability that both one-sided tests reject when the true difference is 0,
# 1 - 2 * pnorm(qnorm(0.95) - margin / se) (floored at 0), at the standard errors of this run.
# Usage: Rscript p32_gateS_compare.R <dir> <out.tsv>
a <- commandArgs(TRUE)
MARGIN_POWER <- 0.01; MARGIN_TYPE1 <- 0.005; MARGIN_EN_REL <- 0.01; ALPHA <- 0.05; THR_POWER <- 0.80; THR_TYPE1 <- 0.10
zc <- qnorm(1 - ALPHA)
fs <- list.files(a[1], "\\.rds$", full.names = TRUE)
runs <- lapply(fs, readRDS)
exp_arm <- function(tab) {
  arms <- unique(tab$Arm_Name)
  if (length(arms) == 1L) arms else if ("Triplet" %in% arms) "Triplet" else stop("unknown arm layout")
}
rows <- do.call(rbind, lapply(runs, function(r) {
  t <- r$tab; arm <- exp_arm(t); t <- t[t$Arm_Name == arm, ]
  t0 <- t[t$scenario == 1L, ]; t1 <- t[t$scenario == 2L, ]
  data.frame(frame = r$frame, engine = r$engine, seed = r$seed, n = r$n, wall = r$wall,
             type1 = t0$Type_I_Error_or_Power, power = t1$Type_I_Error_or_Power,
             EN_null = t0$Exp_N, EN_alt = t1$Exp_N,
             EN_null_var = t0$Var_Exp_N * r$n, EN_alt_var = t1$Var_Exp_N * r$n,   # per-trial variance
             stringsAsFactors = FALSE)
}))
tost <- function(d, se, m) {
  z_lo <- ifelse(se > 0, (d + m) / se, ifelse(d + m > 0, Inf, -Inf)); z_hi <- ifelse(se > 0, (m - d) / se, ifelse(m - d > 0, Inf, -Inf))
  list(z_lo = z_lo, z_hi = z_hi, ok = z_lo > zc & z_hi > zc, pow = ifelse(se > 0, pmax(0, 1 - 2 * pnorm(zc - m / se)), 1))
}
pooled <- function(q, col) sum(q[[col]] * q$n) / sum(q$n)
vtot <- function(q, mcol, vcol) {   # law of total variance over seeds, weighted by n
  m <- pooled(q, mcol); sum(q$n * (q[[vcol]] + (q[[mcol]] - m)^2)) / sum(q$n)
}
out <- list()
for (fr in sort(unique(rows$frame))) {
  r <- rows[rows$frame == fr & rows$engine == "r", ]; x <- rows[rows$frame == fr & rows$engine == "exact", ]
  if (!nrow(r) || !nrow(x)) next
  n1 <- sum(r$n); n2 <- sum(x$n)
  res <- list(frame = fr, R_r = n1, R_exact = n2, seeds_r = paste(r$seed, collapse = ","), seeds_exact = paste(x$seed, collapse = ","),
              wall_r = sum(r$wall), wall_exact = sum(x$wall))
  for (k in c("power", "type1")) {
    q1 <- pooled(r, k); q2 <- pooled(x, k); d <- q2 - q1
    se <- sqrt(q2 * (1 - q2) / n2 + q1 * (1 - q1) / n1)
    m <- if (k == "power") MARGIN_POWER else MARGIN_TYPE1
    t <- tost(d, se, m)
    res[paste0(k, c("_r", "_exact", "_diff", "_se", "_zlo", "_zhi", "_ok", "_test_power"))] <- list(q1, q2, d, se, t$z_lo, t$z_hi, t$ok, t$pow)
  }
  for (k in c("EN_null", "EN_alt")) {
    m1 <- pooled(r, k); m2 <- pooled(x, k); d <- m2 - m1
    v <- vtot(x, k, paste0(k, "_var"))
    se <- sqrt(v / n2 + v / n1); m <- MARGIN_EN_REL * m1
    t <- tost(d, se, m)
    res[paste0(k, c("_r", "_exact", "_diff", "_margin", "_se", "_zlo", "_zhi", "_ok", "_test_power"))] <- list(m1, m2, d, m, se, t$z_lo, t$z_hi, t$ok, t$pow)
  }
  res$class_ok <- ((res$power_exact >= THR_POWER) == (res$power_r >= THR_POWER)) &&
    ((res$type1_exact <= THR_TYPE1) == (res$type1_r <= THR_TYPE1))
  res$frame_pass <- res$power_ok && res$type1_ok && res$EN_null_ok && res$EN_alt_ok && res$class_ok
  out[[fr]] <- as.data.frame(res, stringsAsFactors = FALSE)
}
tab <- do.call(rbind, out)
write.table(tab, a[2], sep = "\t", quote = FALSE, row.names = FALSE)
f <- function(x) sprintf("%.4f", x)
for (i in seq_len(nrow(tab))) with(tab[i, ], cat(sprintf(
  "%-22s R %6d / exact %6d | power %s vs %s (z %.1f, %.1f; test power %.2f) %s | type I %s vs %s (z %.1f, %.1f; %.2f) %s | EN0 %s vs %s (z %.1f, %.1f; %.2f) %s | EN1 %s vs %s (z %.1f, %.1f; %.2f) %s | class %s | %s\n",
  frame, R_r, R_exact, f(power_r), f(power_exact), power_zlo, power_zhi, power_test_power, power_ok,
  f(type1_r), f(type1_exact), type1_zlo, type1_zhi, type1_test_power, type1_ok,
  f(EN_null_r), f(EN_null_exact), EN_null_zlo, EN_null_zhi, EN_null_test_power, EN_null_ok,
  f(EN_alt_r), f(EN_alt_exact), EN_alt_zlo, EN_alt_zhi, EN_alt_test_power, EN_alt_ok, class_ok,
  if (frame_pass) "PASS" else "FAIL")))
cat(sprintf("gate S: %d frames, %d pass\n", nrow(tab), sum(tab$frame_pass)))
