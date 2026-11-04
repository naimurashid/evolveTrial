# p2_compare.R -- gate G5 (dev/plans/2026-11-04-phase2-plan.md Section 3): statistical equivalence of the phase 2 engine
# against dev/perf/evidence/acceptance_set_reference_rows.csv.
#   Rscript p2_compare.R <p2 input> <reference csv> <out dir> [--var-csv <file>] [--subset]
#     <p2 input>  a directory of p2_cert.R outputs (mode p2; the gate input), or a CSV with columns design, seed, R,
#                 power_count, type1_count, EN_null, EN_alt (self-tests only; row checks on mode and engine are skipped
#                 and the report says so).
#     --var-csv   per-design per-trial variances (columns design, EN_null_var, EN_alt_var); required for CSV input,
#                 not allowed for rds input (the rds files carry the variance of the phase 2 run).
#     --subset    do not require all 75 reference rows (smoke tests; the verdict then says SUBSET and never PASS-gate).
# Writes <out dir>/g5_report.tsv (one row per design) and g5_report.md. Exit 0 PASS, 1 FAIL, 2 usage or input error.
#
# Test, exactly as the plan states it, per design (counts pooled over the design's seeds for both engines):
#  1. Power and type I error: two one-sided tests (TOST) on the difference of proportions phase 2 minus reference, unpooled
#     standard error sqrt(p2(1-p2)/n2 + p1(1-p1)/n1), each one-sided test at level 0.05 (z > qnorm(0.95)), margins 0.01 for
#     power (power_count) and 0.005 for type I error (type1_count). Equivalent when both reject.
#  2. EN_null and EN_alt: TOST on the difference of pooled means, level 0.05 each, margin 1% of the pooled reference mean;
#     standard error sqrt(v/n2 + v/n1) with v the per-trial variance of the phase 2 run applied to both samples. v is the
#     pooled variance over the design's seeds by the law of total variance (mean of the within-seed variances plus the
#     variance, divisor S, of the seed means), each within-seed variance being the variance attribute of the wrapper result
#     (variance of the mean) times R.
#  3. Classification: pooled power on the same side of 0.80 (>= 0.80) and pooled type I error on the same side of 0.10
#     (<= 0.10) for both engines.
# A design passes when all of 1 to 3 hold. G5 passes when all 10 designs pass. Margins are constants below; they are not
# arguments, because the plan forbids changing them after a result is seen.
MARGIN_POWER <- 0.01; MARGIN_TYPE1 <- 0.005; MARGIN_EN_REL <- 0.01; ALPHA <- 0.05; THR_POWER <- 0.80; THR_TYPE1 <- 0.10
a <- commandArgs(TRUE)
usage <- "usage: Rscript p2_compare.R <p2 dir|csv> <reference csv> <out dir> [--var-csv <file>] [--subset]"
if (length(a) < 3L) { message(usage); quit(save = "no", status = 2L) }
SUBSET <- "--subset" %in% a; vi <- match("--var-csv", a); VARCSV <- if (!is.na(vi)) a[vi + 1L] else NA_character_
IN <- a[1]; REFF <- a[2]; OUT <- a[3]; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
die <- function(...) { message("[p2_compare] ERROR: ", ...); quit(save = "no", status = 2L) }
need <- c("design", "seed", "R", "power_count", "type1_count", "EN_null", "EN_alt")
ref <- utils::read.csv(REFF, stringsAsFactors = FALSE)
if (!all(need %in% names(ref))) die("reference csv lacks columns")
if (!SUBSET && (nrow(ref) != 75L || length(unique(ref$design)) != 10L || any(ref$R != 10000L))) die("reference is not the 75-row, 10-design, R = 10000 set")
if (anyDuplicated(ref[c("design", "seed")])) die("duplicate reference rows")
isdir <- dir.exists(IN)
if (isdir) {
  if (!is.na(VARCSV)) die("--var-csv is for CSV input only; rds files carry the variance")
  fs <- sort(list.files(IN, "\\.rds$", full.names = TRUE, recursive = TRUE)); if (!length(fs)) die("no rds under ", IN)
  g1 <- function(x, n) { v <- x[[n]]; if (is.null(v) || length(v) != 1L) NA else v }
  p2 <- do.call(rbind, lapply(fs, function(f) { x <- readRDS(f)
    data.frame(file = f, design = g1(x, "design"), seed = as.integer(g1(x, "cert_seed")), mode = g1(x, "mode"), force_r = isTRUE(g1(x, "force_r")),
               engine_ok = isTRUE(g1(x, "engine_ok")), checks_ok = is.logical(x$checks) && length(x$checks) == 8L && all(x$checks),
               R = as.integer(g1(x, "R")), power_count = unname(x$counts[["power_count"]]), type1_count = unname(x$counts[["type1_count"]]),
               EN_null = as.numeric(g1(x, "EN_null")), EN_alt = as.numeric(g1(x, "EN_alt")),
               EN_null_var = as.numeric(g1(x, "EN_null_var")), EN_alt_var = as.numeric(g1(x, "EN_alt_var")),   # variance of the mean
               worktree_sha = g1(x, "worktree_sha"), p2_sha = g1(x, "p2_engine_sha_env"), stringsAsFactors = FALSE) }))
  dup <- duplicated(p2[c("design", "seed")]) | duplicated(p2[c("design", "seed")], fromLast = TRUE)
  p2$rowok <- p2$mode %in% "p2" & !p2$force_r & p2$engine_ok & p2$checks_ok & !dup & p2$R %in% 10000L &
    is.finite(p2$EN_null_var) & is.finite(p2$EN_alt_var)
  p2$vn <- p2$EN_null_var * p2$R; p2$va <- p2$EN_alt_var * p2$R                                    # per-trial variance
  src <- sprintf("%d rds files under %s (mode p2 required; engine_ok, checks, R = 10000, no duplicates, variance present)", nrow(p2), IN)
  shas <- paste(unique(paste0("worktree ", p2$worktree_sha, " P2_ENGINE_SHA ", p2$p2_sha)), collapse = "; ")
} else {
  if (is.na(VARCSV)) die("CSV input needs --var-csv")
  p2 <- utils::read.csv(IN, stringsAsFactors = FALSE); if (!all(need %in% names(p2))) die("input csv lacks columns")
  vt <- utils::read.csv(VARCSV, stringsAsFactors = FALSE); if (!all(c("design", "EN_null_var", "EN_alt_var") %in% names(vt))) die("--var-csv lacks columns")
  dup <- duplicated(p2[c("design", "seed")]) | duplicated(p2[c("design", "seed")], fromLast = TRUE)
  p2$vn <- vt$EN_null_var[match(p2$design, vt$design)]; p2$va <- vt$EN_alt_var[match(p2$design, vt$design)]
  p2$rowok <- !dup & p2$R %in% 10000L & is.finite(p2$vn) & is.finite(p2$va)
  src <- sprintf("CSV %s with per-trial variances from %s (row checks on mode and engine NOT applicable)", IN, VARCSV); shas <- "n/a (csv input)"
}
zc <- qnorm(1 - ALPHA)
tost <- function(d, se, m) { z_lo <- ifelse(se > 0, (d + m) / se, ifelse(d + m > 0, Inf, -Inf)); z_hi <- ifelse(se > 0, (m - d) / se, ifelse(m - d > 0, Inf, -Inf))
  list(z_lo = z_lo, z_hi = z_hi, ok = z_lo > zc & z_hi > zc) }
prop <- function(x2, n2, x1, n1, m) { q2 <- x2 / n2; q1 <- x1 / n1; d <- q2 - q1; se <- sqrt(q2 * (1 - q2) / n2 + q1 * (1 - q1) / n1)
  t <- tost(d, se, m); data.frame(p_ref = q1, p_p2 = q2, diff = d, se = se, margin = m, z_lo = t$z_lo, z_hi = t$z_hi, ok = t$ok) }
meanc <- function(m2, n2, m1, n1, v, rel) { d <- m2 - m1; se <- sqrt(v / n2 + v / n1); m <- rel * m1; t <- tost(d, se, m)
  data.frame(m_ref = m1, m_p2 = m2, diff = d, se = se, margin = m, z_lo = t$z_lo, z_hi = t$z_hi, ok = t$ok) }
ds <- sort(unique(ref$design)); res <- list(); notes <- character(0)
for (dn in ds) {
  r <- ref[ref$design == dn, ]; q <- p2[p2$design == dn, ]
  miss <- setdiff(r$seed, q$seed); extra <- setdiff(q$seed, r$seed)
  q <- q[q$seed %in% r$seed, ]; r <- r[r$seed %in% q$seed, ]
  complete <- !length(miss) && !length(extra) && nrow(q) > 0L
  rowsok <- nrow(q) > 0L && all(q$rowok)
  if (!nrow(q) && SUBSET) next
  if (!nrow(q)) { res[[dn]] <- data.frame(design = dn, n_seeds = 0L, complete = FALSE, rows_ok = FALSE, design_pass = FALSE); next }
  q <- q[match(r$seed, q$seed), ]
  n1 <- sum(r$R); n2 <- sum(q$R); S <- nrow(q)
  P <- prop(sum(q$power_count), n2, sum(r$power_count), n1, MARGIN_POWER)
  T1 <- prop(sum(q$type1_count), n2, sum(r$type1_count), n1, MARGIN_TYPE1)
  vtot <- function(m, v) mean(v) + mean((m - mean(m))^2)                                           # law of total variance, equal R
  E0 <- meanc(mean(q$EN_null), n2, mean(r$EN_null), n1, vtot(q$EN_null, q$vn), MARGIN_EN_REL)
  E1 <- meanc(mean(q$EN_alt), n2, mean(r$EN_alt), n1, vtot(q$EN_alt, q$va), MARGIN_EN_REL)
  cp <- (P$p_p2 >= THR_POWER) == (P$p_ref >= THR_POWER); ct <- (T1$p_p2 <= THR_TYPE1) == (T1$p_ref <= THR_TYPE1)
  tests_ok <- P$ok && T1$ok && E0$ok && E1$ok
  res[[dn]] <- data.frame(design = dn, n_seeds = S, complete = complete, rows_ok = rowsok, R_ref = n1, R_p2 = n2,
    power_ref = P$p_ref, power_p2 = P$p_p2, power_diff = P$diff, power_se = P$se, power_zlo = P$z_lo, power_zhi = P$z_hi, power_ok = P$ok,
    type1_ref = T1$p_ref, type1_p2 = T1$p_p2, type1_diff = T1$diff, type1_se = T1$se, type1_zlo = T1$z_lo, type1_zhi = T1$z_hi, type1_ok = T1$ok,
    ENnull_ref = E0$m_ref, ENnull_p2 = E0$m_p2, ENnull_diff = E0$diff, ENnull_margin = E0$margin, ENnull_se = E0$se, ENnull_zlo = E0$z_lo, ENnull_zhi = E0$z_hi, ENnull_ok = E0$ok,
    ENalt_ref = E1$m_ref, ENalt_p2 = E1$m_p2, ENalt_diff = E1$diff, ENalt_margin = E1$margin, ENalt_se = E1$se, ENalt_zlo = E1$z_lo, ENalt_zhi = E1$z_hi, ENalt_ok = E1$ok,
    power_class_ok = cp, type1_class_ok = ct,
    power_dist_se_from_0.80 = (P$p_ref - THR_POWER) / P$se, type1_dist_se_from_0.10 = (T1$p_ref - THR_TYPE1) / T1$se,
    design_pass = (complete || SUBSET) && rowsok && tests_ok && cp && ct)
  if (length(miss)) notes <- c(notes, sprintf("%s: missing seeds %s", dn, paste(miss, collapse = ",")))
  if (length(extra)) notes <- c(notes, sprintf("%s: seeds not in reference %s", dn, paste(extra, collapse = ",")))
}
allcols <- names(res[[which.max(vapply(res, ncol, 1L))]])
for (dn in names(res)) for (cn in setdiff(allcols, names(res[[dn]]))) res[[dn]][[cn]] <- NA
tab <- do.call(rbind, lapply(res, function(z) z[allcols])); rownames(tab) <- NULL
utils::write.table(tab, file.path(OUT, "g5_report.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
nd <- length(ds); ok <- all(tab$design_pass) && (SUBSET || (nrow(tab) == 10L && all(tab$complete)))
verdict <- if (ok && SUBSET) "SUBSET-PASS (not a gate result)" else if (ok) "PASS" else "FAIL"
f4 <- function(x) ifelse(is.na(x), "NA", sprintf("%.4f", x)); yn <- function(x) ifelse(is.na(x), "-", ifelse(x, "yes", "NO"))
md <- c("# G5 statistical equivalence report", "",
        sprintf("Verdict: **%s** (%d of %d designs pass)", verdict, sum(tab$design_pass), nrow(tab)), "",
        sprintf("Input: %s", src), sprintf("Engine identity: %s", shas), sprintf("Reference: %s", REFF),
        sprintf("Margins: power %.3f, type I %.3f, EN %.0f%% of the pooled reference mean; each one-sided test at level %.2f (z > %.4f).",
                MARGIN_POWER, MARGIN_TYPE1, 100 * MARGIN_EN_REL, ALPHA, zc), "",
        "Each z pair is (z_lo, z_hi) of the two one-sided tests; a test rejects when its z exceeds the critical value.", "",
        "| design | seeds | rows ok | power ref | power p2 | diff | z (lo, hi) | ok | type I ref | type I p2 | diff | z (lo, hi) | ok | EN0 ref | EN0 p2 | z (lo, hi) | ok | EN1 ref | EN1 p2 | z (lo, hi) | ok | class | design |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
zz <- function(l, h) sprintf("%.1f, %.1f", l, h)
for (i in seq_len(nrow(tab))) with(tab[i, ], md <<- c(md, sprintf("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
  design, n_seeds, yn(rows_ok), f4(power_ref), f4(power_p2), sprintf("%+.4f", power_diff), zz(power_zlo, power_zhi), yn(power_ok),
  f4(type1_ref), f4(type1_p2), sprintf("%+.4f", type1_diff), zz(type1_zlo, type1_zhi), yn(type1_ok),
  f4(ENnull_ref), f4(ENnull_p2), zz(ENnull_zlo, ENnull_zhi), yn(ENnull_ok), f4(ENalt_ref), f4(ENalt_p2), zz(ENalt_zlo, ENalt_zhi), yn(ENalt_ok),
  yn(power_class_ok & type1_class_ok), ifelse(design_pass, "PASS", "FAIL"))))
fl <- tab[!tab$design_pass, ]
if (nrow(fl)) md <- c(md, "", "Failing designs: " , paste0("- ", fl$design, ": complete=", fl$complete, " rows_ok=", fl$rows_ok, " power_ok=", fl$power_ok, " type1_ok=", fl$type1_ok,
  " ENnull_ok=", fl$ENnull_ok, " ENalt_ok=", fl$ENalt_ok, " power_class_ok=", fl$power_class_ok, " type1_class_ok=", fl$type1_class_ok))
near <- tab[!is.na(tab$power_dist_se_from_0.80) & (abs(tab$power_dist_se_from_0.80) < 3 | abs(tab$type1_dist_se_from_0.10) < 3), ]
if (nrow(near)) md <- c(md, "", "Classification sensitivity (reference pooled value within 3 standard errors of its threshold): ",
  paste0("- ", near$design, ": power is ", sprintf("%.1f", near$power_dist_se_from_0.80), " SE from 0.80, type I is ", sprintf("%.1f", near$type1_dist_se_from_0.10), " SE from 0.10"))
if (length(notes)) md <- c(md, "", "Input notes:", paste0("- ", notes))
writeLines(md, file.path(OUT, "g5_report.md"))
cat(sprintf("[p2_compare] %s: %d of %d designs pass; report in %s\n", verdict, sum(tab$design_pass), nrow(tab), OUT))
if (!ok) print(tab[!tab$design_pass, c("design", "complete", "rows_ok", "power_ok", "type1_ok", "ENnull_ok", "ENalt_ok", "power_class_ok", "type1_class_ok")], row.names = FALSE)
quit(save = "no", status = if (ok) 0L else 1L)
