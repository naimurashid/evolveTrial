# layer3_compare.R -- layer 3 gate (plan step 6, dev/plans/2026-11-03-step6-layer3.md, decision 2).
#   Rscript layer3_compare.R <dir of layer3_search.R outputs> <E2c results dir> <out.csv> [ref|dev]
# Each output is joined to the stored row (ledger counts, and the run's own rds for the full doubles). A row passes only if
# power_count and type1_count equal the ledger counts as integers, objective (and EN_null and EN_alt for BATON rows) are
# identical() to the rds doubles, its engine check and all its capture checks hold, it ran in the expected mode with
# L3_FORCE_R unset, its input file hashes equal those of the results dir, and it is the only output for its candidate.
# The set must be exactly the three candidates. Pass means 3 of 3. Exit 0 on pass, 1 otherwise.
# For dev, L3_ENGINE_SHA (40 hex) is required and each row's engine commit must equal it.
a <- commandArgs(TRUE)
if (!length(a) %in% 3:4) stop("usage: Rscript layer3_compare.R <out dir> <results dir> <out.csv> [ref|dev]")
MODE <- if (length(a) == 4L) a[4] else "dev"
stopifnot(MODE %in% c("ref", "dev"))
ENGINE <- Sys.getenv("L3_ENGINE_SHA")
if (MODE == "dev" && !grepl("^[0-9a-f]{40}$", ENGINE)) stop("set L3_ENGINE_SHA to the engine commit under test")
CAND <- data.frame(run_id = c("PB-baton-1001", "PB-baton-1003", "PB-ga-1001"), eval_index = c(1L, 95L, 274L), stringsAsFactors = FALSE)
fs <- sort(list.files(a[1], "\\.rds$", full.names = TRUE, recursive = TRUE))
if (!length(fs)) stop("no RDS outputs under ", a[1])
sha256 <- function(f) sub(" .*", "", system2("sha256sum", f, stdout = TRUE))
rows <- do.call(rbind, lapply(fs, function(f) {
  x <- readRDS(f); led <- utils::read.csv(file.path(a[2], x$run_id, "sim_ledger.csv"), stringsAsFactors = FALSE)
  lr <- led[led$purpose == "search" & led$eval_index == x$eval_index, ]
  rdsf <- file.path(a[2], x$run_id, if (isTRUE(x$is_baton)) "f0_out/doublet_triplet_admissible_fut_warmstart_result.rds" else "run_state.rds")
  rd <- readRDS(rdsf)
  e <- if (isTRUE(x$is_baton)) { h <- rd$history_combined; k <- which(h$fidelity == "low" & h$eval_id == x$eval_index)
           list(objective = h$objective[k], EN_null = h$EN_null[k], EN_alt = h$EN_alt[k]) }
       else { h <- rd$history; k <- which(h$eval_index == x$eval_index); list(objective = h$objective[k], EN_null = NA_real_, EN_alt = NA_real_) }
  led_p <- round(lr$power * 3000); led_t <- round(lr$type1 * 3000)
  cnt_ok <- nrow(lr) == 1L && identical(as.integer(x$counts[["power_count"]]), as.integer(led_p)) &&
            identical(as.integer(x$counts[["type1_count"]]), as.integer(led_t))
  dbl_ok <- identical(x$objective, e$objective) && (!isTRUE(x$is_baton) || (identical(x$EN_null, e$EN_null) && identical(x$EN_alt, e$EN_alt)))
  hash_ok <- identical(unname(x$input_sha256[["ledger"]]), sha256(file.path(a[2], x$run_id, "sim_ledger.csv"))) &&
             identical(unname(x$input_sha256[["rds"]]), sha256(rdsf))
  data.frame(file = f, run_id = x$run_id, eval_index = x$eval_index, mode = x$mode, force_r = isTRUE(x$force_r),
             engine_ok = isTRUE(x$engine_ok), checks_ok = is.logical(x$checks) && length(x$checks) == 9L && all(x$checks),
             engine_sha = if (is.null(x$engine_sha) || is.na(x$engine_sha)) NA_character_ else x$engine_sha,
             power_count = as.integer(x$counts[["power_count"]]), power_count_ledger = as.integer(led_p),
             type1_count = as.integer(x$counts[["type1_count"]]), type1_count_ledger = as.integer(led_t),
             counts_match = cnt_ok, objective_identical = identical(x$objective, e$objective),
             EN_null_identical = if (isTRUE(x$is_baton)) identical(x$EN_null, e$EN_null) else NA,
             EN_alt_identical = if (isTRUE(x$is_baton)) identical(x$EN_alt, e$EN_alt) else NA,
             doubles_match = dbl_ok, hash_ok = hash_ok, wall_s = as.numeric(x$wall_s), host = x$host, r_version = x$r_version,
             stringsAsFactors = FALSE)
}))
key <- paste(rows$run_id, rows$eval_index)
rows$candidate <- key %in% paste(CAND$run_id, CAND$eval_index)
dup <- duplicated(key) | duplicated(key, fromLast = TRUE)
rows$pass <- rows$counts_match & rows$doubles_match & rows$hash_ok & rows$engine_ok & rows$checks_ok & rows$candidate & !dup &
  rows$mode %in% MODE & !rows$force_r & (MODE == "ref" | rows$engine_sha %in% ENGINE)
rows <- rows[order(rows$run_id, rows$eval_index), ]
utils::write.csv(rows, a[3], row.names = FALSE)
missing <- sum(!paste(CAND$run_id, CAND$eval_index) %in% key)
cat(sprintf("[layer3] %d outputs, %d candidates, %d counts match, %d doubles identical, %d pass; engine_ok %d, checks_ok %d, mode %s %d, force_r %d, duplicate %d, missing %d\n",
            nrow(rows), sum(rows$candidate), sum(rows$counts_match), sum(rows$doubles_match), sum(rows$pass), sum(rows$engine_ok),
            sum(rows$checks_ok), MODE, sum(rows$mode %in% MODE), sum(rows$force_r), sum(dup), missing))
if (any(!rows$pass)) print(rows[!rows$pass, c("run_id", "eval_index", "mode", "force_r", "engine_ok", "checks_ok", "power_count", "power_count_ledger",
                                              "type1_count", "type1_count_ledger", "objective_identical", "EN_null_identical", "EN_alt_identical", "hash_ok")], row.names = FALSE)
ok <- all(rows$pass) && missing == 0L && nrow(rows) == nrow(CAND) && (MODE == "ref" || identical(unique(rows$engine_sha), ENGINE))
cat(if (ok) "[layer3] PASS" else "[layer3] FAIL", sprintf("(%d of %d)", sum(rows$pass), nrow(CAND)), "\n")
quit(save = "no", status = if (ok) 0L else 1L)
