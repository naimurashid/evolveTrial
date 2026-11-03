# layer2_compare.R -- layer 2 gate (plan step 5, dev/plans/2026-11-03-step5-layer2.md, Section 4).
#   Rscript layer2_compare.R <dir of layer2_cert.R outputs> <acceptance_set_reference_rows.csv> <out.csv> [--subset]
# Joins each output RDS to its reference row on (design, seed) and compares power_count and type1_count as integers and
# EN_null and EN_alt at the stored precision (4 decimals, as in capture check I5). A row passes only if all four match,
# its engine check passed, its capture checks all hold, it ran the dev engine with L2_FORCE_R unset, and it is the only
# output for its (design, seed), and its engine commit equals L2_ENGINE_SHA (required). Without --subset the gate also needs every reference row present (75 of 75).
# Exit 0 on pass, 1 otherwise. Counts come from the objects read, not from logs.
a <- commandArgs(TRUE)
if (!length(a) %in% 3:4) stop("usage: Rscript layer2_compare.R <out dir> <reference csv> <out.csv> [--subset]")
ENGINE <- Sys.getenv("L2_ENGINE_SHA")
if (!grepl("^[0-9a-f]{40}$", ENGINE)) stop("set L2_ENGINE_SHA to the engine commit under test")
SUBSET <- length(a) == 4L && identical(a[4], "--subset")
if (length(a) == 4L && !SUBSET) stop("unknown option ", a[4])
ref <- utils::read.csv(a[2], stringsAsFactors = FALSE)
stopifnot(nrow(ref) == 75L, !anyDuplicated(ref[c("design", "seed")]), all(ref$R == 10000L), all(ref$workers == 19L))
fs <- sort(list.files(a[1], "\\.rds$", full.names = TRUE, recursive = TRUE))
if (!length(fs)) stop("no RDS outputs under ", a[1])
g1 <- function(x, n) { v <- x[[n]]; if (is.null(v) || length(v) != 1L) NA else v }
rows <- do.call(rbind, lapply(fs, function(f) {
  x <- readRDS(f)
  data.frame(file = f, design = g1(x, "design"), seed = as.integer(g1(x, "cert_seed")), mode = g1(x, "mode"),
             force_r = isTRUE(g1(x, "force_r")), engine_ok = isTRUE(g1(x, "engine_ok")),
             checks_ok = is.logical(x$checks) && length(x$checks) == 8L && all(x$checks),
             engine_sha = g1(x, "engine_sha"), R = as.integer(g1(x, "R")),
             power_count = unname(x$counts[["power_count"]]), type1_count = unname(x$counts[["type1_count"]]),
             EN_null = as.numeric(g1(x, "EN_null")), EN_alt = as.numeric(g1(x, "EN_alt")),
             wall_s = as.numeric(g1(x, "wall_s")), host = g1(x, "host"), r_version = g1(x, "r_version"),
             stringsAsFactors = FALSE)
}))
m <- merge(rows, ref[c("design", "seed", "power_count", "type1_count", "EN_null", "EN_alt")], by = c("design", "seed"),
           all.x = TRUE, suffixes = c("", "_ref"))
m$in_ref <- !is.na(m$power_count_ref)
dup <- duplicated(m[c("design", "seed")]) | duplicated(m[c("design", "seed")], fromLast = TRUE)
eq4 <- function(x, y) !is.na(x) & !is.na(y) & abs(round(x, 4) - y) < 1e-9
m$match <- m$in_ref & !is.na(m$power_count) & m$power_count == m$power_count_ref & !is.na(m$type1_count) &
  m$type1_count == m$type1_count_ref & eq4(m$EN_null, m$EN_null_ref) & eq4(m$EN_alt, m$EN_alt_ref)
m$pass <- m$match & m$engine_ok & m$checks_ok & m$mode %in% "dev" & !m$force_r & !dup & m$R %in% 10000L & m$engine_sha %in% ENGINE
m <- m[order(m$design, m$seed), ]
utils::write.csv(m, a[3], row.names = FALSE)
sha <- unique(m$engine_sha)
missing <- if (SUBSET) 0L else nrow(ref) - sum(paste(ref$design, ref$seed) %in% paste(m$design, m$seed))
cat(sprintf("[layer2] %d outputs, %d in reference, %d match, %d pass; engine_ok %d, checks_ok %d, dev %d, force_r %d, duplicate %d, missing %d; engine %s; R %s\n",
            nrow(m), sum(m$in_ref), sum(m$match), sum(m$pass), sum(m$engine_ok), sum(m$checks_ok), sum(m$mode %in% "dev"),
            sum(m$force_r), sum(dup), missing, paste(sha, collapse = "/"), paste(unique(m$r_version), collapse = "/")))
if (any(!m$pass)) print(m[!m$pass, c("design", "seed", "mode", "force_r", "engine_ok", "checks_ok", "power_count", "power_count_ref",
                                     "type1_count", "type1_count_ref", "EN_null", "EN_null_ref", "EN_alt", "EN_alt_ref")], row.names = FALSE)
ok <- all(m$pass) && identical(sha, ENGINE) && missing == 0L && (SUBSET || nrow(m) == nrow(ref))
cat(if (ok) "[layer2] PASS" else "[layer2] FAIL", if (SUBSET) "(subset)" else sprintf("(%d of %d)", sum(m$pass), nrow(ref)), "\n")
quit(save = "no", status = if (ok) 0L else 1L)
