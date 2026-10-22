#!/usr/bin/env Rscript
# run_activedriverwgs.R — cohort-scale ActiveDriverWGS driver.
#
# Why sharding: ActiveDriverWGS parallelises its per-element tests with
# mclapply(). At cohort scale (hundreds of tumours, ~5e5 elements) a small
# number of individual element tests throw inside a worker; mclapply then
# returns error objects mixed into the result list, and the package's final
# aggregation dies with "invalid argument type" — poisoning the ENTIRE run
# for a handful of bad elements. We instead split the element catalogue into
# shards and run ActiveDriverWGS once per shard inside tryCatch(). A failing
# element loses only its shard, not the whole cohort; surviving shards are
# concatenated. Recurrence across donors is unaffected (all mutations are
# passed to every shard; only the tested element set is partitioned).
#
# Usage:
#   run_activedriverwgs.R --mutations <mutations.tsv> --elements <elements.bed> --out <out.tsv>
#                         [--ref-genome hg38] [--cores N] [--shard-size N] [--rlib DIR]
#
# mutations.tsv : chr,pos1,pos2,ref,alt,patient (header)
# elements.bed  : BED (>=4 cols: chr,start,end,id), no header
args <- commandArgs(trailingOnly = TRUE)
usage <- paste("usage: run_activedriverwgs.R --mutations F --elements F --out F",
               "[--ref-genome hg38] [--cores N] [--shard-size N] [--rlib DIR]")
if (length(args) == 0L || length(args) %% 2L != 0L || !all(startsWith(args[c(TRUE, FALSE)], "--"))) stop(usage)
opt <- setNames(as.list(args[c(FALSE, TRUE)]), sub("^--", "", args[c(TRUE, FALSE)]))
known <- c("mutations", "elements", "out", "ref-genome", "cores", "shard-size", "rlib")
if (length(setdiff(names(opt), known))) stop("unknown option(s): ", paste(setdiff(names(opt), known), collapse = ", "), "\n", usage)
if (is.null(opt$mutations) || is.null(opt$elements) || is.null(opt$out)) stop(usage)
mut_path  <- opt$mutations
el_path   <- opt$elements
ref_gen   <- if (is.null(opt[["ref-genome"]])) "hg38" else opt[["ref-genome"]]
out_path  <- opt$out
mc_cores  <- max(1L, as.integer(if (is.null(opt$cores)) 1L else opt$cores))
shard_sz  <- max(1L, as.integer(if (is.null(opt[["shard-size"]])) 5000L else opt[["shard-size"]]))
rlib      <- if (is.null(opt$rlib)) "" else opt$rlib

if (nzchar(rlib)) .libPaths(c(rlib, .libPaths()))
suppressMessages(library(ActiveDriverWGS))

muts <- read.delim(mut_path, stringsAsFactors = FALSE)
# ActiveDriverWGS rejects exact-duplicate (chr,pos1,pos2,ref,alt,patient) rows;
# cohort MAFs contain a few (adjacent indel representations collapsing to one
# coordinate for the same donor). Cross-donor recurrence is preserved.
n0 <- nrow(muts); muts <- muts[!duplicated(muts), ]
message(sprintf("ADWGS: dropped %d exact-duplicate mutation rows (%d -> %d)", n0 - nrow(muts), n0, nrow(muts)))

# THE ROOT CAUSE OF THE SILENT 30% COVERAGE LOSS (audit finding B2).
#
# ActiveDriverWGS validates the ENTIRE mutation set on every call and aborts with
#   "Only autosomal and sex chromosomes may be used in mutation data
#    (24 for human ...). Note that chr23 and chr24 should be formatted as
#    chrX and chrY"
# if ANY mutation sits on a non-primary contig. Because the check is global, the
# error fires no matter which elements are passed -- so an entire 20,000-element
# shard "failed" and was discarded, for a reason that had nothing to do with
# those elements. In the published LIRI-JP/PACA-CA runs that silently removed
# 8 and 9 shards (160,000 / 180,000 elements; 29% and 33% of the catalogue).
#
# The offending rows are a rounding error in the data: 676 of 3,262,697 LIRI-JP
# SNV rows (0.021%) on alt/random/chrUn contigs (chr7_KI270803v1_alt,
# chr4_GL000008v2_random, chrUn_KI270742v1, ...). The element catalogue contains
# ONLY the 24 primary contigs, so these mutations cannot contribute to any
# element test -- dropping them loses no testable signal, while keeping them
# costs a third of the genome.
PRIMARY <- c(paste0("chr", 1:22), "chrX", "chrY")
n1 <- nrow(muts)
off <- !(muts$chr %in% PRIMARY)
if (any(off)) {
  bad_contigs <- sort(unique(muts$chr[off]))
  muts <- muts[!off, ]
  message(sprintf("ADWGS: dropped %d mutation rows (%.3f%%) on %d non-primary contigs: %s",
                  n1 - nrow(muts), 100 * (n1 - nrow(muts)) / n1, length(bad_contigs),
                  paste(utils::head(bad_contigs, 8), collapse = ", ")))
}
if (nrow(muts) == 0L) stop("ADWGS: no mutations remain on primary contigs")

elements <- read.delim(el_path, stringsAsFactors = FALSE, header = FALSE)[, 1:4]
colnames(elements) <- c("chr", "start", "end", "id")
# drop degenerate elements (zero/negative width, NA coords, duplicate ids)
keep <- with(elements, !is.na(start) & !is.na(end) & end > start)
elements <- elements[keep, ]
# DO NOT de-duplicate on id. After CDS subtraction an element is legitimately
# represented by SEVERAL rows (its non-coding fragments), all sharing one id:
# the real catalogue has 595,534 rows for 544,797 ids. ActiveDriverWGS supports
# this natively -- verified empirically: an element supplied as two rows is
# reported ONCE with its fragments pooled ("Tests to do" counts unique ids).
# Dropping duplicate ids here would silently discard 50,737 fragments and
# shrink the tested region, reintroducing the coverage loss this script exists
# to prevent. Only exact duplicate ROWS are dropped.
elements <- elements[!duplicated(elements[, c("chr", "start", "end", "id")]), ]
message(sprintf("ADWGS: %d elements after geometry/id cleaning", nrow(elements)))

# Pre-filter elements whose background window contains NO mutations at all.
# ActiveDriverWGS fits each element against its +/-window_size neighbourhood, so
# an empty window cannot yield a fit and the call errors instead of returning a
# result. Detecting this here costs one interval sweep; discovering it via retry
# costs a full format_muts() pass over the whole cohort per element.
# Measured: 0.42% of the catalogue (2,297 of 544,797 ids) genome-wide, but they
# account for a disproportionate share of failures (27% of those characterised).
# These ids are recorded as untested rather than silently dropped.
failed_elem_ids <- character(0)   # ids we could not test, accumulated below
n_ids_supplied  <- length(unique(elements$id))  # denominator for coverage

ADWGS_WINDOW <- 50000L   # must match ActiveDriverWGS' default window_size
if (nrow(muts) > 0L) {
  mut_by_chr <- split(muts$pos1, muts$chr)
  has_bg <- vapply(seq_len(nrow(elements)), function(i) {
    v <- mut_by_chr[[elements$chr[i]]]
    if (is.null(v)) return(FALSE)
    any(v >= (elements$start[i] - ADWGS_WINDOW) & v <= (elements$end[i] + ADWGS_WINDOW))
  }, logical(1))
  # keep an element if ANY of its fragments has a non-empty window
  ok_ids  <- unique(elements$id[has_bg])
  bad_ids <- setdiff(unique(elements$id), ok_ids)
  if (length(bad_ids)) {
    failed_elem_ids <- c(failed_elem_ids, bad_ids)
    message(sprintf("ADWGS: pre-filtered %d element ids with no mutation within +/-%d bp (untestable background)",
                    length(bad_ids), ADWGS_WINDOW))
    elements <- elements[elements$id %in% ok_ids, , drop = FALSE]
  }
}

n_el   <- nrow(elements)
# Shard by element ID, never by row. A fragmented element (post-CDS-subtraction)
# spans several rows sharing an id; splitting on row index would put its
# fragments in different shards, so each shard would test a PARTIAL element with
# an understated length and the element would be reported more than once.
# Grouping ids keeps every fragment of an element in exactly one shard.
uid        <- unique(elements$id)
grp        <- ceiling(seq_along(uid) / shard_sz)
id_to_grp  <- setNames(grp, uid)
shards     <- split(seq_len(n_el), id_to_grp[elements$id])
message(sprintf("ADWGS: %d shards of <=%d element ids (%d rows / %d ids), mc.cores=%d",
                length(shards), shard_sz, n_el, length(uid), mc_cores))

results <- list()
# NOTE: failed_elem_ids is initialised BEFORE the pre-filter above (which already
# appends to it). Re-initialising here would erase the pre-filtered ids.

# Run one element set, and on failure BISECT and retry each half, recursing down
# to single elements. Rationale: previously a whole 20,000-element shard was
# discarded when any element inside it made ActiveDriverWGS' mclapply aggregation
# throw ("invalid argument type", "numbers of columns of arguments do not match").
# In the published cohort runs that silently dropped 8 shards (160,000 elements,
# 29% of the catalogue) for LIRI-JP and 9 shards (180,000, 33%) for PACA-CA.
# Because an untested element can never be flagged, it is STRUCTURALLY excluded
# from any multi-caller vote -- an unreported coverage bias. Bisection isolates
# the genuinely untestable elements so only they are lost, and logs their ids.
# Recovery strategy: bisect only DOWN TO A SMALL BATCH, then test that batch
# element-by-element in ONE pass.
#
# Why not bisect all the way to singletons: every ActiveDriverWGS() call re-runs
# format_muts() over the FULL cohort mutation set (3.26M rows for LIRI-JP), so
# each retry costs a complete setup pass regardless of how few elements it tests.
# Halving a 20,000-element shard down to individuals is ~15 levels deep and
# hundreds of failures wide; measured on the real cohort that consumed a 48h
# budget without completing even shard 1 of 28 (763 bisection events, depth 14,
# exit 140 = time kill). Bisection is correct but its constant factor is fatal.
#
# Instead: bisect while the batch is large (cheap, few levels), and once a
# failing batch is <= min_batch, run its elements individually in a single loop
# so the cost is linear in the number of BAD elements rather than exponential in
# tree depth. Genuinely untestable elements are still isolated and logged by id.
#
# The underlying cause of these failures, for the record: ActiveDriverWGS fits a
# regression against each element's +/-window_size background, and it errors
# ("object 'fit' not found", "numbers of columns of arguments do not match") when
# that background is too sparse to fit. 27% of the failures we characterised had
# ZERO mutations in their +/-50kb window (vs 0.45% of random elements, a 60x
# enrichment); the remainder merely had sparser-than-average windows (median 72
# vs 94 mutations), which does not by itself explain their failure. NOTE: those
# characterised elements all came from the one shard the run reached before it
# was killed, so they are not a genome-wide sample and no chromosomal pattern
# should be inferred from them. Elements with an empty window are pre-filtered
# below; the rest are caught by the retry path.
# Measured cost on the real LIRI-JP cohort (3.26M mutations, 16 cores), AFTER the
# contig/duplicate cleaning above:
#     100 element ids -> 68.0 s      1,000 ids -> 107.9 s
# Linear fit through those two points: ~0.044 s per element plus a FIXED
# ~63.6 s per call (ActiveDriverWGS re-processes the whole mutation set on every
# invocation). Note the fixed term is the intercept 68 - 100*0.044, NOT the raw
# 108-68 = 40 s delta. Two consequences:
#   * The per-call overhead dominates at small shards: 544,797 ids costs ~16.3 h
#     at shard=1000, ~11.5 h at shard=2000, ~8.6 h at shard=5000. shard 1000-2000
#     trades some wall-clock for cheap failure localisation; the old 20,000
#     amortised overhead best but made each failure enormously expensive to chase.
#   * A failing shard must NOT be bisected to singletons -- at ~64 s per call
#     that is hours per bad element. Hence bisect-then-batch below.
# A 5,000-id shard still hit "invalid argument type" (a genuine per-element
# fault), so recovery logic remains necessary -- just cheap now.
MIN_BATCH <- 32L

run_set <- function(el_sub, depth = 0L) {
  r <- tryCatch(
    ActiveDriverWGS(mutations = muts, elements = el_sub,
                    ref_genome = ref_gen, mc.cores = mc_cores),
    error = function(e) structure(list(msg = conditionMessage(e)), class = "adwgs_err"))
  if (!inherits(r, "adwgs_err")) return(list(r))

  ids <- unique(el_sub$id)
  n_id <- length(ids)

  # Small failing batch: test each element on its own, once.
  if (n_id <= MIN_BATCH) {
    out <- list()
    for (this_id in ids) {
      one <- el_sub[el_sub$id == this_id, , drop = FALSE]
      r1 <- tryCatch(
        ActiveDriverWGS(mutations = muts, elements = one,
                        ref_genome = ref_gen, mc.cores = mc_cores),
        error = function(e) structure(list(msg = conditionMessage(e)), class = "adwgs_err"))
      if (inherits(r1, "adwgs_err")) {
        failed_elem_ids <<- c(failed_elem_ids, this_id)
        message(sprintf("ADWGS: element %s UNTESTABLE: %s", this_id, r1$msg))
      } else {
        out[[length(out) + 1L]] <- r1
      }
    }
    return(out)
  }

  # Large failing batch: split by ID (never mid-element) and recurse.
  mid_id  <- ids[seq_len(n_id %/% 2L)]
  message(sprintf("ADWGS: bisecting %d element ids at depth %d (%s)", n_id, depth, r$msg))
  c(run_set(el_sub[el_sub$id %in% mid_id, , drop = FALSE], depth + 1L),
    run_set(el_sub[!el_sub$id %in% mid_id, , drop = FALSE], depth + 1L))
}

for (si in seq_along(shards)) {
  el_sub <- elements[shards[[si]], , drop = FALSE]
  got <- run_set(el_sub)
  if (length(got)) results <- c(results, got)
  if (si %% 20 == 0) message(sprintf("ADWGS: completed shard %d/%d", si, length(shards)))
}

if (length(results) == 0L) stop("ADWGS: every element failed")
res <- do.call(rbind, results)
# re-derive genome-wide FDR across the pooled per-element p-values (shards were
# FDR-corrected within-shard by the package; recompute over the full set)
if ("pp_element" %in% colnames(res)) res$fdr_element <- p.adjust(res$pp_element, method = "fdr")
if ("pp_site"    %in% colnames(res)) res$fdr_site    <- p.adjust(res$pp_site,    method = "fdr")
res <- res[order(res$pp_element), ]
write.table(res, file = out_path, sep = "\t", quote = FALSE, row.names = FALSE)

# Coverage must be explicit: a silently-missing element is indistinguishable
# downstream from an element that was tested and found unremarkable, which is
# what made the earlier 29-33% loss invisible. Write the untestable ids next to
# the results so ANCHOR_MATCH/CONSENSUS can mark them as not-tested.
# Coverage is per ELEMENT (unique id), not per row: a fragmented element is many
# rows in but one test out, so comparing row counts would understate coverage.
# Denominator must be the ids SUPPLIED, not the post-pre-filter set: using the
# filtered set would report ~100% coverage while silently having dropped the
# pre-filtered elements -- the exact class of overstatement this reporting exists
# to prevent.
n_in  <- n_ids_supplied
n_out <- length(unique(res$id))
cov_path <- sub("\\.tsv$", ".untested.txt", out_path)
writeLines(failed_elem_ids, cov_path)
message(sprintf("ADWGS: wrote %d/%d element results (%.2f%% coverage); %d untestable ids -> %s",
                n_out, n_in, 100 * n_out / n_in, length(failed_elem_ids), basename(cov_path)))
if (n_out < n_in) {
  message(sprintf("ADWGS: WARNING %d elements are absent from the result table", n_in - n_out))
}
