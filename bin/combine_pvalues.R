#!/usr/bin/env Rscript
# COMBINE_PVALUES: combine per-caller element p-values using published packages only.
#   EBM  : EmpiricalBrownsMethod (Poole et al. 2016). Covariance from calculateCovariances() on
#          the callers x elements p matrix of elements tested by every caller (per cohort);
#          each multi-caller element is combined with combinePValues(); one-caller elements
#          pass through unchanged.
#   ACAT : ACAT package (Liu & Xie 2020; github.com/yaowuliu/ACAT @ e737d7d, GPL-3, installed
#          in the image, never vendored here). Inputs with p == 1 are removed before combining
#          (the authors' ACAT-V convention); an element whose inputs are all 1 gets acat_p = 1.
#   RRA  : RobustRankAggreg::rhoScores (Kolde et al. 2012) on per-caller normalised ranks.
#   BH   : stats::p.adjust(method = "BH"), over all elements of the cohort.
#   votes: number of callers with their own q <= --fdr (descriptive only).
# The only in-house logic is argument parsing, column selection and output formatting.
#
# Usage: combine_pvalues.R --input matched.tsv --out consensus.tsv [--method ebm|acat]
#                          [--fdr 0.1] [--skip-acat] [--versions versions.tsv]
# Input: tab-separated, one row per element, an element_id column, one p_<caller> column per
# primary caller and optional q_<caller> columns. Other columns are carried through.

args <- commandArgs(trailingOnly = TRUE)
opt  <- list(input = NULL, out = NULL, method = "ebm", fdr = "0.1", skip_acat = FALSE,
             versions = "combinepvalues.versions.tsv")
i <- 1
while (i <= length(args)) {
    key <- gsub("-", "_", sub("^--", "", args[i]))
    if (key == "skip_acat") { opt[["skip_acat"]] <- TRUE; i <- i + 1; next }
    opt[[key]] <- args[i + 1]
    i <- i + 2
}
fdr <- as.numeric(opt[["fdr"]])
if (is.null(opt[["input"]]) || is.null(opt[["out"]])) stop("--input and --out are required")
if (!opt[["method"]] %in% c("ebm", "acat")) stop("--method must be 'ebm' or 'acat'")
if (opt[["method"]] == "acat" && opt[["skip_acat"]]) stop("--method acat cannot be combined with --skip-acat")

suppressPackageStartupMessages({
    library(EmpiricalBrownsMethod)
    library(RobustRankAggreg)
    if (!opt[["skip_acat"]]) library(ACAT)
})
pkgs <- c("EmpiricalBrownsMethod", "RobustRankAggreg", if (!opt[["skip_acat"]]) "ACAT")
writeLines(c(paste(pkgs, vapply(pkgs, function(p) as.character(packageVersion(p)), ""), sep = "\t"),
             if (opt[["skip_acat"]]) "ACAT\tnot used",
             paste("r-base", paste(R.version[["major"]], R.version[["minor"]], sep = "."), sep = "\t")),
           opt[["versions"]])

d <- read.delim(opt[["input"]], stringsAsFactors = FALSE, check.names = FALSE,
                na.strings = c("NA", ""), colClasses = c(element_id = "character"))
pcols <- grep("^p_", names(d), value = TRUE)
if (length(pcols) == 0) stop("no p_<caller> columns in ", opt[["input"]])
P <- as.matrix(d[, pcols, drop = FALSE])
storage.mode(P) <- "double"
tested <- rowSums(!is.na(P)) > 0
d <- d[tested, , drop = FALSE]
P <- P[tested, , drop = FALSE]
n_tested <- rowSums(!is.na(P))
first_p  <- apply(P, 1, function(v) v[!is.na(v)][1])
complete <- n_tested == ncol(P)
multi    <- n_tested > 1

# --- Empirical Brown ---------------------------------------------------------------------
ebm_p <- first_p
cv <- NULL
if (ncol(P) > 1 && sum(complete) > 2) {
    cv <- EmpiricalBrownsMethod:::calculateCovariances(t(P[complete, , drop = FALSE]))
    ebm_p[multi] <- vapply(which(multi), function(r) {
        k <- which(!is.na(P[r, ]))
        EmpiricalBrownsMethod:::combinePValues(cv[k, k, drop = FALSE], P[r, k])
    }, 0)
}

# --- ACAT (authors' package; p == 1 removed) ---------------------------------------------
acat_p <- rep(NA_real_, nrow(P))
if (!opt[["skip_acat"]]) {
    acat_p <- first_p
    acat_p[multi] <- vapply(which(multi), function(r) {
        v <- P[r, !is.na(P[r, ])]
        v <- v[v < 1]
        if (length(v) == 0) 1 else if (length(v) == 1) v else ACAT(v)
    }, 0)
}

# --- RRA ---------------------------------------------------------------------------------
norm_rank <- function(p) {
    r <- rep(NA_real_, length(p)); k <- !is.na(p)
    r[k] <- rank(p[k], ties.method = "average") / sum(k)
    r
}
R <- apply(P, 2, norm_rank)
if (is.null(dim(R))) R <- matrix(R, ncol = 1)
rra_score <- apply(R, 1, function(r) rhoScores(r[!is.na(r)]))

# --- BH, votes, selection ----------------------------------------------------------------
qcols <- intersect(sub("^p_", "q_", pcols), names(d))
votes <- if (length(qcols)) rowSums(as.matrix(d[, qcols, drop = FALSE]) <= fdr, na.rm = TRUE) else rep(0L, nrow(d))
out <- d
out[["n_tested_callers"]] <- n_tested
out[["ebm_p"]]     <- ebm_p
out[["ebm_q"]]     <- p.adjust(ebm_p, method = "BH")
out[["acat_p"]]    <- acat_p
out[["acat_q"]]    <- if (opt[["skip_acat"]]) NA_real_ else p.adjust(acat_p, method = "BH")
out[["rra_score"]] <- rra_score
out[["rra_q"]]     <- p.adjust(rra_score, method = "BH")
out[["vote_count"]] <- votes
sel_p <- if (opt[["method"]] == "ebm") out[["ebm_p"]] else out[["acat_p"]]
out[["combination_method"]] <- opt[["method"]]
out[["consensus_q"]]    <- if (opt[["method"]] == "ebm") out[["ebm_q"]] else out[["acat_q"]]
out[["consensus_rank"]] <- rank(sel_p, ties.method = "min")
out[["is_consensus"]]   <- ifelse(out[["consensus_q"]] < fdr, "True", "False")
out <- out[order(out[["consensus_rank"]], out[["element_id"]]), , drop = FALSE]
write.table(out, opt[["out"]], sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

# EBM parameters per cohort, for MultiQC and the methods text
if (!is.null(cv)) {
    N <- ncol(P); Var <- 4 * N + 2 * sum(cv[lower.tri(cv)])
    cc <- Var / (4 * N); dfb <- 2 * (2 * N)^2 / Var
    if (dfb > 2 * N) { dfb <- 2 * N; cc <- 1 }
    params <- data.frame(callers = paste(pcols, collapse = ","), elements = nrow(P),
                         complete_elements = sum(complete), cov_offdiag_sum = sum(cv[lower.tri(cv)]),
                         brown_c = cc, brown_df = dfb)
    write.table(params, sub("\\.tsv$", ".ebm_parameters.tsv", opt[["out"]]), sep = "\t",
                quote = FALSE, row.names = FALSE)
}
