#!/usr/bin/env Rscript
#
# MOTIFBREAKR adapter: candidate SNVs of significant elements -> motifbreakR 2.24.0.
#
# Package functions do all of the analysis. universalmotif reads the JASPAR PFMs and turns
# counts into probabilities; motifbreakR scores both alleles (filterp = TRUE, p threshold 1e-4,
# method "default", uniform background) and calculatePvalue() gives allele p-values. These are
# the settings of the 10 Oct 2026 pilot (tests/regression/PROVENANCE.md).
#
# Pseudocount: the pilot added 0.8 to every base count. universalmotif's convert_type() shares
# its `pseudocount` across the four bases, so --pseudocount (per base, default 0.8) is passed to
# it as 4 x 0.8 = 3.2; this reproduces the pilot matrices exactly (max difference 0 over the
# 1,019 JASPAR 2026 CORE vertebrate PFMs).
#
# ANNOTATION ONLY. Nothing in the pipeline filters or ranks an element or a mutation on these
# values. p = 0 means the match is below calculatePvalue()'s resolution (short motifs); it is
# flagged in p_below_resolution and never floored.
#
# Outputs, for --prefix P:
#   P.motifbreakr.tsv         motifbreakR's result table as returned (list columns joined by commas)
#   P.motifbreakr.per_tf.tsv  one row per SNV x TF: the TF's JASPAR matrix with the smallest
#                             best-allele p (TF names compared case-insensitively: JASPAR
#                             vertebrates mixes human and mouse symbols)

usage <- paste("usage: run_motifbreakr.R --candidates <tsv> --motifs <jaspar.pfm> --prefix <name>",
               "[--threshold 1e-4] [--pseudocount 0.8] [--cores 1] [--versions <tsv>]")
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) { cat(usage, "\n", file = stderr()); quit(status = 2) }
opt <- function(name, default = NULL) {
  i <- match(paste0("--", name), args)
  if (is.na(i)) {
    if (is.null(default)) { cat(usage, "\n", file = stderr()); stop("missing --", name, call. = FALSE) }
    return(default)
  }
  args[i + 1]
}
candidates  <- opt("candidates")
motif_file  <- opt("motifs")
prefix      <- opt("prefix")
threshold   <- as.numeric(opt("threshold", "1e-4"))
pseudocount <- as.numeric(opt("pseudocount", "0.8"))
cores       <- as.integer(opt("cores", "1"))
versions    <- opt("versions", "motifbreakr.versions.tsv")

suppressPackageStartupMessages({
  library(motifbreakR); library(universalmotif); library(MotifDb); library(S4Vectors)
  library(BiocParallel); library(GenomicRanges); library(Biostrings)
  library(BSgenome.Hsapiens.UCSC.hg38)
})
genome <- BSgenome.Hsapiens.UCSC.hg38
BKG <- c(A = 0.25, C = 0.25, G = 0.25, T = 0.25)
bp  <- if (cores > 1) MulticoreParam(cores) else SerialParam()

writeLines(c(sprintf("motifbreakR\t%s", packageVersion("motifbreakR")),
             sprintf("universalmotif\t%s", packageVersion("universalmotif")),
             sprintf("BSgenome.Hsapiens.UCSC.hg38\t%s", packageVersion("BSgenome.Hsapiens.UCSC.hg38")),
             sprintf("r-base\t%s.%s", R.version$major, R.version$minor)), versions)

write_tsv <- function(df, path) write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
RAW_COLS <- c("seqnames", "start", "end", "width", "strand", "SNP_id", "REF", "ALT", "varType",
              "geneSymbol", "dataSource", "providerName", "providerId", "seqMatch", "pctRef",
              "pctAlt", "scoreRef", "scoreAlt", "Refpvalue", "Altpvalue", "altPos", "alleleDiff",
              "alleleEffectSize", "effect", "pvalueEffect")
TF_COLS <- c("chrom", "pos", "ref", "alt", "snp_id", "tf", "motif_id", "motif_length", "ref_pvalue",
             "alt_pvalue", "best_pvalue", "direction", "fold_change", "p_below_resolution",
             "n_motifs_for_tf", "motifbreakr_effect")
empty <- function(cols) as.data.frame(setNames(replicate(length(cols), character(), simplify = FALSE), cols))
out_raw <- paste0(prefix, ".motifbreakr.tsv"); out_tf <- paste0(prefix, ".motifbreakr.per_tf.tsv")

# ---- distinct candidate SNVs ------------------------------------------------------------
cand <- read.delim(candidates, colClasses = "character", check.names = FALSE)
need <- c("chrom", "pos", "ref", "alt", "vtype")
miss <- setdiff(need, colnames(cand))
if (length(miss)) stop(candidates, " lacks column(s): ", paste(miss, collapse = ", "), call. = FALSE)
var <- unique(cand[, need, drop = FALSE])
var <- var[order(var$chrom, as.integer(var$pos), var$ref, var$alt), , drop = FALSE]
var$id <- paste(var$chrom, var$pos, var$ref, var$alt, sep = ":")
is_snv <- toupper(var$vtype) == "SNP" & nchar(var$ref) == 1 & nchar(var$alt) == 1 &
  toupper(var$ref) %in% c("A", "C", "G", "T") & toupper(var$alt) %in% c("A", "C", "G", "T")
on_genome <- var$chrom %in% seqnames(genome)
snv <- var[is_snv & on_genome, , drop = FALSE]
cat(sprintf("distinct candidate mutations: %d; SNVs scanned: %d; not SNVs (not assessed): %d; contig not in genome: %d\n",
            nrow(var), nrow(snv), sum(!is_snv), sum(is_snv & !on_genome)))
if (nrow(snv)) {
  got <- as.character(getSeq(genome, GRanges(snv$chrom, IRanges(as.integer(snv$pos), width = 1))))
  bad <- which(toupper(got) != toupper(snv$ref))
  if (length(bad)) stop(sprintf("%d of %d SNVs have a REF that BSgenome.Hsapiens.UCSC.hg38 does not carry (first: %s, genome %s); the input is not chr-named GRCh38",
                                length(bad), nrow(snv), snv$id[bad[1]], got[bad[1]]), call. = FALSE)
}

raw <- NULL
if (nrow(snv)) {
  # ---- motifs: JASPAR PFMs -> PPMs (universalmotif) -> MotifDb MotifList ------------------
  # read_jaspar() returns a bare motif, not a list, for a file with one motif
  um  <- read_jaspar(motif_file)
  if (!is.list(um)) um <- list(um)
  um  <- convert_type(um, "PPM", pseudocount = 4 * pseudocount)
  if (!is.list(um)) um <- list(um)
  ids <- vapply(um, function(m) m@name, "")
  tfs <- vapply(um, function(m) m@altname, "")
  mlen <- setNames(vapply(um, function(m) ncol(m@motif), 0L), ids)
  ml  <- convert_motifs(um, "MotifDb-MotifList")
  names(ml) <- paste0("Hsapiens-JASPAR-", tfs, "-", ids)
  mcols(ml)$providerId <- ids; mcols(ml)$providerName <- tfs; mcols(ml)$geneSymbol <- tfs
  mcols(ml)$dataSource <- "JASPAR"; mcols(ml)$organism <- "Hsapiens"
  cat(sprintf("motifs: %d (%s); pseudocount %.2f per base\n", length(ml), basename(motif_file), pseudocount))

  bed <- tempfile(fileext = ".bed")
  write.table(data.frame(snv$chrom, as.integer(snv$pos) - 1L, as.integer(snv$pos), snv$id, 0, "+"),
              bed, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  snps <- snps.from.file(bed, search.genome = genome, format = "bed")
  lost <- setdiff(snv$id, names(snps))
  if (length(lost)) warning(length(lost), " SNV(s) not read by motifbreakR: ", paste(head(lost, 5), collapse = ", "), call. = FALSE)

  res <- motifbreakR(snpList = snps, filterp = TRUE, pwmList = ml, threshold = threshold,
                     method = "default", bkg = BKG, BPPARAM = bp)
  if (length(res)) {
    res <- calculatePvalue(res, background = BKG, BPPARAM = bp)
    raw <- as.data.frame(res, row.names = NULL)
    for (j in which(vapply(raw, is.list, logical(1))))
      raw[[j]] <- vapply(raw[[j]], function(x) paste(as.character(unlist(x)), collapse = ","), "")
  }
}

if (is.null(raw) || !nrow(raw)) {
  write_tsv(empty(RAW_COLS), out_raw); write_tsv(empty(TF_COLS), out_tf)
  cat("motifbreakR effects: 0\n"); quit(status = 0)
}
write_tsv(raw, out_raw)

# ---- one row per SNV x TF ---------------------------------------------------------------
pr <- as.numeric(raw$Refpvalue); pa <- as.numeric(raw$Altpvalue)
d <- data.frame(snp_id = raw$SNP_id, tf = raw$geneSymbol, motif_id = raw$providerId,
                ref_pvalue = pr, alt_pvalue = pa, best_pvalue = suppressWarnings(pmin(pr, pa, na.rm = TRUE)),
                motifbreakr_effect = raw$effect, stringsAsFactors = FALSE)
d$tf_key <- toupper(d$tf)
d <- d[order(d$snp_id, d$tf_key, d$best_pvalue, d$motif_id), , drop = FALSE]
grp <- paste(d$snp_id, d$tf_key, sep = "\r")
n_mot <- tapply(d$motif_id, grp, function(x) length(unique(x)))
d <- d[!duplicated(grp), , drop = FALSE]
d$n_motifs_for_tf <- as.integer(n_mot[paste(d$snp_id, d$tf_key, sep = "\r")])
d$p_below_resolution <- (!is.na(d$ref_pvalue) & d$ref_pvalue == 0) | (!is.na(d$alt_pvalue) & d$alt_pvalue == 0)
d$direction <- ifelse(is.na(d$ref_pvalue) | is.na(d$alt_pvalue), NA,
                      ifelse(d$alt_pvalue < d$ref_pvalue, "gain", ifelse(d$ref_pvalue < d$alt_pvalue, "loss", "none")))
d$fold_change <- ifelse(d$p_below_resolution | is.na(d$direction), NA,
                        pmax(d$ref_pvalue, d$alt_pvalue) / pmin(d$ref_pvalue, d$alt_pvalue))
d$motif_length <- as.integer(mlen[d$motif_id])
v <- var[match(d$snp_id, var$id), ]
d$chrom <- v$chrom; d$pos <- v$pos; d$ref <- v$ref; d$alt <- v$alt
d <- d[order(d$chrom, as.integer(d$pos), d$ref, d$alt, d$best_pvalue, d$tf_key), TF_COLS, drop = FALSE]
write_tsv(d, out_tf)
cat(sprintf("motifbreakR effects: %d rows; per SNV x TF: %d rows at %d SNVs; p below resolution: %d\n",
            nrow(raw), nrow(d), length(unique(d$snp_id)), sum(d$p_below_resolution)))
