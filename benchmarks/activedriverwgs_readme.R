#!/usr/bin/env Rscript
#
# Reproduce the ActiveDriverWGS README example (ATM, MYD88 and NOTCH1 in the package's CLL data) through
# bin/run_activedriverwgs.R, the script of the ACTIVEDRIVERWGS module, and compare it with the values the
# package README prints (ActiveDriverWGS 1.2.1, hg19, package defaults).
#
# Manual check, not run by CI: the example is hg19, so it needs BSgenome.Hsapiens.UCSC.hg19 (about 710 MB),
# which the pipeline image does not carry. Install it into a scratch library and run inside the image:
#
#   R CMD INSTALL -l $RLIB BSgenome.Hsapiens.UCSC.hg19_1.4.3.tar.gz
#   Rscript benchmarks/activedriverwgs_readme.R --rlib $RLIB --wrapper bin/run_activedriverwgs.R --outdir out
#
# Checks: the direct ActiveDriverWGS() call and the wrapper at shard sizes 5000 and 1 match the README
# (relative 1e-6: the README prints 7-8 significant digits), and the wrapper tables equal the direct call.
# Also reports, without a pass criterion, the result on SNVs only (what the module passes to ActiveDriverWGS).
# Exit status 1 if a check fails.
#
args <- commandArgs(trailingOnly = TRUE)
opt  <- setNames(as.list(args[c(FALSE, TRUE)]), sub("^--", "", args[c(TRUE, FALSE)]))
if (is.null(opt$rlib)) stop("usage: activedriverwgs_readme.R --rlib DIR [--wrapper FILE] [--outdir DIR]")
wrapper <- if (is.null(opt$wrapper)) "bin/run_activedriverwgs.R" else opt$wrapper
outdir  <- if (is.null(opt$outdir)) "activedriverwgs_readme" else opt$outdir
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(opt$rlib, .libPaths()))
suppressMessages(library(ActiveDriverWGS))

# README values (ActiveDriverWGS 1.2.1)
expected <- data.frame(
    id               = c("ATM", "MYD88", "NOTCH1"),
    pp_element       = c(0.002092488, 0.003309630, 0.004964277),
    element_muts_obs = c(2L, 1L, 2L),
    element_muts_exp = c(0.11059736, 0.01337728, 0.16562752),
    fdr_element      = c(0.004964277, 0.004964277, 0.004964277)
)

data("cll_mutations"); data("cancer_genes")
el <- cancer_genes[cancer_genes$id %in% expected$id, ]
f_mut <- file.path(outdir, "readme_mutations.tsv"); f_el <- file.path(outdir, "readme_elements.bed")
f_snv <- file.path(outdir, "readme_mutations_snv.tsv")
write.table(cll_mutations, f_mut, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(el, f_el, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
acgt <- c("A", "C", "G", "T")
write.table(cll_mutations[cll_mutations$ref %in% acgt & cll_mutations$alt %in% acgt, ], f_snv,
            sep = "\t", quote = FALSE, row.names = FALSE)

# the direct call goes through the same write/read as the wrapper output, so the tables compare exactly
f_direct <- file.path(outdir, "direct.tsv")
write.table(ActiveDriverWGS(mutations = cll_mutations, elements = el), f_direct, sep = "\t", quote = FALSE,
            row.names = FALSE)
runs <- list(direct = read.delim(f_direct, stringsAsFactors = FALSE))
for (s in c(5000L, 1L)) {
    out <- file.path(outdir, sprintf("wrapper_shard%d.tsv", s))
    rc  <- system2("Rscript", c(wrapper, "--mutations", f_mut, "--elements", f_el, "--out", out,
                                "--ref-genome", "hg19", "--rlib", opt$rlib, "--shard-size", s))
    if (rc != 0) stop("wrapper failed at shard size ", s)
    runs[[sprintf("wrapper_shard%d", s)]] <- read.delim(out, stringsAsFactors = FALSE)
}
keys <- c("pp_element", "element_muts_obs", "element_muts_exp", "fdr_element")
ok <- TRUE
for (nm in names(runs)) {
    r <- runs[[nm]][match(expected$id, runs[[nm]]$id), ]
    rel <- max(sapply(c("pp_element", "element_muts_exp", "fdr_element"),
                      function(k) max(abs(r[[k]] - expected[[k]]) / expected[[k]])))
    pass <- rel <= 1e-6 && all(r$element_muts_obs == expected$element_muts_obs)
    if (nm != "direct") {
        d <- runs$direct[match(expected$id, runs$direct$id), ]
        pass <- pass && all(sapply(keys, function(k) identical(r[[k]], d[[k]])))
    }
    cat(sprintf("%-16s max relative difference from README %.2e  %s\n", nm, rel, if (pass) "PASS" else "FAIL"))
    ok <- ok && pass
}
out <- file.path(outdir, "wrapper_snv_only.tsv")
system2("Rscript", c(wrapper, "--mutations", f_snv, "--elements", f_el, "--out", out,
                     "--ref-genome", "hg19", "--rlib", opt$rlib))
snv <- read.delim(out, stringsAsFactors = FALSE)
cat("SNVs only (no pass criterion):\n"); print(snv[, c("id", keys)], row.names = FALSE)
quit(status = if (ok) 0L else 1L)
