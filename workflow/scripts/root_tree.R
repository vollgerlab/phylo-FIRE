#!/usr/bin/env Rscript
# Midpoint-root an ML tree ONCE into a canonical rooted Newick.
#
# Both phylo_mm.R and the plot read this rooted tree and neither re-roots, so
# the clade each branch represents is identical everywhere — fixing the
# scattered leader lines / detached annotations / mirror branches that came
# from phylo_mm and the plot rooting independently. Midpoint rooting is the
# justifiable default for these no-outgroup intraspecific haplotype trees
# (it assumes roughly clock-like divergence over a short locus). Rooted with
# ape (not Bio.Phylo) because phylo_mm's downstream chronos()/inverseA() need
# an ape-friendly rooting; the plot's Bio.Phylo reads the same rooted Newick.
#
# Usage: Rscript root_tree.R <in.treefile> <out.rooted.nwk>
suppressPackageStartupMessages(library(ape))

args <- commandArgs(trailingOnly = TRUE)

midpoint_root <- function(tr) {
    if (length(tr$tip.label) < 3L) return(tr)
    D  <- cophenetic(tr)
    mx <- which(D == max(D), arr.ind = TRUE)[1, ]
    half <- D[mx[1], mx[2]] / 2
    np <- ape::nodepath(tr, mx[1], mx[2])
    el <- setNames(tr$edge.length, paste(tr$edge[, 1], tr$edge[, 2]))
    cum <- 0
    for (i in seq_len(length(np) - 1L)) {
        a <- np[i]; b <- np[i + 1L]
        e <- el[paste(a, b)]; if (is.na(e)) e <- el[paste(b, a)]
        if (cum + e >= half) return(root(tr, node = b, resolve.root = TRUE))
        cum <- cum + e
    }
    tr
}

tr <- read.tree(args[1])
if (!is.rooted(tr)) {
    tr <- tryCatch(
        midpoint_root(tr),
        error = function(e) root(tr, outgroup = tr$tip.label[1], resolve.root = TRUE)
    )
}
write.tree(tr, args[2])
