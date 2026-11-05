#!/usr/bin/env Rscript
# Root an ML tree ONCE into a canonical rooted Newick.
#
# phylo_mm.R and the plot both read this tree and never re-root, so clades match.
#
# ROOTING = ANCESTRAL PROXY. We root on the haplotype(s) closest to the MSA
# consensus (the modal allele ≈ ancestral in young intraspecific variation).
# This polarizes the tree ancestral→derived without midpoint's molecular-clock
# assumption, so the reference allele lands at the base on every peak. When the
# min-distance set isn't monophyletic we fall back to midpoint, and if even
# that fails to the first tip. Rooted with ape (not Bio.Phylo) because
# phylo_mm's downstream chronos()/inverseA() need an ape-friendly rooting; the
# plot's Bio.Phylo reads the same rooted Newick.
#
# Usage: Rscript root_tree.R <in.treefile> <alignment.fa> <out.rooted.nwk>
suppressPackageStartupMessages(library(ape))

args <- commandArgs(trailingOnly = TRUE)
treefile <- args[1]
alnfile  <- args[2]
outfile  <- args[3]

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

# Minimal FASTA reader → named character matrix (one row per leaf, one col per
# alignment position). Avoids ape::read.dna quirks on mixed-case / gap chars.
read_aln_matrix <- function(path) {
    lines <- readLines(path)
    nm <- NULL; seqs <- list(); buf <- character(0)
    for (ln in lines) {
        if (startsWith(ln, ">")) {
            if (!is.null(nm)) seqs[[nm]] <- paste(buf, collapse = "")
            nm <- sub("^>", "", strsplit(ln, "\\s+")[[1]][1])
            buf <- character(0)
        } else {
            buf <- c(buf, ln)
        }
    }
    if (!is.null(nm)) seqs[[nm]] <- paste(buf, collapse = "")
    do.call(rbind, lapply(seqs, function(s) strsplit(toupper(s), "")[[1]]))
}

# ── pick the ancestral-proxy outgroup = leaves closest to the MSA consensus ──
ancestral_proxy_tips <- function(tr, mat) {
    strip <- function(x) sub("^_R_", "", x)
    rownames(mat) <- strip(rownames(mat))
    tip_strip <- strip(tr$tip.label)
    keep <- rownames(mat) %in% tip_strip
    mat <- mat[keep, , drop = FALSE]
    if (nrow(mat) < 3L) return(character(0))
    bases <- c("A", "C", "G", "T")
    cons <- apply(mat, 2, function(col) {
        col <- col[col %in% bases]
        if (!length(col)) return("-")
        names(sort(table(col), decreasing = TRUE))[1]
    })
    # Hamming distance to consensus over scored (A/C/G/T) positions only.
    dist <- apply(mat, 1, function(row) {
        sum(row %in% bases & cons %in% bases & row != cons)
    })
    proxy_strip <- rownames(mat)[dist == min(dist)]
    # Map back to the tree's actual (possibly _R_-prefixed) tip labels.
    tr$tip.label[tip_strip %in% proxy_strip]
}

tr  <- read.tree(treefile)
mat <- tryCatch(read_aln_matrix(alnfile), error = function(e) NULL)

rooted <- NULL
if (!is.null(mat) && length(tr$tip.label) >= 3L) {
    proxy <- ancestral_proxy_tips(tr, mat)
    # Need a proper outgroup: non-empty and not all tips. root() errors if the
    # set isn't monophyletic, so try it and fall through on failure.
    if (length(proxy) >= 1L && length(proxy) < length(tr$tip.label)) {
        rooted <- tryCatch(
            root(tr, outgroup = proxy, resolve.root = TRUE),
            error = function(e) NULL
        )
    }
}
if (is.null(rooted)) {
    rooted <- tryCatch(
        midpoint_root(tr),
        error = function(e) root(tr, outgroup = tr$tip.label[1], resolve.root = TRUE)
    )
}
write.tree(rooted, outfile)
