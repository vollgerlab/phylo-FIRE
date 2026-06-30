#!/usr/bin/env Rscript
# Per-peak phylogenetic analysis: heritability sensitivity + joint LASSO
# branch selection + per-branch diagnostic scan.
#
# Three outputs are produced per peak, all summarized in a single TSV +
# diagnostic PDF:
#
# 1. **Phylogenetic heritability H²_p — multi-prior sensitivity.**
#    The intercept-only MCMCglmm is fit three times, each with a different
#    parameter-expansion scale on the phylogenetic random effect:
#      - alpha.V = 1000  (very wide, current default)
#      - alpha.V = 1     (Half-Cauchy(0,1), moderate)
#      - alpha.V = 0.25  (tight, says "V_phylo is probably small")
#    All three posteriors are reported. If they agree, H²_p is data-driven
#    and trustworthy; if they disagree, H²_p is prior-dominated for this
#    peak (typical of nearly-clonal short-haplotype trees — see
#    Hadfield 2015 MEE for the underlying identifiability issue).
#
# 2. **Per-split sister-clade contrast — headline inference + annotation.**
#    Each internal node defines a split into two daughter subtrees. For the
#    SMALLER daughter (the split's representative, `is_split`) we run a
#    coverage-weighted WLS contrast of actuation against its SISTER subtree
#    (the other daughter), NOT against the rest of the tree. This is the local
#    phylogenetic contrast — "does this lineage differ from the one it just
#    split from?" — so a clade that looks significant only because the global
#    mean lumps in a distant high/low clade is correctly NOT significant. BH
#    FDR is over the split representatives (one test per node). p_fdr is the
#    headline p; Δ = mean_in − mean_out (clade − sister) is the headline
#    effect. The plot draws splits with p_fdr < 0.05.
#
# 3. **LASSO joint branch selection — secondary effect size.**
#    For each clean bipartition (a split no identical-sequence cluster
#    straddles), cv.glmnet (alpha=1, family=gaussian, lambda.min) +
#    selectiveInference::fixedLassoInf give a joint post-selection β + CI,
#    retained in the TSV (`beta`/`ci_*`/`selected_by_lasso`) and the plot
#    hover as a secondary effect-size view. This NO LONGER drives the
#    comparison or which branches are annotated (the sister contrast in 2
#    does), because the joint β over correlated nested clades is unstable
#    (can sign-flip vs the marginal Δ). The linked-variants column lists MSA
#    positions whose allele pattern matches a bipartition exactly.
#
# Usage:
#   Rscript phylo_mm.R \
#       --tree      tree.treefile \
#       --actuation peak.tsv \
#       --alignment aligned.fa \
#       --tsv       out.tsv \
#       --pdf       out.pdf \
#       --peak-id   peak_id

suppressPackageStartupMessages({
    library(MCMCglmm)
    library(ape)
    library(glmnet)
    library(selectiveInference)   # post-LASSO valid p-values + CIs
})

# ── arg parsing ───────────────────────────────────────────────────────────────
parse_args <- function(args) {
    out <- list()
    i <- 1L
    while (i <= length(args)) {
        out[[sub("^--", "", args[i])]] <- args[i + 1L]
        i <- i + 2L
    }
    out
}
opt <- parse_args(commandArgs(trailingOnly = TRUE))
stopifnot(!is.null(opt$tree), !is.null(opt$actuation),
          !is.null(opt$alignment),
          !is.null(opt$tsv))
# --pdf is optional: the per-peak diagnostic PDF is no longer produced in the
# pipeline (everything it showed is in the TSV header + the main category plot,
# which only renders for the plotted subset). Pass --pdf for ad-hoc debugging.
peak_id <- if (is.null(opt[["peak-id"]])) tools::file_path_sans_ext(basename(opt$tree)) else opt[["peak-id"]]

# Reproducibility: MCMCglmm (MCMC) and cv.glmnet (random CV folds) are both
# stochastic. Seed from the peak_id so each peak is deterministic across reruns
# but peaks don't share an identical fold structure.
set.seed(sum(utf8ToInt(peak_id)) %% .Machine$integer.max)

# Minimum effective sample size for the headline (alpha.V=1) H2_p posterior to
# be called converged. Reported in the TSV header as h2_converged.
ESS_MIN <- 200

# Minimum per-leaf coverage (fibers) for a leaf to enter the per-branch
# clade-vs-rest test. Actuation from very few fibers is a noisy proportion;
# the H2_p model already downweights it via `mev`, but the branch scan needs
# an explicit floor. Leaves are also coverage-WEIGHTED within the test (WLS).
MIN_COV <- as.numeric(if (is.null(opt[["min-cov"]])) 4 else opt[["min-cov"]])

# bp of flank that add_flanks padded each side of the RE (0 = no flank). When
# >0, site-based association (linked variants, LASSO anchoring) is restricted to
# the regulatory-element columns; the flanks only contribute to the tree.
FLANK <- as.integer(if (is.null(opt$flank)) 0 else opt$flank)

# Heritability gate for variant attribution (see config h2_hpd_*). A H²ₚ fit is
# "stable + significant" — i.e. actuation is a heritable function of sequence —
# when its 95% HPD is no wider than H2_WIDTH_MAX (identified) AND its lower bound
# is at least H2_LOWER_MIN (confidently > 0). RE-only H² is tried first; if it
# fails, the whole-region H² is used (and attribution widens to flank variants).
H2_WIDTH_MAX <- as.numeric(if (is.null(opt[["h2-width-max"]])) 0.5 else opt[["h2-width-max"]])
H2_LOWER_MIN <- as.numeric(if (is.null(opt[["h2-lower-min"]])) 0.1 else opt[["h2-lower-min"]])
# RE-only tree + alignment (for the RE-only H²). Default to the whole-region
# tree/alignment when not supplied (e.g. flank=0, where the RE *is* the whole).
RE_TREE_PATH <- if (is.null(opt[["re-tree"]])) opt$tree else opt[["re-tree"]]
RE_ALN_PATH  <- if (is.null(opt[["re-alignment"]])) opt$alignment else opt[["re-alignment"]]

# A variant must be carried by >= MIN_CARRIERS haplotypes on each side to be a
# candidate. A singleton (1 carrier, e.g. a private insertion) trivially
# "explains" that one haplotype's actuation — overfitting, not a driver.
MIN_CARRIERS <- as.integer(if (is.null(opt[["min-carriers"]])) 2 else opt[["min-carriers"]])

NITT   <- as.integer(if (is.null(opt$nitt))   50000 else opt$nitt)
BURNIN <- as.integer(if (is.null(opt$burnin)) 10000 else opt$burnin)
THIN   <- as.integer(if (is.null(opt$thin))      25 else opt$thin)

# ── load tree + actuation TSV ─────────────────────────────────────────────────
# The tree arrives pre-midpoint-rooted (rule root_tree, via Bio.Phylo) so
# phylo_mm and the plot share ONE rooting — the clade each branch represents is
# then identical everywhere. `midpoint_root` is kept only as a fallback if an
# unrooted tree ever slips in.
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
tree <- read.tree(opt$tree)
tree$tip.label <- gsub("^_R_", "", tree$tip.label)
if (!is.rooted(tree)) {
    tree <- tryCatch(
        midpoint_root(tree),
        error = function(e) root(tree, outgroup = tree$tip.label[1], resolve.root = TRUE)
    )
}
zero_edges <- tree$edge.length <= 0
if (any(zero_edges)) tree$edge.length[zero_edges] <- 1e-6
tree$node.label <- NULL

df <- read.table(opt$actuation, header = TRUE, sep = "\t",
                 stringsAsFactors = FALSE, comment.char = "", quote = "")
keep <- !is.na(df$actuation) & df$actuation != "NA"
trait <- suppressWarnings(as.numeric(df$actuation[keep]))
names(trait) <- df$sample[keep]
trait <- trait[!is.na(trait)]

# Coverage (number of fibers) per sample — used to compute logit-scale
# measurement error variance for MCMCglmm (Hadfield 2010 meta-analysis
# trick: per-observation precision in the residual likelihood).
cov_vec <- suppressWarnings(as.numeric(df$coverage[keep]))
names(cov_vec) <- df$sample[keep]
cov_vec <- cov_vec[!is.na(cov_vec) & cov_vec > 0]

common <- intersect(tree$tip.label, names(trait))

empty_branch_cols <- c(
    "branch_id", "n_below", "descendants",
    "mean_in", "mean_out", "t_stat", "p_raw", "p_fdr", "is_split", "linked_sites",
    "selected_by_lasso", "lasso_coef",
    "beta", "se", "ci_lower", "ci_upper", "p"
)
write_empty_header <- function(reason) {
    header <- c(
        paste0("# peak_id=", peak_id),
        paste0("# n_tips=", length(common)),
        paste0("# n_unique_seqs=NA"),
        paste0("# h2_alpha1000_mean=NA"),
        paste0("# h2_alpha1000_lower=NA"),
        paste0("# h2_alpha1000_upper=NA"),
        paste0("# h2_alpha1_mean=NA"),
        paste0("# h2_alpha1_lower=NA"),
        paste0("# h2_alpha1_upper=NA"),
        paste0("# h2_alpha025_mean=NA"),
        paste0("# h2_alpha025_lower=NA"),
        paste0("# h2_alpha025_upper=NA"),
        paste0("# h2_alpha1_ess=NA"),
        paste0("# h2_converged=NA"),
        paste0("# h2_re_mean=NA"),
        paste0("# h2_re_lower=NA"),
        paste0("# h2_re_upper=NA"),
        paste0("# h2_re_ok=NA"),
        paste0("# h2_region_mean=NA"),
        paste0("# h2_region_lower=NA"),
        paste0("# h2_region_upper=NA"),
        paste0("# h2_region_ok=NA"),
        paste0("# attribution_mode=none"),
        paste0("# min_cov=", MIN_COV),
        paste0("# flank=", FLANK),
        paste0("# re_first_col=NA"),
        paste0("# re_last_col=NA"),
        paste0("# n_variants_tested=0"),
        paste0("# n_variants_selected=0"),
        paste0("# n_branches_selected=0"),
        paste0("# inference_method=none"),
        paste0("# note=", reason)
    )
    writeLines(c(header, paste(empty_branch_cols, collapse = "\t")), opt$tsv)
    if (!is.null(opt$pdf)) {
        pdf(opt$pdf, width = 6, height = 4)
        plot.new(); title(main = paste(peak_id, "—", reason))
        dev.off()
    }
}

if (length(common) < 4) {
    write_empty_header(sprintf("only %d tips with actuation data", length(common)))
    quit(status = 0)
}

trait <- trait[common]
tree  <- keep.tip(tree, common)

# ── alignment for linked-variant detection ────────────────────────────────────
read_fasta <- function(path) {
    lines <- readLines(path)
    out <- list(); nm <- NULL; buf <- character(0)
    for (ln in lines) {
        if (startsWith(ln, ">")) {
            if (!is.null(nm)) out[[nm]] <- paste(buf, collapse = "")
            nm <- sub("^_R_", "", sub("^>", "", strsplit(ln, " ")[[1]][1]))
            buf <- character(0)
        } else buf <- c(buf, ln)
    }
    if (!is.null(nm)) out[[nm]] <- paste(buf, collapse = "")
    out
}
aln <- read_fasta(opt$alignment)
aln <- aln[common]; aln <- aln[!sapply(aln, is.null)]
seq_mat <- if (length(aln) >= 4) {
    do.call(rbind, strsplit(toupper(unlist(aln)), ""))
} else {
    NULL
}
if (!is.null(seq_mat)) rownames(seq_mat) <- names(aln)

# Regulatory-element column band [re_lo, re_hi] (1-based alignment columns).
# When the alignment carries flanks, the RE is the interior: skip `FLANK`
# non-gap characters from each end of every row, union the interior spans.
# flank=0 (or no alignment) -> whole alignment is the RE.
re_band <- function(seqs, flank) {
    if (flank <= 0L || length(seqs) == 0L) return(NULL)
    lo <- Inf; hi <- -Inf
    for (s in seqs) {
        chars  <- strsplit(s, "")[[1]]
        nongap <- which(chars != "-")
        if (length(nongap) <= 2L * flank) next
        lo <- min(lo, nongap[flank + 1L])
        hi <- max(hi, nongap[length(nongap) - flank])
    }
    if (is.finite(lo) && hi >= lo) c(lo, hi) else NULL
}
re_cols <- re_band(toupper(unlist(aln)), FLANK)
re_lo <- if (is.null(re_cols)) 1L else re_cols[1]
re_hi <- if (is.null(re_cols)) (if (is.null(seq_mat)) 0L else ncol(seq_mat)) else re_cols[2]

# Variable column: ≥2 distinct characters at this position, treating gap
# and N as their own allele states (so a few-bp deletion that segregates
# cleanly between clades can anchor a branch, not just A/T/G/C substitutions).
# BUT skip gap/N-DOMINATED columns (>50% of sequences missing a real base):
# those are an indel/extraction-window edge, not a callable variant, and they
# produce spurious weak associations anchored on a near-empty column.
# Restricted to the RE band so neutral flank variation never drives the
# site-based association — the flanks resolve the tree, the RE drives the test.
var_cols_all <- if (is.null(seq_mat)) integer(0) else which(apply(seq_mat, 2, function(col) {
    real <- col %in% c("A", "C", "G", "T")
    length(unique(col)) >= 2L && mean(real) >= 0.5
}))
var_cols_re <- var_cols_all[var_cols_all >= re_lo & var_cols_all <= re_hi]
# `var_cols` (the candidates for attribution) is set AFTER the heritability gate
# below: RE-band variants if the RE-only H² is heritable, the whole-alignment
# variants if we expand to the region, none otherwise.
cat(sprintf("[phylo_mm] %s RE band: cols %d-%d (flank=%d); %d variable RE sites (of %d total)\n",
            peak_id, re_lo, re_hi, FLANK, length(var_cols_re), length(var_cols_all)),
    file = stderr())

# ── 1. Heritability — RE-only first, expand to whole region ───────────────────
# H²ₚ is fit TWICE: on the RE-only tree/alignment (does actuation track the RE
# sequence?) and on the whole flanked tree (does it track the broader region?).
# Variant attribution is gated on these: RE-only first, expand to region only if
# the RE fails but the region works, else no drivers. Recipe per fit (Halliwell/
# de Villemereuil/Hadfield): logit-transform the [0,1] trait with a
# coverage-derived offset; pass per-obs measurement-error variance `mev`
# (binomial SE on the logit scale = 1/(n·p·(1-p))); fit on an ultrametric tree
# (chronos). The shared trait prep:
cov_used <- cov_vec[common]
if (any(is.na(cov_used))) cov_used[is.na(cov_used)] <- median(cov_used, na.rm = TRUE)
trait_eps  <- max(1e-4, 0.5 / max(cov_vec[common], na.rm = TRUE))
trait_clip <- pmin(pmax(trait[common], trait_eps), 1 - trait_eps)
y_logit    <- setNames(qlogis(trait_clip), common)
mev_logit  <- setNames(1 / (cov_used * trait_clip * (1 - trait_clip)), common)

summarise_h2 <- function(fit) {
    if (is.null(fit)) return(list(median = NA_real_, lo = NA_real_, hi = NA_real_,
                                  ess = NA_real_, sigma_p = NA_real_, sigma_r = NA_real_))
    s_p <- fit$VCV[, "animal"]
    s_r <- fit$VCV[, "units"]
    h2  <- s_p / (s_p + s_r)
    ci  <- HPDinterval(h2, prob = 0.95)
    list(median = median(h2), lo = ci[1], hi = ci[2],
         ess = unname(effectiveSize(h2)), sigma_p = median(s_p), sigma_r = median(s_r))
}

# Fit the H² sensitivity (3 alphas + tight-R) on a (tree, alignment) pair. The
# alignment defines the identical-sequence CLUSTERS (within-cluster spread → the
# residual σ²_R that, when large, makes H² unidentified); the tree gives the
# phylogenetic covariance. `ok` = the alpha.V=1 headline is identified (HPD width
# <= H2_WIDTH_MAX) AND significant (HPD lower >= H2_LOWER_MIN).
compute_h2 <- function(tree_in, aln_in) {
    # Drop IQ-TREE's bootstrap internal-node labels: they repeat (e.g. two nodes
    # labelled "100"), and inverseA rejects non-unique tip/node labels. (Local
    # copy — the plot's tree keeps its support labels.)
    tree_in$node.label <- NULL
    NA3 <- list(median = NA_real_, lo = NA_real_, hi = NA_real_, ess = NA_real_,
                sigma_p = NA_real_, sigma_r = NA_real_)
    fail <- list(n_unique = NA_integer_, ok = FALSE,
                 a1000 = NA3, a1 = NA3, a025 = NA3, tightR = NA3, fits = list())
    # Restrict to tips present in BOTH this scope's tree and alignment. The
    # RE-scoped alignment/tree can be missing tips the main scope has; without
    # this, aln_in[common] re-introduces NULLs and vapply() errors hard (which
    # would fail the whole peak instead of degrading to ok=FALSE).
    common_h <- intersect(intersect(common, tree_in$tip.label), names(aln_in))
    if (length(common_h) < 4L) return(fail)
    strings <- vapply(aln_in[common_h], function(s) s, character(1))
    cidx <- match(strings, unique(strings))
    n_unique <- max(cidx)
    if (n_unique < 2L) { fail$n_unique <- n_unique; return(fail) }
    reps <- common_h[!duplicated(cidx)]
    tip_to_rep <- reps[cidx]
    # Ultrametrify for inverseA(scale=TRUE); if chronos can't (fragile on
    # near-zero RE branches) fall back to the raw tree with scale=FALSE, and if
    # inverseA still fails treat this scope's H² as unavailable (→ ok=FALSE).
    tree_ultra <- tryCatch(
        suppressWarnings(chronos(tree_in, model = "correlated",
                                 control = chronos.control(tol = 1e-8))),
        error = function(e) NULL)
    ultra_ok <- !is.null(tree_ultra)
    if (!ultra_ok) tree_ultra <- tree_in
    tcol <- ape::keep.tip(tree_ultra, reps)
    invA <- tryCatch(inverseA(tcol, scale = ultra_ok)$Ainv,
                     error = function(e) tryCatch(inverseA(tcol, scale = FALSE)$Ainv,
                                                  error = function(e2) NULL))
    if (is.null(invA)) { fail$n_unique <- n_unique; return(fail) }
    df <- data.frame(animal = tip_to_rep, actuation = unname(y_logit[common_h]),
                     mev = unname(mev_logit[common_h]), stringsAsFactors = FALSE)
    df$animal <- factor(df$animal, levels = tcol$tip.label)
    fit_h2 <- function(alpha_V, r_nu = 0.002, r_V = 1) {
        prior <- list(R = list(V = r_V, nu = r_nu),
                      G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = alpha_V)))
        tryCatch(
            MCMCglmm(actuation ~ 1, random = ~animal, ginverse = list(animal = invA),
                     mev = df$mev, data = df, prior = prior,
                     nitt = NITT, burnin = BURNIN, thin = THIN, verbose = FALSE),
            error = function(e) NULL)
    }
    f1000 <- fit_h2(1000); f1 <- fit_h2(1); f025 <- fit_h2(0.25)
    ftR   <- fit_h2(1, r_nu = 1, r_V = 1)
    a1 <- summarise_h2(f1)
    ok <- !is.na(a1$median) && (a1$hi - a1$lo) <= H2_WIDTH_MAX && a1$lo >= H2_LOWER_MIN
    list(n_unique = n_unique, ok = ok, a1000 = summarise_h2(f1000), a1 = a1,
         a025 = summarise_h2(f025), tightR = summarise_h2(ftR),
         fits = list(a1000 = f1000, a1 = f1, a025 = f025))
}

# RE-only tree + alignment (the whole alignment/tree when flank=0).
aln_re <- read_fasta(RE_ALN_PATH)
aln_re <- aln_re[common]; aln_re <- aln_re[!sapply(aln_re, is.null)]
re_tree <- read.tree(RE_TREE_PATH)
re_tree$tip.label <- gsub("^_R_", "", re_tree$tip.label)
re_tree <- tryCatch(keep.tip(re_tree, common), error = function(e) re_tree)

h2_region <- compute_h2(tree, aln)
h2_re     <- compute_h2(re_tree, aln_re)

n_unique_seqs <- h2_region$n_unique
if (n_unique_seqs < 2L) {
    write_empty_header(sprintf(
        "monomorphic alignment — only %d unique sequence across %d tips",
        n_unique_seqs, length(common)))
    quit(status = 0)
}

# ── Attribution mode: RE-only first, expand to whole region, else none ────────
if (isTRUE(h2_re$ok)) {
    attribution_mode <- "RE";     h2_head <- h2_re;     var_cols <- var_cols_re
} else if (isTRUE(h2_region$ok)) {
    attribution_mode <- "region"; h2_head <- h2_region; var_cols <- var_cols_all
} else {
    attribution_mode <- "none"
    h2_head  <- if (!is.na(h2_re$a1$median)) h2_re else h2_region
    var_cols <- integer(0)
}
# Headline H² (TSV h2_alpha1_* etc.) = whichever mode won.
h2_1000 <- h2_head$a1000; h2_1 <- h2_head$a1
h2_025  <- h2_head$a025;  h2_tightR <- h2_head$tightR
cat(sprintf(
    "[phylo_mm] %s H² gate: RE %s (%.2f [%.2f,%.2f]) | region %s (%.2f [%.2f,%.2f]) -> mode=%s\n",
    peak_id,
    if (isTRUE(h2_re$ok)) "OK" else "fail", h2_re$a1$median, h2_re$a1$lo, h2_re$a1$hi,
    if (isTRUE(h2_region$ok)) "OK" else "fail",
    h2_region$a1$median, h2_region$a1$lo, h2_region$a1$hi,
    attribution_mode), file = stderr())

# ── 2. Per-branch coverage-weighted clade-vs-rest test (headline) ─────────────
nTips <- length(tree$tip.label)
root_node <- nTips + 1L
internal_nodes <- setdiff((nTips + 1L):(nTips + tree$Nnode), root_node)

y <- trait[tree$tip.label]
# Coverage = WLS weight; cov_used (median-filled) is named by `common` == tips.
w_cov  <- cov_used[tree$tip.label]
cov_ok <- !is.na(w_cov) & w_cov >= MIN_COV
cat(sprintf("[phylo_mm] %s coverage filter: %d/%d tips have coverage >= %g\n",
            peak_id, sum(cov_ok), length(cov_ok), MIN_COV), file = stderr())

indicator_for_site <- function(col_idx) {
    # Treat every character (A/T/G/C/gap/N) as a valid allele state. The
    # indicator marks the second-most-common state — so a gap polymorphism
    # against a majority base produces a proper 0/1 indicator vector with
    # no missing entries.
    col <- seq_mat[, col_idx]
    if (length(unique(col)) < 2L) return(NULL)
    counts <- sort(table(col), decreasing = TRUE)
    minor <- names(counts)[2]
    col == minor
}
site_indicators <- lapply(var_cols, indicator_for_site)
names(site_indicators) <- as.character(var_cols)

# ── Variant candidates: one block per polymorphic site (correlated sites merged)
# The goal is variant attribution — "which sequence variants contribute to this
# heritable trait" — so the candidates are the VARIABLE SITES themselves (RE-
# scoped via var_cols), not tree bipartitions. Each site is a 0/1 minor-allele
# indicator over the tips. Sites carrying the SAME partition (perfectly
# correlated / complementary) are statistically indistinguishable, so we merge
# them into one block and report all member positions together — "multiple
# variants on a single branch".
tip_order <- match(tree$tip.label, rownames(seq_mat))
canon_key <- function(col) {
    if (length(col) && !is.na(col[1]) && col[1] == 1L) col <- 1L - col
    paste(ifelse(is.na(col), "?", col), collapse = "")
}
blocks <- list()
block_of <- list()
for (cc in names(site_indicators)) {
    ind <- as.integer(site_indicators[[cc]][tip_order])
    if (all(is.na(ind)) || length(unique(ind[!is.na(ind)])) < 2L) next
    # Drop singleton / near-private variants: both alleles need >= MIN_CARRIERS
    # haplotypes, else the "driver" is just one or two samples.
    if (min(sum(ind == 1L, na.rm = TRUE), sum(ind == 0L, na.rm = TRUE)) < MIN_CARRIERS) next
    k <- canon_key(ind)
    if (!is.null(block_of[[k]])) {
        bi <- block_of[[k]]
        blocks[[bi]]$cols <- c(blocks[[bi]]$cols, cc)
    } else {
        blocks[[length(blocks) + 1L]] <- list(cols = cc, ind = ind)
        block_of[[k]] <- length(blocks)
    }
}
# NB: `blocks` may be empty when flanks resolve the tree but the RE itself is
# monomorphic (epigenetic peaks). We do NOT bail here — H² was already fit on
# the flanked tree above and must be reported; the variant table just comes out
# empty (0 contributing variants), handled downstream.

# ── Forward-selection joint regression: which variants jointly contribute ─────
# Coverage-weighted WLS of actuation on the variant indicators. Start empty; at
# each step add the candidate that explains the most ADDITIONAL variance
# (largest gain in weighted R²), and stop when the best gain falls below
# `R2_MIN` or stops being significant. We select on variance explained, NOT on
# p-value: at N≈78 even a β≈0.08 variant that adds <2% of the variance clears
# p<0.05, so a significance gate over-selects weak passengers (the chr19 case:
# pos 77 alone explains 85%, two others add ~1% each yet are "significant").
# The JOINT (conditional) coefficient is the variant's contribution after
# accounting for those already in — so a passenger that merely tracks a stronger
# variant adds ~0 R² and is dropped, while a NESTED variant that switches
# actuation back off among a stronger site's carriers explains real additional
# variance and is kept (with a negative β). Few candidates (2-10 site-blocks),
# so this is exact and stable — no LASSO / selectiveInference needed.
R2_MIN <- 0.05    # a variant must explain >=5% additional actuation variance
if (length(blocks) > 0L) {
    V <- do.call(cbind, lapply(blocks, function(b) b$ind))
    colnames(V) <- paste0("v", seq_along(blocks))
    ok <- cov_ok & !apply(is.na(V), 1L, any)    # rows usable in the joint model
} else {
    V  <- matrix(numeric(0), nrow = length(tree$tip.label), ncol = 0)
    ok <- cov_ok
}
ss_tot <- if (sum(ok) >= 4L) {
    sum(w_cov[ok] * (y[ok] - weighted.mean(y[ok], w_cov[ok]))^2)
} else 0
# Heritability gate: only attribute variants to a HERITABLE trait. When H²ₚ is
# near zero, actuation does not track sequence at all, so any "contributing
# variant" is a spurious weak association (and often anchored on a gap/N column)
# — the premise "which variants drive this heritable trait" doesn't hold. Below
# H2_MIN we leave the candidate variants in the TSV with their marginal stats
# but select NONE (n_variants_selected=0, nothing drawn). H² itself is still
# reported so the reader sees WHY nothing was attributed.
H2_MIN <- 0.1
h2_point <- if (is.na(h2_1$median)) 0 else h2_1$median
selected  <- integer(0)
remaining <- if (h2_point >= H2_MIN) seq_along(blocks) else integer(0)
r2_cur    <- 0
while (length(remaining) > 0L && sum(ok) >= 4L && ss_tot > 0) {
    best <- NA_integer_; best_gain <- -Inf; best_p <- NA_real_; best_r2 <- r2_cur
    for (j in remaining) {
        dfj  <- data.frame(y = y[ok], V[ok, c(selected, j), drop = FALSE])
        fitj <- tryCatch(lm(y ~ ., data = dfj, weights = w_cov[ok]),
                         error = function(e) NULL)
        if (is.null(fitj)) next
        sm <- summary(fitj)$coefficients
        nm <- paste0("v", j)
        if (!(nm %in% rownames(sm))) next         # aliased (collinear) → skip
        r2   <- 1 - sum(w_cov[ok] * residuals(fitj)^2) / ss_tot
        gain <- r2 - r2_cur
        if (gain > best_gain) {
            best_gain <- gain; best <- j
            best_p <- sm[nm, "Pr(>|t|)"]; best_r2 <- r2
        }
    }
    if (is.na(best) || best_gain < R2_MIN ||
        !is.finite(best_p) || best_p >= 0.05) break
    selected  <- c(selected, best)
    remaining <- setdiff(remaining, best)
    r2_cur    <- best_r2
}

# Final joint fit on the selected variants → conditional β, SE, CI, p per block.
joint <- list()
if (length(selected) > 0L && sum(ok) >= 4L) {
    dfs  <- data.frame(y = y[ok], V[ok, selected, drop = FALSE])
    fits <- lm(y ~ ., data = dfs, weights = w_cov[ok])
    sm   <- summary(fits)$coefficients
    ci   <- suppressWarnings(confint(fits))
    for (j in selected) {
        nm <- paste0("v", j)
        if (nm %in% rownames(sm)) {
            joint[[nm]] <- list(
                beta = sm[nm, "Estimate"], se = sm[nm, "Std. Error"],
                ci_lower = ci[nm, 1], ci_upper = ci[nm, 2], p = sm[nm, "Pr(>|t|)"])
        }
    }
}
# Aliases kept for the TSV header / PDF panel below.
selected_branch_ids <- selected
inference_method    <- "forward_wls"
lasso_lambda_used   <- NA_real_
cat(sprintf("[phylo_mm] %s variant attribution: %d candidate site-block(s), %d selected (forward WLS)\n",
            peak_id, length(blocks), length(selected)),
    file = stderr())

# ── Per-variant rows: marginal carriers-vs-REST contrast + joint contribution.
# mean_in/mean_out are the raw "haplotypes WITH this allele vs the rest"
# (marginal); `beta` is the JOINT (conditional) contribution from the forward
# model — the headline effect, since a nested site's marginal effect can be
# muddied by the stronger site it sits under. `is_split` flags the selected
# (contributing) variants; the plot draws those. `descendants` = the allele
# carriers (so the plot highlights them); `linked_sites` = the position(s).
rows <- list()
for (j in seq_along(blocks)) {
    b   <- blocks[[j]]
    ind <- b$ind
    carriers  <- tree$tip.label[which(ind == 1L)]
    in_clade  <- !is.na(ind) & ind == 1L
    out_clade <- !is.na(ind) & ind == 0L
    in_ok  <- in_clade  & cov_ok
    out_ok <- out_clade & cov_ok
    mean_in <- mean_out <- t_stat <- p_raw <- NA_real_
    if (sum(in_ok) >= 2L && sum(out_ok) >= 2L) {
        sub <- in_ok | out_ok
        fit <- tryCatch(lm(y[sub] ~ in_clade[sub], weights = w_cov[sub]),
                        error = function(e) NULL)
        if (!is.null(fit)) {
            sm <- summary(fit)$coefficients
            if (nrow(sm) >= 2L) {
                mean_in  <- weighted.mean(y[in_ok],  w_cov[in_ok])
                mean_out <- weighted.mean(y[out_ok], w_cov[out_ok])
                t_stat   <- unname(sm[2, "t value"])
                p_raw    <- sm[2, "Pr(>|t|)"]
            }
        }
    }
    jt  <- joint[[paste0("v", j)]]
    sel <- j %in% selected
    rows[[length(rows) + 1L]] <- list(
        branch_id         = j,
        n_below           = length(carriers),
        descendants       = paste(carriers, collapse = ","),
        mean_in           = mean_in,
        mean_out          = mean_out,
        t_stat            = t_stat,
        p_raw             = p_raw,
        is_split          = sel,
        linked_sites      = paste(b$cols, collapse = ","),
        selected_by_lasso = sel,
        lasso_coef        = NA_real_,
        beta     = if (!is.null(jt)) jt$beta     else NA_real_,
        se       = if (!is.null(jt)) jt$se       else NA_real_,
        ci_lower = if (!is.null(jt)) jt$ci_lower else NA_real_,
        ci_upper = if (!is.null(jt)) jt$ci_upper else NA_real_,
        p        = if (!is.null(jt)) jt$p        else NA_real_
    )
}
if (length(rows) > 0L) {
    branch_df <- do.call(rbind, lapply(rows, function(r) as.data.frame(r, stringsAsFactors = FALSE)))
    # BH-FDR over the marginal carriers-vs-rest p of all candidate variants.
    branch_df$p_fdr <- p.adjust(branch_df$p_raw, method = "BH")
} else {
    # No RE variants (e.g. epigenetic peak with flanks): empty variant table,
    # but the H² header below is still the real fit on the flanked tree.
    branch_df <- as.data.frame(
        matrix(nrow = 0, ncol = length(empty_branch_cols),
               dimnames = list(NULL, empty_branch_cols)))
}

# Sort: selected (contributing) variants first by joint p, then the rest by
# marginal p.
ord <- order(!branch_df$is_split,
             ifelse(branch_df$is_split, branch_df$p, Inf),
             branch_df$p_raw)
branch_df <- branch_df[ord, ]
# Re-order columns into the canonical layout.
branch_df <- branch_df[, empty_branch_cols]

# ── Write TSV ─────────────────────────────────────────────────────────────────
header <- c(
    paste0("# peak_id=", peak_id),
    paste0("# n_tips=", length(common)),
    paste0("# n_unique_seqs=", n_unique_seqs),
    paste0("# h2_alpha1000_mean=",  round(h2_1000$median, 4)),
    paste0("# h2_alpha1000_lower=", round(h2_1000$lo,     4)),
    paste0("# h2_alpha1000_upper=", round(h2_1000$hi,     4)),
    paste0("# h2_alpha1_mean=",     round(h2_1$median,    4)),
    paste0("# h2_alpha1_lower=",    round(h2_1$lo,        4)),
    paste0("# h2_alpha1_upper=",    round(h2_1$hi,        4)),
    paste0("# h2_alpha025_mean=",   round(h2_025$median,  4)),
    paste0("# h2_alpha025_lower=",  round(h2_025$lo,      4)),
    paste0("# h2_alpha025_upper=",  round(h2_025$hi,      4)),
    paste0("# h2_alpha1_ess=",      if (is.na(h2_1$ess)) "NA" else round(h2_1$ess, 1)),
    paste0("# h2_converged=",
           if (!is.na(h2_1$ess) && h2_1$ess >= ESS_MIN) "TRUE" else "FALSE"),
    # Both H² fits (alpha.V=1) + which one gated attribution. The headline
    # h2_alpha1_* above is whichever mode won (h2_re if RE-driven, else region).
    paste0("# h2_re_mean=",      round(h2_re$a1$median, 4)),
    paste0("# h2_re_lower=",     round(h2_re$a1$lo,     4)),
    paste0("# h2_re_upper=",     round(h2_re$a1$hi,     4)),
    paste0("# h2_re_ok=",        if (isTRUE(h2_re$ok)) "TRUE" else "FALSE"),
    paste0("# h2_region_mean=",  round(h2_region$a1$median, 4)),
    paste0("# h2_region_lower=", round(h2_region$a1$lo,     4)),
    paste0("# h2_region_upper=", round(h2_region$a1$hi,     4)),
    paste0("# h2_region_ok=",    if (isTRUE(h2_region$ok)) "TRUE" else "FALSE"),
    paste0("# attribution_mode=", attribution_mode),
    paste0("# min_cov=", MIN_COV),
    paste0("# flank=", FLANK),
    paste0("# re_first_col=", if (FLANK > 0L) re_lo else "NA"),
    paste0("# re_last_col=",  if (FLANK > 0L) re_hi else "NA"),
    paste0("# n_variants_tested=",  sum(!is.na(branch_df$p_raw))),
    paste0("# n_variants_selected=", length(selected_branch_ids)),
    paste0("# n_branches_selected=", length(selected_branch_ids)),
    paste0("# inference_method=",   inference_method),
    paste0("# nitt=", NITT, " burnin=", BURNIN, " thin=", THIN)
)
con <- file(opt$tsv, "w")
writeLines(header, con)
write.table(branch_df, con, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
close(con)

# ── PDF diagnostics (ad-hoc only; the pipeline does NOT pass --pdf) ────────────
if (!is.null(opt$pdf)) {
pdf(opt$pdf, width = 9, height = 7)
op <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

# Posterior densities of H² under three priors
plot_h2 <- function(fit, label, color) {
    if (is.null(fit)) {
        plot.new(); title(main = paste(label, "(failed)"))
        return(invisible())
    }
    h2 <- fit$VCV[, "animal"] / (fit$VCV[, "animal"] + fit$VCV[, "units"])
    plot(density(h2, from = 0, to = 1), col = color, lwd = 2,
         xlim = c(0, 1), xlab = "H²_p", ylab = "density",
         main = paste0(label, "  median=", sprintf("%.2f", median(h2)),
                       "  HPD=[", sprintf("%.2f", HPDinterval(h2)[1]),
                       ", ", sprintf("%.2f", HPDinterval(h2)[2]), "]"))
    abline(v = median(h2), lty = 2, col = color)
}
# Posteriors for the WINNING H² mode (RE if it worked, else region) — per the
# spec we only render the RE-only panels when the RE-only H² is the one used.
.pdf_fits <- h2_head$fits
.pdf_tag  <- sprintf(" [%s]", attribution_mode)
plot_h2(.pdf_fits$a1000, paste0("alpha.V=1000 (wide)", .pdf_tag),   "#1f77b4")
plot_h2(.pdf_fits$a1,    paste0("alpha.V=1 (moderate)", .pdf_tag),  "#2ca02c")
plot_h2(.pdf_fits$a025,  paste0("alpha.V=0.25 (tight)", .pdf_tag),  "#d62728")

# Bar plot of the contributing variants' joint β with CIs
if (nrow(branch_df) > 0L && any(branch_df$selected_by_lasso)) {
    sel <- branch_df[branch_df$selected_by_lasso, ]
    sel <- sel[order(abs(sel$beta), decreasing = TRUE), ]
    sel <- head(sel, 12)
    nm <- sprintf("pos %s (n=%d)", sel$linked_sites, sel$n_below)
    plot_ci_lo <- pmax(sel$ci_lower, -2, na.rm = FALSE)
    plot_ci_hi <- pmin(sel$ci_upper,  2, na.rm = FALSE)
    yvals <- c(plot_ci_lo, plot_ci_hi, 0)
    yvals <- yvals[is.finite(yvals)]
    ylim_use <- if (length(yvals) >= 2L) range(yvals) else c(-1, 1)
    bp <- barplot(sel$beta, names.arg = nm, las = 2,
                  col = ifelse(sel$p < 0.05, "#d62728", "grey60"),
                  border = NA, ylab = "joint β",
                  ylim = ylim_use,
                  main = "Contributing variants (forward WLS)",
                  cex.names = 0.7)
    arrows(bp, plot_ci_lo, bp, plot_ci_hi,
           angle = 90, code = 3, length = 0.04, col = "black")
    abline(h = 0, col = "black", lty = 1, lwd = 0.6)
} else {
    plot.new(); title(main = "No contributing variants selected")
}
par(op)
dev.off()
}
