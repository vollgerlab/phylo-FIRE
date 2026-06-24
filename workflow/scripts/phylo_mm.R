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
# 2. **LASSO joint branch selection — branch annotation.**
#    For each sequence-anchored bipartition (one mirrored by >=1 polymorphic
#    MSA column) we build a binary indicator over tips. cv.glmnet (alpha=1,
#    family=gaussian, lambda.min) picks the minimal set of branches that
#    *jointly* explain actuation variance, collapsing correlated nested
#    branches to one representative each. selectiveInference::fixedLassoInf
#    then gives valid post-selection beta + CI for the selected set (OLS
#    fallback when the polyhedral check fails). LASSO+SI supply the *effect
#    size* shown on the plot; they are NOT the headline p-value.
#
# 3. **Per-branch coverage-weighted clade-vs-rest scan — headline inference.**
#    For each non-trivial bipartition: a WLS clade-vs-rest test of actuation
#    (weights = per-leaf coverage; leaves below `--min-cov` fibers are dropped
#    from the test), then BH FDR. The bipartitions are pre-specified by the
#    (sequence-derived) tree topology, so this is multiple-testing over fixed
#    hypotheses, not double-dipping; p_fdr is the headline p the plot
#    annotates. The linked-variants column lists MSA positions whose allele
#    pattern matches that bipartition exactly.
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
          !is.null(opt$tsv), !is.null(opt$pdf))
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

NITT   <- as.integer(if (is.null(opt$nitt))   50000 else opt$nitt)
BURNIN <- as.integer(if (is.null(opt$burnin)) 10000 else opt$burnin)
THIN   <- as.integer(if (is.null(opt$thin))      25 else opt$thin)

# ── load tree + actuation TSV ─────────────────────────────────────────────────
# Midpoint-root to MATCH the plot's Bio.Phylo root_at_midpoint(). The rooting
# decides which side of each edge is the reported "clade"; with an arbitrary
# outgroup root that side can be the non-contiguous complement of the displayed
# clade, which scatters the coloured leader lines and mis-places branch
# annotations. (phangorn::midpoint isn't installed — root by hand at the
# midpoint of the longest tip-to-tip path, the same definition Bio.Phylo uses.)
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
tree <- tryCatch(midpoint_root(tree), error = function(e) {
    if (!is.rooted(tree)) {
        root(tree, outgroup = tree$tip.label[1], resolve.root = TRUE)
    } else {
        tree
    }
})
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
    "mean_in", "mean_out", "t_stat", "p_raw", "p_fdr", "linked_sites",
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
        paste0("# min_cov=", MIN_COV),
        paste0("# flank=", FLANK),
        paste0("# re_first_col=NA"),
        paste0("# re_last_col=NA"),
        paste0("# n_branches_tested=0"),
        paste0("# n_branches_selected=0"),
        paste0("# lasso_lambda=NA"),
        paste0("# inference_method=none"),
        paste0("# note=", reason)
    )
    writeLines(c(header, paste(empty_branch_cols, collapse = "\t")), opt$tsv)
    pdf(opt$pdf, width = 6, height = 4)
    plot.new(); title(main = paste(peak_id, "—", reason))
    dev.off()
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
# Restricted to the RE band so neutral flank variation never drives the
# site-based association — the flanks resolve the tree, the RE drives the test.
var_cols <- if (is.null(seq_mat)) integer(0) else which(apply(seq_mat, 2, function(col) {
    length(unique(col)) >= 2L
}))
var_cols <- var_cols[var_cols >= re_lo & var_cols <= re_hi]
cat(sprintf("[phylo_mm] %s RE band: cols %d-%d (flank=%d); %d variable RE sites\n",
            peak_id, re_lo, re_hi, FLANK, length(var_cols)), file = stderr())

# ── 1. Heritability sensitivity (three alpha priors) ──────────────────────────
# Primer recipe applied here (per the Halliwell/de Villemereuil/Hadfield notes):
#   - Logit-transform the bounded [0,1] trait with a coverage-derived offset
#     so the Gaussian likelihood is on an unconstrained scale (`y_logit`).
#   - Pass per-observation measurement-error variance (`mev`) so MCMCglmm
#     properly downweights low-coverage observations rather than treating
#     a 5-fiber proportion the same as a 50-fiber one. On the logit scale
#     the binomial-sampling SE squared transforms by the delta method to
#     1 / (n * p * (1-p)).
#   - Run the model on an ultrametric tree (chronos) so the random-effect
#     covariance matrix has interpretable depth and the prior on σ²_phylo
#     applies to a consistent time scale across peaks.
# These three changes don't affect the per-branch / LASSO sections (those
# stay on the raw 0-1 actuation scale, which is what makes β values
# directly readable as "actuation difference between clade and rest").

# Ultrametrify the tree (used ONLY for the H² model). chronos can be
# numerically fragile on near-zero branches — fall back to the original
# phylogram if it errors.
tree_ultra <- tryCatch(
    suppressWarnings(chronos(tree, model = "correlated",
                             control = chronos.control(tol = 1e-8))),
    error = function(e) tree
)

# ── Identical-haplotype collapse (Option A) ────────────────────────────────
# Cluster tips by EXACT aligned-sequence string (gaps count — only collapse
# truly identical sequences, never gap-vs-base or differing-gap patterns).
# Each cluster gets one rep tip in the collapsed tree; all members of the
# cluster share that rep's name as their `animal` factor level in the MCMCglmm
# data frame. Multiple observations of the same cluster contribute as
# replicates: variation among them is absorbed by the Gaussian residual
# (σ²_R), keeping σ²_animal as the strictly between-sequence variance.
#
# This fixes three problems on nearly-clonal haplotype trees: (1) inverse-A
# was nearly singular for same-sequence tips, (2) σ²_animal was absorbing
# within-cluster noise that should sit in σ²_R, and (3) the effective N
# we report for H² should be N_unique not N_tips.
aln_full_strings <- sapply(aln[common], function(s) s)   # exact, including gaps
cluster_idx      <- match(aln_full_strings, unique(aln_full_strings))
n_unique_seqs    <- max(cluster_idx)
rep_per_cluster  <- common[!duplicated(cluster_idx)]      # one rep per cluster
tip_to_rep       <- rep_per_cluster[cluster_idx]          # per-leaf rep name

if (n_unique_seqs < 2L) {
    write_empty_header(sprintf(
        "monomorphic alignment — only %d unique sequence across %d tips",
        n_unique_seqs, length(common)))
    quit(status = 0)
}

tree_collapsed <- ape::keep.tip(tree_ultra, rep_per_cluster)
inv_A <- inverseA(tree_collapsed, scale = TRUE)$Ainv

# Logit trait + mev. eps offset: half a fiber's worth of probability,
# floored by 1/(2 * max_cov) so we never collapse to numerical infinity.
trait_eps <- max(1e-4, 0.5 / max(cov_vec[common], na.rm = TRUE))
trait_clip <- pmin(pmax(trait[common], trait_eps), 1 - trait_eps)
y_logit    <- qlogis(trait_clip)
# Logit-scale measurement-error variance via the delta method.
# var_logit ≈ var_raw / (p*(1-p))^2 = (p(1-p)/n) / (p(1-p))^2 = 1/(n*p*(1-p))
cov_used <- cov_vec[common]
if (any(is.na(cov_used))) {
    cov_used[is.na(cov_used)] <- median(cov_used, na.rm = TRUE)
}
mev_logit <- 1 / (cov_used * trait_clip * (1 - trait_clip))

mm_df <- data.frame(
    animal    = tip_to_rep,                  # rep name (shared across cluster)
    actuation = unname(y_logit),
    mev       = unname(mev_logit),
    stringsAsFactors = FALSE
)
mm_df$animal <- factor(mm_df$animal, levels = tree_collapsed$tip.label)

fit_h2 <- function(alpha_V, r_nu = 0.002, r_V = 1) {
    prior <- list(
        R = list(V = r_V, nu = r_nu),
        G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = alpha_V))
    )
    tryCatch(
        MCMCglmm(actuation ~ 1, random = ~animal,
                 ginverse = list(animal = inv_A),
                 mev      = mm_df$mev,
                 data     = mm_df, prior = prior,
                 nitt = NITT, burnin = BURNIN, thin = THIN, verbose = FALSE),
        error = function(e) NULL
    )
}

# Posterior summary: median + 95% HPD plus the underlying σ² medians so we
# can diagnose whether H² inflation is being driven by a shrinking σ²_R
# (denominator collapse) vs. a genuinely large σ²_phylo (numerator growth).
summarise_h2 <- function(fit) {
    if (is.null(fit)) return(list(median = NA_real_, lo = NA_real_, hi = NA_real_,
                                  ess = NA_real_, sigma_p = NA_real_, sigma_r = NA_real_))
    s_p <- fit$VCV[, "animal"]
    s_r <- fit$VCV[, "units"]
    h2  <- s_p / (s_p + s_r)
    ci  <- HPDinterval(h2, prob = 0.95)
    list(
        median  = median(h2),
        lo      = ci[1],
        hi      = ci[2],
        ess     = unname(effectiveSize(h2)),
        sigma_p = median(s_p),
        sigma_r = median(s_r)
    )
}

# ── Diagnostic: is mev overstated relative to within-cluster spread? ───────
# If MEDIAN(mev) >> MEDIAN(within-cluster var(y_logit)), the model is being
# told most residual noise is already "explained" → σ²_R shrinks → H² inflates.
mev_q <- quantile(mm_df$mev, c(0.25, 0.5, 0.75), na.rm = TRUE)
within_var <- tapply(mm_df$actuation, mm_df$animal,
                     function(v) if (length(v) >= 2L) var(v) else NA_real_)
within_var <- within_var[!is.na(within_var)]
cat(sprintf(
    "[phylo_mm] %s mev: q25=%.3f median=%.3f q75=%.3f  | within-cluster var(y_logit): n_rep_clusters=%d %s\n",
    peak_id, mev_q[1], mev_q[2], mev_q[3], length(within_var),
    if (length(within_var) > 0)
        sprintf("median=%.3f range=[%.3f, %.3f]",
                median(within_var), min(within_var), max(within_var))
    else "(no replicated clusters — within-cluster σ²_R is unidentified)"
), file = stderr())

fit_a1000 <- fit_h2(1000)
fit_a1    <- fit_h2(1)
fit_a025  <- fit_h2(0.25)
# Tight-R sensitivity: same G prior as the default (alpha.V=1) but with an
# informative R prior pinning σ²_R near 1 logit². If H² drops materially vs.
# the default, R was being squeezed to zero by the weakly-informative nu=0.002.
fit_a1_tightR <- fit_h2(1, r_nu = 1, r_V = 1)
h2_1000   <- summarise_h2(fit_a1000)
h2_1      <- summarise_h2(fit_a1)
h2_025    <- summarise_h2(fit_a025)
h2_tightR <- summarise_h2(fit_a1_tightR)
cat(sprintf(
    "[phylo_mm] %s H²: alpha=1000 -> %.3f [%.3f, %.3f]  alpha=1 -> %.3f [%.3f, %.3f]  alpha=0.25 -> %.3f [%.3f, %.3f]  alpha=1+tightR -> %.3f [%.3f, %.3f]\n",
    peak_id,
    h2_1000$median,   h2_1000$lo,   h2_1000$hi,
    h2_1$median,      h2_1$lo,      h2_1$hi,
    h2_025$median,    h2_025$lo,    h2_025$hi,
    h2_tightR$median, h2_tightR$lo, h2_tightR$hi
), file = stderr())
cat(sprintf(
    "[phylo_mm] %s σ²: alpha=1000 (σ²p=%.3f σ²r=%.3f)  alpha=1 (σ²p=%.3f σ²r=%.3f)  alpha=0.25 (σ²p=%.3f σ²r=%.3f)  alpha=1+tightR (σ²p=%.3f σ²r=%.3f)\n",
    peak_id,
    h2_1000$sigma_p,   h2_1000$sigma_r,
    h2_1$sigma_p,      h2_1$sigma_r,
    h2_025$sigma_p,    h2_025$sigma_r,
    h2_tightR$sigma_p, h2_tightR$sigma_r
), file = stderr())
cat(sprintf(
    "[phylo_mm] %s convergence: H²(alpha=1) ESS=%s -> %s (threshold %d)\n",
    peak_id,
    if (is.na(h2_1$ess)) "NA" else sprintf("%.0f", h2_1$ess),
    if (!is.na(h2_1$ess) && h2_1$ess >= ESS_MIN) "CONVERGED" else "NOT CONVERGED",
    ESS_MIN
), file = stderr())

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

# Build a single design matrix of all non-trivial bipartition indicators
# (re-used by both per-branch and LASSO).
X_cols <- list()
br_meta <- list()
for (node in internal_nodes) {
    desc <- extract.clade(tree, node)$tip.label
    n_below <- length(desc)
    if (n_below < 2L || n_below > (nTips - 2L)) next
    ind <- as.integer(tree$tip.label %in% desc)
    X_cols[[length(X_cols) + 1L]] <- ind
    br_meta[[length(br_meta) + 1L]] <- list(node = node, desc = desc, n_below = n_below)
}
n_branches_raw <- length(br_meta)
if (n_branches_raw == 0L) {
    write_empty_header("no non-trivial bipartitions")
    quit(status = 0)
}

X <- do.call(cbind, X_cols)
rownames(X) <- tree$tip.label

# ── RE-anchoring: which RE variants exactly mirror each bipartition ──────────
# RE-anchored = `linked` non-empty. This drives the LASSO candidate set AND
# which branches get a p-value. The COMPARISON GROUP for each branch's test is
# decided after the LASSO (below), conditional on the significance hierarchy —
# so we only need the linked sites here.
linked_list <- vector("list", length(br_meta))
for (i in seq_along(br_meta)) {
    indicator <- X[, i] == 1L
    linked <- character(0)
    for (col_chr in names(site_indicators)) {
        ind_site <- site_indicators[[col_chr]]
        if (is.null(ind_site)) next
        ind_site_t <- ind_site[match(tree$tip.label, rownames(seq_mat))]
        if (anyNA(ind_site_t)) {
            mask <- !is.na(ind_site_t)
            if (sum(mask) < 4) next
            match_eq  <- all(ind_site_t[mask] == indicator[mask])
            match_neg <- all(ind_site_t[mask] != indicator[mask])
        } else {
            match_eq  <- all(ind_site_t == indicator)
            match_neg <- all(ind_site_t != indicator)
        }
        if (match_eq || match_neg) linked <- c(linked, col_chr)
    }
    linked_list[[i]] <- linked
}

# ── 3. LASSO joint branch selection ───────────────────────────────────────────
# Restrict the design matrix to **sequence-anchored bipartitions**: branches
# whose in/out partition is mirrored by at least one polymorphic MSA column
# (i.e., `linked_sites` non-empty in the per-branch scan). This prevents
# LASSO from picking spurious within-identical-cluster bipartitions that
# IQ-TREE has to arbitrarily resolve when sequences are tied. On chr17
# (only 1 variable site) this collapses the candidate set to exactly the
# one biologically-meaningful bipartition; on richer peaks the filter is
# permissive (any branch backed by ≥1 perfect-match SNP qualifies).
# After filtering, we deduplicate identical columns (perfectly nested
# bipartitions that share the same indicator vector).
linked_per_branch <- vapply(linked_list, function(l) length(l) > 0L, logical(1))
seq_anchored <- linked_per_branch
n_filtered_out <- sum(!seq_anchored)
if (any(seq_anchored)) {
    X_seq <- X[, seq_anchored, drop = FALSE]
    br_meta_seq <- br_meta[seq_anchored]
} else {
    # No sequence-anchored bipartition — leave the design empty so LASSO
    # selects 0 branches (will hit the n_branches_selected=0 branch below).
    X_seq <- X[, integer(0), drop = FALSE]
    br_meta_seq <- list()
}
col_keys <- if (ncol(X_seq) > 0L) apply(X_seq, 2, paste, collapse = "") else character(0)
unique_idx <- !duplicated(col_keys)
X_uniq <- X_seq[, unique_idx, drop = FALSE]
br_meta_uniq <- br_meta_seq[unique_idx]
if (ncol(X_uniq) > 0L) {
    colnames(X_uniq) <- paste0("br", sapply(br_meta_uniq, function(m) m$node))
}
cat(sprintf("[phylo_mm] %s LASSO design: %d branches passed sequence-anchor filter (%d filtered out), %d unique cols\n",
            peak_id, sum(seq_anchored), n_filtered_out, ncol(X_uniq)),
    file = stderr())

selected_branch_ids <- integer(0)
lasso_lambda_used   <- NA_real_
inference_table     <- NULL
inference_method    <- "none"   # selectiveInference | ols_fallback | univariate_ols | none

if (ncol(X_uniq) == 1L && length(y) >= 4L) {
    # Single sequence-anchored candidate — no LASSO needed (the candidate
    # set has size 1, there's no selection event, so no double-dipping).
    # Fit a univariate OLS directly. This recovers the headline branch
    # for peaks with one variant (chr17) or fully-correlated variants
    # collapsing to one bipartition (chr1's sites 21/38).
    x1 <- as.numeric(X_uniq[, 1L])
    ols_fit <- lm(y ~ x1)
    sm <- summary(ols_fit)$coefficients
    if (nrow(sm) >= 2L) {
        ci_ols <- confint(ols_fit)[-1, , drop = FALSE]
        bid <- br_meta_uniq[[1L]]$node
        selected_branch_ids <- bid
        lasso_lambda_used   <- 0   # univariate path, no LASSO penalty
        inference_table <- data.frame(
            branch_id  = bid,
            lasso_coef = sm[2, "Estimate"],
            beta       = sm[2, "Estimate"],
            se         = sm[2, "Std. Error"],
            ci_lower   = ci_ols[1, 1],
            ci_upper   = ci_ols[1, 2],
            p          = sm[2, "Pr(>|t|)"],
            stringsAsFactors = FALSE
        )
        inference_method <- "univariate_ols"
    }
} else if (ncol(X_uniq) >= 2L && length(y) >= 4L) {
    # standardize=FALSE matches the selectiveInference vignette
    # (Tibshirani et al.) — our design matrix is already 0/1 branch
    # indicators on a common scale, so standardisation adds nothing
    # and creates a scale mismatch that fails the polyhedral check.
    # thresh=1e-12 tightens KKT for the post-selection check.
    cv <- tryCatch(
        cv.glmnet(X_uniq, y, alpha = 1, family = "gaussian",
                  nfolds = min(10L, length(y)),
                  standardize = FALSE,
                  thresh = 1e-12),
        error = function(e) NULL
    )
    if (!is.null(cv)) {
        # lambda.min (CV-optimal) rather than lambda.1se — the latter is
        # routinely too aggressive with small effective sample sizes (e.g.
        # peaks with only 2 unique phenotype levels) and picks 0 branches.
        lasso_lambda_used <- cv$lambda.min
        # exact=TRUE refits glmnet at this exact lambda (instead of
        # interpolating between gridpoints) — required by fixedLassoInf
        # to satisfy its polyhedral-constraint check.
        coefs_lasso <- as.numeric(coef(cv, s = "lambda.min",
                                       exact = TRUE,
                                       x = X_uniq, y = y))[-1]
        sel_cols <- which(coefs_lasso != 0)
        if (length(sel_cols) > 0L) {
            selected_branch_ids <- sapply(br_meta_uniq[sel_cols], function(m) m$node)
            n_obs <- length(y)
            # Post-selection inference via Lee/Sun/Sun/Taylor 2016. Lambda
            # is rescaled to selectiveInference's convention: glmnet's
            # loss is (1/2N)·RSS + λ·|β|₁, while fixedLassoInf expects
            # (1/2)·RSS + λ·|β|₁, so λ_si = N · λ_glmnet.
            # Suppress only errors — fixedLassoInf emits informative
            # warnings (e.g. truncated CIs) that we want surfaced rather
            # than treated as failures.
            si <- tryCatch(
                suppressWarnings(
                    fixedLassoInf(X_uniq, y,
                                  beta   = coefs_lasso,
                                  lambda = lasso_lambda_used * n_obs,
                                  alpha  = 0.05)   # 95% CIs
                ),
                error = function(e) {
                    message("[phylo_mm] fixedLassoInf failed for ", peak_id, ": ", conditionMessage(e))
                    NULL
                }
            )
            if (!is.null(si) && length(si$vars) > 0L &&
                length(si$pv) == length(si$vars)) {
                # Map selectiveInference's active-set indices back to branch IDs.
                sel_meta <- br_meta_uniq[si$vars]
                inference_table <- data.frame(
                    branch_id  = sapply(sel_meta, function(m) m$node),
                    lasso_coef = coefs_lasso[si$vars],
                    beta       = si$coef0,
                    se         = si$sd,
                    ci_lower   = si$ci[, 1],
                    ci_upper   = si$ci[, 2],
                    p          = si$pv,
                    stringsAsFactors = FALSE
                )
                selected_branch_ids <- inference_table$branch_id
                inference_method <- "selectiveInference"
            } else {
                # Fallback to naive OLS so we still produce *some* point
                # estimates + CIs. p-values are flagged as naive in the
                # header (`inference_method=ols_fallback`); consumers
                # should treat them as descriptive, not confirmatory.
                X_sel <- X_uniq[, sel_cols, drop = FALSE]
                ols_df <- data.frame(y = y, X_sel)
                ols_fit <- lm(y ~ ., data = ols_df)
                sm <- summary(ols_fit)$coefficients[-1, , drop = FALSE]
                ci_ols <- confint(ols_fit)[-1, , drop = FALSE]
                inference_table <- data.frame(
                    branch_id  = selected_branch_ids,
                    lasso_coef = coefs_lasso[sel_cols],
                    beta       = sm[, "Estimate"],
                    se         = sm[, "Std. Error"],
                    ci_lower   = ci_ols[, 1],
                    ci_upper   = ci_ols[, 2],
                    p          = sm[, "Pr(>|t|)"],
                    stringsAsFactors = FALSE
                )
                inference_method <- "ols_fallback"
            }
        }
    }
}

cat(sprintf("[phylo_mm] %s LASSO selected %d branches (lambda.min=%.4g, %d unique bipartitions, inference=%s)\n",
            peak_id, length(selected_branch_ids),
            lasso_lambda_used, ncol(X_uniq), inference_method),
    file = stderr())

# ── Per-branch test, conditional on the LASSO significance hierarchy ──────────
# Comparison group ("out") for a branch B with tip set B_tips:
#   - if a LASSO-selected clade STRICTLY contains B, the innermost such clade
#     minus B_tips — B is tested against the rest of the significant clade it
#     sits in (its *incremental* effect within the signal); else
#   - the rest of the tree — B is a top-level signal (no significant clade
#     above it), tested for its overall effect.
# Falls back to vs-rest when the conditional remainder has <2 testable tips.
# Using the LASSO selection as the hierarchy (rather than the p-values
# themselves) avoids circularity, and makes the marginal Δ a CONDITIONAL effect
# that aligns with the LASSO joint β instead of depending on an arbitrary
# tree-sibling. Coverage-weighted WLS, min-coverage floor, RE-anchored only.
node_to_desc <- setNames(
    lapply(br_meta, function(m) m$desc),
    vapply(br_meta, function(m) as.character(m$node), character(1)))
sel_desc <- node_to_desc[as.character(selected_branch_ids)]
comparison_tips <- function(B_tips) {
    containing <- Filter(
        function(s) length(s) > length(B_tips) && all(B_tips %in% s), sel_desc)
    if (length(containing) > 0L) {
        setdiff(containing[[which.min(vapply(containing, length, integer(1)))]],
                B_tips)
    } else {
        setdiff(tree$tip.label, B_tips)
    }
}

rows <- list()
for (i in seq_along(br_meta)) {
    meta   <- br_meta[[i]]
    linked <- linked_list[[i]]
    in_clade <- tree$tip.label %in% meta$desc
    in_out   <- tree$tip.label %in% comparison_tips(meta$desc)
    if (sum(in_out & cov_ok) < 2L) in_out <- !in_clade   # fall back to vs-rest
    in_ok  <- in_clade & cov_ok
    out_ok <- in_out   & cov_ok
    mean_in <- mean_out <- t_stat <- p_raw <- NA_real_
    if (length(linked) > 0L && sum(in_ok) >= 2L && sum(out_ok) >= 2L) {
        sub <- in_ok | out_ok                # restrict the test to clade ∪ out
        yk  <- y[sub]; ind_k <- in_clade[sub]; wk <- w_cov[sub]
        fit <- tryCatch(lm(yk ~ ind_k, weights = wk), error = function(e) NULL)
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
    rows[[length(rows) + 1L]] <- list(
        branch_id    = meta$node,
        n_below      = meta$n_below,
        descendants  = paste(meta$desc, collapse = ","),
        mean_in      = mean_in,
        mean_out     = mean_out,
        t_stat       = t_stat,
        p_raw        = p_raw,
        linked_sites = paste(linked, collapse = ",")
    )
}
branch_df <- do.call(rbind, lapply(rows, as.data.frame))
branch_df$p_fdr <- p.adjust(branch_df$p_raw, method = "BH")

# ── Annotate branch_df with LASSO selection + inference columns ───────────────
branch_df$selected_by_lasso <- branch_df$branch_id %in% selected_branch_ids
for (col in c("lasso_coef", "beta", "se", "ci_lower", "ci_upper", "p")) {
    branch_df[[col]] <- NA_real_
}
if (!is.null(inference_table)) {
    for (i in seq_len(nrow(inference_table))) {
        bid <- inference_table$branch_id[i]
        rmask <- branch_df$branch_id == bid
        branch_df$lasso_coef[rmask] <- inference_table$lasso_coef[i]
        branch_df$beta[rmask]       <- inference_table$beta[i]
        branch_df$se[rmask]         <- inference_table$se[i]
        branch_df$ci_lower[rmask]   <- inference_table$ci_lower[i]
        branch_df$ci_upper[rmask]   <- inference_table$ci_upper[i]
        branch_df$p[rmask]          <- inference_table$p[i]
    }
}

# Sort: LASSO-selected first (by post-selection p), then remainder by p_fdr.
ord <- order(!branch_df$selected_by_lasso,
             ifelse(branch_df$selected_by_lasso,
                    branch_df$p, branch_df$p_fdr),
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
    paste0("# min_cov=", MIN_COV),
    paste0("# flank=", FLANK),
    paste0("# re_first_col=", if (FLANK > 0L) re_lo else "NA"),
    paste0("# re_last_col=",  if (FLANK > 0L) re_hi else "NA"),
    paste0("# n_branches_tested=",  sum(!is.na(branch_df$p_raw))),
    paste0("# n_branches_selected=", length(selected_branch_ids)),
    paste0("# lasso_lambda=",       if (is.na(lasso_lambda_used)) "NA" else format(lasso_lambda_used, digits = 6)),
    paste0("# inference_method=",   inference_method),
    paste0("# nitt=", NITT, " burnin=", BURNIN, " thin=", THIN)
)
con <- file(opt$tsv, "w")
writeLines(header, con)
write.table(branch_df, con, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
close(con)

# ── PDF diagnostics: H² posteriors + LASSO selection plot ─────────────────────
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
plot_h2(fit_a1000, "alpha.V=1000 (wide)",   "#1f77b4")
plot_h2(fit_a1,    "alpha.V=1   (moderate)", "#2ca02c")
plot_h2(fit_a025,  "alpha.V=0.25 (tight)",   "#d62728")

# LASSO bar plot of selected branches' OLS β with CIs
if (nrow(branch_df) > 0L && any(branch_df$selected_by_lasso)) {
    sel <- branch_df[branch_df$selected_by_lasso, ]
    sel <- sel[order(abs(sel$beta), decreasing = TRUE), ]
    sel <- head(sel, 12)
    nm <- sprintf("br%d (n=%d)", sel$branch_id, sel$n_below)
    # selectiveInference can return ±Inf when the truncated Gaussian is
    # degenerate. Clip CIs to a sensible plotting range so plot.window
    # doesn't choke on infinite limits.
    plot_ci_lo <- pmax(sel$ci_lower, -2, na.rm = FALSE)
    plot_ci_hi <- pmin(sel$ci_upper,  2, na.rm = FALSE)
    yvals <- c(plot_ci_lo, plot_ci_hi, 0)
    yvals <- yvals[is.finite(yvals)]
    ylim_use <- if (length(yvals) >= 2L) range(yvals) else c(-1, 1)
    bp <- barplot(sel$beta, names.arg = nm, las = 2,
                  col = ifelse(sel$p < 0.05, "#d62728", "grey60"),
                  border = NA, ylab = "β",
                  ylim = ylim_use,
                  main = sprintf("LASSO-selected branches (%s)", inference_method),
                  cex.names = 0.7)
    arrows(bp, plot_ci_lo, bp, plot_ci_hi,
           angle = 90, code = 3, length = 0.04, col = "black")
    abline(h = 0, col = "black", lty = 1, lwd = 0.6)
} else {
    plot.new(); title(main = "LASSO selected 0 branches")
}
par(op)
dev.off()
