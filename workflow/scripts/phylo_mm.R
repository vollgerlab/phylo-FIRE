#!/usr/bin/env Rscript
# Per-peak phylogenetic heritability (MCMCglmm) + variant attribution.
#
# 1. H²ₚ: intercept-only MCMCglmm on the RE-only tree and on the flanked region
#    tree, each under alpha.V = 1000 / 1 / 0.25. alpha.V=1 is the headline; if
#    the priors disagree, H²ₚ is prior-dominated (Hadfield 2015).
# 2. Attribution runs only if an H² gate passes (RE first, then region). Forward
#    selection over variable-site blocks adds the block with the largest
#    coverage-weighted R² gain and stops when the gain is < R2_MIN (5%) or its
#    p >= 0.05. The joint fit gives beta/se/ci/p; `is_split` marks those rows.
# 3. Per-site carriers-vs-rest WLS (mean_in/out, t_stat, p_raw, BH p_fdr) is
#    descriptive only. selected_by_lasso/lasso_coef are legacy column names.
#
# Usage: Rscript phylo_mm.R --tree T --actuation A --alignment F --tsv OUT
#   [--re-tree --re-alignment --peak-id --flank --min-cov --min-carriers
#    --h2-width-max --h2-lower-min --nitt --burnin --thin]

suppressPackageStartupMessages({
    library(MCMCglmm)
    library(ape)
})

# Minimum effective sample size for the headline (alpha.V=1) H2_p posterior to
# be called converged. Reported in the TSV header as h2_converged.
ESS_MIN <- 200

# A variant must explain >= R2_MIN additional actuation variance to be selected.
R2_MIN <- 0.05

# Below H2_MIN nothing is selected (candidates and H² are still reported). This only
# binds if H2_LOWER_MIN is set below H2_MIN; the H² gate already implies it.
H2_MIN <- 0.1

BRANCH_COLS <- c(
    "branch_id", "n_below", "descendants",
    "mean_in", "mean_out", "t_stat", "p_raw", "p_fdr", "is_split", "linked_sites",
    "selected_by_lasso", "lasso_coef",
    "beta", "se", "ci_lower", "ci_upper", "p"
)

# ── arguments ─────────────────────────────────────────────────────────────────
parse_args <- function(args) {
    out <- list()
    i <- 1L
    while (i <= length(args)) {
        out[[sub("^--", "", args[i])]] <- args[i + 1L]
        i <- i + 2L
    }
    out
}

opt_or <- function(opt, key, default) if (is.null(opt[[key]])) default else opt[[key]]

# Run settings from the command line (defaults match config.yaml).
load_config <- function(opt) {
    list(
        # Leaves with < min_cov fibers are dropped from the (coverage-weighted) WLS.
        min_cov = as.numeric(opt_or(opt, "min-cov", 4)),
        # bp of flank add_flanks padded each side of the RE (0 = none).
        flank = as.integer(opt_or(opt, "flank", 0)),
        # Attribution gate: alpha.V=1 HPD width <= h2_width_max and lower >= h2_lower_min.
        h2_width_max = as.numeric(opt_or(opt, "h2-width-max", 0.5)),
        h2_lower_min = as.numeric(opt_or(opt, "h2-lower-min", 0.1)),
        # RE-only tree + alignment; the whole region when not supplied (flank=0).
        re_tree = opt_or(opt, "re-tree", opt$tree),
        re_alignment = opt_or(opt, "re-alignment", opt$alignment),
        # Each allele needs >= min_carriers haplotypes: a singleton trivially
        # "explains" its own actuation.
        min_carriers = as.integer(opt_or(opt, "min-carriers", 2)),
        nitt = as.integer(opt_or(opt, "nitt", 50000)),
        burnin = as.integer(opt_or(opt, "burnin", 10000)),
        thin = as.integer(opt_or(opt, "thin", 25))
    )
}

# ── inputs ────────────────────────────────────────────────────────────────────
# root_tree.R already rooted this tree (ancestral proxy); midpoint_root is a fallback.
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

load_tree <- function(path) {
    tree <- read.tree(path)
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
    tree
}

# Actuation per sample, and coverage (number of fibers) per sample. Coverage gives
# the logit-scale measurement error variance for MCMCglmm (Hadfield 2010
# meta-analysis trick: per-observation precision in the residual likelihood).
load_actuation <- function(path) {
    df <- read.table(path, header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE, comment.char = "", quote = "")
    keep <- !is.na(df$actuation) & df$actuation != "NA"
    trait <- suppressWarnings(as.numeric(df$actuation[keep]))
    names(trait) <- df$sample[keep]
    trait <- trait[!is.na(trait)]
    cov_vec <- suppressWarnings(as.numeric(df$coverage[keep]))
    names(cov_vec) <- df$sample[keep]
    cov_vec <- cov_vec[!is.na(cov_vec) & cov_vec > 0]
    list(trait = trait, cov_vec = cov_vec)
}

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

# Alignment restricted to the tips in `common` (missing tips dropped).
read_alignment <- function(path, common) {
    aln <- read_fasta(path)
    aln <- aln[common]
    aln[!sapply(aln, is.null)]
}

# RE column band [re_lo, re_hi] (1-based): skip FLANK non-gap chars at each row's ends
# and union the interior spans. NULL (whole alignment) when flank=0.
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

# Variable column: >=2 states (gap/N count as alleles so segregating indels qualify)
# and >=50% real bases (skips indel / extraction-window edge columns).
variable_columns <- function(seq_mat) {
    if (is.null(seq_mat)) return(integer(0))
    which(apply(seq_mat, 2, function(col) {
        real <- col %in% c("A", "C", "G", "T")
        length(unique(col)) >= 2L && mean(real) >= 0.5
    }))
}

# ── heritability ──────────────────────────────────────────────────────────────
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

# Recipe (de Villemereuil/Hadfield): logit trait with a coverage offset, binomial
# mev = 1/(n·p·(1-p)), ultrametric (chronos) tree. Shared by both H² scopes.
h2_trait <- function(trait, cov_vec, common) {
    cov_used <- cov_vec[common]
    if (any(is.na(cov_used))) cov_used[is.na(cov_used)] <- median(cov_used, na.rm = TRUE)
    trait_eps  <- max(1e-4, 0.5 / max(cov_vec[common], na.rm = TRUE))
    trait_clip <- pmin(pmax(trait[common], trait_eps), 1 - trait_eps)
    list(common = common, cov_used = cov_used,
         y_logit = setNames(qlogis(trait_clip), common),
         mev_logit = setNames(1 / (cov_used * trait_clip * (1 - trait_clip)), common))
}

# H² sensitivity (3 alpha.V priors + tight-R) on one (tree, alignment) scope. Identical
# sequences collapse to one tip; within-cluster spread is the residual. `ok` = the
# alpha.V=1 HPD width <= h2_width_max and lower >= h2_lower_min.
compute_h2 <- function(tree_in, aln_in, ht, cfg) {
    # Drop IQ-TREE support labels: inverseA rejects duplicates (e.g. two "100" nodes).
    tree_in$node.label <- NULL
    NA3 <- list(median = NA_real_, lo = NA_real_, hi = NA_real_, ess = NA_real_,
                sigma_p = NA_real_, sigma_r = NA_real_)
    fail <- list(n_unique = NA_integer_, ok = FALSE,
                 a1000 = NA3, a1 = NA3, a025 = NA3, tightR = NA3)
    # Intersect tips per scope (the RE scope can lack some); otherwise vapply() errors
    # and fails the whole peak instead of returning ok=FALSE.
    common_h <- intersect(intersect(ht$common, tree_in$tip.label), names(aln_in))
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
    df <- data.frame(animal = tip_to_rep, actuation = unname(ht$y_logit[common_h]),
                     mev = unname(ht$mev_logit[common_h]), stringsAsFactors = FALSE)
    df$animal <- factor(df$animal, levels = tcol$tip.label)
    fit_h2 <- function(alpha_V, r_nu = 0.002, r_V = 1) {
        prior <- list(R = list(V = r_V, nu = r_nu),
                      G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = alpha_V)))
        tryCatch(
            MCMCglmm(actuation ~ 1, random = ~animal, ginverse = list(animal = invA),
                     mev = df$mev, data = df, prior = prior,
                     nitt = cfg$nitt, burnin = cfg$burnin, thin = cfg$thin,
                     verbose = FALSE),
            error = function(e) NULL)
    }
    f1000 <- fit_h2(1000); f1 <- fit_h2(1); f025 <- fit_h2(0.25)
    ftR   <- fit_h2(1, r_nu = 1, r_V = 1)
    a1 <- summarise_h2(f1)
    ok <- !is.na(a1$median) && (a1$hi - a1$lo) <= cfg$h2_width_max &&
        a1$lo >= cfg$h2_lower_min
    list(n_unique = n_unique, ok = ok, a1000 = summarise_h2(f1000), a1 = a1,
         a025 = summarise_h2(f025), tightR = summarise_h2(ftR))
}

# Attribution mode: RE-only first, expand to whole region, else none. The headline
# H² (TSV h2_alpha1_* etc.) is whichever mode won.
choose_attribution <- function(h2_re, h2_region, var_cols_re, var_cols_all) {
    if (isTRUE(h2_re$ok)) {
        list(mode = "RE", h2 = h2_re, var_cols = var_cols_re)
    } else if (isTRUE(h2_region$ok)) {
        list(mode = "region", h2 = h2_region, var_cols = var_cols_all)
    } else {
        list(mode = "none",
             h2 = if (!is.na(h2_re$a1$median)) h2_re else h2_region,
             var_cols = integer(0))
    }
}

# ── variant candidates ────────────────────────────────────────────────────────
indicator_for_site <- function(seq_mat, col_idx) {
    # Every character (A/C/G/T/gap/N) is an allele state; the indicator marks the
    # second-most-common one, so gap polymorphisms give a complete 0/1 vector.
    col <- seq_mat[, col_idx]
    if (length(unique(col)) < 2L) return(NULL)
    counts <- sort(table(col), decreasing = TRUE)
    minor <- names(counts)[2]
    col == minor
}

# Variable sites as 0/1 minor-allele indicators; identical or complementary
# partitions are indistinguishable, so they merge into one block.
build_blocks <- function(seq_mat, var_cols, tip_labels, min_carriers) {
    site_indicators <- lapply(var_cols, function(cc) indicator_for_site(seq_mat, cc))
    names(site_indicators) <- as.character(var_cols)
    tip_order <- match(tip_labels, rownames(seq_mat))
    canon_key <- function(col) {
        if (length(col) && !is.na(col[1]) && col[1] == 1L) col <- 1L - col
        paste(ifelse(is.na(col), "?", col), collapse = "")
    }
    blocks <- list()
    block_of <- list()
    for (cc in names(site_indicators)) {
        ind <- as.integer(site_indicators[[cc]][tip_order])
        if (all(is.na(ind)) || length(unique(ind[!is.na(ind)])) < 2L) next
        # Drop singleton / near-private variants: both alleles need >= min_carriers
        # haplotypes, else the "driver" is just one or two samples.
        if (min(sum(ind == 1L, na.rm = TRUE), sum(ind == 0L, na.rm = TRUE)) < min_carriers) next
        k <- canon_key(ind)
        if (!is.null(block_of[[k]])) {
            bi <- block_of[[k]]
            blocks[[bi]]$cols <- c(blocks[[bi]]$cols, cc)
        } else {
            blocks[[length(blocks) + 1L]] <- list(cols = cc, ind = ind)
            block_of[[k]] <- length(blocks)
        }
    }
    blocks
}

# ── forward-selection joint regression ────────────────────────────────────────
# Coverage-weighted WLS: add the block with the largest R² gain until the gain is
# < R2_MIN or p >= 0.05. Select on R² gain, not p: at N~78 a variant adding <2% of
# the variance still has p<0.05. The joint fit drops passengers that track a
# stronger site and keeps nested opposite-sign variants (negative beta).
forward_select <- function(V, y, w_cov, ok, n_blocks, h2_point) {
    ss_tot <- if (sum(ok) >= 4L) {
        sum(w_cov[ok] * (y[ok] - weighted.mean(y[ok], w_cov[ok]))^2)
    } else 0
    selected  <- integer(0)
    remaining <- if (h2_point >= H2_MIN) seq_len(n_blocks) else integer(0)
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
    selected
}

# Final joint fit on the selected variants → conditional β, SE, CI, p per block.
joint_fit <- function(V, y, w_cov, ok, selected) {
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
    joint
}

# ── per-variant rows ──────────────────────────────────────────────────────────
# mean_in/mean_out = marginal carriers vs rest; beta = joint (conditional) effect,
# the headline. descendants = carriers; linked_sites = alignment columns.
variant_rows <- function(blocks, tip_labels, y, w_cov, cov_ok, joint, selected) {
    rows <- list()
    for (j in seq_along(blocks)) {
        b   <- blocks[[j]]
        ind <- b$ind
        carriers  <- tip_labels[which(ind == 1L)]
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
        # but the H² header is still the real fit on the flanked tree.
        branch_df <- as.data.frame(
            matrix(nrow = 0, ncol = length(BRANCH_COLS),
                   dimnames = list(NULL, BRANCH_COLS)))
    }
    # Selected (contributing) variants first by joint p, then the rest by marginal p.
    ord <- order(!branch_df$is_split,
                 ifelse(branch_df$is_split, branch_df$p, Inf),
                 branch_df$p_raw)
    branch_df[ord, BRANCH_COLS]
}

# ── output ────────────────────────────────────────────────────────────────────
# Header values in output order. Every key defaults to NA (or none / 0) so a peak
# that stops early writes the same keys; `...` overrides values and appends keys.
header_vals <- function(peak_id, n_tips, cfg, ...) {
    na_keys <- c("n_unique_seqs",
                 paste0("h2_alpha", rep(c("1000", "1", "025"), each = 3), "_",
                        c("mean", "lower", "upper")),
                 "h2_alpha1_ess", "h2_converged",
                 paste0("h2_re_", c("mean", "lower", "upper", "ok")),
                 paste0("h2_region_", c("mean", "lower", "upper", "ok")))
    vals <- c(list(peak_id = peak_id, n_tips = n_tips),
              setNames(as.list(rep("NA", length(na_keys))), na_keys),
              list(attribution_mode = "none", min_cov = cfg$min_cov, flank = cfg$flank,
                   re_first_col = "NA", re_last_col = "NA",
                   n_variants_tested = 0, n_variants_selected = 0,
                   n_branches_selected = 0, inference_method = "none"))
    over <- list(...)
    vals[names(over)] <- over
    vals
}

header_lines <- function(vals) {
    paste0("# ", names(vals), "=", vapply(vals, function(v) paste0(v), character(1)))
}

write_tsv <- function(path, vals, branch_df) {
    con <- file(path, "w")
    writeLines(header_lines(vals), con)
    write.table(branch_df, con, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
    close(con)
}

write_empty <- function(path, peak_id, n_tips, cfg, reason) {
    vals <- header_vals(peak_id, n_tips, cfg, note = reason)
    writeLines(c(header_lines(vals), paste(BRANCH_COLS, collapse = "\t")), path)
}

h2_vals <- function(prefix, a) {
    setNames(list(round(a$median, 4), round(a$lo, 4), round(a$hi, 4)),
             paste0(prefix, c("_mean", "_lower", "_upper")))
}

# ── main ──────────────────────────────────────────────────────────────────────
main <- function(args = commandArgs(trailingOnly = TRUE)) {
    opt <- parse_args(args)
    stopifnot(!is.null(opt$tree), !is.null(opt$actuation),
              !is.null(opt$alignment),
              !is.null(opt$tsv))
    peak_id <- if (is.null(opt[["peak-id"]])) tools::file_path_sans_ext(basename(opt$tree)) else opt[["peak-id"]]

    # Reproducibility: MCMCglmm is stochastic. Seed from the peak_id so each peak is
    # deterministic across reruns. The RNG order (region H², then RE H²) is part
    # of the output: reordering the fits changes every peak's draws.
    set.seed(sum(utf8ToInt(peak_id)) %% .Machine$integer.max)
    cfg <- load_config(opt)

    tree <- load_tree(opt$tree)
    act <- load_actuation(opt$actuation)
    common <- intersect(tree$tip.label, names(act$trait))
    if (length(common) < 4) {
        write_empty(opt$tsv, peak_id, length(common), cfg,
                    sprintf("only %d tips with actuation data", length(common)))
        return(invisible(NULL))
    }
    trait <- act$trait[common]
    tree  <- keep.tip(tree, common)

    # ── alignment + RE band ──
    aln <- read_alignment(opt$alignment, common)
    seq_mat <- if (length(aln) >= 4) {
        do.call(rbind, strsplit(toupper(unlist(aln)), ""))
    } else {
        NULL
    }
    if (!is.null(seq_mat)) rownames(seq_mat) <- names(aln)
    re_cols <- re_band(toupper(unlist(aln)), cfg$flank)
    re_lo <- if (is.null(re_cols)) 1L else re_cols[1]
    re_hi <- if (is.null(re_cols)) (if (is.null(seq_mat)) 0L else ncol(seq_mat)) else re_cols[2]
    var_cols_all <- variable_columns(seq_mat)
    var_cols_re <- var_cols_all[var_cols_all >= re_lo & var_cols_all <= re_hi]
    cat(sprintf("[phylo_mm] %s RE band: cols %d-%d (flank=%d); %d variable RE sites (of %d total)\n",
                peak_id, re_lo, re_hi, cfg$flank, length(var_cols_re), length(var_cols_all)),
        file = stderr())

    # ── 1. Heritability: RE-only first, expand to whole region ──
    ht <- h2_trait(trait, act$cov_vec, common)
    aln_re <- read_alignment(cfg$re_alignment, common)
    re_tree <- read.tree(cfg$re_tree)
    re_tree$tip.label <- gsub("^_R_", "", re_tree$tip.label)
    re_tree <- tryCatch(keep.tip(re_tree, common), error = function(e) re_tree)

    h2_region <- compute_h2(tree, aln, ht, cfg)
    h2_re     <- compute_h2(re_tree, aln_re, ht, cfg)

    n_unique_seqs <- h2_region$n_unique
    if (n_unique_seqs < 2L) {
        write_empty(opt$tsv, peak_id, length(common), cfg, sprintf(
            "monomorphic alignment — only %d unique sequence across %d tips",
            n_unique_seqs, length(common)))
        return(invisible(NULL))
    }

    attrib <- choose_attribution(h2_re, h2_region, var_cols_re, var_cols_all)
    h2_head <- attrib$h2
    cat(sprintf(
        "[phylo_mm] %s H² gate: RE %s (%.2f [%.2f,%.2f]) | region %s (%.2f [%.2f,%.2f]) -> mode=%s\n",
        peak_id,
        if (isTRUE(h2_re$ok)) "OK" else "fail", h2_re$a1$median, h2_re$a1$lo, h2_re$a1$hi,
        if (isTRUE(h2_region$ok)) "OK" else "fail",
        h2_region$a1$median, h2_region$a1$lo, h2_region$a1$hi,
        attrib$mode), file = stderr())

    # ── 2. Coverage filter + candidate blocks ──
    y <- trait[tree$tip.label]
    # Coverage = WLS weight; cov_used (median-filled) is named by `common` == tips.
    w_cov  <- ht$cov_used[tree$tip.label]
    cov_ok <- !is.na(w_cov) & w_cov >= cfg$min_cov
    cat(sprintf("[phylo_mm] %s coverage filter: %d/%d tips have coverage >= %g\n",
                peak_id, sum(cov_ok), length(cov_ok), cfg$min_cov), file = stderr())
    # `blocks` may be empty (monomorphic RE with flanks); H² is still reported.
    blocks <- build_blocks(seq_mat, attrib$var_cols, tree$tip.label, cfg$min_carriers)

    # ── 3. Forward selection + joint fit ──
    if (length(blocks) > 0L) {
        V <- do.call(cbind, lapply(blocks, function(b) b$ind))
        colnames(V) <- paste0("v", seq_along(blocks))
        ok <- cov_ok & !apply(is.na(V), 1L, any)    # rows usable in the joint model
    } else {
        V  <- matrix(numeric(0), nrow = length(tree$tip.label), ncol = 0)
        ok <- cov_ok
    }
    h2_point <- if (is.na(h2_head$a1$median)) 0 else h2_head$a1$median
    selected <- forward_select(V, y, w_cov, ok, length(blocks), h2_point)
    joint <- joint_fit(V, y, w_cov, ok, selected)
    cat(sprintf("[phylo_mm] %s variant attribution: %d candidate site-block(s), %d selected (forward WLS)\n",
                peak_id, length(blocks), length(selected)),
        file = stderr())

    branch_df <- variant_rows(blocks, tree$tip.label, y, w_cov, cov_ok, joint, selected)

    # ── write ──
    h2_1 <- h2_head$a1
    vals <- do.call(header_vals, c(
        list(peak_id, length(common), cfg, n_unique_seqs = n_unique_seqs),
        h2_vals("h2_alpha1000", h2_head$a1000),
        h2_vals("h2_alpha1", h2_1),
        h2_vals("h2_alpha025", h2_head$a025),
        list(h2_alpha1_ess = if (is.na(h2_1$ess)) "NA" else round(h2_1$ess, 1),
             h2_converged = if (!is.na(h2_1$ess) && h2_1$ess >= ESS_MIN) "TRUE" else "FALSE"),
        # Both H² fits (alpha.V=1) + which one gated attribution. The headline
        # h2_alpha1_* above is whichever mode won (h2_re if RE-driven, else region).
        h2_vals("h2_re", h2_re$a1),
        list(h2_re_ok = if (isTRUE(h2_re$ok)) "TRUE" else "FALSE"),
        h2_vals("h2_region", h2_region$a1),
        list(h2_region_ok = if (isTRUE(h2_region$ok)) "TRUE" else "FALSE",
             attribution_mode = attrib$mode,
             re_first_col = if (cfg$flank > 0L) re_lo else "NA",
             re_last_col  = if (cfg$flank > 0L) re_hi else "NA",
             n_variants_tested = sum(!is.na(branch_df$p_raw)),
             n_variants_selected = length(selected),
             n_branches_selected = length(selected),
             inference_method = "forward_wls",
             nitt = paste0(cfg$nitt, " burnin=", cfg$burnin, " thin=", cfg$thin))
    ))
    write_tsv(opt$tsv, vals, branch_df)
}

main()
