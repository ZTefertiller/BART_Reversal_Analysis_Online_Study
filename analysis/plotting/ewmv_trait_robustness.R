# analysis/plotting/ewmv_trait_robustness.R
# ===========================================================================
# Robustness / stress-test of the EWMV trait-parameter -> questionnaire
# Pearson correlations reported in the paper (the cut-down heatmap built by
# analysis/plotting/ewmv_cutdown_heatmap.R::ewmv_cutdown_heatmap()).
#
# This script is ADDITIVE: it reuses the published pipeline's exact definitions
# and helpers and never modifies them. It does NOT recompute the model; it uses
# the full EWMV fits (mcmc/model_fits_list.rds) and the intermediate
# modeling/ewmv/ewmv_by_block.csv that the heatmap reads.
#
# Pipeline reused verbatim:
#   - parameters .ECH_PARAMS, questionnaires .ECH_Q, prefixes pre_/post_/con_
#     (from ewmv_cutdown_heatmap.R)
#   - Pearson .cor_test (complete.cases, exact = FALSE) (from param_q_heatmap.R)
#   - Holm correction per (param, q) across the 3 conditions (m = 3)
#   - participant set / ordering = the CSV rows (== fit draws column order)
#
# Two checks, for every correlation significant (Holm) in the FULL sample:
#   1) LEAVE-ONE-PARTICIPANT-OUT jackknife. Drop each participant in turn,
#      recompute the focal correlation on n-1, and -- to test whether the
#      reported significance is fragile -- re-run the SAME Holm correction
#      (all 3 conditions for that param x q, recomputed on n-1) and read off the
#      focal cell's Holm-adjusted p. Report full-sample r, min/max LOO r,
#      whether Holm significance ever flips on any single removal, and the most
#      influential participant (largest |r_full - r_LOO|).
#   2) POSTERIOR PROPAGATION. Using the subject-level posterior draws from the
#      cmdstanr fits, recompute the correlation at EACH draw (the participant
#      parameter vector at that draw vs the fixed trait scores) and report the
#      median and 95% interval of the posterior r.
#
# A correlation "survives" iff it stays Holm-significant across ALL LOO drops
# AND its posterior-r 95% interval excludes 0; otherwise it is "fragile".
#
# Public entry point: ewmv_trait_robustness(...) -> list(table, summary, ...).
# ===========================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(tibble)
})

source(here::here("analysis", "plotting", "param_q_heatmap.R"))      # .cor_test, .param_cols
source(here::here("analysis", "plotting", "ewmv_cutdown_heatmap.R")) # .ECH_PARAMS/.ECH_Q/.ech_cor_long
source(here::here("analysis", "plotting", "ewmv_trait_robustness_verify.R")) # ewmv_trait_long

# Map condition label -> CSV column prefix and -> fit name in model_fits_list.
.ROB_COND <- tibble::tribble(
  ~cond,            ~prefix,  ~fit_pat,
  "pre reversal",   "pre_",   "^ewmv.*blue",
  "post reversal",  "post_",  "^ewmv.*orange",
  "control",        "con_",   "^ewmv.*pink"
)

# ---------------------------------------------------------------------------
# Helper: recompute the per-(param,q) Holm-adjusted p for ONE focal cell from a
# data.frame `d` that carries pre_/post_/con_ <param> + the questionnaire cols.
# Returns the focal condition's Holm-adjusted p (the exact quantity .ech_plot
# uses to decide significance). Used inside the LOO loop on n-1 rows.
.holm_p_for_cell <- function(d, param, q, focal_prefix, method = "pearson") {
  ps <- vapply(.ROB_COND$prefix, function(pf) {
    .cor_test(d[[q]], d[[paste0(pf, param)]], method = method)$p
  }, numeric(1))
  padj <- p.adjust(ps, method = "holm")          # m = 3, per (param, q)
  unname(padj[match(focal_prefix, .ROB_COND$prefix)])
}

# ---------------------------------------------------------------------------
# 1) Leave-one-participant-out jackknife for a single focal cell.
.loo_one <- function(d, param, q, cond, method = "pearson", alpha = 0.05) {
  pf <- .ROB_COND$prefix[.ROB_COND$cond == cond]
  pc <- paste0(pf, param)

  x_all <- d[[q]]; y_all <- d[[pc]]
  ok    <- complete.cases(x_all, y_all)          # participants used by this cell
  full  <- .cor_test(x_all, y_all, method = method)

  idx <- which(ok)
  n   <- length(idx)
  loo <- vapply(idx, function(i) {
    keep <- setdiff(seq_len(nrow(d)), i)
    .cor_test(d[[q]][keep], d[[pc]][keep], method = method)$r
  }, numeric(1))

  # Holm-adjusted p for the focal cell after each single removal, recomputed
  # over the same 3-condition family the paper uses.
  loo_padj <- vapply(idx, function(i) {
    .holm_p_for_cell(d[setdiff(seq_len(nrow(d)), i), , drop = FALSE],
                     param, q, pf, method = method)
  }, numeric(1))

  most_inf <- idx[which.max(abs(loo - full$r))]

  tibble(
    param = param, q = q, cond = cond,
    n = n,
    r_full = full$r,
    p_holm_full = .holm_p_for_cell(d, param, q, pf, method = method),
    loo_r_min = min(loo), loo_r_max = max(loo),
    # significance "flips" if it was significant in full sample but any single
    # removal pushes the Holm-adjusted p to >= alpha (or vice-versa).
    loo_padj_max = max(loo_padj),
    loo_sig_flips = any(loo_padj >= alpha),
    most_influential_pid = d$participant_id[most_inf],
    most_influential_delta_r = (loo - full$r)[which.max(abs(loo - full$r))],
    .loo_r = list(setNames(loo, d$participant_id[idx]))
  )
}

# ---------------------------------------------------------------------------
# Pull subject-level posterior draws (n_draws x n_subj) for one parameter from
# a CmdStanMCMC fit, columns ordered by subject index (== CSV row order).
.subject_draws <- function(fit, param) {
  if (inherits(fit, "CmdStanMCMC") || inherits(fit, "CmdStanGQ")) {
    draws <- fit$draws(variables = param, format = "df")
  } else if (is.data.frame(fit)) {
    draws <- fit
  } else if (is.list(fit) && !is.null(fit$draws_df)) {
    draws <- fit$draws_df
  } else stop("Unrecognised fit object for .subject_draws()")
  cols <- .param_cols(draws, param)              # ordered phi[1], phi[2], ...
  as.matrix(draws[, cols, drop = FALSE])
}

# ---------------------------------------------------------------------------
# 2) Posterior propagation for a single focal cell. For each posterior draw,
# correlate the subject parameter vector (that draw) with the fixed trait
# scores; summarise the resulting posterior distribution of r.
.post_one <- function(draws_mat, trait, param, q, cond, method = "pearson",
                      probs = c(0.025, 0.5, 0.975)) {
  ok <- is.finite(trait)
  tr <- trait[ok]
  dm <- draws_mat[, ok, drop = FALSE]
  # constant-trait guard (matches .cor_test's degenerate-input handling)
  rvec <- if (length(unique(tr)) < 2 || sd(tr) == 0) {
    rep(NA_real_, nrow(dm))
  } else {
    apply(dm, 1L, function(pr) suppressWarnings(stats::cor(pr, tr, method = method)))
  }
  qs <- stats::quantile(rvec, probs = probs, na.rm = TRUE)
  tibble(
    param = param, q = q, cond = cond,
    post_r_median = unname(qs[2]),
    post_r_lo = unname(qs[1]), post_r_hi = unname(qs[3]),
    post_r_excludes_0 = (qs[1] > 0 & qs[3] > 0) | (qs[1] < 0 & qs[3] < 0),
    n_draws = sum(is.finite(rvec))
  )
}

# ===========================================================================
# MAIN
# ===========================================================================
ewmv_trait_robustness <- function(
    by_block_csv = here::here("data_published", "bart_rl_online_ewmv_params.csv"),
    fits_rds     = here::here("mcmc", "model_fits_list.rds"),
    draws_rds    = here::here("mcmc", "ewmv_subject_draws.rds"),
    method = "pearson", alpha = 0.05, verbose = TRUE) {

  # ---- reproduce the published matrix + identify the significant set -------
  long <- ewmv_trait_long(by_block_csv, method)
  sig <- long %>%
    filter(!is.na(p_adj), p_adj < alpha) %>%
    transmute(param, q, cond,
              prefix = .ROB_COND$prefix[match(cond, .ROB_COND$cond)]) %>%
    arrange(cond, q, param)

  if (verbose)
    cat(sprintf("Significant (Holm p<%.2f) full-sample correlations: %d\n",
                alpha, nrow(sig)))

  d <- utils::read.csv(by_block_csv)
  if (!"participant_id" %in% names(d) && "sub_id" %in% names(d)) d$participant_id <- d$sub_id

  # ---- (1) LOO jackknife for every significant cell ------------------------
  loo_tbl <- pmap_dfr(list(sig$param, sig$q, sig$cond),
                      ~ .loo_one(d, ..1, ..2, ..3, method = method, alpha = alpha))

  # ---- (2) posterior propagation -------------------------------------------
  # Subject-level draws come from the compact per-condition draws file when it
  # is present (mcmc/ewmv_subject_draws.rds, ~47 MB, shipped with the release);
  # otherwise they are pulled from the full 3.1 GB combined fits
  # (mcmc/model_fits_list.rds). Both yield identical subject-ordered draws.
  if (file.exists(draws_rds)) {
    if (verbose) cat("Loading compact subject draws for posterior propagation...\n")
    cdraws <- readRDS(draws_rds)                    # named list keyed by condition
    fit_for_cond <- function(cond) cdraws[[cond]]
  } else {
    if (verbose) cat("Loading full fits for posterior propagation...\n")
    fits <- readRDS(fits_rds)
    get_fit <- function(pat) {
      i <- grep(pat, names(fits), ignore.case = TRUE)
      if (length(i) != 1) stop("fit pattern '", pat, "' matched ", length(i))
      fits[[i]]
    }
    fit_for_cond <- function(cond) get_fit(.ROB_COND$fit_pat[.ROB_COND$cond == cond])
  }
  # cache subject-draw matrices per (cond, param) so we only pull each once
  needed <- distinct(sig, cond, param)
  draw_cache <- new.env(parent = emptyenv())
  for (k in seq_len(nrow(needed))) {
    cd <- needed$cond[k]; pm <- needed$param[k]
    assign(paste(cd, pm, sep = "||"),
           .subject_draws(fit_for_cond(cd), pm), envir = draw_cache)
  }

  post_tbl <- pmap_dfr(list(sig$param, sig$q, sig$cond), function(pm, qq, cd) {
    dm <- get(paste(cd, pm, sep = "||"), envir = draw_cache)
    # fixed trait scores aligned to CSV/draw column order
    trait <- d[[qq]]
    stopifnot(length(trait) == ncol(dm))
    .post_one(dm, trait, pm, qq, cd, method = method)
  })

  # ---- combine into one table per correlation ------------------------------
  tbl <- loo_tbl %>%
    select(-.loo_r) %>%
    left_join(post_tbl, by = c("param", "q", "cond")) %>%
    mutate(
      param_label = unname(.ECH_PARAMS[param]),
      q_label     = unname(.ECH_Q[q]),
      survives    = (!loo_sig_flips) & post_r_excludes_0
    ) %>%
    relocate(cond, q_label, param_label) %>%
    arrange(cond, q_label, param_label)

  # tuck the raw LOO r vectors alongside for any later drill-down
  loo_r_vectors <- setNames(loo_tbl$.loo_r,
                            paste(loo_tbl$cond, loo_tbl$q, loo_tbl$param,
                                  sep = " | "))

  summary_txt <- .robustness_summary(tbl)
  if (verbose) {
    cat("\n", strrep("=", 70), "\n", sep = "")
    cat("ROBUSTNESS TABLE (LOO jackknife + posterior propagation)\n")
    cat(strrep("=", 70), "\n", sep = "")
    print(.robustness_display(tbl), row.names = FALSE)
    cat("\n", summary_txt, "\n", sep = "")
  }

  invisible(list(table = tbl, long = long, sig = sig,
                 loo_r_vectors = loo_r_vectors, summary = summary_txt))
}

# Compact display version of the table (rounded, readable column names).
.robustness_display <- function(tbl) {
  tbl %>%
    transmute(
      Condition = cond,
      Trait = q_label,
      Parameter = param_label,
      n,
      `r (full)` = sprintf("%+.3f", r_full),
      `p_holm` = signif(p_holm_full, 3),
      `LOO r min` = sprintf("%+.3f", loo_r_min),
      `LOO r max` = sprintf("%+.3f", loo_r_max),
      `LOO max p_holm` = signif(loo_padj_max, 3),
      `Sig flips?` = ifelse(loo_sig_flips, "YES", "no"),
      `Most infl. (Δr)` = sprintf("%s (%+.3f)", most_influential_pid,
                                  most_influential_delta_r),
      `Post r [95%]` = sprintf("%+.3f [%+.3f, %+.3f]",
                               post_r_median, post_r_lo, post_r_hi),
      `Excl 0?` = ifelse(post_r_excludes_0, "yes", "NO"),
      Survives = ifelse(survives, "SURVIVES", "fragile")
    ) %>%
    as.data.frame()
}

# Short text summary of which survive vs which are fragile (+ why).
.robustness_summary <- function(tbl) {
  surv <- tbl %>% filter(survives)
  frag <- tbl %>% filter(!survives)
  lab <- function(r) sprintf("%s x %s (%s)", r$q_label, r$param_label, r$cond)
  out <- c(sprintf("SUMMARY: %d of %d significant correlations survive both checks.",
                   nrow(surv), nrow(tbl)))
  if (nrow(surv))
    out <- c(out, "SURVIVE (Holm-sig across all LOO drops AND posterior-r 95% excludes 0):",
             paste0("  - ", purrr::map_chr(seq_len(nrow(surv)),
                                           ~ lab(surv[.x, ]))))
  if (nrow(frag)) {
    reasons <- purrr::map_chr(seq_len(nrow(frag)), function(i) {
      r <- frag[i, ]
      why <- c(if (r$loo_sig_flips) "LOO flips Holm-sig" else NULL,
               if (!r$post_r_excludes_0) "posterior-r 95% includes 0" else NULL)
      sprintf("  - %s: %s", lab(r), paste(why, collapse = "; "))
    })
    out <- c(out, "FRAGILE:", reasons)
  }
  paste(out, collapse = "\n")
}

# ---------------------------------------------------------------------------
# CONCISE OUTPUTS for the end of the report.
# ---------------------------------------------------------------------------

# Compact at-a-glance table: one row per significant correlation with just the
# full r, the LOO r range, the posterior-r 95% interval, and the verdict.
ewmv_robustness_compact <- function(tbl) {
  tbl %>%
    arrange(desc(survives), cond, q_label, param_label) %>%
    transmute(
      Correlation = sprintf("%s × %s", q_label, param_label),
      Condition   = sub(" reversal", "", cond),           # pre / post / control
      `r`         = sprintf("%+.2f", r_full),
      `LOO r range` = sprintf("%+.2f to %+.2f", loo_r_min, loo_r_max),
      `Posterior r [95%]` = sprintf("%+.2f [%+.2f, %+.2f]",
                                    post_r_median, post_r_lo, post_r_hi),
      Verdict = ifelse(survives, "Survives", "Fragile")
    ) %>%
    as.data.frame()
}

# Single summary forest figure: each significant correlation on its own row,
# point = full-sample r, thick segment = leave-one-out r range, thin segment =
# posterior-r 95% interval; colour = survive (solid) vs fragile (faded).
# Rows grouped survivors-first, ordered by |r|. Returns a ggplot.
ewmv_robustness_forest <- function(tbl, base_size = 13, family = "Arial") {
  d <- tbl %>%
    mutate(
      lab = sprintf("%s × %s  (%s)", q_label, param_label,
                    sub(" reversal", "", cond)),
      verdict = ifelse(survives, "Survives both checks", "Fragile")
    ) %>%
    arrange(survives, abs(r_full)) %>%        # bottom -> top: fragile then survivors
    mutate(lab = factor(lab, levels = lab))

  pal <- c("Survives both checks" = "#1B7837", "Fragile" = "#B2ABD2")

  ggplot(d, aes(y = lab, colour = verdict)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    # posterior 95% interval (thin)
    geom_linerange(aes(xmin = post_r_lo, xmax = post_r_hi), linewidth = 0.9,
                   alpha = 0.55) +
    # leave-one-out r range (thick)
    geom_linerange(aes(xmin = loo_r_min, xmax = loo_r_max), linewidth = 2.6,
                   alpha = 0.9) +
    # full-sample r (point)
    geom_point(aes(x = r_full), size = 3.4) +
    scale_colour_manual(values = pal, name = NULL,
                        breaks = c("Survives both checks", "Fragile")) +
    scale_x_continuous(limits = c(-0.4, 0.5), breaks = seq(-0.4, 0.4, 0.2)) +
    labs(x = "Pearson r (point = full sample; thick = LOO range; thin = posterior 95%)",
         y = NULL) +
    theme_minimal(base_size = base_size, base_family = family) +
    theme(
      text = element_text(family = family, colour = "black"),
      axis.text.y = element_text(size = base_size, colour = "black"),
      axis.text.x = element_text(size = base_size, colour = "black"),
      axis.title.x = element_text(size = base_size, face = "bold"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.text = element_text(size = base_size)
    )
}

if (sys.nframe() == 0) {
  res <- ewmv_trait_robustness()
}
