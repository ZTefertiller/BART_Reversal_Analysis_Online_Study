# analysis/plotting/ewmv_pre_post_spearman.R
# Spearman analogues of the Steiger/Williams pre- vs post-reversal trait tests,
# built from the lightweight modeling/ewmv/ewmv_by_block.csv (no Stan objects).
#
# Two pieces:
#  (1) ewmv_steiger_auto():  Steiger/Williams dependent-correlation dumbbell +
#      table on the pre-specified significant EWMV param x trait pairs PLUS the
#      updating-exponent (xi) associations. Pre- vs post-reversal only. Same look
#      as ewmv_steiger_byblock.
#
#  (2) ewmv_bootstrap_prepost():  the Spearman equivalent of the Steiger test.
#      Resample participants (BCa bootstrap), recompute the pre and post Spearman
#      correlations AND their difference (post - pre) per replicate, and form
#      BCa CIs for r_pre, r_post, and the difference. The difference CI excluding
#      0 is the Spearman analogue of a significant dependent-correlation test.
#      Returns a dumbbell (matched style) + a stats table.
#
# Both run on the same fixed pair set by default (.PPS_SIG_SET).

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(cocor); library(boot)
})

# .ech_cor_long / .ECH_PARAMS / .ECH_Q (and .cor_test) drive the data-driven pair
# selection below so it matches the cut-down heatmap exactly.
source(here::here("analysis", "plotting", "ewmv_cutdown_heatmap.R"))

.PPS_PARAMS <- c(phi = "Prior Weight (ψ)", eta = "Updating Exponent (ξ)",
                 rho = "Risk Preference (ρ)", tau = "Inverse Temperature (τ)",
                 lambda = "Loss Aversion (λ)")
.PPS_TRAITS <- c(caps_total = "CAPS Total", pdi_total = "PDI Total",
                 spq_total  = "SPQ Total",  spq_cogper = "SPQ Cog-Per",
                 spq_interp = "SPQ Interp", spq_disorg = "SPQ Disorg")

# The pre-specified significant pairs (tau x interp / PDI / SPQ, lambda x SPQ)
# plus the updating-exponent (xi) associations. This is the fixed set shown in
# the "all significant associations" figures/tables.
.PPS_SIG_SET <- tibble::tribble(
  ~param,   ~trait,
  "eta",    "pdi_total",
  "eta",    "spq_total",
  "eta",    "spq_cogper",
  "eta",    "spq_interp",
  "tau",    "spq_interp",
  "tau",    "pdi_total",
  "tau",    "spq_total",
  "lambda", "spq_total"
)

# ---- data-driven pair selection --------------------------------------------
# Returns the (param, trait) pairs whose association is significant in pre OR
# post reversal in modeling/ewmv/ewmv_by_block.csv, using the SAME parameter set,
# questionnaire set and Holm-per-(param, q)-across-conditions scheme as the
# cut-down heatmap. adjust = "holm" (default; matches the heatmap) or "none"
# (nominal p). Pairs are selected from the current fits.
ewmv_prepost_sig_pairs <- function(
    by_block_csv = here::here("data_published", "bart_rl_online_ewmv_params.csv"),
    adjust = c("holm", "none"), alpha = 0.05, method = "pearson") {
  adjust <- match.arg(adjust)
  d <- utils::read.csv(by_block_csv)
  if (!"participant_id" %in% names(d) && "sub_id" %in% names(d)) d$participant_id <- d$sub_id
  long <- dplyr::bind_rows(
    .ech_cor_long(d, "pre_",  "pre",     method),
    .ech_cor_long(d, "post_", "post",    method),
    .ech_cor_long(d, "con_",  "control", method)) %>%
    dplyr::group_by(param, q) %>%
    dplyr::mutate(p_use = if (adjust == "holm") p.adjust(p, "holm") else p) %>%
    dplyr::ungroup()
  sel <- long %>%
    dplyr::filter(cond %in% c("pre", "post"), p_use < alpha) %>%
    dplyr::distinct(param, q) %>%
    dplyr::transmute(param, trait = q)
  if (nrow(sel) == 0L) return(sel)
  sel %>% dplyr::arrange(match(param, names(.PPS_PARAMS)),
                         match(trait, names(.PPS_TRAITS)))
}

.pps_pair_factor <- function(df) {
  lv <- df %>% distinct(param, trait) %>%
    mutate(ord = match(param, names(.PPS_PARAMS)) * 10 +
                 match(trait, names(.PPS_TRAITS))) %>%
    arrange(ord) %>%
    mutate(lab = sprintf("%s × %s", .PPS_PARAMS[param], .PPS_TRAITS[trait]))
  factor(sprintf("%s × %s", .PPS_PARAMS[df$param], .PPS_TRAITS[df$trait]),
         levels = rev(lv$lab))
}

.pps_dumbbell <- function(rows, base_size, family) {
  long <- rows %>%
    select(pair_label, r_pre, r_post) %>%
    pivot_longer(c(r_pre, r_post), names_to = "phase", values_to = "r") %>%
    mutate(phase = factor(recode(phase, r_pre = "Pre reversal",
                                 r_post = "Post reversal"),
                          levels = c("Pre reversal", "Post reversal")))
  star <- rows %>% filter(diff_sig) %>% mutate(r_star = pmax(r_pre, r_post) + 0.045)

  # ALL text — axis title, parameter pair labels, r tick labels, legend — at one
  # uniform black bold Arial size (`txt`) for maximum readability.
  txt <- base_size
  ggplot(rows, aes(y = pair_label)) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.5) +
    geom_segment(aes(x = r_pre, xend = r_post, yend = pair_label),
                 color = "grey55", linewidth = 1.6) +
    geom_point(data = long, aes(x = r, color = phase), size = 6) +
    geom_text(data = star, aes(x = r_star, label = "*"), size = 11,
              vjust = 0.75, color = "black") +
    scale_color_manual(values = c("Pre reversal" = "#2630F5",
                                  "Post reversal" = "#E68D33"), name = NULL) +
    labs(x = "Correlation with trait (r)", y = NULL) +
    theme_minimal(base_size = base_size, base_family = family) +
    theme(legend.position   = "bottom",
          legend.text       = element_text(size = txt, colour = "black",
                                            family = family),
          panel.grid.major.y = element_blank(),
          panel.grid.minor   = element_blank(),
          axis.text.y  = element_text(size = txt, colour = "black",
                                      family = family),
          axis.text.x  = element_text(size = txt, colour = "black",
                                      family = family),
          axis.title.x = element_text(size = txt, face = "bold", colour = "black",
                                      family = family),
          axis.title.y = element_text(size = txt, face = "bold", colour = "black",
                                      family = family))
}

# ---- (1) Steiger/Williams on the fixed significant pair set (Pearson) -------
ewmv_steiger_auto <- function(
    by_block_csv = here::here("data_published", "bart_rl_online_ewmv_params.csv"),
    base_size = 20, family = "Arial",
    pairs = ewmv_prepost_sig_pairs(by_block_csv)) {
  d <- utils::read.csv(by_block_csv)
  if (!"participant_id" %in% names(d) && "sub_id" %in% names(d)) d$participant_id <- d$sub_id
  if (nrow(pairs) == 0L)
    return(list(plot = NULL, table = NULL, results = pairs[0, ], pairs = pairs))

  rows <- purrr::pmap_dfr(pairs, function(param, trait) {
    xp <- d[[paste0("pre_",  param)]]; xq <- d[[paste0("post_", param)]]
    y  <- d[[trait]]
    ok <- is.finite(xp) & is.finite(xq) & is.finite(y)
    xp <- xp[ok]; xq <- xq[ok]; y <- y[ok]; n <- length(y)
    r_pre <- cor(y, xp); r_post <- cor(y, xq); r_pp <- cor(xp, xq)
    ct <- suppressWarnings(cocor::cocor.dep.groups.overlap(
      r.jk = r_pre, r.jh = r_post, r.kh = r_pp, n = n,
      alternative = "two.sided", test = c("williams1959", "steiger1980")))
    tibble(param = param, trait = trait, n = n,
           r_pre = r_pre, r_post = r_post, r_pp = r_pp,
           p_williams = ct@williams1959$p.value,
           p_steiger  = ct@steiger1980$p.value)
  }) %>%
    mutate(diff_sig = p_steiger < 0.05,
           pair_label = .pps_pair_factor(.))

  p <- .pps_dumbbell(rows, base_size, family)

  tbl <- rows %>%
    transmute(Parameter = unname(.PPS_PARAMS[param]),
              Trait     = unname(.PPS_TRAITS[trait]),
              n,
              `r pre`       = round(r_pre, 2),
              `r post`      = round(r_post, 2),
              `r(pre,post)` = round(r_pp, 2),
              `Williams p`  = signif(p_williams, 2),
              `Steiger p`   = signif(p_steiger, 2))

  list(plot = p, table = tbl, results = rows, pairs = pairs)
}

# ---- (2) Bootstrap Spearman pre/post + difference (BCa) ---------------------
# pairs: data.frame(param, trait); defaults to the fixed significant pair set.
ewmv_bootstrap_prepost <- function(
    by_block_csv = here::here("data_published", "bart_rl_online_ewmv_params.csv"),
    pairs = .PPS_SIG_SET, R = 5000, base_size = 20, family = "Arial",
    seed = 42) {
  d <- utils::read.csv(by_block_csv)
  if (!"participant_id" %in% names(d) && "sub_id" %in% names(d)) d$participant_id <- d$sub_id

  set.seed(seed)
  rows <- purrr::pmap_dfr(pairs, function(param, trait) {
    xp <- d[[paste0("pre_",  param)]]; xq <- d[[paste0("post_", param)]]
    y  <- d[[trait]]
    ok <- is.finite(xp) & is.finite(xq) & is.finite(y)
    sub <- data.frame(xp = xp[ok], xq = xq[ok], y = y[ok])
    n <- nrow(sub)

    stat_fn <- function(dat, idx) {
      s <- dat[idx, ]
      r_pre  <- suppressWarnings(cor(s$y, s$xp, method = "spearman"))
      r_post <- suppressWarnings(cor(s$y, s$xq, method = "spearman"))
      c(r_pre = r_pre, r_post = r_post, diff = r_post - r_pre)
    }
    bo <- boot::boot(sub, stat_fn, R = R)

    bca <- function(i) {
      ci <- tryCatch(boot::boot.ci(bo, type = "bca", index = i),
                     error = function(e) NULL)
      if (is.null(ci) || is.null(ci$bca)) c(NA_real_, NA_real_)
      else ci$bca[4:5]
    }
    pre_ci  <- bca(1); post_ci <- bca(2); diff_ci <- bca(3)

    tibble(param = param, trait = trait, n = n,
           r_pre = bo$t0[["r_pre"]], r_post = bo$t0[["r_post"]],
           diff = bo$t0[["diff"]],
           pre_lo = pre_ci[1],  pre_hi = pre_ci[2],
           post_lo = post_ci[1], post_hi = post_ci[2],
           diff_lo = diff_ci[1], diff_hi = diff_ci[2])
  }) %>%
    mutate(
      # "significant" pre/post difference = BCa CI of (post - pre) excludes 0
      diff_sig = is.finite(diff_lo) & is.finite(diff_hi) &
                 (diff_lo > 0 | diff_hi < 0),
      pair_label = .pps_pair_factor(.))

  p <- .pps_dumbbell(rows, base_size, family)

  tbl <- rows %>%
    transmute(Parameter = unname(.PPS_PARAMS[param]),
              Trait     = unname(.PPS_TRAITS[trait]),
              n,
              `rho pre`  = sprintf("%.2f [%.2f, %.2f]", r_pre, pre_lo, pre_hi),
              `rho post` = sprintf("%.2f [%.2f, %.2f]", r_post, post_lo, post_hi),
              `Δ (post − pre)` = sprintf("%.2f [%.2f, %.2f]", diff,
                                         diff_lo, diff_hi),
              `CI excl. 0` = ifelse(diff_sig, "yes", "no"))

  list(plot = p, table = tbl, results = rows, pairs = pairs)
}
