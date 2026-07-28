# ================================================================
# analysis/plotting/bayesian_block_contrasts.R
#
# Bayesian group-level contrasts for block comparisons
# (Blue vs Orange, Blue vs Pink, Orange vs Pink)
#
# RATIONALE: Comparing individual-level posterior means with t-tests
# inflates Type I error (Böhm et al.) because it discards posterior
# uncertainty. The correct approach is to contrast the group-level
# hyperparameter posteriors (mu_* from the Stan hierarchy) directly.
# Each MCMC draw gives one plausible value for the group mean; the
# distribution of draw-wise differences IS the posterior of the
# contrast — no frequentist test required.
#
# Since blue/orange/pink are fit in separate Stan runs their posteriors
# are independent, so pairing draws from each fit gives the exact
# posterior of (condition_B - condition_A) under independence.
# ================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggtext)
  library(patchwork)
})

# ----------------------------------------------------------------
# Core contrast engine
# ----------------------------------------------------------------

#' Extract full posterior draws for group-level parameters from three
#' block fits and compute draw-wise differences for every pair.
#'
#' @param fit_blue,fit_orange,fit_pink  CmdStanFit / stanfit objects
#' @param group_params  character vector of variable names to extract
#'   (must be scalars present in each fit, e.g. "mu_phi")
#' @param prob  credible-interval width (default 0.95)
#' @param param_labels  optional named character vector mapping
#'   group_params → human-readable labels (HTML OK)
#'
#' @return tibble with one row per parameter × contrast
.bayesian_group_contrasts <- function(fit_blue, fit_orange, fit_pink,
                                      group_params,
                                      prob         = 0.95,
                                      param_labels = NULL) {

  stopifnot(length(group_params) >= 1)

  # ---- extract posterior draws (format="df" drops chain/iter cols) ----
  draws_b <- fit_blue$draws(variables   = group_params, format = "df")
  draws_o <- fit_orange$draws(variables = group_params, format = "df")
  draws_p <- fit_pink$draws(variables   = group_params, format = "df")

  # align draw counts (trim to the smallest chain × iter product)
  n <- min(nrow(draws_b), nrow(draws_o), nrow(draws_p))
  draws_b <- draws_b[seq_len(n), , drop = FALSE]
  draws_o <- draws_o[seq_len(n), , drop = FALSE]
  draws_p <- draws_p[seq_len(n), , drop = FALSE]

  alpha <- (1 - prob) / 2

  # Experimentally meaningful condition labels. The task ran in this order:
  # pre-reversal (blue) first, post-reversal (orange) after the contingency swap,
  # control (pink, single fixed contingency) last. Every contrast is the LATER
  # condition minus the EARLIER one, so positive = parameter higher later in the task.
  cond_lab <- c(Blue = "Pre-reversal", Orange = "Post-reversal", Pink = "Control")

  pairs <- list(
    list(a = "Blue",   fit_a = draws_b, b = "Orange", fit_b = draws_o),
    list(a = "Blue",   fit_a = draws_b, b = "Pink",   fit_b = draws_p),
    list(a = "Orange", fit_a = draws_o, b = "Pink",   fit_b = draws_p)
  )

  purrr::map_dfr(group_params, function(gp) {
    # human-readable label (strip HTML/entities for console; keep raw for plots)
    label_html <- if (!is.null(param_labels) && !is.na(param_labels[gp]))
      unname(param_labels[gp]) else gp
    label_plain <- gsub("<[^>]+>",  "",  label_html)
    label_plain <- gsub("&[^;]+;",  "",  label_plain)
    label_plain <- trimws(label_plain)

    purrr::map_dfr(pairs, function(pr) {
      a_vec <- pr$fit_a[[gp]]
      b_vec <- pr$fit_b[[gp]]

      if (is.null(a_vec) || is.null(b_vec)) return(NULL)

      delta <- b_vec - a_vec   # posterior of (B - A); positive = B > A

      tibble::tibble(
        parameter      = gp,
        label_plain    = label_plain,
        label_html     = label_html,
        contrast       = paste0(cond_lab[pr$b], " \u2212 ", cond_lab[pr$a]),  # later − earlier
        cond_a         = pr$a,
        cond_b         = pr$b,
        n_draws        = n,
        post_mean      = mean(delta),
        post_median    = median(delta),
        ci_lower       = quantile(delta, alpha,     names = FALSE),
        ci_upper       = quantile(delta, 1 - alpha, names = FALSE),
        p_b_gt_a       = mean(delta > 0),   # P(B > A | data)
        p_direction    = pmax(mean(delta > 0), mean(delta < 0))
      )
    })
  })
}

# ----------------------------------------------------------------
# Credible-interval annotation helper (replaces sig bars)
# ----------------------------------------------------------------

#' Build a small annotation data frame for adding Bayesian credible-
#' interval text to a spaghetti/violin plot.
#'
#' @param contrasts  output of .bayesian_group_contrasts()
#' @param prob       CI width matching the contrasts table
.build_ci_annotations <- function(contrasts, prob = 0.95) {
  prob_pct <- round(prob * 100)
  contrasts %>%
    mutate(
      direction_label = case_when(
        p_direction >= 0.99 ~ "***",
        p_direction >= 0.95 ~ "**",
        p_direction >= 0.90 ~ "*",
        TRUE ~ ""
      ),
      annotation = sprintf(
        "%s: \u0394=%.3g [%.3g, %.3g] P\u209c=%.2f %s",
        contrast,
        post_mean, ci_lower, ci_upper,
        p_direction, direction_label
      )
    )
}

# ----------------------------------------------------------------
# Forest plot of contrasts for one model
# ----------------------------------------------------------------

#' @param contrasts  tibble from .bayesian_group_contrasts()
#' @param model_tag  character label for title, e.g. "EWMV"
#' @param prob       CI width (for axis label)
.plot_contrast_forest <- function(contrasts, model_tag, prob = 0.95) {

  prob_pct <- round(prob * 100)
  # y-axis top->bottom: reversal effect first, then control comparisons
  contrast_levels <- c("Control \u2212 Post-reversal",
                       "Control \u2212 Pre-reversal",
                       "Post-reversal \u2212 Pre-reversal")

  pd <- contrasts %>%
    mutate(
      contrast   = factor(contrast, levels = contrast_levels),
      label_html = factor(label_html, levels = unique(label_html))
    )

  ggplot(pd, aes(x = post_mean, y = contrast)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(
      aes(xmin = ci_lower, xmax = ci_upper),
      height = 0.25, linewidth = 0.8, colour = "black"
    ) +
    geom_point(size = 2.5, colour = "black") +
    facet_wrap(~ label_html, scales = "free_x") +
    labs(
      title = sprintf(
        "%s \u2014 Bayesian group-level contrasts",
        model_tag
      ),
      x = sprintf("Posterior mean difference (%d%% credible interval)", prob_pct),
      y = NULL
    ) +
    theme_classic() +
    theme(
      plot.title    = ggtext::element_markdown(face = "bold"),
      strip.text    = ggtext::element_markdown(),
      axis.text.y   = element_text(size = 9)
    )
}

# ----------------------------------------------------------------
# Console summary printer
# ----------------------------------------------------------------

.print_contrasts_table <- function(contrasts, model_tag, prob = 0.95) {
  prob_pct <- round(prob * 100)
  cat("\n================================================\n")
  cat(sprintf("Bayesian group-level contrasts | Model: %s\n", toupper(model_tag)))
  cat(sprintf("%d%% credible intervals | P_direction = P(contrast in dominant direction)\n", prob_pct))
  cat("(*** p_dir >= 0.99 | ** >= 0.95 | * >= 0.90)\n")
  cat("================================================\n\n")

  out <- contrasts %>%
    mutate(
      stars = case_when(
        p_direction >= 0.99 ~ "***",
        p_direction >= 0.95 ~ "**",
        p_direction >= 0.90 ~ "*",
        TRUE ~ ""
      )
    ) %>%
    select(
      Parameter  = label_plain,
      Contrast   = contrast,
      Mean_diff  = post_mean,
      CI_lower   = ci_lower,
      CI_upper   = ci_upper,
      P_B_gt_A   = p_b_gt_a,
      P_direction= p_direction,
      Evidence   = stars
    )

  print(as.data.frame(out), row.names = FALSE, digits = 3)
  cat("\n")
  invisible(contrasts)
}

# ----------------------------------------------------------------
# Model-specific wrappers
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# Table formatter (kableExtra) for use in Rmd
# ----------------------------------------------------------------

#' Format a combined contrasts tibble as a publication-ready kableExtra table.
#'
#' @param contrasts_list  named list of tibbles from run_bayesian_contrasts_*()
#' @param prob            CI width (for column header)
#' @param caption         table caption string
.format_contrasts_kable <- function(contrasts_list, prob = 0.95, caption = NULL) {
  prob_pct <- round(prob * 100)

  combined <- dplyr::bind_rows(contrasts_list, .id = "Model") %>%
    dplyr::mutate(
      Model = toupper(Model),
      `P direction` = sprintf("%.3f%s", p_direction,
                              dplyr::case_when(
                                p_direction >= 0.99 ~ "***",
                                p_direction >= 0.95 ~ "**",
                                p_direction >= 0.90 ~ "*",
                                TRUE ~ ""
                              )),
      CrI = sprintf("[%.3g, %.3g]", ci_lower, ci_upper)
    ) %>%
    dplyr::select(
      Model,
      Parameter   = label_plain,
      Contrast    = contrast,
      `Post. mean` = post_mean,
      !!sprintf("%d%% CrI", prob_pct) := CrI,
      `P direction`
    )

  if (is.null(caption))
    caption <- sprintf(
      paste0(
        "Bayesian group-level block contrasts (%d%% credible intervals). ",
        "P direction = posterior probability that the contrast is in the ",
        "reported direction (*** \u2265 0.99, ** \u2265 0.95, * \u2265 0.90)."
      ),
      prob_pct
    )

  combined %>%
    kableExtra::kbl(
      digits  = 3,
      caption = caption,
      align   = c("l", "l", "l", "r", "c", "r")
    ) %>%
    kableExtra::kable_styling(
      bootstrap_options = c("striped", "hover", "condensed"),
      full_width = TRUE
    ) %>%
    kableExtra::collapse_rows(columns = 1:2, valign = "top") %>%
    kableExtra::footnote(
      general = paste0(
        "Contrasts derived from draw-wise differences of group-level ",
        "hyperparameter posteriors (\u03bc\u03b8) across independently fitted ",
        "hierarchical Stan models per block condition."
      ),
      general_title = "Note: ",
      footnote_as_chunk = TRUE
    )
}

# ----------------------------------------------------------------
#' Run Bayesian group-level contrasts for the EWMV model.
#' Group-level parameters are in generated quantities: mu_phi, mu_eta,
#' mu_rho, mu_tau, mu_lambda.
#'
#' @param ewmv_blue,ewmv_orange,ewmv_pink  CmdStanFit objects
#' @param out_dir  directory for saved plots (created if needed)
#' @param prob  credible-interval width
#' @return invisibly, the contrasts tibble
run_bayesian_contrasts_ewmv <- function(ewmv_blue, ewmv_orange, ewmv_pink,
                                        out_dir = here::here("reporting", "figures",
                                                             "param_evolution",
                                                             "bayesian_contrasts", "ewmv"),
                                        prob = 0.95) {

  group_params <- c("mu_phi", "mu_eta", "mu_rho", "mu_tau", "mu_lambda")

  param_labels <- c(
    mu_phi    = "Prior Weight (&psi;)",
    mu_eta    = "Updating Exponent (&xi;)",
    mu_rho    = "Risk Preference (&rho;)",
    mu_tau    = "Inverse Temperature (&tau;)",
    mu_lambda = "Loss Aversion (&lambda;)"
  )

  contrasts <- .bayesian_group_contrasts(
    ewmv_blue, ewmv_orange, ewmv_pink,
    group_params  = group_params,
    prob          = prob,
    param_labels  = param_labels
  )

  .print_contrasts_table(contrasts, "EWMV", prob)

  p <- .plot_contrast_forest(contrasts, "EWMV", prob)
  print(p)

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(out_dir, "ewmv_bayesian_contrasts.png"), p,
         width = 12, height = 6, dpi = 600, bg = "white")
  ggsave(file.path(out_dir, "ewmv_bayesian_contrasts.svg"), p,
         width = 12, height = 6, bg = "white", device = "svg")

  invisible(contrasts)
}


#' Run Bayesian group-level contrasts for the STL model.
#' Group-level parameters live directly in parameters{}: mu_vwin,
#' mu_vloss, mu_beta, mu_omegaone.
#'
#' @param stl_blue,stl_orange,stl_pink  CmdStanFit objects
#' @param out_dir  directory for saved plots (created if needed)
#' @param prob  credible-interval width
#' @return invisibly, the contrasts tibble
run_bayesian_contrasts_stl <- function(stl_blue, stl_orange, stl_pink,
                                       out_dir = here::here("reporting", "figures",
                                                            "param_evolution",
                                                            "bayesian_contrasts", "stl"),
                                       prob = 0.95) {

  group_params <- c("mu_vwin", "mu_vloss", "mu_beta", "mu_omegaone")

  param_labels <- c(
    mu_vwin     = "Reward Learning Rate (v<sub>win</sub>)",
    mu_vloss    = "Loss Learning Rate (v<sub>loss</sub>)",
    mu_beta     = "Behavioural Consistency (&beta;)",
    mu_omegaone = "Initial Pump Target (&omega;<sub>1</sub>)"
  )

  contrasts <- .bayesian_group_contrasts(
    stl_blue, stl_orange, stl_pink,
    group_params  = group_params,
    prob          = prob,
    param_labels  = param_labels
  )

  .print_contrasts_table(contrasts, "STL", prob)

  p <- .plot_contrast_forest(contrasts, "STL", prob)
  print(p)

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(out_dir, "stl_bayesian_contrasts.png"), p,
         width = 11, height = 5, dpi = 600, bg = "white")
  ggsave(file.path(out_dir, "stl_bayesian_contrasts.svg"), p,
         width = 11, height = 5, bg = "white", device = "svg")

  invisible(contrasts)
}


#' Run Bayesian group-level contrasts for the four-parameter model.
#' Group-level parameters are in generated quantities: mu_phi, mu_eta,
#' mu_gam, mu_tau.
#'
#' @param fourpar_blue,fourpar_orange,fourpar_pink  CmdStanFit objects
#' @param out_dir  directory for saved plots (created if needed)
#' @param prob  credible-interval width
#' @return invisibly, the contrasts tibble
run_bayesian_contrasts_fourpar <- function(fourpar_blue, fourpar_orange, fourpar_pink,
                                           out_dir = here::here("reporting", "figures",
                                                                "param_evolution",
                                                                "bayesian_contrasts", "fourpar"),
                                           prob = 0.95) {

  group_params <- c("mu_phi", "mu_eta", "mu_gam", "mu_tau")

  param_labels <- c(
    mu_phi = "Prior Weight (&phi;)",
    mu_eta = "Updating Exponent (&eta;)",
    mu_gam = "Risk/Curvature (&gamma;)",
    mu_tau = "Inverse Temperature (&tau;)"
  )

  contrasts <- .bayesian_group_contrasts(
    fourpar_blue, fourpar_orange, fourpar_pink,
    group_params  = group_params,
    prob          = prob,
    param_labels  = param_labels
  )

  .print_contrasts_table(contrasts, "4PAR", prob)

  p <- .plot_contrast_forest(contrasts, "4PAR", prob)
  print(p)

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(out_dir, "fourpar_bayesian_contrasts.png"), p,
         width = 11, height = 5, dpi = 600, bg = "white")
  ggsave(file.path(out_dir, "fourpar_bayesian_contrasts.svg"), p,
         width = 11, height = 5, bg = "white", device = "svg")

  invisible(contrasts)
}
