# analysis/plotting/bayesian_contrasts_variants.R
# Layout variants of the Bayesian group-level contrast forest from
# bayesian_block_contrasts.R (.plot_contrast_forest), in the analysis-report
# style: theme_classic, dashed 0-line, horizontal 95% CI bars, point at the
# posterior mean, x-axis "Posterior mean difference (95% credible interval)".
#
# Input is the contrasts tibble returned by run_bayesian_contrasts_ewmv() (one
# row per parameter x contrast, columns label_html / contrast / post_mean /
# ci_lower / ci_upper). No new headers / footnotes / titles are added.
#
#   contrasts_forest_tworow(): faceted, wraps to two rows (matches the report).
#   contrasts_forest_onerow(): the same panels forced onto a single row.
#   contrasts_boxes_nolabels(): every parameter as its own panel, NO contrast
#       labels on the y-axis (clean boxes to drag-and-drop).
#   contrasts_labels_only(): the three contrast labels as a standalone strip,
#       same vertical order/spacing as the boxes, to drop beside them.

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(ggtext); library(patchwork)
})

# y-axis order top->bottom (same as .plot_contrast_forest)
.BCV_CLEVELS <- c("Control − Post-reversal",
                  "Control − Pre-reversal",
                  "Post-reversal − Pre-reversal")

.bcv_prep <- function(contrasts) {
  contrasts %>%
    mutate(
      contrast   = factor(contrast, levels = .BCV_CLEVELS),
      label_html = factor(label_html, levels = unique(label_html))
    )
}

.bcv_xlab <- function(prob = 0.95)
  sprintf("Posterior mean difference (%d%% credible interval)", round(prob * 100))

.bcv_axis_theme <- function(base_size, family) {
  theme(
    text         = element_text(family = family),
    axis.title.x = element_text(face = "bold", family = family, size = 12),
    axis.title.y = element_text(face = "bold", family = family, size = 12)
  )
}

# ---- (1) two-row faceted forest (matches the analysis report) ---------------
contrasts_forest_tworow <- function(contrasts, prob = 0.95,
                                    base_size = 16, family = "Arial",
                                    nrow = 2) {
  pd <- .bcv_prep(contrasts)
  ggplot(pd, aes(x = post_mean, y = contrast)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper),
                   height = 0.25, linewidth = 0.9, colour = "black") +
    geom_point(size = 3, colour = "black") +
    facet_wrap(~ label_html, scales = "free_x", nrow = nrow) +
    labs(x = .bcv_xlab(prob), y = NULL) +
    theme_classic(base_size = base_size, base_family = family) +
    .bcv_axis_theme(base_size, family) +
    theme(strip.text  = ggtext::element_markdown(size = base_size,
                                                 face = "bold"),
          strip.background = element_blank(),
          axis.text.y = element_text(size = base_size - 2, family = family),
          axis.text.x = element_text(size = base_size - 3, family = family))
}

# ---- (2) one-row faceted forest --------------------------------------------
contrasts_forest_onerow <- function(contrasts, prob = 0.95,
                                    base_size = 16, family = "Arial") {
  pd <- .bcv_prep(contrasts)
  ggplot(pd, aes(x = post_mean, y = contrast)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper),
                   height = 0.25, linewidth = 0.9, colour = "black") +
    geom_point(size = 3, colour = "black") +
    facet_wrap(~ label_html, scales = "free_x", nrow = 1) +
    labs(x = .bcv_xlab(prob), y = NULL) +
    theme_classic(base_size = base_size, base_family = family) +
    .bcv_axis_theme(base_size, family) +
    theme(strip.text  = ggtext::element_markdown(size = base_size,
                                                 face = "bold"),
          strip.background = element_blank(),
          axis.text.y = element_text(size = base_size - 2, family = family),
          axis.text.x = element_text(size = base_size - 3, family = family))
}

# Parameter (group_param) display order for the individual boxes, with the
# inverse temperature (mu_tau) LAST.
.BCV_BOX_ORDER <- c("mu_phi", "mu_eta", "mu_rho", "mu_lambda", "mu_tau")

# A "nice" step (1, 2, 2.5 or 5 x 10^k) just below `target`.
.bcv_nice_step <- function(target) {
  if (!is.finite(target) || target <= 0) return(1)
  p    <- 10^floor(log10(target))
  cand <- c(1, 2, 2.5, 5, 10) * p
  cand <- cand[cand <= target]
  if (!length(cand)) p else max(cand)
}

# Symmetric x-limits centred on zero (so each box is visually balanced about 0).
.bcv_box_xlim <- function(pd_one) {
  half <- max(abs(c(pd_one$ci_lower, pd_one$ci_upper, 0)), na.rm = TRUE)
  if (!is.finite(half) || half == 0) half <- 1
  c(-half, half) * 1.12
}

# Exactly THREE x-axis ticks: one below zero, zero, one above (symmetric). The
# single non-zero step is a nice value that fits inside the data half-extent.
.bcv_box_breaks <- function(pd_one) {
  half <- max(abs(c(pd_one$ci_lower, pd_one$ci_upper, 0)), na.rm = TRUE)
  if (!is.finite(half) || half == 0) return(c(-1, 0, 1))
  # pick the largest nice step that still sits inside the data range
  step <- .bcv_nice_step(half)
  for (s in c(step, step * 0.5, step * 0.4, step * 0.2, step * 0.1)) {
    if (s <= half) { step <- s; break }
  }
  c(-step, 0, step)
}

# trim trailing zeros so labels read e.g. "0.01", "0", "-0.01"
.bcv_fmt <- function(x) {
  out <- formatC(x, format = "f", digits = 6)
  out <- sub("0+$", "", out)
  out <- sub("\\.$", "", out)
  out[as.numeric(x) == 0] <- "0"
  out
}

# ---- (3) individual boxes, NO y-axis labels, FULL border per panel ----------
# One self-contained box PER parameter (so each can be saved/dragged on its
# own): full rectangular outline, parameter name on top (strip), the three
# contrasts with the y-axis labels stripped, dashed 0-line. No x-axis title
# (that lives on the contrast-label strip from contrasts_labels_only()).
.bcv_one_box <- function(pd_one, base_size, family, strip = TRUE) {
  p <- ggplot(pd_one, aes(x = post_mean, y = contrast)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50",
               linewidth = 0.9) +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper),
                   height = 0.25, linewidth = 1.6, colour = "black") +
    geom_point(size = 4.4, colour = "black") +
    scale_x_continuous(breaks = .bcv_box_breaks(pd_one), labels = .bcv_fmt) +
    coord_cartesian(xlim = .bcv_box_xlim(pd_one)) +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = base_size, base_family = family) +
    .bcv_axis_theme(base_size, family) +
    theme(
      strip.background = element_blank(),
      axis.line    = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.4),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x  = element_text(size = base_size, family = family,
                                  face = "bold"),
      plot.margin  = margin(8, 10, 8, 10)
    )
  # strip = TRUE keeps the parameter name on top (default standalone box);
  # strip = FALSE drops it (for the combined grid where the name sits on the
  # paired line graph above instead, so it isn't duplicated).
  if (strip) {
    p + facet_wrap(~ label_html) +
      theme(strip.text = ggtext::element_markdown(size = base_size + 3,
                                                  face = "bold"))
  } else {
    p
  }
}

# Returns a NAMED list of individual per-parameter bordered boxes (names are the
# group_params, e.g. mu_phi), each a standalone drag-and-drop image. Inverse
# temperature (mu_tau) is ordered last. strip = FALSE omits the parameter-name
# strip (used by the combined sweep+contrast grid).
contrasts_boxes_nolabels <- function(contrasts, base_size = 20,
                                     family = "Arial", strip = TRUE) {
  pd <- .bcv_prep(contrasts)
  params <- intersect(.BCV_BOX_ORDER, unique(as.character(pd$parameter)))
  params <- c(params, setdiff(unique(as.character(pd$parameter)), params))
  setNames(lapply(params, function(gp) {
    .bcv_one_box(dplyr::filter(pd, parameter == gp), base_size, family, strip)
  }), params)
}

# Convenience: lay the individual boxes out in one row for inline display.
contrasts_boxes_row_display <- function(boxes) {
  patchwork::wrap_plots(boxes, nrow = 1)
}

# ---- (4) standalone contrast-label strip (drag beside the boxes) -----------
# The three contrast names in the SAME top->bottom order/spacing as the box
# panels, plus the "Posterior mean difference (95% credible interval)" x-axis
# title at the SAME size as the contrast labels (so it can serve as the shared
# axis title when dropped next to the boxes).
contrasts_labels_only <- function(prob = 0.95,
                                   base_size = 20, family = "Arial") {
  lab_df <- data.frame(
    contrast = factor(.BCV_CLEVELS, levels = .BCV_CLEVELS),
    x = 0)
  xlab_df <- data.frame(x = 0, label = .bcv_xlab(prob))
  # left-anchor all text so the longest string (the axis title) never clips.
  ggplot(lab_df, aes(x = x, y = contrast)) +
    geom_text(aes(label = contrast), hjust = 0, family = family,
              size = base_size / .pt) +
    # the x-axis title, same size as the contrast labels, beneath them
    geom_text(data = xlab_df,
              aes(x = x, y = 0.35, label = label),
              hjust = 0, family = family, fontface = "bold",
              size = base_size / .pt, inherit.aes = FALSE) +
    scale_x_continuous(limits = c(0, 1), expand = expansion(add = c(0.02, 0))) +
    scale_y_discrete(expand = expansion(add = c(1.0, 0.6))) +
    labs(x = NULL, y = NULL) +
    theme_void(base_size = base_size, base_family = family) +
    theme(text = element_text(family = family),
          plot.margin = margin(8, 4, 8, 4))
}
