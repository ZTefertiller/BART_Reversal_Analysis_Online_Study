# analysis/plotting/breakpoint_combined.R
# Combined balloon-breakpoint scatter: the reversal paradigm (session 1) and the
# control (session 2) side by side, sharing a common y-axis, with a value-tier
# HEADER stacked directly above each breakpoint panel.
#
# Header (consistent Arial fonts, centred over each region):
#   * top:    bold "Pre-reversal" / "Post-reversal" labels (one per half);
#             "Control" over the control panel.
#   * below:  one entry per value tier — a balloon-coloured swatch, the tier name
#             ("<Low/Medium/High>-value"), and (mode = "avg" only) the average
#             breakpoint. Tiers ordered Low -> Medium -> High, side by side,
#             centred over their half (pre = trials 1-90, post = 91-180). Control
#             shows a single High-value entry.
#
# Two modes:
#   mode = "avg"   (default) header shows swatch + tier name + "Avg breakpoint: N".
#   mode = "lines" header shows swatch + tier name only; each colour/condition's
#                  average breakpoint is drawn as a light, colour-coded dotted
#                  horizontal line on the panel (3 pre, 3 post, 1 control).
#
# A dotted vertical line marks trial 91 on the reversal panel.
# Returns a patchwork object.

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(patchwork)
})

# Value-tier definitions per phase (which balloon colour is which tier). The
# reversal swaps blue<->orange between halves; yellow stays medium throughout.
.BP_PRE_TIERS  <- c(low = "o", medium = "y", high = "b")   # pre-reversal
.BP_POST_TIERS <- c(low = "b", medium = "y", high = "o")   # post-reversal

# Build one header panel from a set of tier entries spanning x in [x0, x1].
# show_avg: include the "Avg breakpoint: N" line under the tier name.
.bp_header_panel <- function(entries, top_labs, xlim, family, base_size,
                             show_avg = TRUE, text_size = NULL) {
  sw   <- base_size * 0.48          # swatch point size
  # All header text (phase heading AND tier/value labels) at one bold size.
  ts_h <- if (is.null(text_size)) base_size * 0.52 else text_size
  ts_t <- ts_h                       # tier/value labels match the heading size

  # vertical layout (arbitrary canvas; heading near top, tiers below).
  if (show_avg) {
    y_head <- 9.4; y_sw <- 7.0; y_tier <- 5.7; y_avg <- 4.2
    ylo <- 3.4
  } else {
    y_head <- 9.4; y_sw <- 6.4; y_tier <- 4.9; y_avg <- NA
    ylo <- 4.0
  }

  df_sw  <- do.call(rbind, lapply(entries, function(e)
    data.frame(x = e$x_centre, y = y_sw, col = e$swatch_col)))
  df_tier <- do.call(rbind, lapply(entries, function(e)
    data.frame(x = e$x_centre, y = y_tier, lab = e$tier_lab)))

  p <- ggplot() +
    geom_text(data = top_labs, aes(x = x, y = y_head, label = label),
              family = family, fontface = "bold", size = ts_h, hjust = 0.5) +
    geom_point(data = df_sw, aes(x = x, y = y, colour = col), size = sw,
               shape = 16, show.legend = FALSE) +
    scale_colour_identity() +
    geom_text(data = df_tier, aes(x = x, y = y, label = lab),
              family = family, fontface = "bold", size = ts_t, hjust = 0.5)

  if (show_avg) {
    df_avg <- do.call(rbind, lapply(entries, function(e)
      data.frame(x = e$x_centre, y = y_avg, lab = e$avg)))
    p <- p + geom_text(data = df_avg, aes(x = x, y = y, label = lab),
                       family = family, fontface = "bold", size = ts_t,
                       hjust = 0.5)
  }

  p +
    scale_x_continuous(limits = xlim, expand = expansion(add = 1)) +
    scale_y_continuous(limits = c(ylo, 10.2)) +
    labs(x = NULL, y = NULL) +
    theme_void(base_family = family) +
    theme(plot.margin = margin(2, 2, 0, 2))
}

# bal_cols: named vector with Blue/Orange/Yellow/Pink hex codes (BAL_COLS).
# mode: "avg" (tier name + avg breakpoint in header) or "lines" (tier name only;
#       avg breakpoints drawn as colour-coded dotted lines on the panels).
breakpoint_combined <- function(data, bal_cols,
                                family = "Arial", base_size = 13,
                                mode = c("avg", "lines")) {
  mode <- match.arg(mode)
  col_hex <- c(b = unname(bal_cols["Blue"]), o = unname(bal_cols["Orange"]),
               y = unname(bal_cols["Yellow"]), p = unname(bal_cols["Pink"]))

  bp_s1 <- data %>%
    filter(session == 1) %>%
    distinct(trial_number, .keep_all = TRUE) %>%
    select(trial_number, balloon_color, optimal_inflations) %>%
    arrange(trial_number)

  bp_s2 <- data %>%
    filter(session == 2) %>%
    distinct(trial_number, .keep_all = TRUE) %>%
    select(trial_number, balloon_color, optimal_inflations) %>%
    arrange(trial_number)

  # ---- per-tier average breakpoint (rows = phase-filtered slice of bp_s1) ----
  .tier_avg <- function(rows, colour_code) {
    v <- rows$optimal_inflations[rows$balloon_color == colour_code]
    mean(v[is.finite(v)])
  }
  tier_label <- c(low = "Low-value", medium = "Medium-value", high = "High-value")

  .make_entries <- function(rows, tiers, centres) {
    lapply(seq_along(tiers), function(i) {
      cc <- tiers[[i]]
      a  <- .tier_avg(rows, cc)
      list(x_centre   = centres[i],
           swatch_col = col_hex[[cc]],
           colour_code = cc,
           avg_val    = a,
           tier_lab   = tier_label[[names(tiers)[i]]],
           avg        = sprintf("Avg breakpoint: %d", round(a)))
    })
  }

  pre_centres  <- c(14, 44, 73)
  post_centres <- c(104, 133, 161)
  pre_rows  <- bp_s1[bp_s1$trial_number < 91, ]
  post_rows <- bp_s1[bp_s1$trial_number > 90, ]
  pre_entries  <- .make_entries(pre_rows,  .BP_PRE_TIERS,  pre_centres)
  post_entries <- .make_entries(post_rows, .BP_POST_TIERS, post_centres)

  # ---- breakpoint panels ----------------------------------------------------
  rev_shading <- list(
    annotate("rect", xmin=0.5,   xmax=30.5,  ymin=-Inf, ymax=Inf, fill="grey50",  alpha=0.12),
    annotate("rect", xmin=30.5,  xmax=50.5,  ymin=-Inf, ymax=Inf, fill="#2596be", alpha=0.09),
    annotate("rect", xmin=50.5,  xmax=70.5,  ymin=-Inf, ymax=Inf, fill="#E07B30", alpha=0.09),
    annotate("rect", xmin=70.5,  xmax=90.5,  ymin=-Inf, ymax=Inf, fill="#C8A800", alpha=0.09),
    annotate("rect", xmin=90.5,  xmax=120.5, ymin=-Inf, ymax=Inf, fill="grey50",  alpha=0.12),
    annotate("rect", xmin=120.5, xmax=140.5, ymin=-Inf, ymax=Inf, fill="#2596be", alpha=0.09),
    annotate("rect", xmin=140.5, xmax=160.5, ymin=-Inf, ymax=Inf, fill="#E07B30", alpha=0.09),
    annotate("rect", xmin=160.5, xmax=180.5, ymin=-Inf, ymax=Inf, fill="#C8A800", alpha=0.09)
  )

  yr   <- range(c(bp_s1$optimal_inflations, bp_s2$optimal_inflations), na.rm = TRUE)
  ylim <- yr + diff(yr) * c(-0.02, 0.05)
  rev_xlim  <- c(1, 180)
  ctrl_xlim <- c(1, 30)

  # Header text size (geom_text units); axis titles + tick labels use the same
  # visual size (converted to pt) and bold, so ALL text matches the headings.
  hdr_size <- base_size * 0.52
  ax_size  <- hdr_size * ggplot2::.pt
  base_theme <- theme(
    text         = element_text(family = family),
    axis.title.x = element_text(face = "bold", family = family, size = ax_size),
    axis.title.y = element_text(face = "bold", family = family, size = ax_size),
    axis.text.x  = element_text(face = "bold", family = family, size = ax_size),
    axis.text.y  = element_text(face = "bold", family = family, size = ax_size))

  # average-breakpoint dotted lines (mode = "lines"): one per tier per half,
  # segment spanning that half; colour-coded by balloon colour (shares the
  # points' scale_color_manual via the balloon_color code).
  # scatter-dot opacity: the "lines" figure uses slightly more translucent dots
  # so the colour-coded average lines read clearly through them.
  pt_alpha <- if (mode == "lines") 0.7 else 1.0

  avg_segs <- NULL
  if (mode == "lines") {
    seg_rows <- function(entries, x0, x1)
      do.call(rbind, lapply(entries, function(e)
        data.frame(x = x0, xend = x1, y = e$avg_val,
                   balloon_color = e$colour_code)))
    avg_df <- rbind(seg_rows(pre_entries, 0.5, 90.5),
                    seg_rows(post_entries, 90.5, 180.5))
    avg_segs <- geom_segment(
      data = avg_df,
      aes(x = x, xend = xend, y = y, yend = y, colour = balloon_color),
      linetype = "dotted", linewidth = 2.2, alpha = 0.7, inherit.aes = FALSE)
  }

  p_rev <- ggplot(bp_s1, aes(x = trial_number, y = optimal_inflations)) +
    rev_shading +
    geom_vline(xintercept = 91, linetype = "dotted", color = "grey20",
               linewidth = 1.8) +
    avg_segs +
    geom_point(aes(color = balloon_color), size = 5.5, alpha = pt_alpha, shape = 16) +
    scale_color_manual(values = col_hex) +
    scale_x_continuous(breaks = c(1, 30, 50, 70, 90, 120, 140, 160, 180),
                       expand = expansion(add = 1)) +
    coord_cartesian(xlim = rev_xlim, ylim = ylim) +
    labs(x = "Trial", y = "Breakpoint") +
    theme_minimal(base_size = base_size, base_family = family) +
    base_theme +
    theme(legend.position = "none")

  ctrl_avg <- mean(bp_s2$optimal_inflations[is.finite(bp_s2$optimal_inflations)])
  ctrl_avg_seg <- if (mode == "lines")
    list(geom_segment(data = data.frame(x = 0.5, xend = 30.5, y = ctrl_avg),
                      aes(x = x, xend = xend, y = y, yend = y),
                      colour = col_hex[["p"]], linetype = "dotted",
                      linewidth = 2.2, alpha = 0.7, inherit.aes = FALSE))
    else NULL

  p_ctrl <- ggplot(bp_s2, aes(x = trial_number, y = optimal_inflations)) +
    annotate("rect", xmin=0.5, xmax=30.5, ymin=-Inf, ymax=Inf,
             fill = unname(bal_cols["Pink"]), alpha = 0.08) +
    ctrl_avg_seg +
    geom_point(color = unname(bal_cols["Pink"]), size = 5.5, alpha = pt_alpha,
               shape = 16) +
    scale_x_continuous(breaks = c(1, 10, 20, 30), expand = expansion(add = 0.5)) +
    coord_cartesian(xlim = ctrl_xlim, ylim = ylim) +
    labs(x = "Trial", y = NULL) +
    theme_minimal(base_size = base_size, base_family = family) +
    base_theme +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

  # ---- headers --------------------------------------------------------------
  show_avg <- (mode == "avg")
  rev_top_labs <- data.frame(
    x = c(mean(c(1, 90)), mean(c(91, 180))),
    label = c("Pre-reversal", "Post-reversal"))
  rev_header <- .bp_header_panel(c(pre_entries, post_entries), rev_top_labs,
                                 rev_xlim, family, base_size, show_avg = show_avg,
                                 text_size = hdr_size)

  ctrl_entries <- list(list(
    x_centre = mean(ctrl_xlim), swatch_col = col_hex[["p"]],
    tier_lab = "High-value",
    avg = sprintf("Avg breakpoint: %d", round(ctrl_avg))))
  ctrl_top_labs <- data.frame(x = mean(ctrl_xlim), label = "Control")
  ctrl_header <- .bp_header_panel(ctrl_entries, ctrl_top_labs, ctrl_xlim,
                                  family, base_size, show_avg = show_avg,
                                  text_size = hdr_size)

  # ---- assemble: header over panel, per column, then columns side by side ---
  hr <- if (show_avg) 1 else 0.8
  rev_col  <- rev_header  / p_rev  + plot_layout(heights = c(hr, 2.3))
  ctrl_col <- ctrl_header / p_ctrl + plot_layout(heights = c(hr, 2.3))

  (rev_col | ctrl_col) + plot_layout(widths = c(14, 4.5))
}
