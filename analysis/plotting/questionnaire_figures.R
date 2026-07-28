# analysis/plotting/questionnaire_figures.R
# Reusable questionnaire density plots and summary-stat cards, driven by
# analysis/plotting/questionnaire_meta.R (consistent colours/labels everywhere).

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(ggrepel)
})
source(here::here("analysis", "plotting", "questionnaire_meta.R"))

# Per-participant long table of normalised scores for the chosen vars.
.q_long <- function(per_id, vars, norm = c("theoretical", "observed")) {
  norm <- match.arg(norm)
  meta <- Q_INFO[match(vars, Q_INFO$var), ]
  per_id %>%
    select(all_of(vars)) %>%
    mutate(.row = dplyr::row_number()) %>%
    pivot_longer(all_of(vars), names_to = "var", values_to = "value") %>%
    left_join(meta, by = "var") %>%
    group_by(var) %>%
    mutate(
      norm = if (norm == "theoretical" && all(is.finite(tmax))) {
        (value - tmin) / pmax(1e-9, tmax - tmin)
      } else {
        (value - min(value, na.rm = TRUE)) /
          pmax(1e-9, max(value, na.rm = TRUE) - min(value, na.rm = TRUE))
      },
      is_factor = label %in% SPQ_FACTOR_LABELS_DENS
    ) %>%
    ungroup() %>%
    filter(is.finite(norm), norm >= 0, norm <= 1) %>%
    mutate(label = factor(label, levels = meta$label))
}

# Density overlay. label_pos = "top" or "bottom" (labels in the margin on that
# side, with leader lines to each curve's peak). caption = NULL to omit.
q_density <- function(per_id, vars, norm = "theoretical", label_pos = "top",
                      caption = NA, base_size = 14, line_size = 2.0,
                      family = "Arial", label_size = 4.6) {
  ql <- .q_long(per_id, vars, norm)
  pk <- ql %>% group_by(label) %>%
    summarise(
      d_x = { d <- density(norm, from = 0, to = 1, n = 512); d$x[which.max(d$y)] },
      d_y = { d <- density(norm, from = 0, to = 1, n = 512); max(d$y) },
      is_factor = first(is_factor), .groups = "drop")
  ymax <- max(pk$d_y, na.rm = TRUE)

  if (identical(label_pos, "bottom")) {
    repel_ylim <- c(NA, -0.04 * ymax)
    coord_ylim <- c(-0.45 * ymax, ymax * 1.04)
  } else {
    repel_ylim <- c(ymax * 1.05, NA)
    coord_ylim <- c(0, ymax * 1.5)
  }

  p <- ggplot(ql, aes(x = norm, colour = label, linetype = is_factor)) +
    geom_density(linewidth = line_size) +
    ggrepel::geom_label_repel(
      data = pk, aes(x = d_x, y = d_y, label = label, colour = label),
      fill = "white", size = label_size, fontface = "bold", family = family,
      label.padding = unit(0.26, "lines"), label.size = 0.45,
      box.padding = unit(0.4, "lines"), point.padding = unit(0.2, "lines"),
      min.segment.length = 0, segment.size = 0.55, segment.colour = "grey45",
      ylim = repel_ylim, direction = "both", force = 3, max.overlaps = Inf,
      show.legend = FALSE, seed = 42, inherit.aes = FALSE) +
    scale_colour_manual(values = Q_PAL) +
    scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "longdash"),
                          guide = "none") +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                       labels = scales::percent_format(accuracy = 1)) +
    coord_cartesian(ylim = coord_ylim, clip = "off") +
    labs(x = "Score (proportion of scale range)", y = "Density",
         caption = if (length(caption) == 1 && is.na(caption)) NULL else caption) +
    theme_minimal(base_size = base_size, base_family = family) +
    theme(legend.position = "none", panel.grid.minor = element_blank(),
          plot.caption = element_text(hjust = 0, colour = "grey40"))
  p
}

# Colour-coded stats card. nested=TRUE indents SPQ subscales under their factor.
# family/base_text control the font; column headers (Mean (SD)/Observed/Possible/
# Scale) and their cells underneath are all black, only the scale name (left
# column) keeps its family colour.
q_stats_card <- function(per_id, vars, nested = FALSE, base_text = 4.4,
                         family = "Arial", row_gap = 1) {
  meta <- Q_INFO[match(vars, Q_INFO$var), ]
  rows <- lapply(seq_along(vars), function(i) {
    x <- per_id[[vars[i]]]
    pos <- if (is.finite(meta$tmax[i])) sprintf("%g–%g", meta$tmin[i], meta$tmax[i]) else "—"
    data.frame(label = meta$label[i], color = meta$color[i],
               mean_sd = sprintf("%.0f (%.0f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE)),
               obs = sprintf("%g–%g", min(x, na.rm = TRUE), max(x, na.rm = TRUE)),
               poss = pos, scl = meta$scale[i],
               indent = if (nested && grepl("^SPQ ", meta$label[i]) &&
                            !meta$label[i] %in% c("SPQ", "SPQ-B")) 0.18 else 0,
               stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, rows)
  # row_gap > 1 spreads the rows further apart vertically (more breathing room).
  df$y <- rev(seq_len(nrow(df))) * row_gap
  ny <- nrow(df) * row_gap
  # Column x-positions spread proportionally to base_text so larger fonts don't
  # collide. The label column starts at 0; each data column is pushed right by
  # `sp` (which grows with the font) past a wide label gutter.
  sp  <- base_text / 4.4               # spread factor relative to the old default
  # Wider label gutter (mean_sd pushed right) so long labels like "SPQ Cog-Per"
  # clear the first data column, and wider inter-column gaps so the italic
  # headers ("Mean (SD)", "Observed", ...) don't touch at large fonts.
  xs  <- c(mean_sd = 2.9, obs = 4.5, poss = 5.95, scl = 7.55) * sp
  x0  <- -0.25 * sp                    # left edge of the bounding box / labels
  xr  <- xs["scl"] + 1.05 * sp         # right edge of the bounding box
  ggplot(df) +
    annotate("rect", xmin = x0, xmax = xr, ymin = 0.4, ymax = ny + 1.3,
             fill = "white", colour = "grey70", linewidth = 0.5) +
    annotate("text", x = xs["mean_sd"], y = ny + 0.85, label = "Mean (SD)",
             fontface = "italic", size = base_text - 1, colour = "black", family = family) +
    annotate("text", x = xs["obs"], y = ny + 0.85, label = "Observed",
             fontface = "italic", size = base_text - 1, colour = "black", family = family) +
    annotate("text", x = xs["poss"], y = ny + 0.85, label = "Possible",
             fontface = "italic", size = base_text - 1, colour = "black", family = family) +
    annotate("text", x = xs["scl"], y = ny + 0.85, label = "Scale",
             fontface = "italic", size = base_text - 1, colour = "black", family = family) +
    geom_text(aes(x = x0 + 0.1 * sp + indent, y = y, label = label, colour = label),
              hjust = 0, fontface = "bold", size = base_text, family = family) +
    geom_text(aes(x = xs["mean_sd"], y = y, label = mean_sd), hjust = 0.5,
              size = base_text - 0.3, colour = "black", family = family) +
    geom_text(aes(x = xs["obs"], y = y, label = obs), hjust = 0.5,
              size = base_text - 0.3, colour = "black", family = family) +
    geom_text(aes(x = xs["poss"], y = y, label = poss), hjust = 0.5,
              size = base_text - 0.3, colour = "black", family = family) +
    geom_text(aes(x = xs["scl"], y = y, label = scl), hjust = 0.5,
              size = base_text - 0.6, colour = "black", family = family) +
    scale_colour_manual(values = Q_PAL, guide = "none") +
    coord_cartesian(xlim = c(x0 - 0.05 * sp, xr + 0.05 * sp),
                    ylim = c(0.3, ny + 1.5), clip = "off") +
    theme_void(base_family = family)
}
