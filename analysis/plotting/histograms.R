# analysis/plotting/histograms.R
# Purpose: Generate histograms of the average inflations per color 
# (Blue, Orange, Pink) for each participant. Horizontal layout, no labels, no title.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

generate_color_avg_histograms <- function(
    data,
    dv_mode      = c("adjusted", "all"), # "adjusted" uses *_avg_adj, "all" uses *_avg
    save_name    = "color_avg_histograms", # Base name for saved plots
    base_size    = 12,                      # base font size for the figure theme
    family       = "Arial",                 # font family (kept consistent everywhere)
    show_legend  = FALSE                    # show the condition legend on this panel
) {
  dv_mode <- match.arg(dv_mode)
  
  # --- 1. Sanity Checks and Column Selection ---
  
  id_col <- "sub_id"
  if (!(id_col %in% names(data))) {
    stop("Expected column 'sub_id' in `data`.")
  }
  
  if (dv_mode == "adjusted") {
    avg_cols <- c("b1_avg_adj", "b2_avg_adj", "b3_avg_adj",
                  "o4_avg_adj", "o5_avg_adj", "o6_avg_adj",
                  "p1_avg_adj", "p2_avg_adj", "p3_avg_adj")
    y_lab <- "Adjusted pumps"
  } else {
    avg_cols <- c("b1_avg", "b2_avg", "b3_avg",
                  "o4_avg", "o5_avg", "o6_avg",
                  "p1_avg", "p2_avg", "p3_avg")
    y_lab <- "Pumps"
  }
  
  miss <- setdiff(avg_cols, names(data))
  if (length(miss)) {
    stop(paste("Missing required columns:", paste(miss, collapse = ", ")))
  }
  
  # Select only required columns
  sel <- data %>%
    dplyr::distinct(!!sym(id_col), .keep_all = TRUE) %>%
    dplyr::select(
      !!sym(id_col),
      all_of(avg_cols)
    )
  
  # --- 2. Calculate Average Inflation per Color per Participant ---
  
  long_df <- sel %>%
    tidyr::pivot_longer(
      cols = all_of(avg_cols),
      names_to  = "cell",
      values_to = "block_mean_infl"
    ) %>%
    dplyr::mutate(
      Color = factor(dplyr::case_when(
        startsWith(cell, "b") ~ "Blue",
        startsWith(cell, "o") ~ "Orange",
        startsWith(cell, "p") ~ "Pink",
        TRUE ~ NA_character_
      ), levels = c("Blue", "Orange", "Pink"))
    ) %>%
    dplyr::filter(!is.na(Color))
  
  color_avg_df <- long_df %>%
    dplyr::group_by(!!sym(id_col), Color) %>%
    dplyr::summarise(
      avg_infl_per_color = mean(block_mean_infl, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(is.finite(avg_infl_per_color)) 
  
  n_participants <- dplyr::n_distinct(color_avg_df[[id_col]])
  
  if (!n_participants) {
    stop("No valid data points found to generate histograms.")
  }
  
  # --- 3. Generate Histograms (Plotting) ---
  
  # Condition colours match the rest of the deck (BAL_COLS): the orange is the
  # deeper #E68D33 so it never reads as yellow.
  fill_map  <- c(Blue = "#2630F5", Orange = "#E68D33", Pink = "#CF3160")
  cond_labs <- c(Blue = "pre-reversal", Orange = "post-reversal",
                 Pink = "control")

  p_hist <- ggplot(color_avg_df, aes(x = avg_infl_per_color, fill = Color)) +
    geom_histogram(
      aes(y = after_stat(density)),
      bins = 15,
      alpha = 0.7,
      color = "black",
      na.rm = TRUE
    ) +
    geom_density(
      aes(color = Color),
      linewidth = 1,
      alpha = 0.5,
      na.rm = TRUE
    ) +
    # Shared y-scale (the three densities span nearly the same range), so the
    # y-axis ticks/label are drawn ONLY on the leftmost panel.
    facet_wrap(~ Color, scales = "free_x", ncol = 3) +
    scale_fill_manual(values = fill_map, labels = cond_labs, name = NULL) +
    scale_color_manual(values = fill_map, labels = cond_labs, name = NULL,
                       guide = "none") +
    labs(
      title = NULL,
      x     = paste("Average", y_lab, "per Participant"),
      y     = "Density"
    ) +
    theme_minimal(base_size = base_size, base_family = family) +
    theme(
      text = element_text(family = family),
      plot.title = element_blank(),
      legend.position = if (show_legend) "bottom" else "none",
      legend.text = element_text(size = base_size + 2, family = family,
                                 colour = "black"),
      strip.text.x = element_blank(),
      axis.title = element_text(size = base_size + 2, family = family,
                                colour = "black"),
      # uniform black axis numbers/ticks to match the rest of the figure
      axis.text  = element_text(size = base_size, family = family,
                                face = "bold", colour = "black")
    )
  
  print(p_hist)
  
  
  # --- 4. Save Plot ---
  fig_base <- file.path("output", "figures", save_name)
  dir.create(dirname(fig_base), recursive = TRUE, showWarnings = FALSE)
  
  ggsave(
    paste0(fig_base, ".png"), p_hist,
    width = 3500, # Increased width for a slightly wider plot
    height = 1200, 
    dpi = 300, 
    units = "px", bg = "white", limitsize = FALSE
  )
  
  invisible(list(
    data = color_avg_df,
    plot = p_hist
  ))
}