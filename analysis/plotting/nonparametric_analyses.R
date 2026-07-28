# analysis/plotting/nonparametric_analyses.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(grid)
  library(WRS2)
  library(purrr)
})

generate_binned_anova_wrs2_yuend <- function(
    data,
    dv_mode         = c("adjusted", "all", "explosions"),   # "adjusted" = *_avg_adj, "all" = *_avg, "explosions" = popped counts
    save_name       = "binned_rm_wrs2_yuend",
    tr              = 0.2,                    # trim level
    p_adjust_method = c("none", "holm"),      # adjustment for yuend p-values
    y_lab_override  = NULL,                    # override the y-axis label
    y_limits        = NULL,                    # override y-axis limits (shared scales)
    bar_gap         = 1.0,                     # multiplier on spacing between stacked sig bars
    bar_base        = 1.0,                     # multiplier on the first bar's offset above the box
    base_size       = 20,                      # base font size for the figure theme
    sig_size        = 9,                       # size of the significance stars/markers
    sig_linewidth   = 0.9,                     # thickness of the significance comparison bars
    # To align significance bars ACROSS panels (so e.g. First-10 sits at the
    # same y in both), supply a shared per-block first-bar y and a shared delta.
    # Named numeric over c("First 10","Middle 10","Last 10"); overrides the
    # per-panel auto-placement when given. bar_delta_override sets the spacing.
    bar_y_override     = NULL,
    bar_delta_override = NULL
) {
  dv_mode         <- match.arg(dv_mode)
  p_adjust_method <- match.arg(p_adjust_method)
  
  # ---------- required columns ----------
  need_adj <- c(
    "participant_id",
    "b1_avg_adj","b2_avg_adj","b3_avg_adj",
    "o4_avg_adj","o5_avg_adj","o6_avg_adj",
    "p1_avg_adj","p2_avg_adj","p3_avg_adj"
  )
  need_all <- c(
    "participant_id",
    "b1_avg","b2_avg","b3_avg",
    "o4_avg","o5_avg","o6_avg",
    "p1_avg","p2_avg","p3_avg"
  )
  need <- switch(dv_mode,
                 adjusted   = need_adj,
                 all        = need_all,
                 explosions = c("participant_id", "balloon_color",
                                "trial_number", "popped"))
  miss <- setdiff(need, names(data))
  if (length(miss)) stop("Missing required columns: ", paste(miss, collapse = ", "))

  # ---------- select wide per-participant means ----------
  if (dv_mode == "adjusted") {
    sel <- data %>%
      distinct(participant_id, .keep_all = TRUE) %>%
      select(
        participant_id,
        b1 = b1_avg_adj, b2 = b2_avg_adj, b3 = b3_avg_adj,
        o4 = o4_avg_adj, o5 = o5_avg_adj, o6 = o6_avg_adj,
        p1 = p1_avg_adj, p2 = p2_avg_adj, p3 = p3_avg_adj
      )
    y_lab <- "Adjusted pumps"
  } else if (dv_mode == "all") {
    sel <- data %>%
      distinct(participant_id, .keep_all = TRUE) %>%
      select(
        participant_id,
        b1 = b1_avg, b2 = b2_avg, b3 = b3_avg,
        o4 = o4_avg, o5 = o5_avg, o6 = o6_avg,
        p1 = p1_avg, p2 = p2_avg, p3 = p3_avg
      )
    y_lab <- "Pumps"
  } else {
    # ---- explosions: count popped trials per participant x color x block-of-10 ----
    # Matches the b1-b3 / o4-o6 / p1-p3 cells used for the pump analyses:
    # blocks are 10-trial bins within each colour, ordered by trial number.
    is_pop <- function(x) tolower(as.character(x)) %in% c("1", "true", "t", "yes", "y")
    sel <- data %>%
      filter(balloon_color %in% c("b", "o", "p")) %>%
      mutate(.pop = is_pop(popped)) %>%
      arrange(participant_id, balloon_color, trial_number) %>%
      group_by(participant_id, balloon_color) %>%
      mutate(.blk = ceiling(row_number() / 10)) %>%
      ungroup() %>%
      mutate(cell = paste0(balloon_color, .blk)) %>%
      filter(cell %in% c("b1","b2","b3","o4","o5","o6","p1","p2","p3")) %>%
      group_by(participant_id, cell) %>%
      summarise(n_pop = sum(.pop, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = cell, values_from = n_pop) %>%
      select(participant_id, b1, b2, b3, o4, o5, o6, p1, p2, p3)
    y_lab <- "Explosions"
  }
  
  # ---------- long format with Color / Block ----------
  long <- sel %>%
    pivot_longer(-participant_id, names_to = "cell", values_to = "mean_infl") %>%
    mutate(
      Color = case_when(
        startsWith(cell, "b") ~ "Blue",
        startsWith(cell, "o") ~ "Orange",
        startsWith(cell, "p") ~ "Pink",
        TRUE ~ NA_character_
      ),
      Block = case_when(
        cell %in% c("b1","o4","p1") ~ "First 10",
        cell %in% c("b2","o5","p2") ~ "Middle 10",
        cell %in% c("b3","o6","p3") ~ "Last 10",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(Color), !is.na(Block))
  
  # ---------- keep only complete 9-cell participants ----------
  wide_check <- long %>%
    unite(Cell, Color, Block, sep = ":", remove = FALSE) %>%
    select(participant_id, Cell, mean_infl) %>%
    pivot_wider(names_from = Cell, values_from = mean_infl)
  
  complete_ids <- wide_check %>% tidyr::drop_na() %>% pull(participant_id)
  
  anova_df <- long %>%
    filter(participant_id %in% complete_ids) %>%
    mutate(
      participant_id = factor(participant_id),
      Color = factor(Color, levels = c("Blue","Orange","Pink")),
      Block = factor(Block, levels = c("First 10","Middle 10","Last 10")),
      cell  = factor(cell,
                     levels = c("b1","o4","p1",
                                "b2","o5","p2",
                                "b3","o6","p3"))
    )
  
  if (!nrow(anova_df)) stop("No complete participants (all 9 cells) for WRS2 ANOVA + yuend.")
  
  # Wide again, but only complete cases, for yuend
  sel_complete <- sel %>%
    filter(participant_id %in% complete_ids)
  
  # ==========================================================
  # 1) ONE-WAY trimmed-means RM ANOVA over 9 cells
  # ==========================================================
  cat("\n--- WRS2 trimmed means repeated-measures ANOVA (one-way over 9 cells) ---\n")
  cat("Trim level (tr): ", tr, "\n")
  
  rm_fit <- WRS2::rmanova(
    y      = anova_df$mean_infl,
    groups = anova_df$cell,
    blocks = anova_df$participant_id,
    tr     = tr
  )
  
  print(rm_fit)
  
  Fv  <- as.numeric(rm_fit$test)
  df1 <- as.numeric(rm_fit$df1)
  df2 <- as.numeric(rm_fit$df2)
  pv  <- as.numeric(rm_fit$p.value)
  
  anova_line <- sprintf(
    "One-way rmanova over 9 cells: F(%.2f, %.2f) = %.2f, p = %s",
    df1, df2, Fv, format.pval(pv, digits = 3, eps = .001)
  )
  cat("\nANOVA summary:\n", anova_line, "\n", sep = "")
  
  # ==========================================================
  # 2) Yuen’s paired trimmed-mean tests (yuend) for requested contrasts
  # ==========================================================
  pair_list <- list(
    list(cell1 = "b1", cell2 = "o4"),
    list(cell1 = "b1", cell2 = "p1"),
    list(cell1 = "o4", cell2 = "p1"),
    list(cell1 = "b2", cell2 = "o5"),
    list(cell1 = "b2", cell2 = "p2"),
    list(cell1 = "o5", cell2 = "p2"),
    list(cell1 = "b3", cell2 = "o6"),
    list(cell1 = "b3", cell2 = "p3"),
    list(cell1 = "o6", cell2 = "p3")
  )
  
  posthoc_res <- purrr::map_dfr(pair_list, function(p) {
    c1 <- p$cell1
    c2 <- p$cell2
    x  <- sel_complete[[c1]]
    y  <- sel_complete[[c2]]
    
    ok <- !is.na(x) & !is.na(y)
    x  <- x[ok]; y <- y[ok]
    
    if (!length(x) || !length(y) || length(x) != length(y)) {
      return(tibble::tibble(
        cell1 = c1,
        cell2 = c2,
        n     = length(x),
        test  = NA_real_,
        df    = NA_real_,
        p_raw = NA_real_,
        diff  = NA_real_
      ))
    }
    
    yu <- WRS2::yuend(x, y, tr = tr)
    
    tibble::tibble(
      cell1 = c1,
      cell2 = c2,
      n     = length(x),
      test  = as.numeric(yu$test),
      df    = as.numeric(yu$df),
      p_raw = as.numeric(yu$p.value),
      diff  = as.numeric(yu$diff)
    )
  })
  
  posthoc_res <- posthoc_res %>%
    mutate(
      p_adj = stats::p.adjust(p_raw, method = p_adjust_method),
      pair  = paste0(cell1, " vs ", cell2)
    )
  
  # stars based on chosen adjusted p
  posthoc_res <- posthoc_res %>%
    mutate(
      p_used = if (p_adjust_method == "none") p_raw else p_adj,
      sig    = case_when(
        is.na(p_used)        ~ "n.s.",
        p_used < 0.001       ~ "***",
        p_used < 0.01        ~ "**",
        p_used < 0.05        ~ "*",
        TRUE                 ~ "n.s."
      )
    )
  
  cat("\n--- Yuen's trimmed-mean paired tests (yuend) ---\n")
  cat("Trim level (tr): ", tr, "\n")
  cat("p-value adjustment used in reporting: ", p_adjust_method, "\n\n")
  print(
    posthoc_res %>%
      select(pair, n, test, df, p_raw, p_adj, sig)
  )
  
  # ==========================================================
  # 3) Build geometry info for significance bars
  # ==========================================================
  posthoc_res <- posthoc_res %>%
    mutate(
      Block = case_when(
        cell1 %in% c("b1","o4","p1") ~ "First 10",
        cell1 %in% c("b2","o5","p2") ~ "Middle 10",
        cell1 %in% c("b3","o6","p3") ~ "Last 10",
        TRUE ~ NA_character_
      )
    )

  posthoc_res <- posthoc_res %>%
    mutate(
      block_num = as.numeric(factor(Block,
                                    levels = c("First 10","Middle 10","Last 10")))
    )
  
  color_offset <- function(cell) {
    first_letter <- substr(cell, 1, 1)
    offs <- c(b = -0.25, o = 0, p = 0.25)
    offs[first_letter]
  }
  
  posthoc_res <- posthoc_res %>%
    mutate(
      x1 = block_num + vapply(cell1, color_offset, numeric(1)),
      x2 = block_num + vapply(cell2, color_offset, numeric(1))
    )
  
  # Anchor sig bars on the VISIBLE box top = the boxplot upper whisker (largest
  # value within Q3 + 1.5*IQR), per colour, then the max across colours in the
  # block. Outliers aren't drawn (outlier.shape = NA), so using the raw data max
  # would float the bars above empty space and inflate the y-axis.
  .upper_whisker <- function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    qs <- stats::quantile(x, c(0.25, 0.75), names = FALSE)
    cap <- qs[2] + 1.5 * (qs[2] - qs[1])
    max(x[x <= cap], na.rm = TRUE)
  }
  max_y_by_block <- anova_df %>%
    group_by(Block, Color) %>%
    summarise(w = .upper_whisker(mean_infl), .groups = "drop_last") %>%
    summarise(ymax = max(w, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      block_num = as.numeric(factor(Block,
                                    levels = c("First 10","Middle 10","Last 10")))
    )
  
  posthoc_res <- posthoc_res %>%
    left_join(max_y_by_block, by = c("Block", "block_num"))
  
  y_span <- max(anova_df$mean_infl, na.rm = TRUE) -
    min(anova_df$mean_infl, na.rm = TRUE)
  delta <- if (!is.null(bar_delta_override)) {
    bar_delta_override
  } else if (is.finite(y_span) && y_span > 0) {
    y_span * 0.08
  } else {
    0.2
  }

  posthoc_res <- posthoc_res %>%
    group_by(Block) %>%
    arrange(pair, .by_group = TRUE) %>%
    mutate(
      rank_in_block = row_number(),
      # First-bar base: shared per-block override (aligns bars across panels) or
      # this panel's own box top + bar_base offset. bar_gap spaces the stack.
      bar_y0        = if (!is.null(bar_y_override))
                        unname(bar_y_override[as.character(Block[1])])
                      else ymax + bar_base * delta,
      y             = bar_y0 + (rank_in_block - 1) * delta * bar_gap,
      y_text        = y + 0.3 * delta
    ) %>%
    ungroup()
  
  # ==========================================================
  # 4) Plot: Blue / Orange / Pink by Block (no text overlay, no outlier dots, just bars + stars)
  # ==========================================================
  # Match the rest of the deck (BAL_COLS): deeper orange so it never reads yellow.
  fill_map  <- c(Blue = "#2630F5", Orange = "#E68D33", Pink = "#CF3160")
  label_map <- c(Blue = "pre-reversal", Orange = "post-reversal", Pink = "control")
  
  pd <- position_dodge(width = 0.8)

  # Title and N intentionally omitted from the figure (added later in Word).
  box_df <- anova_df

  # Explosion counts are small integers, so plain boxplots collapse to thin
  # bars. For that mode, lighten the boxes and overlay jittered points (with a
  # little vertical jitter to separate tied integer counts) so the spread shows.
  is_expl   <- dv_mode == "explosions"
  box_alpha <- if (is_expl) 0.45 else 1
  expl_jitter <- if (is_expl) {
    geom_point(
      aes(fill = Color),
      position = position_jitterdodge(jitter.width = 0.18,
                                      jitter.height = 0.25,
                                      dodge.width = 0.8),
      shape = 21, size = 1.3, stroke = 0.2, colour = "grey25",
      alpha = 0.5, show.legend = FALSE, na.rm = TRUE
    )
  } else {
    NULL
  }

  p_box <- ggplot(box_df, aes(x = Block, y = mean_infl)) +
    geom_boxplot(
      aes(fill = Color),
      position = pd, width = 0.7,
      colour = "black",
      outlier.shape = NA,   # <- no outlier points
      alpha = box_alpha,
      na.rm = TRUE, show.legend = TRUE
    ) +
    expl_jitter +
    geom_segment(
      data = posthoc_res,
      inherit.aes = FALSE,
      aes(x = x1, xend = x2, y = y, yend = y),
      linewidth = sig_linewidth
    ) +
    geom_text(
      data = posthoc_res,
      inherit.aes = FALSE,
      aes(x = (x1 + x2) / 2, y = y_text, label = sig),
      vjust = 0, size = sig_size, fontface = "bold"
    ) +
    scale_fill_manual(values = fill_map, labels = label_map) +
    scale_color_manual(values = fill_map, labels = label_map) +
    guides(fill = guide_legend(title = NULL), colour = "none") +
    labs(
      y = if (!is.null(y_lab_override)) y_lab_override else y_lab,
      x = NULL
    ) +
    theme_minimal(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(face = "bold"),
      legend.position = "bottom",
      legend.text  = element_text(size = base_size, family = "Arial", face = "bold",
                                  colour = "black"),
      plot.title = element_text(hjust = 0.5),
      # Uniform axis text: same Arial size, same weight, all BLACK — for both the
      # First/Middle/Last 10 labels (x) and the numeric ticks (y).
      axis.text.x  = element_text(size = base_size, family = "Arial",
                                  face = "bold", colour = "black"),
      axis.text.y  = element_text(size = base_size, family = "Arial",
                                  face = "bold", colour = "black"),
      axis.title.y = element_text(size = base_size + 2, family = "Arial",
                                  face = "bold", colour = "black")
    )

  if (!is.null(y_limits)) {
    p_box <- p_box + coord_cartesian(ylim = y_limits)
  }
  
  print(p_box)
  
  base_path <- file.path("output", "figures", save_name)
  dir.create(dirname(base_path), recursive = TRUE, showWarnings = FALSE)
  
  ggsave(paste0(base_path, ".png"), p_box,
         width = 2000, height = 1200, dpi = 300,
         units = "px", bg = "white", limitsize = FALSE)
  ggsave(paste0(base_path, ".svg"), p_box,
         width = 2000, height = 1200, dpi = 300,
         units = "px", bg = "white", limitsize = FALSE, device = "svg")
  
  # per-block box tops (max whisker top per block), for callers computing shared
  # significance-bar anchors across panels.
  block_ymax <- max_y_by_block %>%
    select(Block, ymax) %>%
    { setNames(.$ymax, as.character(.$Block)) }
  n_bars_per_block <- max(posthoc_res$rank_in_block, na.rm = TRUE)

  invisible(list(
    plot        = p_box,
    rmanova     = rm_fit,
    posthoc_tbl = posthoc_res %>% select(pair, n, test, df, p_raw, p_adj, sig),
    # geometry for callers that want an exact y-axis ceiling: the highest star
    # anchor and the per-bar spacing unit (delta), so a star-height cushion can
    # be added without guessing.
    top_text    = max(posthoc_res$y_text, na.rm = TRUE),
    delta       = delta,
    block_ymax  = block_ymax,           # named over First/Middle/Last 10
    n_bars      = n_bars_per_block       # stacked bars per block (3 here)
  ))
}


