# analysis/plotting/inflation_regression.R
# OLS regressions: per-participant mean adjusted inflations (per color) ~ SPQ/PDI/CAPS.
# Returns the regression table and the merged data frame used for scatterplotting.
#
# Usage:
#   source(here::here("analysis", "plotting", "inflation_regression.R"))
#   res <- run_inflation_questionnaire_regressions(data)
#   res$table   # tidy stats table
#   res$data    # long df ready for ggplot (one row per participant × color)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(broom)
  library(purrr)
})

run_inflation_questionnaire_regressions <- function(
  data,
  questionnaires = c(
    spq_total  = "SPQ Total",
    pdi_total  = "PDI Total",
    caps_total = "CAPS Total"
  ),
  colors = c(
    b = "Blue (pre-reversal)",
    o = "Orange (post-reversal)",
    p = "Pink (control)"
  ),
  out_dir = NULL
) {
  is_popped <- function(x) tolower(as.character(x)) %in% c("1", "true", "t", "yes", "y")

  # Per-participant mean adjusted inflation per color.
  # Condition = the HIGH-VALUE phase of each colour (matches the EWMV fits):
  #   Blue   = pre-reversal  (trial_number < 91)
  #   Orange = post-reversal (trial_number > 90)
  #   Pink   = control (all trials). Other colours: all trials.
  adj_by_color <- data %>%
    filter(
      balloon_color %in% names(colors),
      balloon_color != "b" | trial_number < 91,
      balloon_color != "o" | trial_number > 90
    ) %>%
    mutate(.pop = is_popped(popped)) %>%
    group_by(sub_id, balloon_color) %>%
    summarise(
      adj_infl = mean(inflations[!.pop], na.rm = TRUE),
      .groups  = "drop"
    ) %>%
    mutate(
      Color       = recode(balloon_color, !!!colors),
      color_code  = balloon_color
    )

  # One questionnaire row per participant
  q_df <- data %>%
    distinct(sub_id, .keep_all = TRUE) %>%
    select(sub_id, all_of(names(questionnaires)))

  merged <- adj_by_color %>%
    left_join(q_df, by = "sub_id")

  # Regression: adj_infl ~ questionnaire_score, per (color × questionnaire)
  results <- map_dfr(names(colors), function(col_code) {
    d_color <- filter(merged, color_code == col_code)
    imap_dfr(questionnaires, function(q_label, q_col) {
      df_pair <- d_color %>%
        transmute(y = adj_infl, x = .data[[q_col]]) %>%
        filter(is.finite(y), is.finite(x))

      if (nrow(df_pair) < 3) return(tibble())

      fit <- lm(y ~ x, data = df_pair)
      td  <- tidy(fit, conf.int = TRUE) %>% filter(term == "x")
      gl  <- glance(fit)

      tibble(
        Condition     = colors[[col_code]],
        Questionnaire = q_label,
        n             = nrow(df_pair),
        b             = td$estimate,
        SE            = td$std.error,
        t             = td$statistic,
        p             = td$p.value,
        `CI lo`       = td$conf.low,
        `CI hi`       = td$conf.high,
        R2            = gl$r.squared
      )
    })
  })

  # Add significance stars
  results <- results %>%
    mutate(
      stars = case_when(
        p < 0.001 ~ "***",
        p < 0.01  ~ "**",
        p < 0.05  ~ "*",
        TRUE      ~ ""
      ),
      p_fmt = paste0(format.pval(p, digits = 2, eps = 0.001), stars)
    )

  if (!is.null(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    readr::write_csv(results, file.path(out_dir, "inflation_questionnaire_regressions.csv"))
  }

  invisible(list(table = results, data = merged))
}

# OLS regressions: per-participant EXPLOSION COUNT (per condition) ~ questionnaire.
# Same phase-filtering / structure as the adjusted-pumps version above.
run_explosion_questionnaire_regressions <- function(
  data,
  questionnaires = c(spq_total = "SPQ Total", pdi_total = "PDI Total",
                     caps_total = "CAPS Total"),
  colors = c(b = "Blue (pre-reversal)", o = "Orange (post-reversal)",
             p = "Pink (control)")
) {
  is_popped <- function(x) tolower(as.character(x)) %in% c("1","true","t","yes","y")

  expl_by_color <- data %>%
    filter(
      balloon_color %in% names(colors),
      balloon_color != "b" | trial_number < 91,
      balloon_color != "o" | trial_number > 90
    ) %>%
    mutate(.pop = is_popped(popped)) %>%
    group_by(sub_id, balloon_color) %>%
    summarise(n_expl = sum(.pop, na.rm = TRUE), .groups = "drop") %>%
    mutate(Color = recode(balloon_color, !!!colors), color_code = balloon_color)

  q_df <- data %>% distinct(sub_id, .keep_all = TRUE) %>%
    select(sub_id, all_of(names(questionnaires)))
  merged <- expl_by_color %>% left_join(q_df, by = "sub_id")

  results <- map_dfr(names(colors), function(col_code) {
    d_color <- filter(merged, color_code == col_code)
    imap_dfr(questionnaires, function(q_label, q_col) {
      df_pair <- d_color %>% transmute(y = n_expl, x = .data[[q_col]]) %>%
        filter(is.finite(y), is.finite(x))
      if (nrow(df_pair) < 3) return(tibble())
      fit <- lm(y ~ x, data = df_pair)
      td  <- tidy(fit, conf.int = TRUE) %>% filter(term == "x")
      gl  <- glance(fit)
      tibble(Condition = colors[[col_code]], Questionnaire = q_label,
             n = nrow(df_pair), b = td$estimate, SE = td$std.error,
             t = td$statistic, p = td$p.value, R2 = gl$r.squared)
    })
  })
  results <- results %>%
    mutate(stars = case_when(p < .001 ~ "***", p < .01 ~ "**", p < .05 ~ "*", TRUE ~ ""),
           p_fmt = paste0(format.pval(p, digits = 2, eps = .001), stars))
  invisible(list(table = results, data = merged))
}
