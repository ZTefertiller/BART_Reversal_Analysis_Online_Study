# This file contains common helper functions that can be
# shared by multiple analysis scripts.

mean_or_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

# Drop popped trials ONLY for averaging (do NOT use this to assign bins)
drop_popped <- function(df) {
  flag <- intersect(names(df), c("exploded","popped","pop","burst","explosion","exploded_flag"))
  if (length(flag) == 0) return(df)  # assume 'inflations' already excludes pops
  f <- flag[1]
  df %>% filter(!(!!as.name(f) %in% c(TRUE, 1, "yes", "Yes", "YES")))
}

# Assign fixed 10-trial bins from trial order using ALL trials (popped + non-popped)
assign_fixed_bins <- function(df, n_keep = 30L) {
  df %>%
    arrange(sub_id, trial_number) %>%
    group_by(sub_id) %>%
    mutate(rank30 = row_number()) %>%
    filter(rank30 <= n_keep) %>%
    mutate(
      BlockNum = ((rank30 - 1) %/% 10) + 1L,  # 1, 2, 3
      Block = factor(BlockNum, levels = c(1,2,3),
                     labels = c("First 10","Second 10","Last 10"))
    ) %>%
    ungroup()
}

# Summarise means within fixed bins by LEFT-JOINING ONLY the successful 'inflations'
fixedbin_mean_success <- function(all_trials_slice, success_slice) {
  all_trials_slice %>%
    select(-any_of("inflations")) %>%  # avoid .x/.y suffixing
    left_join(
      success_slice %>% select(sub_id, trial_number, inflations),
      by = c("sub_id", "trial_number")
    ) %>%
    group_by(sub_id, Block) %>%
    summarise(mean_infl = mean_or_na(inflations), .groups = "drop")
}


# --- utils.R additions ---

# Ensure adjusted_inflations exists (0 if popped, else inflations)
ensure_adjusted <- function(df) {
  stopifnot(all(c("inflations","popped") %in% names(df)))
  if (!"adjusted_inflations" %in% names(df)) {
    df <- df %>%
      mutate(
        inflations = as.numeric(inflations),
        .popped_flag = tolower(as.character(popped)) %in% c("1","true","yes","y","t"),
        adjusted_inflations = if_else(.popped_flag, 0, inflations, missing = inflations)
      ) %>%
      relocate(adjusted_inflations, .after = popped) %>%
      select(-.popped_flag)
  }
  df
}

# Compute blocks by trial order WITHIN COLOR (groups of 10)
add_blocks_per_color <- function(df, block_size = 10L) {
  stopifnot(all(c("sub_id","balloon_color","trial_number") %in% names(df)))
  if (!"block" %in% names(df)) {
    df <- df %>%
      arrange(sub_id, balloon_color, trial_number) %>%
      group_by(sub_id, balloon_color) %>%
      mutate(block = ceiling(row_number() / block_size)) %>%
      ungroup()
  }
  df
}

# Get per-block means with a toggle:
# mean_mode = "adjusted" (default; uses adjusted_inflations)
#           = "all" (uses inflations as-is, popped included)
block_means <- function(df, mean_mode = c("adjusted","all"),
                        block_size = 10L, keep_blocks = 1:3) {
  mean_mode <- match.arg(mean_mode)
  df <- ensure_adjusted(df)
  df <- add_blocks_per_color(df, block_size)
  
  value_col <- if (mean_mode == "adjusted") "adjusted_inflations" else "inflations"
  dv_label  <- if (mean_mode == "adjusted") "Average adjusted pumps" else "Average pumps"
  
  out <- df %>%
    filter(block %in% keep_blocks) %>%
    group_by(sub_id, balloon_color, block) %>%
    summarise(mean_value = mean(.data[[value_col]], na.rm = TRUE), .groups = "drop") %>%
    mutate(dv_label = dv_label)
  
  out
}



# age and sex 
summarize_recruitment <- function(df, recruitment_name) {
  df %>%
    distinct(sub_id, .keep_all = TRUE) %>%
    summarise(
      recruitment = recruitment_name,
      n = n(),
      mean_age = mean(age, na.rm = TRUE),
      age_min = min(age, na.rm = TRUE),
      age_max = max(age, na.rm = TRUE),
      n_female = sum(sex %in% c("female", "f", 1), na.rm = TRUE),
      n_male   = sum(sex %in% c("male", "m", 0), na.rm = TRUE)
    )
}
