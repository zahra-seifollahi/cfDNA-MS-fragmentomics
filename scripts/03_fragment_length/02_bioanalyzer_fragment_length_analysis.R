# ============================================================
# Bioanalyzer fragment-size distribution analysis
#
# Purpose:
#   Analyze raw Bioanalyzer fragment-size interval data.
#
# Important:
#   This script uses all Bioanalyzer samples.
#   It does NOT apply BAM-derived QC filtering.
#
# Input:
#   Excel file with sheet: Intervals_50_300
#
# Required columns:
#   Sample ID
#   Range
#   %
#
# Optional columns:
#   Raw %
#   Avg. bp
#
# Usage:
#   Rscript scripts/03_fragment_length/02_bioanalyzer_fragment_length_analysis.R /path/to/bioanalyzer_bins_50_300_intervals.xlsx
#
# Output:
#   results/tables/fragment_length/
#   results/figures/fragment_length/
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(mgcv)
  library(ggrepel)
  library(writexl)
})

graphics.off()

# ============================================================
# 1. Settings
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

bioanalyzer_file <- ifelse(
  length(args) >= 1,
  args[1],
  "data/example/bioanalyzer_bins_50_300_intervals.xlsx"
)

sheet_name <- "Intervals_50_300"

table_dir <- "results/tables/fragment_length"
figure_dir <- "results/figures/fragment_length"

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c("Healthy", "Remission", "Relapse")

group_cols <- c(
  "Healthy" = "darkgreen",
  "Remission" = "deeppink3",
  "Relapse" = "blue3"
)

group_fills <- c(
  "Healthy" = "darkseagreen3",
  "Remission" = "lightpink",
  "Relapse" = "lightblue"
)

# Use generic serif font.
# This avoids the Mac error: invalid font type / Times New Roman not found.
plot_family <- "serif"

theme_bio <- function(base_size = 14) {
  theme_classic(base_size = base_size, base_family = plot_family) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, family = plot_family),
      plot.subtitle = element_text(hjust = 0.5, family = plot_family),
      axis.title = element_text(color = "black", family = plot_family),
      axis.text = element_text(color = "black", family = plot_family),
      legend.title = element_text(family = plot_family),
      legend.text = element_text(family = plot_family),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.line = element_blank()
    )
}

# ============================================================
# 2. Helper functions
# ============================================================

assign_group_from_sample <- function(sample_id) {
  sample_num <- as.integer(str_extract(as.character(sample_id), "\\d+"))

  case_when(
    sample_num >= 1  & sample_num <= 30 ~ "Healthy",
    sample_num >= 31 & sample_num <= 60 ~ "Remission",
    sample_num >= 61 & sample_num <= 84 ~ "Relapse",
    TRUE ~ NA_character_
  )
}

safe_kruskal_p <- function(value, group) {
  ok <- !is.na(value) & !is.na(group)

  value <- value[ok]
  group <- droplevels(group[ok])

  if (length(value) < 2) {
    return(NA_real_)
  }

  if (length(unique(group)) < 2) {
    return(NA_real_)
  }

  out <- tryCatch(
    kruskal.test(value ~ group)$p.value,
    error = function(e) NA_real_
  )

  as.numeric(out)
}

safe_pairwise_wilcox <- function(value, group) {
  ok <- !is.na(value) & !is.na(group)

  value <- value[ok]
  group <- droplevels(group[ok])

  if (length(value) < 2) {
    return(NULL)
  }

  if (length(unique(group)) < 2) {
    return(NULL)
  }

  tryCatch(
    pairwise.wilcox.test(
      x = value,
      g = group,
      p.adjust.method = "BH",
      exact = FALSE
    ),
    error = function(e) NULL
  )
}

get_pw <- function(pw, row, col) {
  if (is.null(pw)) {
    return(NA_real_)
  }

  value <- tryCatch(
    pw$p.value[row, col],
    error = function(e) NA_real_
  )

  as.numeric(value)
}

# ============================================================
# 3. Load Bioanalyzer data
# ============================================================

if (!file.exists(bioanalyzer_file)) {
  stop("Bioanalyzer file not found: ", bioanalyzer_file)
}

raw_df <- read_excel(
  path = bioanalyzer_file,
  sheet = sheet_name
)

required_cols <- c("Sample ID", "Range", "%")
missing_cols <- setdiff(required_cols, colnames(raw_df))

if (length(missing_cols) > 0) {
  stop(
    "Bioanalyzer file is missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

# Parse range correctly.
# Example: "100 - 105" -> range_start = 100, range_end = 105, range_mid = 102.5
range_parts <- str_match(
  as.character(raw_df$Range),
  "^\\s*(\\d+(?:\\.\\d+)?)\\s*-\\s*(\\d+(?:\\.\\d+)?)\\s*$"
)

df_clean <- raw_df %>%
  transmute(
    sample_id = as.character(`Sample ID`),
    range = as.character(`Range`),
    range_start = as.numeric(range_parts[, 2]),
    range_end = as.numeric(range_parts[, 3]),
    range_mid = (range_start + range_end) / 2,
    percent = readr::parse_number(as.character(`%`)),
    raw_percent = if ("Raw %" %in% colnames(raw_df)) {
      readr::parse_number(as.character(`Raw %`))
    } else {
      NA_real_
    },
    avg_bp_original = if ("Avg. bp" %in% colnames(raw_df)) {
      readr::parse_number(as.character(`Avg. bp`))
    } else {
      NA_real_
    }
  ) %>%
  filter(
    !is.na(sample_id),
    !is.na(range),
    !is.na(range_start),
    !is.na(range_end),
    !is.na(range_mid),
    !is.na(percent)
  ) %>%
  mutate(
    group = assign_group_from_sample(sample_id),
    group = factor(group, levels = group_levels)
  ) %>%
  filter(!is.na(group)) %>%
  arrange(sample_id, range_mid)

cat("\nBioanalyzer samples loaded:\n")
print(df_clean %>% distinct(sample_id, group) %>% count(group))

cat("\nRange check:\n")
print(
  df_clean %>%
    distinct(range, range_start, range_end, range_mid) %>%
    arrange(range_start) %>%
    head(10)
)

write_tsv(
  df_clean,
  file.path(table_dir, "bioanalyzer_clean_interval_data.tsv")
)

# ============================================================
# 4. Group-level mean per interval
# ============================================================

group_interval_mean <- df_clean %>%
  group_by(group, range, range_start, range_end, range_mid) %>%
  summarise(
    n_samples = n_distinct(sample_id),

    mean_percent = mean(percent, na.rm = TRUE),
    sd_percent = sd(percent, na.rm = TRUE),
    se_percent = sd_percent / sqrt(n_samples),

    mean_raw_percent = mean(raw_percent, na.rm = TRUE),
    sd_raw_percent = sd(raw_percent, na.rm = TRUE),
    se_raw_percent = sd_raw_percent / sqrt(n_samples),

    .groups = "drop"
  ) %>%
  arrange(group, range_start)

write_tsv(
  group_interval_mean,
  file.path(table_dir, "bioanalyzer_group_interval_mean.tsv")
)

# ============================================================
# 5. Fit GAM curve and identify fitted-line peak
# ============================================================

x_grid <- seq(50, 300, by = 0.5)

gam_curve_df <- group_interval_mean %>%
  group_by(group) %>%
  group_modify(~ {
    fit <- mgcv::gam(mean_percent ~ s(range_mid, k = 10), data = .x)

    pred <- predict(
      fit,
      newdata = data.frame(range_mid = x_grid),
      se.fit = TRUE
    )

    tibble(
      range_mid = x_grid,
      fitted_percent = as.numeric(pred$fit),
      fitted_se = as.numeric(pred$se.fit)
    )
  }) %>%
  ungroup()

fitted_peak_df <- gam_curve_df %>%
  group_by(group) %>%
  slice_max(fitted_percent, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    peak_type = "Fitted GAM peak",
    peak_label = paste0(
      group, "\n",
      round(range_mid, 1), " bp\n",
      round(fitted_percent, 2), "%"
    )
  )

raw_peak_df <- group_interval_mean %>%
  group_by(group) %>%
  slice_max(mean_percent, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    group,
    peak_type = "Raw mean-bin peak",
    peak_range = range,
    peak_mid_bp = range_mid,
    peak_percent = mean_percent,
    n_samples
  )

fitted_peak_summary <- fitted_peak_df %>%
  transmute(
    group,
    peak_type,
    peak_range = NA_character_,
    peak_mid_bp = range_mid,
    peak_percent = fitted_percent,
    n_samples = NA_integer_
  )

peak_summary <- bind_rows(
  raw_peak_df,
  fitted_peak_summary
) %>%
  arrange(group, peak_type)

write_tsv(
  gam_curve_df,
  file.path(table_dir, "bioanalyzer_gam_curve.tsv")
)

write_tsv(
  peak_summary,
  file.path(table_dir, "bioanalyzer_peak_summary.tsv")
)

# ============================================================
# 6. Clean combined plot with peak point on fitted line
# ============================================================

fitted_peak_df <- fitted_peak_df %>%
  mutate(
    nudge_x = case_when(
      group == "Healthy" ~ -28,
      group == "Remission" ~ 0,
      group == "Relapse" ~ 28,
      TRUE ~ 0
    ),
    nudge_y = case_when(
      group == "Healthy" ~ 0.55,
      group == "Remission" ~ 0.95,
      group == "Relapse" ~ 0.55,
      TRUE ~ 0.55
    )
  )

p_combined_clean <- ggplot() +
  geom_col(
    data = group_interval_mean,
    aes(
      x = range_mid,
      y = mean_percent,
      fill = group,
      color = group
    ),
    position = position_identity(),
    width = 4.2,
    alpha = 0.22,
    linewidth = 0.25
  ) +
  geom_line(
    data = gam_curve_df,
    aes(
      x = range_mid,
      y = fitted_percent,
      color = group
    ),
    linewidth = 1.3
  ) +
  geom_point(
    data = fitted_peak_df,
    aes(
      x = range_mid,
      y = fitted_percent,
      color = group
    ),
    size = 3.5
  ) +
  geom_label_repel(
    data = fitted_peak_df,
    aes(
      x = range_mid,
      y = fitted_percent,
      label = peak_label,
      color = group
    ),
    fill = "white",
    size = 3.4,
    family = plot_family,
    box.padding = 0.35,
    point.padding = 0.35,
    segment.color = "grey30",
    segment.size = 0.4,
    min.segment.length = 0,
    seed = 123,
    show.legend = FALSE,
    nudge_x = fitted_peak_df$nudge_x,
    nudge_y = fitted_peak_df$nudge_y
  ) +
  scale_fill_manual(values = group_fills) +
  scale_color_manual(values = group_cols) +
  scale_x_continuous(
    breaks = seq(50, 300, by = 25),
    limits = c(48, 302)
  ) +
  labs(
    x = "Fragment size bin midpoint (bp)",
    y = "Mean fragment percentage (%)",
    title = "Bioanalyzer fragment-size distribution by group",
    subtitle = "Bars show group mean per 5-bp interval; points indicate fitted GAM peak values",
    fill = "Group",
    color = "Group"
  ) +
  theme_bio(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  file.path(figure_dir, "bioanalyzer_combined_percent_GAM_peak_clean.png"),
  p_combined_clean,
  width = 10,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(figure_dir, "bioanalyzer_combined_percent_GAM_peak_clean.pdf"),
  p_combined_clean,
  width = 10,
  height = 7,
  device = "pdf"
)

# ============================================================
# 7. Per-bin Kruskal-Wallis tests
# ============================================================

bin_tests <- df_clean %>%
  group_by(range, range_start, range_end, range_mid) %>%
  summarise(
    n_groups = n_distinct(group),
    n_samples = n_distinct(sample_id),
    p_kruskal_percent = safe_kruskal_p(percent, group),
    p_kruskal_raw_percent = safe_kruskal_p(raw_percent, group),
    .groups = "drop"
  ) %>%
  mutate(
    padj_kruskal_percent = p.adjust(p_kruskal_percent, method = "BH"),
    padj_kruskal_raw_percent = p.adjust(p_kruskal_raw_percent, method = "BH"),
    percent_significant_BH_0.05 = padj_kruskal_percent < 0.05,
    raw_percent_significant_BH_0.05 = padj_kruskal_raw_percent < 0.05
  ) %>%
  arrange(padj_kruskal_percent)

write_tsv(
  bin_tests,
  file.path(table_dir, "bioanalyzer_per_bin_kruskal_tests.tsv")
)

# ============================================================
# 8. Pairwise Wilcoxon tests per bin
# ============================================================

pairwise_bin_tests <- df_clean %>%
  group_by(range, range_start, range_end, range_mid) %>%
  group_modify(~ {
    pw <- safe_pairwise_wilcox(.x$percent, .x$group)

    tibble(
      comparison = c(
        "Remission_vs_Healthy",
        "Relapse_vs_Healthy",
        "Relapse_vs_Remission"
      ),
      p_value = c(
        get_pw(pw, "Remission", "Healthy"),
        get_pw(pw, "Relapse", "Healthy"),
        get_pw(pw, "Relapse", "Remission")
      )
    )
  }) %>%
  ungroup() %>%
  mutate(
    padj_global_BH = p.adjust(p_value, method = "BH"),
    significant_global_BH_0.05 = padj_global_BH < 0.05
  ) %>%
  arrange(padj_global_BH)

write_tsv(
  pairwise_bin_tests,
  file.path(table_dir, "bioanalyzer_pairwise_bin_wilcoxon_percent.tsv")
)

# ============================================================
# 9. Sample-level whole-profile metrics
# ============================================================

sample_profile <- df_clean %>%
  select(sample_id, group, range_mid, percent) %>%
  arrange(sample_id, range_mid) %>%
  group_by(sample_id, group) %>%
  mutate(
    percent_norm = percent / sum(percent, na.rm = TRUE)
  ) %>%
  ungroup()

healthy_reference <- sample_profile %>%
  filter(group == "Healthy") %>%
  group_by(range_mid) %>%
  summarise(
    healthy_mean_percent_norm = mean(percent_norm, na.rm = TRUE),
    .groups = "drop"
  )

profile_distance <- sample_profile %>%
  left_join(healthy_reference, by = "range_mid") %>%
  group_by(sample_id, group) %>%
  summarise(
    L1_distance_from_healthy = sum(
      abs(percent_norm - healthy_mean_percent_norm),
      na.rm = TRUE
    ),
    L2_distance_from_healthy = sqrt(
      sum((percent_norm - healthy_mean_percent_norm)^2, na.rm = TRUE)
    ),
    KS_like_distance_from_healthy = max(
      abs(cumsum(percent_norm) - cumsum(healthy_mean_percent_norm)),
      na.rm = TRUE
    ),
    peak_bp = range_mid[which.max(percent)],
    peak_percent = max(percent, na.rm = TRUE),
    .groups = "drop"
  )

write_tsv(
  profile_distance,
  file.path(table_dir, "bioanalyzer_sample_profile_distance_from_healthy.tsv")
)

# ============================================================
# 10. Global statistical tests
# ============================================================

kruskal_L1 <- kruskal.test(L1_distance_from_healthy ~ group, data = profile_distance)
kruskal_L2 <- kruskal.test(L2_distance_from_healthy ~ group, data = profile_distance)
kruskal_KS <- kruskal.test(KS_like_distance_from_healthy ~ group, data = profile_distance)
kruskal_peak_bp <- kruskal.test(peak_bp ~ group, data = profile_distance)
kruskal_peak_percent <- kruskal.test(peak_percent ~ group, data = profile_distance)

pairwise_L1 <- safe_pairwise_wilcox(
  profile_distance$L1_distance_from_healthy,
  profile_distance$group
)

pairwise_KS <- safe_pairwise_wilcox(
  profile_distance$KS_like_distance_from_healthy,
  profile_distance$group
)

pairwise_peak_bp <- safe_pairwise_wilcox(
  profile_distance$peak_bp,
  profile_distance$group
)

pairwise_peak_percent <- safe_pairwise_wilcox(
  profile_distance$peak_percent,
  profile_distance$group
)

# ============================================================
# 11. Group-level summary
# ============================================================

group_profile_summary <- profile_distance %>%
  group_by(group) %>%
  summarise(
    n_samples = n(),

    mean_L1_distance = mean(L1_distance_from_healthy, na.rm = TRUE),
    median_L1_distance = median(L1_distance_from_healthy, na.rm = TRUE),

    mean_L2_distance = mean(L2_distance_from_healthy, na.rm = TRUE),
    median_L2_distance = median(L2_distance_from_healthy, na.rm = TRUE),

    mean_KS_like_distance = mean(KS_like_distance_from_healthy, na.rm = TRUE),
    median_KS_like_distance = median(KS_like_distance_from_healthy, na.rm = TRUE),

    mean_peak_bp = mean(peak_bp, na.rm = TRUE),
    median_peak_bp = median(peak_bp, na.rm = TRUE),

    mean_peak_percent = mean(peak_percent, na.rm = TRUE),
    median_peak_percent = median(peak_percent, na.rm = TRUE),

    .groups = "drop"
  )

write_tsv(
  group_profile_summary,
  file.path(table_dir, "bioanalyzer_group_profile_summary.tsv")
)

# ============================================================
# 12. Single-sheet statistical summary
# ============================================================

overall_tests_summary <- tibble(
  section = "Overall Kruskal-Wallis tests",
  metric = c(
    "L1 distance from Healthy profile",
    "L2 distance from Healthy profile",
    "KS-like cumulative distance from Healthy profile",
    "Sample-level peak bp",
    "Sample-level peak percent"
  ),
  comparison = "Healthy vs Remission vs Relapse",
  statistic = c(
    as.numeric(kruskal_L1$statistic),
    as.numeric(kruskal_L2$statistic),
    as.numeric(kruskal_KS$statistic),
    as.numeric(kruskal_peak_bp$statistic),
    as.numeric(kruskal_peak_percent$statistic)
  ),
  p_value = c(
    kruskal_L1$p.value,
    kruskal_L2$p.value,
    kruskal_KS$p.value,
    kruskal_peak_bp$p.value,
    kruskal_peak_percent$p.value
  ),
  p_adjusted = NA_real_,
  significant = p_value < 0.05,
  interpretation = case_when(
    metric == "Sample-level peak bp" & p_value < 0.05 ~
      "Significant: fragment peak position differs between groups.",
    metric == "Sample-level peak percent" & p_value >= 0.05 ~
      "Not significant: peak height/intensity is not different between groups.",
    str_detect(metric, "distance") & p_value >= 0.05 ~
      "Not significant: no strong evidence that the whole profile differs globally.",
    TRUE ~ "Check result."
  )
)

pairwise_tests_summary <- tibble(
  section = "Pairwise Wilcoxon tests",
  metric = c(
    "L1 distance from Healthy profile",
    "L1 distance from Healthy profile",
    "L1 distance from Healthy profile",
    "KS-like cumulative distance from Healthy profile",
    "KS-like cumulative distance from Healthy profile",
    "KS-like cumulative distance from Healthy profile",
    "Sample-level peak bp",
    "Sample-level peak bp",
    "Sample-level peak bp",
    "Sample-level peak percent",
    "Sample-level peak percent",
    "Sample-level peak percent"
  ),
  comparison = c(
    "Remission vs Healthy",
    "Relapse vs Healthy",
    "Relapse vs Remission",
    "Remission vs Healthy",
    "Relapse vs Healthy",
    "Relapse vs Remission",
    "Remission vs Healthy",
    "Relapse vs Healthy",
    "Relapse vs Remission",
    "Remission vs Healthy",
    "Relapse vs Healthy",
    "Relapse vs Remission"
  ),
  statistic = NA_real_,
  p_value = c(
    get_pw(pairwise_L1, "Remission", "Healthy"),
    get_pw(pairwise_L1, "Relapse", "Healthy"),
    get_pw(pairwise_L1, "Relapse", "Remission"),
    get_pw(pairwise_KS, "Remission", "Healthy"),
    get_pw(pairwise_KS, "Relapse", "Healthy"),
    get_pw(pairwise_KS, "Relapse", "Remission"),
    get_pw(pairwise_peak_bp, "Remission", "Healthy"),
    get_pw(pairwise_peak_bp, "Relapse", "Healthy"),
    get_pw(pairwise_peak_bp, "Relapse", "Remission"),
    get_pw(pairwise_peak_percent, "Remission", "Healthy"),
    get_pw(pairwise_peak_percent, "Relapse", "Healthy"),
    get_pw(pairwise_peak_percent, "Relapse", "Remission")
  ),
  p_adjusted = p_value,
  significant = p_adjusted < 0.05,
  interpretation = case_when(
    metric == "Sample-level peak bp" & comparison == "Relapse vs Healthy" & significant ~
      "Significant: Relapse has a different peak position from Healthy.",
    metric == "Sample-level peak bp" & comparison == "Relapse vs Remission" & significant ~
      "Significant: Relapse has a different peak position from Remission.",
    metric == "Sample-level peak bp" & comparison == "Remission vs Healthy" & !significant ~
      "Not significant: Remission peak position is not clearly different from Healthy.",
    !significant ~ "Not significant.",
    TRUE ~ "Significant."
  )
)

bin_summary_rows <- bin_tests %>%
  summarise(
    n_bins_tested = n(),
    n_significant_percent_bins_BH_0.05 = sum(padj_kruskal_percent < 0.05, na.rm = TRUE),
    significant_percent_bins = paste(range[padj_kruskal_percent < 0.05], collapse = "; "),
    n_significant_raw_percent_bins_BH_0.05 = sum(padj_kruskal_raw_percent < 0.05, na.rm = TRUE),
    significant_raw_percent_bins = paste(range[padj_kruskal_raw_percent < 0.05], collapse = "; ")
  ) %>%
  transmute(
    section = "Per-bin Kruskal-Wallis summary",
    metric = "Percent and raw percent per 5-bp bin",
    comparison = "Healthy vs Remission vs Relapse",
    statistic = NA_real_,
    p_value = NA_real_,
    p_adjusted = NA_real_,
    significant = n_significant_percent_bins_BH_0.05 > 0,
    interpretation = paste0(
      n_significant_percent_bins_BH_0.05,
      " percent-based bins were significant after BH correction: ",
      ifelse(significant_percent_bins == "", "none", significant_percent_bins),
      ". Raw-percent significant bins: ",
      ifelse(significant_raw_percent_bins == "", "none", significant_raw_percent_bins),
      "."
    )
  )

group_profile_summary_rows <- group_profile_summary %>%
  transmute(
    section = "Group descriptive summary",
    metric = paste0("Group = ", group),
    comparison = NA_character_,
    statistic = NA_real_,
    p_value = NA_real_,
    p_adjusted = NA_real_,
    significant = NA,
    interpretation = paste0(
      "n = ", n_samples,
      "; mean peak bp = ", round(mean_peak_bp, 2),
      "; median peak bp = ", round(median_peak_bp, 2),
      "; mean peak percent = ", round(mean_peak_percent, 2),
      "; mean L1 distance = ", round(mean_L1_distance, 3),
      "; mean KS-like distance = ", round(mean_KS_like_distance, 3)
    )
  )

fitted_peak_rows <- fitted_peak_df %>%
  transmute(
    section = "Fitted GAM peak summary",
    metric = paste0("Group = ", group),
    comparison = NA_character_,
    statistic = NA_real_,
    p_value = NA_real_,
    p_adjusted = NA_real_,
    significant = NA,
    interpretation = paste0(
      "Fitted peak at ",
      round(range_mid, 1),
      " bp with fitted mean percentage ",
      round(fitted_percent, 2),
      "%."
    )
  )

statistical_summary_single_sheet <- bind_rows(
  group_profile_summary_rows,
  fitted_peak_rows,
  overall_tests_summary,
  pairwise_tests_summary,
  bin_summary_rows
)

write_tsv(
  statistical_summary_single_sheet,
  file.path(table_dir, "bioanalyzer_statistical_summary_single_sheet.tsv")
)

write_xlsx(
  statistical_summary_single_sheet,
  file.path(table_dir, "bioanalyzer_statistical_summary_single_sheet.xlsx")
)

# ============================================================
# 13. Plot sample-level profile distances
# ============================================================

p_distance <- profile_distance %>%
  pivot_longer(
    cols = c(
      L1_distance_from_healthy,
      L2_distance_from_healthy,
      KS_like_distance_from_healthy
    ),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(
      metric,
      L1_distance_from_healthy = "L1 distance",
      L2_distance_from_healthy = "L2 distance",
      KS_like_distance_from_healthy = "KS-like distance"
    )
  ) %>%
  ggplot(aes(x = group, y = value, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.65, color = "black") +
  geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white", color = "black") +
  geom_jitter(aes(color = group), width = 0.08, size = 2, alpha = 0.8) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_fill_manual(values = group_fills) +
  scale_color_manual(values = group_cols) +
  labs(
    x = "Group",
    y = "Distance from Healthy reference",
    title = "Bioanalyzer profile distance from Healthy reference"
  ) +
  theme_bio(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", family = plot_family)
  )

ggsave(
  file.path(figure_dir, "bioanalyzer_profile_distance_from_healthy.png"),
  p_distance,
  width = 9,
  height = 5.5,
  dpi = 300
)

# ============================================================
# 14. Final report
# ============================================================

cat("\nBioanalyzer fragment-size analysis completed.\n")

cat("\nInput file:\n")
cat(bioanalyzer_file, "\n")

cat("\nSamples used:\n")
print(df_clean %>% distinct(sample_id, group) %>% count(group))

cat("\nPeak summary:\n")
print(peak_summary)

cat("\nOverall statistical summary:\n")
print(overall_tests_summary)

cat("\nOutputs saved to:\n")
cat(table_dir, "\n")
cat(figure_dir, "\n")

cat("\nDone.\n")
