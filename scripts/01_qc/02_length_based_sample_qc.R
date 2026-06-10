# ============================================================
# Length-based sample QC for cfDNA fragmentomics
#
# Purpose:
#   1. Read FinaleToolkit fragment-length output files
#   2. Calculate fragment-length QC metrics per sample
#   3. Assign clinical groups based on Cap sample IDs
#   4. Apply degradation-based QC filtering
#   5. Save QC tables and QC figures
#
# Input:
#   results/intermediate/fragment_length/*.frag_length.tsv
#
# Output:
#   results/tables/qc/
#   results/figures/qc/
#
# Notes:
#   Samples are assigned to groups using sample numbers:
#   Cap01-Cap30 = Healthy
#   Cap31-Cap60 = Remission
#   Cap61-Cap84 = Relapse
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(tidyr)
  library(forcats)
  library(ggrepel)
})

# ============================================================
# 1. User settings
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

input_dir <- ifelse(
  length(args) >= 1,
  args[1],
  "results/intermediate/fragment_length"
)
table_dir  <- "results/tables/qc"
figure_dir <- "results/figures/qc"

dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

# QC thresholds
median_length_cutoff <- 125
mono_fraction_cutoff <- 0.50
short_fraction_cutoff <- 0.50

# ============================================================
# 2. Plot theme and colors
# ============================================================

theme_thesis <- function(base_size = 13) {
  theme_classic(base_size = base_size, base_family = "Times New Roman") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 3),
      plot.subtitle = element_text(hjust = 0.5, size = base_size - 1),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 2),
      legend.title = element_text(size = base_size - 1),
      legend.text = element_text(size = base_size - 2),
      strip.text = element_text(face = "bold", size = base_size - 2)
    )
}

qc_colors <- c(
  "Kept" = "#2C7BB6",
  "Removed" = "#D7191C"
)

group_colors <- c(
  "Healthy" = "darkgreen",
  "Remission" = "deeppink3",
  "Relapse" = "blue3",
  "Unknown" = "grey40"
)

# ============================================================
# 3. Functions
# ============================================================

read_length_file <- function(file_path) {
  df <- read.delim(
    file_path,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
  )

  colnames(df) <- tolower(colnames(df))

  required_cols <- c("min", "max", "count")
  missing_cols <- setdiff(required_cols, colnames(df))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in ",
      basename(file_path),
      ": ",
      paste(missing_cols, collapse = ", ")
    )
  }

  df %>%
    mutate(
      min = as.numeric(min),
      max = as.numeric(max),
      count = as.numeric(count),
      length = (min + max) / 2
    ) %>%
    filter(
      !is.na(length),
      !is.na(count),
      count > 0,
      length < 300
    )
}

weighted_median <- function(length, count) {
  df <- data.frame(length = length, count = count) %>%
    arrange(length) %>%
    mutate(cum_count = cumsum(count))

  total <- sum(df$count, na.rm = TRUE)

  if (total == 0 || nrow(df) == 0) {
    return(NA_real_)
  }

  df$length[which(df$cum_count >= total / 2)[1]]
}

assign_group_from_sample <- function(sample_name) {
  sample_number <- as.integer(str_extract(sample_name, "\\d+"))

  case_when(
    sample_number >= 1  & sample_number <= 30 ~ "Healthy",
    sample_number >= 31 & sample_number <= 60 ~ "Remission",
    sample_number >= 61 & sample_number <= 84 ~ "Relapse",
    TRUE ~ "Unknown"
  )
}

calculate_length_qc <- function(file_path) {
  sample_name <- basename(file_path) %>%
    str_remove("\\.frag_length\\.tsv$") %>%
    str_remove("\\.tsv$") %>%
    str_remove("\\.dedup$")

  df <- read_length_file(file_path)

  total_fragments_below300 <- sum(df$count, na.rm = TRUE)

  if (total_fragments_below300 == 0) {
    warning("No valid fragments below 300 bp for sample: ", sample_name)

    return(tibble(
      sample = sample_name,
      total_fragments_below300 = 0,
      mean_length_below300 = NA_real_,
      median_length_below300 = NA_real_,
      mode_length_below300 = NA_real_,
      fraction_below100 = NA_real_,
      fraction_below120 = NA_real_,
      fraction_below130 = NA_real_,
      fraction_below150 = NA_real_,
      fraction_lt50_below300 = NA_real_,
      fraction_50_120_below300 = NA_real_,
      fraction_100_150_below300 = NA_real_,
      fraction_120_220_below300 = NA_real_,
      fraction_151_220_below300 = NA_real_,
      fraction_180_220_below300 = NA_real_,
      fraction_220_300_below300 = NA_real_,
      short_to_mono_ratio_50_120_vs_120_220 = NA_real_,
      delfi_like_ratio_100_150_vs_151_220 = NA_real_
    ))
  }

  mean_length_below300 <- sum(df$length * df$count, na.rm = TRUE) / total_fragments_below300
  median_length_below300 <- weighted_median(df$length, df$count)
  mode_length_below300 <- df$length[which.max(df$count)]

  n_below100 <- sum(df$count[df$length < 100], na.rm = TRUE)
  n_below120 <- sum(df$count[df$length < 120], na.rm = TRUE)
  n_below130 <- sum(df$count[df$length < 130], na.rm = TRUE)
  n_below150 <- sum(df$count[df$length < 150], na.rm = TRUE)

  n_lt50 <- sum(df$count[df$length < 50], na.rm = TRUE)
  n_50_120 <- sum(df$count[df$length >= 50 & df$length < 120], na.rm = TRUE)
  n_100_150 <- sum(df$count[df$length >= 100 & df$length < 150], na.rm = TRUE)
  n_120_220 <- sum(df$count[df$length >= 120 & df$length < 220], na.rm = TRUE)
  n_151_220 <- sum(df$count[df$length >= 151 & df$length < 220], na.rm = TRUE)
  n_180_220 <- sum(df$count[df$length >= 180 & df$length < 220], na.rm = TRUE)
  n_220_300 <- sum(df$count[df$length >= 220 & df$length < 300], na.rm = TRUE)

  tibble(
    sample = sample_name,
    total_fragments_below300 = total_fragments_below300,

    mean_length_below300 = mean_length_below300,
    median_length_below300 = median_length_below300,
    mode_length_below300 = mode_length_below300,

    fraction_below100 = n_below100 / total_fragments_below300,
    fraction_below120 = n_below120 / total_fragments_below300,
    fraction_below130 = n_below130 / total_fragments_below300,
    fraction_below150 = n_below150 / total_fragments_below300,

    fraction_lt50_below300 = n_lt50 / total_fragments_below300,
    fraction_50_120_below300 = n_50_120 / total_fragments_below300,
    fraction_100_150_below300 = n_100_150 / total_fragments_below300,
    fraction_120_220_below300 = n_120_220 / total_fragments_below300,
    fraction_151_220_below300 = n_151_220 / total_fragments_below300,
    fraction_180_220_below300 = n_180_220 / total_fragments_below300,
    fraction_220_300_below300 = n_220_300 / total_fragments_below300,

    short_to_mono_ratio_50_120_vs_120_220 =
      ifelse(n_120_220 == 0, NA_real_, n_50_120 / n_120_220),

    delfi_like_ratio_100_150_vs_151_220 =
      ifelse(n_151_220 == 0, NA_real_, n_100_150 / n_151_220)
  )
}

save_plot <- function(plot_object, file_name, width = 8, height = 6) {
  ggsave(
    filename = file.path(figure_dir, file_name),
    plot = plot_object,
    width = width,
    height = height,
    dpi = 300
  )
}

# ============================================================
# 4. Read fragment-length files
# ============================================================

files <- list.files(
  input_dir,
  pattern = "\\.frag_length\\.tsv$",
  full.names = TRUE
)

if (length(files) == 0) {
  stop(
    "No .frag_length.tsv files were found in: ",
    input_dir,
    "\nPlace FinaleToolkit fragment-length outputs there or change input_dir."
  )
}

cat("Number of fragment-length files found:", length(files), "\n")

length_qc <- map_dfr(files, calculate_length_qc)

# ============================================================
# 5. Assign groups
# ============================================================

length_qc <- length_qc %>%
  mutate(
    sample_id = as.integer(str_extract(sample, "\\d+")),
    group = assign_group_from_sample(sample),
    group = factor(group, levels = c("Healthy", "Remission", "Relapse", "Unknown"))
  )

cat("\nGroup assignment:\n")
print(length_qc %>% count(group))

unknown_samples <- length_qc %>%
  filter(group == "Unknown") %>%
  select(sample, sample_id)

if (nrow(unknown_samples) > 0) {
  warning("Some samples were assigned to Unknown group.")
  print(unknown_samples)
}

# ============================================================
# 6. Apply degradation-based QC filter
# ============================================================

length_qc <- length_qc %>%
  mutate(
    reason_median_low = median_length_below300 < median_length_cutoff,
    reason_low_mono_fraction = fraction_120_220_below300 < mono_fraction_cutoff,
    reason_high_short_fraction = fraction_below130 > short_fraction_cutoff,

    exclude = reason_median_low |
      reason_low_mono_fraction |
      reason_high_short_fraction,

    qc_status = ifelse(exclude, "Removed", "Kept"),
    qc_status = factor(qc_status, levels = c("Kept", "Removed")),

    removal_reason = case_when(
      reason_median_low & reason_low_mono_fraction & reason_high_short_fraction ~
        "median<125; mono<50%; short<130>50%",
      reason_median_low & reason_low_mono_fraction ~
        "median<125; mono<50%",
      reason_median_low & reason_high_short_fraction ~
        "median<125; short<130>50%",
      reason_low_mono_fraction & reason_high_short_fraction ~
        "mono<50%; short<130>50%",
      reason_median_low ~
        "median<125",
      reason_low_mono_fraction ~
        "mono<50%",
      reason_high_short_fraction ~
        "short<130>50%",
      TRUE ~
        "passed"
    ),

    include_analysis = ifelse(exclude, "no", "yes")
  )

clean_samples <- length_qc %>% filter(!exclude)
removed_samples <- length_qc %>% filter(exclude)

# ============================================================
# 7. Save QC and metadata tables
# ============================================================

final_sample_metadata <- length_qc %>%
  transmute(
    sample,
    sample_id,
    group = as.character(group),
    include_analysis,
    exclusion_reason = ifelse(include_analysis == "yes", NA_character_, removal_reason)
  ) %>%
  arrange(sample_id)

qc_summary_by_group <- length_qc %>%
  group_by(group, qc_status) %>%
  summarise(
    n = n(),
    median_of_median_length = median(median_length_below300, na.rm = TRUE),
    median_fraction_below130 = median(fraction_below130, na.rm = TRUE),
    median_fraction_120_220 = median(fraction_120_220_below300, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(length_qc, file.path(table_dir, "length_QC_all_samples_with_flags.csv"))
write_csv(clean_samples, file.path(table_dir, "length_QC_clean_samples.csv"))
write_csv(removed_samples, file.path(table_dir, "length_QC_removed_samples.csv"))
write_csv(final_sample_metadata, file.path(table_dir, "final_sample_metadata.csv"))
write_csv(qc_summary_by_group, file.path(table_dir, "QC_summary_by_group.csv"))

cat("\n=== QC Filter Summary ===\n")
cat("Total samples: ", nrow(length_qc), "\n")
cat("Clean samples: ", nrow(clean_samples), "\n")
cat("Removed samples: ", nrow(removed_samples), "\n\n")

cat("Breakdown of exclusion reasons:\n")
cat("Median length < 125 bp: ", sum(length_qc$reason_median_low, na.rm = TRUE), "\n")
cat("Mononucleosomal fraction 120-220 bp < 50%: ", sum(length_qc$reason_low_mono_fraction, na.rm = TRUE), "\n")
cat("Short fragment fraction <130 bp > 50%: ", sum(length_qc$reason_high_short_fraction, na.rm = TRUE), "\n\n")

cat("Removed samples:\n")
print(
  removed_samples %>%
    select(
      sample,
      sample_id,
      group,
      median_length_below300,
      fraction_below130,
      fraction_120_220_below300,
      short_to_mono_ratio_50_120_vs_120_220,
      removal_reason
    ) %>%
    arrange(median_length_below300)
)

# ============================================================
# 8. Build long fragment-length distribution table
# ============================================================

length_distribution_long <- map_dfr(files, function(file_path) {
  sample_name <- basename(file_path) %>%
    str_remove("\\.frag_length\\.tsv$") %>%
    str_remove("\\.tsv$") %>%
    str_remove("\\.dedup$")

  read_length_file(file_path) %>%
    mutate(sample = sample_name)
}) %>%
  left_join(
    length_qc %>%
      select(sample, sample_id, group, qc_status, exclude),
    by = "sample"
  )

write_csv(
  length_distribution_long,
  file.path(table_dir, "length_distribution_long_below300.csv")
)

# ============================================================
# 9. QC plots
# ============================================================

p_density_qc <- length_distribution_long %>%
  ggplot(aes(
    x = length,
    weight = count,
    color = qc_status,
    linetype = qc_status
  )) +
  geom_density(linewidth = 1.1, adjust = 1.0, na.rm = TRUE) +
  geom_vline(
    xintercept = c(120, 125, 130, 220),
    linetype = "dashed",
    color = "grey35",
    linewidth = 0.45
  ) +
  scale_color_manual(values = qc_colors) +
  scale_linetype_manual(values = c("Kept" = "solid", "Removed" = "dashed")) +
  labs(
    title = "Fragment Length Distribution by QC Status",
    subtitle = "Dashed vertical lines indicate 120, 125, 130, and 220 bp thresholds",
    x = "Fragment length (bp)",
    y = "Weighted density",
    color = "QC status",
    linetype = "QC status"
  ) +
  theme_thesis(base_size = 13)

save_plot(p_density_qc, "fragment_length_density_by_qc_status.png", width = 8, height = 6)

p_density_group <- length_distribution_long %>%
  ggplot(aes(
    x = length,
    weight = count,
    color = group,
    linetype = group
  )) +
  geom_density(linewidth = 1.0, adjust = 1.0, na.rm = TRUE) +
  geom_vline(
    xintercept = c(120, 220),
    linetype = "dashed",
    color = "grey35",
    linewidth = 0.45
  ) +
  scale_color_manual(values = group_colors) +
  labs(
    title = "Fragment Length Distribution by Clinical Group",
    subtitle = "Dashed vertical lines indicate the 120-220 bp mononucleosomal window",
    x = "Fragment length (bp)",
    y = "Weighted density",
    color = "Group",
    linetype = "Group"
  ) +
  theme_thesis(base_size = 13)

save_plot(p_density_group, "fragment_length_density_by_group.png", width = 8, height = 6)

length_qc_plot <- length_qc %>%
  mutate(sample_label = as.character(sample_id))

p_qc_scatter <- length_qc_plot %>%
  ggplot(aes(
    x = median_length_below300,
    y = fraction_120_220_below300,
    color = qc_status,
    shape = qc_status,
    size = fraction_below130
  )) +
  geom_point(alpha = 0.85) +
  geom_vline(
    xintercept = median_length_cutoff,
    linetype = "dashed",
    color = "grey35",
    linewidth = 0.4
  ) +
  geom_hline(
    yintercept = mono_fraction_cutoff,
    linetype = "dashed",
    color = "grey35",
    linewidth = 0.45
  ) +
  geom_text_repel(
    data = subset(length_qc_plot, qc_status == "Removed"),
    aes(label = sample_label),
    size = 3.5,
    family = "Times New Roman",
    color = "black",
    box.padding = 0.3,
    point.padding = 0.25,
    segment.color = "grey40",
    show.legend = FALSE,
    max.overlaps = Inf
  ) +
  scale_color_manual(values = qc_colors) +
  labs(
    title = "Length-Based Degradation QC",
    subtitle = "Point size represents the fraction of fragments shorter than 130 bp",
    x = "Weighted median fragment length below 300 bp",
    y = "Fraction of fragments in 120-220 bp",
    color = "QC status",
    shape = "QC status",
    size = "Fraction <130 bp"
  ) +
  theme_thesis(base_size = 13)

save_plot(p_qc_scatter, "length_based_qc_scatter.png", width = 8, height = 6)

# ============================================================
# 10. PCA of QC metrics
# ============================================================

pca_features <- c(
  "mean_length_below300",
  "median_length_below300",
  "mode_length_below300",
  "fraction_below100",
  "fraction_below120",
  "fraction_below130",
  "fraction_below150",
  "fraction_50_120_below300",
  "fraction_100_150_below300",
  "fraction_120_220_below300",
  "fraction_151_220_below300",
  "fraction_180_220_below300",
  "fraction_220_300_below300",
  "short_to_mono_ratio_50_120_vs_120_220",
  "delfi_like_ratio_100_150_vs_151_220"
)

pca_input_all <- length_qc %>%
  select(sample, sample_id, group, qc_status, all_of(pca_features)) %>%
  drop_na()

if (nrow(pca_input_all) >= 3) {
  pca_matrix_all <- pca_input_all %>%
    select(all_of(pca_features)) %>%
    as.matrix()

  pca_all <- prcomp(pca_matrix_all, center = TRUE, scale. = TRUE)

  pca_scores_all <- as.data.frame(pca_all$x[, 1:2]) %>%
    mutate(
      sample = pca_input_all$sample,
      sample_id = pca_input_all$sample_id,
      group = pca_input_all$group,
      qc_status = pca_input_all$qc_status,
      sample_label = as.character(sample_id)
    )

  var_explained_all <- summary(pca_all)$importance[2, 1:2] * 100

  p_pca_all <- pca_scores_all %>%
    ggplot(aes(x = PC1, y = PC2, color = qc_status, shape = qc_status)) +
    geom_point(size = 3, alpha = 0.85) +
    geom_text_repel(
      data = subset(pca_scores_all, qc_status == "Removed"),
      aes(label = sample_label),
      size = 3.5,
      family = "Times New Roman",
      color = "black",
      box.padding = 0.35,
      point.padding = 0.25,
      segment.color = "grey45",
      show.legend = FALSE,
      max.overlaps = Inf
    ) +
    scale_color_manual(values = qc_colors) +
    labs(
      title = "PCA of Fragment Length QC Metrics",
      subtitle = "All samples before degradation-based filtering",
      x = paste0("PC1 (", round(var_explained_all[1], 1), "%)"),
      y = paste0("PC2 (", round(var_explained_all[2], 1), "%)"),
      color = "QC status",
      shape = "QC status"
    ) +
    theme_thesis(base_size = 13)

  save_plot(p_pca_all, "pca_length_metrics_before_qc.png", width = 8, height = 6)

  pca_loadings_all <- as.data.frame(pca_all$rotation[, 1:2]) %>%
    rownames_to_column("feature") %>%
    arrange(desc(abs(PC1)))

  write_csv(
    pca_loadings_all,
    file.path(table_dir, "PCA_length_metrics_loadings_all_samples.csv")
  )
}

pca_input_clean <- length_qc %>%
  filter(!exclude) %>%
  select(sample, sample_id, group, qc_status, all_of(pca_features)) %>%
  drop_na()

if (nrow(pca_input_clean) >= 3) {
  pca_matrix_clean <- pca_input_clean %>%
    select(all_of(pca_features)) %>%
    as.matrix()

  pca_clean <- prcomp(pca_matrix_clean, center = TRUE, scale. = TRUE)

  pca_scores_clean <- as.data.frame(pca_clean$x[, 1:2]) %>%
    mutate(
      sample = pca_input_clean$sample,
      sample_id = pca_input_clean$sample_id,
      group = pca_input_clean$group,
      qc_status = pca_input_clean$qc_status
    )

  var_explained_clean <- summary(pca_clean)$importance[2, 1:2] * 100

  p_pca_clean <- pca_scores_clean %>%
    ggplot(aes(x = PC1, y = PC2, color = group, shape = group)) +
    geom_point(size = 3, alpha = 0.85) +
    scale_color_manual(values = group_colors) +
    labs(
      title = "PCA of Fragment Length QC Metrics",
      subtitle = "QC-passed samples only",
      x = paste0("PC1 (", round(var_explained_clean[1], 1), "%)"),
      y = paste0("PC2 (", round(var_explained_clean[2], 1), "%)"),
      color = "Group",
      shape = "Group"
    ) +
    theme_thesis(base_size = 13)

  save_plot(p_pca_clean, "pca_length_metrics_after_qc.png", width = 8, height = 6)
}

# ============================================================
# 11. Boxplot of key degradation metrics
# ============================================================

metric_labels <- c(
  "fraction_120_220_below300" = "Mononucleosomal fraction\n120-220 bp",
  "fraction_below130" = "Short-fragment fraction\n<130 bp",
  "median_length_below300" = "Weighted median\nfragment length",
  "short_to_mono_ratio_50_120_vs_120_220" = "Short-to-mono ratio\n50-120 / 120-220 bp"
)

qc_metrics_long <- length_qc %>%
  select(
    sample,
    sample_id,
    group,
    qc_status,
    median_length_below300,
    fraction_below130,
    fraction_120_220_below300,
    short_to_mono_ratio_50_120_vs_120_220
  ) %>%
  pivot_longer(
    cols = c(
      median_length_below300,
      fraction_below130,
      fraction_120_220_below300,
      short_to_mono_ratio_50_120_vs_120_220
    ),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric_label = recode(metric, !!!metric_labels),
    metric_label = factor(metric_label, levels = unname(metric_labels))
  )

p_qc_box <- qc_metrics_long %>%
  ggplot(aes(x = qc_status, y = value, color = qc_status)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.7) +
  geom_jitter(width = 0.12, alpha = 0.75, size = 1.8) +
  facet_wrap(~ metric_label, scales = "free_y") +
  scale_color_manual(values = qc_colors) +
  labs(
    title = "Key Degradation Metrics by QC Status",
    x = "QC status",
    y = "Metric value",
    color = "QC status"
  ) +
  theme_thesis(base_size = 13)

save_plot(p_qc_box, "key_degradation_metrics_by_qc_status.png", width = 9, height = 6)

# ============================================================
# 12. Sensitivity analysis for QC thresholds
# ============================================================

threshold_grid <- expand_grid(
  median_cutoff = c(120, 125, 130),
  mono_fraction_cutoff = c(0.45, 0.50, 0.55),
  short_fraction_cutoff = c(0.45, 0.50, 0.55)
)

sensitivity_results <- threshold_grid %>%
  mutate(
    n_removed = pmap_int(
      list(median_cutoff, mono_fraction_cutoff, short_fraction_cutoff),
      function(median_cutoff, mono_fraction_cutoff, short_fraction_cutoff) {
        sum(
          length_qc$median_length_below300 < median_cutoff |
            length_qc$fraction_120_220_below300 < mono_fraction_cutoff |
            length_qc$fraction_below130 > short_fraction_cutoff,
          na.rm = TRUE
        )
      }
    ),
    n_kept = nrow(length_qc) - n_removed,
    percent_removed = 100 * n_removed / nrow(length_qc)
  )

write_csv(
  sensitivity_results,
  file.path(table_dir, "QC_threshold_sensitivity_analysis.csv")
)

p_sensitivity <- sensitivity_results %>%
  mutate(
    threshold_set = paste0(
      "median<", median_cutoff,
      "; mono<", mono_fraction_cutoff,
      "; short>", short_fraction_cutoff
    )
  ) %>%
  arrange(percent_removed) %>%
  mutate(threshold_set = fct_inorder(threshold_set)) %>%
  ggplot(aes(x = threshold_set, y = percent_removed)) +
  geom_col(fill = "lavenderblush3") +
  coord_flip() +
  labs(
    title = "Sensitivity Analysis of Degradation-Based QC Thresholds",
    x = "Threshold combination",
    y = "Removed samples (%)"
  ) +
  theme_thesis(base_size = 11)

save_plot(p_sensitivity, "qc_threshold_sensitivity.png", width = 9, height = 7)

# ============================================================
# 13. Final message
# ============================================================

cat("\nAll QC tables were saved to:\n")
cat(table_dir, "\n")

cat("\nAll QC figures were saved to:\n")
cat(figure_dir, "\n")
