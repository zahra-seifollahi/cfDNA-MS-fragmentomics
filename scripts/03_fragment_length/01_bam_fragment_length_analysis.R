# ============================================================
# BAM-derived fragment length distribution analysis
#
# Purpose:
#   Analyze FinaleToolkit fragment-length output files derived
#   from BAM/fragment files.
#
# Important:
#   This script uses QC-passed samples only.
#
# Input:
#   .frag_length.tsv files
#   results/tables/qc/final_sample_metadata.csv
#
# Usage:
#   Rscript scripts/03_fragment_length/01_bam_fragment_length_analysis.R /path/to/frag_length_tsv
#
# If no input path is provided, default is:
#   results/intermediate/fragment_length
#
# Output:
#   results/tables/fragment_length/
#   results/figures/fragment_length/
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
})

graphics.off()

# ============================================================
# 1. Settings
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

frag_path <- ifelse(
  length(args) >= 1,
  args[1],
  "results/intermediate/fragment_length"
)

qc_metadata_file <- "results/tables/qc/final_sample_metadata.csv"

table_dir <- "results/tables/fragment_length"
figure_dir <- "results/figures/fragment_length"

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c("Healthy", "Remission", "Relapse")

group_colors <- c(
  Healthy = "darkgreen",
  Remission = "deeppink3",
  Relapse = "blue3"
)

box_cols <- c(
  Healthy = "darkseagreen3",
  Remission = "lightpink",
  Relapse = "lightblue"
)

point_cols <- c(
  Healthy = "darkgreen",
  Remission = "deeppink3",
  Relapse = "blue3"
)

ms_colors <- c(
  Healthy = "cornflowerblue",
  MS = "tomato"
)

plot_family <- "serif"

base_theme <- theme_classic(base_size = 16, base_family = plot_family) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, family = plot_family),
    plot.subtitle = element_text(hjust = 0.5, family = plot_family),
    axis.text = element_text(color = "black", family = plot_family),
    axis.title = element_text(color = "black", family = plot_family),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(family = plot_family)
  )

# ============================================================
# 2. Helper functions
# ============================================================

clean_sample_name <- function(file_path) {
  basename(file_path) %>%
    str_remove("\\.frag_length\\.tsv$") %>%
    str_remove("\\.tsv$") %>%
    str_remove("\\.dedup$")
}

assign_group_from_sample <- function(sample_name) {
  sample_num <- as.integer(str_extract(sample_name, "\\d+"))

  case_when(
    sample_num >= 1  & sample_num <= 30 ~ "Healthy",
    sample_num >= 31 & sample_num <= 60 ~ "Remission",
    sample_num >= 61 & sample_num <= 84 ~ "Relapse",
    TRUE ~ NA_character_
  )
}

smooth7 <- function(x) {
  s <- stats::filter(x, rep(1 / 7, 7), sides = 2)
  s[is.na(s)] <- x[is.na(s)]
  as.numeric(s)
}

weighted_median <- function(x, w) {
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  x[which(cumsum(w) / sum(w) >= 0.5)[1]]
}

peak_in <- function(df, lo, hi) {
  df %>%
    filter(length >= lo, length <= hi) %>%
    group_by(group) %>%
    slice_max(mean_freq_smooth, n = 1, with_ties = FALSE) %>%
    ungroup()
}

make_labels <- function(peak_df, x_offsets, y_offsets) {
  peak_df %>%
    mutate(
      label = paste0(length, " bp"),
      x_lab = length + x_offsets[match(as.character(group), group_levels)],
      y_lab = mean_freq_smooth + y_offsets[match(as.character(group), group_levels)]
    )
}

label_layers <- function(lab_df) {
  geom_text(
    data = lab_df,
    aes(x = x_lab, y = y_lab, label = label),
    size = 4,
    family = plot_family,
    show.legend = FALSE
  )
}

format_p <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p < 0.001, "p < 0.001", paste0("p = ", signif(p, 2))))
}

p_to_star <- function(p) {
  case_when(
    is.na(p) ~ "NA",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

get_pw <- function(pw, row, col) {
  value <- tryCatch(
    pw$p.value[row, col],
    error = function(e) NA_real_
  )

  as.numeric(value)
}

pairwise_to_table <- function(pairwise_result, variable_name) {
  as.data.frame(as.table(pairwise_result$p.value)) %>%
    rename(
      group1 = Var1,
      group2 = Var2,
      p_value_BH = Freq
    ) %>%
    filter(!is.na(p_value_BH)) %>%
    mutate(
      variable = variable_name,
      test = "Pairwise Wilcoxon rank-sum test",
      p_adjust_method = "BH",
      significance = p_to_star(p_value_BH)
    ) %>%
    select(variable, test, group1, group2, p_value_BH, p_adjust_method, significance)
}

add_sig_bar <- function(x1, x2, y, h = 1, label, cex = 1) {
  segments(x1, y, x2, y, lwd = 1.5)
  segments(c(x1, x2), y, c(x1, x2), y - h, lwd = 1.5)
  text((x1 + x2) / 2, y + 3, labels = label, cex = cex, family = plot_family)
}

# ============================================================
# 3. Load QC metadata
# ============================================================

if (!file.exists(qc_metadata_file)) {
  stop(
    "QC metadata file not found: ",
    qc_metadata_file,
    "\nRun scripts/01_qc/02_length_based_sample_qc.R first."
  )
}

qc_metadata <- read_csv(qc_metadata_file, show_col_types = FALSE) %>%
  mutate(sample = as.character(sample))

required_qc_cols <- c("sample", "include_analysis")
missing_qc_cols <- setdiff(required_qc_cols, colnames(qc_metadata))

if (length(missing_qc_cols) > 0) {
  stop("QC metadata is missing columns: ", paste(missing_qc_cols, collapse = ", "))
}

included_samples <- qc_metadata %>%
  filter(include_analysis == "yes") %>%
  pull(sample)

cat("\nQC-passed samples:", length(included_samples), "\n")

# ============================================================
# 4. Load FinaleToolkit fragment-length files
# ============================================================

frag_files <- list.files(
  frag_path,
  pattern = "\\.frag_length\\.tsv$",
  full.names = TRUE
)

if (length(frag_files) == 0) {
  stop("No .frag_length.tsv files found in: ", frag_path)
}

frag_df <- frag_files %>%
  lapply(function(f) {
    sample_name <- clean_sample_name(f)

    read_tsv(f, comment = "#", show_col_types = FALSE) %>%
      mutate(sample = sample_name)
  }) %>%
  bind_rows()

required_cols <- c("min", "max", "count")
missing_cols <- setdiff(required_cols, colnames(frag_df))

if (length(missing_cols) > 0) {
  stop("Fragment-length files are missing columns: ", paste(missing_cols, collapse = ", "))
}

frag_df <- frag_df %>%
  mutate(
    length = ifelse(min == max, min, (min + max) / 2),
    sample_num = as.integer(str_extract(sample, "\\d+")),
    group = assign_group_from_sample(sample),
    group = factor(group, levels = group_levels)
  ) %>%
  filter(!is.na(group)) %>%
  filter(sample %in% included_samples) %>%
  group_by(sample) %>%
  mutate(freq = count / sum(count, na.rm = TRUE)) %>%
  ungroup()

cat("\nSamples loaded after QC filtering:\n")
print(frag_df %>% distinct(sample, group) %>% count(group))

write_tsv(
  frag_df,
  file.path(table_dir, "bam_fragment_length_long_data.tsv")
)

# ============================================================
# 5. Group-mean smoothed distributions
# ============================================================

group_freq_smooth <- frag_df %>%
  group_by(group, length) %>%
  summarise(mean_freq = mean(freq, na.rm = TRUE), .groups = "drop") %>%
  group_by(group) %>%
  arrange(length, .by_group = TRUE) %>%
  mutate(mean_freq_smooth = smooth7(mean_freq)) %>%
  ungroup()

write_tsv(
  group_freq_smooth,
  file.path(table_dir, "bam_fragment_length_group_mean_smooth.tsv")
)

# ============================================================
# 6. Landmark peaks
# ============================================================

shoulder_df <- peak_in(group_freq_smooth, 80, 105)
peak1_df <- peak_in(group_freq_smooth, 132, 145)
peak2_df <- peak_in(group_freq_smooth, 155, 170)

shoulder_lab <- make_labels(
  shoulder_df,
  x_offsets = c(0, -3, -3),
  y_offsets = c(3e-4, 4e-4, -4e-4)
)

peak1_lab <- make_labels(
  peak1_df,
  x_offsets = c(-8.5, 2, -2),
  y_offsets = c(2e-4, 7e-4, 4e-4)
)

peak2_lab <- make_labels(
  peak2_df,
  x_offsets = c(17.5, 8, 6),
  y_offsets = c(4e-4, 6e-4, 7.5e-4)
)

peak_summary <- bind_rows(
  shoulder_df %>% mutate(peak_region = "80-105 bp shoulder"),
  peak1_df %>% mutate(peak_region = "132-145 bp peak"),
  peak2_df %>% mutate(peak_region = "155-170 bp peak")
) %>%
  select(peak_region, group, length, mean_freq_smooth)

write_tsv(
  peak_summary,
  file.path(table_dir, "bam_fragment_length_landmark_peaks.tsv")
)

# ============================================================
# 7. Plot A: group mean curves with landmark labels
# ============================================================

plot_group_mean <- ggplot(
  group_freq_smooth,
  aes(x = length, y = mean_freq_smooth, color = group)
) +
  geom_line(linewidth = 1.3) +
  geom_point(data = shoulder_df, size = 2.3, show.legend = FALSE) +
  geom_point(data = peak1_df, size = 2.3, show.legend = FALSE) +
  geom_point(data = peak2_df, size = 2.3, show.legend = FALSE) +
  label_layers(shoulder_lab) +
  label_layers(peak1_lab) +
  label_layers(peak2_lab) +
  scale_color_manual(values = group_colors) +
  coord_cartesian(xlim = c(35, 300)) +
  labs(
    x = "Fragment length (bp)",
    y = "Mean relative frequency",
    title = "Normalized BAM-derived fragment length distribution by group"
  ) +
  base_theme

ggsave(
  filename = file.path(figure_dir, "bam_fragment_length_group_mean.png"),
  plot = plot_group_mean,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  filename = file.path(figure_dir, "bam_fragment_length_group_mean.pdf"),
  plot = plot_group_mean,
  width = 8,
  height = 5,
  device = "pdf"
)

# ============================================================
# 8. Plot B: per-sample lines
# ============================================================

plot_per_sample <- ggplot(
  frag_df,
  aes(x = length, y = freq, group = sample, color = group)
) +
  geom_line(alpha = 0.35, linewidth = 0.5) +
  scale_color_manual(values = group_colors) +
  coord_cartesian(xlim = c(35, 300)) +
  labs(
    x = "Fragment length (bp)",
    y = "Relative frequency",
    title = "BAM-derived fragment length distribution per sample"
  ) +
  base_theme

ggsave(
  filename = file.path(figure_dir, "bam_fragment_length_per_sample_lines.png"),
  plot = plot_per_sample,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  filename = file.path(figure_dir, "bam_fragment_length_per_sample_lines.pdf"),
  plot = plot_per_sample,
  width = 8,
  height = 5,
  device = "pdf"
)

# ============================================================
# 9. Plot C: faceted per-sample plot
# ============================================================

plot_facet <- ggplot(
  frag_df,
  aes(x = length, y = freq, color = group)
) +
  geom_line(linewidth = 0.6, show.legend = FALSE) +
  facet_wrap(~ sample, scales = "free_y") +
  scale_color_manual(values = group_colors) +
  coord_cartesian(xlim = c(35, 400)) +
  labs(
    x = "Fragment length (bp)",
    y = "Relative frequency",
    title = "Per-sample BAM-derived fragment length distributions"
  ) +
  theme_classic(base_size = 12, base_family = plot_family) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, family = plot_family),
    axis.text = element_text(color = "black", family = plot_family),
    axis.title = element_text(color = "black", family = plot_family),
    strip.text = element_text(family = plot_family)
  )

ggsave(
  filename = file.path(figure_dir, "bam_fragment_length_per_sample_facets.png"),
  plot = plot_facet,
  width = 12,
  height = 10,
  dpi = 300
)

ggsave(
  filename = file.path(figure_dir, "bam_fragment_length_per_sample_facets.pdf"),
  plot = plot_facet,
  width = 12,
  height = 10,
  device = "pdf"
)

# ============================================================
# 10. Plot D: Healthy vs MS area plot
# ============================================================

ms_smooth <- frag_df %>%
  mutate(
    group2 = ifelse(group == "Healthy", "Healthy", "MS"),
    group2 = factor(group2, levels = c("Healthy", "MS"))
  ) %>%
  group_by(group2, length) %>%
  summarise(mean_freq = mean(freq, na.rm = TRUE), .groups = "drop") %>%
  group_by(group2) %>%
  arrange(length, .by_group = TRUE) %>%
  mutate(mean_freq_smooth = smooth7(mean_freq)) %>%
  ungroup()

plot_ms <- ggplot(
  ms_smooth,
  aes(x = length, y = mean_freq_smooth, fill = group2, color = group2)
) +
  geom_area(alpha = 0.20, position = "identity", linewidth = 0) +
  geom_line(linewidth = 1) +
  scale_fill_manual(values = ms_colors) +
  scale_color_manual(values = ms_colors) +
  scale_x_continuous(limits = c(20, 300), expand = c(0, 0)) +
  labs(
    x = "Fragment length (bp)",
    y = NULL,
    title = "BAM-derived fragment length distribution: Healthy vs MS"
  ) +
  theme_classic(base_size = 15, base_family = plot_family) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, family = plot_family),
    legend.position = c(0.18, 0.88),
    legend.title = element_blank(),
    legend.background = element_blank(),
    legend.text = element_text(family = plot_family),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.text = element_text(color = "black", family = plot_family),
    axis.title.x = element_text(color = "black", family = plot_family)
  )

ggsave(
  filename = file.path(figure_dir, "bam_fragment_length_healthy_vs_ms_area.png"),
  plot = plot_ms,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  filename = file.path(figure_dir, "bam_fragment_length_healthy_vs_ms_area.pdf"),
  plot = plot_ms,
  width = 7,
  height = 5,
  device = "pdf"
)

# ============================================================
# 11. Sample-level fragment statistics
# ============================================================

sample_stats <- frag_df %>%
  group_by(sample, group) %>%
  summarise(
    mean_length = weighted.mean(length, count),
    median_length = weighted_median(length, count),
    short_fraction_lt150 = sum(count[length < 150], na.rm = TRUE) / sum(count, na.rm = TRUE),
    fraction_120_220 = sum(count[length >= 120 & length <= 220], na.rm = TRUE) / sum(count, na.rm = TRUE),
    fraction_lt130 = sum(count[length < 130], na.rm = TRUE) / sum(count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    group = factor(group, levels = group_levels)
  )

write_tsv(
  sample_stats,
  file.path(table_dir, "bam_fragment_length_sample_statistics.tsv")
)

# ============================================================
# 12. Statistical tests
# ============================================================

kruskal_short <- kruskal.test(short_fraction_lt150 ~ group, data = sample_stats)
kruskal_mean <- kruskal.test(mean_length ~ group, data = sample_stats)
kruskal_median <- kruskal.test(median_length ~ group, data = sample_stats)

pairwise_short <- pairwise.wilcox.test(
  sample_stats$short_fraction_lt150,
  sample_stats$group,
  p.adjust.method = "BH",
  exact = FALSE
)

pairwise_mean <- pairwise.wilcox.test(
  sample_stats$mean_length,
  sample_stats$group,
  p.adjust.method = "BH",
  exact = FALSE
)

pairwise_median <- pairwise.wilcox.test(
  sample_stats$median_length,
  sample_stats$group,
  p.adjust.method = "BH",
  exact = FALSE
)

kruskal_results <- tibble(
  variable = c("short_fraction_lt150", "mean_length", "median_length"),
  test = "Kruskal-Wallis rank-sum test",
  statistic = c(
    as.numeric(kruskal_short$statistic),
    as.numeric(kruskal_mean$statistic),
    as.numeric(kruskal_median$statistic)
  ),
  df = c(
    as.numeric(kruskal_short$parameter),
    as.numeric(kruskal_mean$parameter),
    as.numeric(kruskal_median$parameter)
  ),
  p_value = c(
    kruskal_short$p.value,
    kruskal_mean$p.value,
    kruskal_median$p.value
  ),
  significance = p_to_star(p_value)
)

pairwise_all <- bind_rows(
  pairwise_to_table(pairwise_short, "short_fraction_lt150"),
  pairwise_to_table(pairwise_mean, "mean_length"),
  pairwise_to_table(pairwise_median, "median_length")
)

write_tsv(
  kruskal_results,
  file.path(table_dir, "bam_fragment_length_kruskal_results.tsv")
)

write_tsv(
  pairwise_all,
  file.path(table_dir, "bam_fragment_length_pairwise_wilcoxon_results.tsv")
)

# ============================================================
# 13. Base R boxplot: mean fragment length
# ============================================================

p_H_R <- get_pw(pairwise_mean, "Remission", "Healthy")
p_Rel_Rem <- get_pw(pairwise_mean, "Relapse", "Remission")
p_H_Rel <- get_pw(pairwise_mean, "Relapse", "Healthy")

ymax <- max(sample_stats$mean_length, na.rm = TRUE)
ymin <- min(sample_stats$mean_length, na.rm = TRUE)

png(
  filename = file.path(figure_dir, "bam_mean_fragment_length_boxplot.png"),
  width = 7,
  height = 5,
  units = "in",
  res = 300
)

par(mar = c(5, 5, 3, 2), family = plot_family)

boxplot(
  mean_length ~ group,
  data = sample_stats,
  col = box_cols,
  border = "black",
  outline = FALSE,
  xaxt = "n",
  yaxt = "n",
  xlab = "Group",
  ylab = "Mean fragment length (bp)",
  cex.lab = 1.4,
  lwd = 2,
  ylim = c(ymin - 5, ymax + 25)
)

axis(2, las = 1, cex.axis = 1.1, lwd = 1.2)

axis(
  1,
  at = 1:3,
  labels = levels(sample_stats$group),
  cex.axis = 1.1,
  lwd = 1.2
)

for (i in seq_along(levels(sample_stats$group))) {
  grp <- levels(sample_stats$group)[i]
  y <- sample_stats$mean_length[sample_stats$group == grp]

  points(
    jitter(rep(i, length(y)), amount = 0.12),
    y,
    pch = 16,
    cex = 1,
    col = point_cols[grp]
  )
}

add_sig_bar(1, 2, ymax + 8, label = format_p(p_H_R))
add_sig_bar(2, 3, ymax + 9, label = format_p(p_Rel_Rem))
add_sig_bar(1, 3, ymax + 15, label = format_p(p_H_Rel))

dev.off()

pdf(
  file = file.path(figure_dir, "bam_mean_fragment_length_boxplot.pdf"),
  width = 7,
  height = 5,
  family = plot_family
)

par(mar = c(5, 5, 3, 2), family = plot_family)

boxplot(
  mean_length ~ group,
  data = sample_stats,
  col = box_cols,
  border = "black",
  outline = FALSE,
  xaxt = "n",
  yaxt = "n",
  xlab = "Group",
  ylab = "Mean fragment length (bp)",
  cex.lab = 1.4,
  lwd = 2,
  ylim = c(ymin - 5, ymax + 25)
)

axis(2, las = 1, cex.axis = 1.1, lwd = 1.2)

axis(
  1,
  at = 1:3,
  labels = levels(sample_stats$group),
  cex.axis = 1.1,
  lwd = 1.2
)

for (i in seq_along(levels(sample_stats$group))) {
  grp <- levels(sample_stats$group)[i]
  y <- sample_stats$mean_length[sample_stats$group == grp]

  points(
    jitter(rep(i, length(y)), amount = 0.12),
    y,
    pch = 16,
    cex = 1,
    col = point_cols[grp]
  )
}

add_sig_bar(1, 2, ymax + 8, label = format_p(p_H_R))
add_sig_bar(2, 3, ymax + 9, label = format_p(p_Rel_Rem))
add_sig_bar(1, 3, ymax + 15, label = format_p(p_H_Rel))

dev.off()

# ============================================================
# 14. Final report
# ============================================================

cat("\nBAM-derived fragment length analysis completed.\n")

cat("\nInput folder:\n")
cat(frag_path, "\n")

cat("\nSamples used after QC:\n")
print(sample_stats %>% count(group))

cat("\nKruskal-Wallis results:\n")
print(kruskal_results)

cat("\nOutputs saved to:\n")
cat(table_dir, "\n")
cat(figure_dir, "\n")
cat("\nDone.\n")
