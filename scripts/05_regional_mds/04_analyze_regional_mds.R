# ============================================================
# Regional MDS downstream analysis
#
# Purpose:
#   Analyze regional MDS matrix produced by:
#   scripts/05_regional_mds/03_calculate_regional_mds.R
#
# Input:
#   results/tables/regional_mds/regional_MDS_matrix_samples_by_regions.tsv
#
# Main outputs:
#   Sample-level mean regional MDS summary
#   Kruskal-Wallis and pairwise Wilcoxon tests
#   PCA / t-SNE / UMAP
#   Jonckheere trend test per region
#   Heatmap of top altered regions
#
# Usage:
#   Rscript scripts/05_regional_mds/04_analyze_regional_mds.R
#
# Optional:
#   Rscript scripts/05_regional_mds/04_analyze_regional_mds.R /path/to/regional_MDS_matrix_samples_by_regions.tsv
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(clinfun)
  library(pheatmap)
})

graphics.off()

# Optional packages
has_Rtsne <- requireNamespace("Rtsne", quietly = TRUE)
has_uwot <- requireNamespace("uwot", quietly = TRUE)

# ============================================================
# 1. Paths and settings
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

matrix_file <- ifelse(
  length(args) >= 1,
  args[1],
  "results/tables/regional_mds/regional_MDS_matrix_samples_by_regions.tsv"
)

table_dir <- "results/tables/regional_mds"
figure_dir <- "results/figures/regional_mds"
sanity_dir <- file.path(table_dir, "sanity_checks")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sanity_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c("Healthy", "Remission", "Relapse")

base_font <- "serif"

group_colors <- c(
  Healthy = "darkgreen",
  Remission = "deeppink3",
  Relapse = "blue3"
)

group_fill <- c(
  Healthy = "darkseagreen3",
  Remission = "lightpink",
  Relapse = "lightblue"
)

theme_set(
  theme_classic(base_size = 16, base_family = base_font) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, family = base_font),
      axis.text = element_text(color = "black", family = base_font),
      axis.title = element_text(color = "black", family = base_font),
      legend.title = element_blank(),
      legend.text = element_text(family = base_font)
    )
)

# ============================================================
# 2. Helper functions
# ============================================================

assign_group_from_sample <- function(sample_id) {
  sample_number <- as.integer(str_extract(sample_id, "\\d+"))

  case_when(
    sample_number >= 1  & sample_number <= 30 ~ "Healthy",
    sample_number >= 31 & sample_number <= 60 ~ "Remission",
    sample_number >= 61 & sample_number <= 84 ~ "Relapse",
    TRUE ~ NA_character_
  )
}

p_to_label <- function(p) {
  case_when(
    is.na(p) ~ "NA",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

safe_pairwise_p <- function(pairwise_result, row_name, col_name) {
  mat <- pairwise_result$p.value

  if (is.null(mat)) return(NA_real_)
  if (!(row_name %in% rownames(mat))) return(NA_real_)
  if (!(col_name %in% colnames(mat))) return(NA_real_)

  as.numeric(mat[row_name, col_name])
}

add_violin_sig <- function(p, df, y_col, pairwise_result) {
  y_max <- max(df[[y_col]], na.rm = TRUE)
  y_min <- min(df[[y_col]], na.rm = TRUE)
  y_range <- y_max - y_min

  if (!is.finite(y_range) || y_range == 0) {
    y_range <- abs(y_max) * 0.1 + 1e-6
  }

  p_H_Rem <- safe_pairwise_p(pairwise_result, "Remission", "Healthy")
  p_H_Rel <- safe_pairwise_p(pairwise_result, "Relapse", "Healthy")
  p_Rem_Rel <- safe_pairwise_p(pairwise_result, "Relapse", "Remission")

  lab_H_Rem <- p_to_label(p_H_Rem)
  lab_H_Rel <- p_to_label(p_H_Rel)
  lab_Rem_Rel <- p_to_label(p_Rem_Rel)

  y1 <- y_max + 0.15 * y_range
  y2 <- y_max + 0.22 * y_range
  y3 <- y_max + 0.29 * y_range

  tick <- 0.025 * y_range
  text_offset <- 0.025 * y_range

  p +
    annotate("segment", x = 1, xend = 2, y = y1, yend = y1, linewidth = 0.7) +
    annotate("segment", x = 1, xend = 1, y = y1, yend = y1 - tick, linewidth = 0.7) +
    annotate("segment", x = 2, xend = 2, y = y1, yend = y1 - tick, linewidth = 0.7) +
    annotate("text", x = 1.5, y = y1 + text_offset, label = lab_H_Rem, size = 5, family = base_font) +

    annotate("segment", x = 2, xend = 3, y = y2, yend = y2, linewidth = 0.7) +
    annotate("segment", x = 2, xend = 2, y = y2, yend = y2 - tick, linewidth = 0.7) +
    annotate("segment", x = 3, xend = 3, y = y2, yend = y2 - tick, linewidth = 0.7) +
    annotate("text", x = 2.5, y = y2 + text_offset, label = lab_Rem_Rel, size = 5, family = base_font) +

    annotate("segment", x = 1, xend = 3, y = y3, yend = y3, linewidth = 0.7) +
    annotate("segment", x = 1, xend = 1, y = y3, yend = y3 - tick, linewidth = 0.7) +
    annotate("segment", x = 3, xend = 3, y = y3, yend = y3 - tick, linewidth = 0.7) +
    annotate("text", x = 2, y = y3 + text_offset, label = lab_H_Rel, size = 5, family = base_font) +

    coord_cartesian(
      ylim = c(y_min - 0.15 * y_range, y_max + 0.38 * y_range),
      clip = "off"
    )
}

median_impute_matrix <- function(mat) {
  out <- as.data.frame(mat)

  for (j in seq_len(ncol(out))) {
    med_j <- median(out[[j]], na.rm = TRUE)
    if (!is.finite(med_j)) med_j <- 0
    out[[j]][is.na(out[[j]])] <- med_j
  }

  out
}

# ============================================================
# 3. Read regional MDS matrix
# ============================================================

if (!file.exists(matrix_file)) {
  stop("Regional MDS matrix file not found: ", matrix_file)
}

mds_matrix_df <- read_tsv(
  matrix_file,
  show_col_types = FALSE
)

if (!"sample_id" %in% colnames(mds_matrix_df)) {
  stop("The matrix must contain a column named sample_id.")
}

mds_matrix <- as.data.frame(mds_matrix_df)
rownames(mds_matrix) <- mds_matrix$sample_id
mds_matrix$sample_id <- NULL

mds_matrix <- mds_matrix %>%
  mutate(across(everything(), as.numeric)) %>%
  as.data.frame()

cat("\nInput MDS matrix dimensions:\n")
print(dim(mds_matrix))

cat("\nFirst sample IDs:\n")
print(head(rownames(mds_matrix)))

# ============================================================
# 4. Assign groups
# ============================================================

sample_group_df <- data.frame(
  sample_id = rownames(mds_matrix),
  stringsAsFactors = FALSE
) %>%
  mutate(
    sample_number = as.integer(str_extract(sample_id, "\\d+")),
    group = assign_group_from_sample(sample_id),
    group = factor(group, levels = group_levels)
  )

cat("\nSample group assignment:\n")
print(table(sample_group_df$group, useNA = "ifany"))

if (any(is.na(sample_group_df$group))) {
  print(sample_group_df %>% filter(is.na(group)))
  stop("Some samples could not be assigned to a group. Check sample_id names.")
}

group <- sample_group_df$group
names(group) <- sample_group_df$sample_id

write_tsv(
  sample_group_df,
  file.path(sanity_dir, "sample_group_assignment.tsv")
)

# ============================================================
# 5. Matrix sanity checks
# ============================================================

matrix_summary <- data.frame(
  total_samples = nrow(mds_matrix),
  total_regions = ncol(mds_matrix),
  total_missing_values = sum(is.na(mds_matrix)),
  min_MDS = min(as.matrix(mds_matrix), na.rm = TRUE),
  median_MDS = median(as.matrix(mds_matrix), na.rm = TRUE),
  mean_MDS = mean(as.matrix(mds_matrix), na.rm = TRUE),
  max_MDS = max(as.matrix(mds_matrix), na.rm = TRUE)
)

print(matrix_summary)

write_tsv(
  matrix_summary,
  file.path(table_dir, "regional_mds_matrix_summary.tsv")
)

missing_per_sample <- data.frame(
  sample_id = rownames(mds_matrix),
  n_missing_regions = rowSums(is.na(mds_matrix)),
  percent_missing_regions = 100 * rowMeans(is.na(mds_matrix))
) %>%
  mutate(
    group = group[match(sample_id, names(group))]
  )

missing_per_region <- data.frame(
  region = colnames(mds_matrix),
  n_missing_samples = colSums(is.na(mds_matrix)),
  percent_missing_samples = 100 * colMeans(is.na(mds_matrix))
)

write_tsv(
  missing_per_sample,
  file.path(sanity_dir, "sanity_missing_per_sample.tsv")
)

write_tsv(
  missing_per_region,
  file.path(sanity_dir, "sanity_missing_per_region.tsv")
)

cat("\nMissing per sample summary:\n")
print(summary(missing_per_sample$n_missing_regions))

cat("\nMissing per region summary:\n")
print(summary(missing_per_region$n_missing_samples))

# Keep regions with <= 8 missing values
keep_regions <- colSums(is.na(mds_matrix)) <= 8
mds_matrix_filt <- mds_matrix[, keep_regions, drop = FALSE]

mds_matrix_complete <- mds_matrix_filt[
  ,
  colSums(is.na(mds_matrix_filt)) == 0,
  drop = FALSE
]

cat("\nFiltered matrix dimensions:\n")
print(dim(mds_matrix_filt))

cat("\nComplete matrix dimensions:\n")
print(dim(mds_matrix_complete))

filtering_summary <- data.frame(
  filtering_step = c("raw", "missing_le_8", "complete_only"),
  n_samples = c(
    nrow(mds_matrix),
    nrow(mds_matrix_filt),
    nrow(mds_matrix_complete)
  ),
  n_regions = c(
    ncol(mds_matrix),
    ncol(mds_matrix_filt),
    ncol(mds_matrix_complete)
  )
)

write_tsv(
  filtering_summary,
  file.path(table_dir, "regional_mds_filtering_summary.tsv")
)

# ============================================================
# 6. Long-format table
# ============================================================

all_mds <- mds_matrix_df %>%
  pivot_longer(
    cols = -sample_id,
    names_to = "region",
    values_to = "MDS"
  ) %>%
  mutate(
    group = group[match(sample_id, names(group))],
    group = factor(group, levels = group_levels)
  )

write_tsv(
  all_mds,
  file.path(table_dir, "regional_mds_long_from_matrix.tsv")
)

# ============================================================
# 7. Sample-level regional MDS summary
# ============================================================

sample_summary <- all_mds %>%
  group_by(sample_id, group) %>%
  summarise(
    mean_MDS = mean(MDS, na.rm = TRUE),
    median_MDS = median(MDS, na.rm = TRUE),
    sd_MDS = sd(MDS, na.rm = TRUE),
    n_regions = sum(!is.na(MDS)),
    n_missing = sum(is.na(MDS)),
    .groups = "drop"
  )

write_tsv(
  sample_summary,
  file.path(table_dir, "sample_level_regional_mds_summary.tsv")
)

regional_mds_group_summary <- sample_summary %>%
  group_by(group) %>%
  summarise(
    n_samples = n(),
    mean_MDS = mean(mean_MDS, na.rm = TRUE),
    median_MDS = median(mean_MDS, na.rm = TRUE),
    sd_MDS = sd(mean_MDS, na.rm = TRUE),
    mean_n_regions = mean(n_regions, na.rm = TRUE),
    median_n_regions = median(n_regions, na.rm = TRUE),
    min_n_regions = min(n_regions, na.rm = TRUE),
    max_n_regions = max(n_regions, na.rm = TRUE),
    .groups = "drop"
  )

print(regional_mds_group_summary)

write_tsv(
  regional_mds_group_summary,
  file.path(table_dir, "regional_mds_group_summary.tsv")
)

kw_mean_mds <- kruskal.test(mean_MDS ~ group, data = sample_summary)

pairwise_mean_mds <- pairwise.wilcox.test(
  sample_summary$mean_MDS,
  sample_summary$group,
  p.adjust.method = "BH",
  exact = FALSE
)

sample_level_tests <- tibble(
  variable = "mean_regional_MDS",
  comparison = "Healthy vs Remission vs Relapse",
  test = "Kruskal-Wallis rank-sum test",
  statistic = as.numeric(kw_mean_mds$statistic),
  df = as.numeric(kw_mean_mds$parameter),
  p_value = kw_mean_mds$p.value
)

pairwise_mean_mds_table <- as.data.frame(as.table(pairwise_mean_mds$p.value)) %>%
  rename(
    group2 = Var1,
    group1 = Var2,
    p_adjusted = Freq
  ) %>%
  filter(!is.na(p_adjusted)) %>%
  mutate(
    variable = "mean_regional_MDS",
    test = "Pairwise Wilcoxon rank-sum test",
    p_value = p_adjusted,
    p_adjust_method = "BH",
    significance = p_to_label(p_adjusted)
  ) %>%
  select(variable, test, group1, group2, p_value, p_adjusted, p_adjust_method, significance)

write_tsv(
  sample_level_tests,
  file.path(table_dir, "sample_level_mean_regional_mds_kruskal.tsv")
)

write_tsv(
  pairwise_mean_mds_table,
  file.path(table_dir, "sample_level_mean_regional_mds_pairwise_wilcoxon.tsv")
)

# Violin plot
mean_mds_plot <- ggplot(
  sample_summary,
  aes(x = group, y = mean_MDS, fill = group)
) +
  geom_violin(trim = FALSE, color = "black", alpha = 0.6) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    fill = "white",
    color = "black"
  ) +
  geom_jitter(
    aes(color = group),
    width = 0.1,
    size = 2.2,
    alpha = 0.9
  ) +
  scale_fill_manual(values = group_fill) +
  scale_color_manual(values = group_colors) +
  labs(
    x = "Group",
    y = "Mean regional MDS",
    title = "Sample-level mean regional MDS"
  ) +
  theme(
    legend.position = "none",
    plot.margin = ggplot2::margin(10, 20, 20, 20)
  )

mean_mds_plot <- add_violin_sig(
  p = mean_mds_plot,
  df = sample_summary,
  y_col = "mean_MDS",
  pairwise_result = pairwise_mean_mds
)

ggsave(
  file.path(figure_dir, "sample_level_mean_regional_mds.png"),
  plot = mean_mds_plot,
  width = 7,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(figure_dir, "sample_level_mean_regional_mds.pdf"),
  plot = mean_mds_plot,
  width = 7,
  height = 6,
  device = "pdf"
)

# ============================================================
# 8. Sanity check: mean MDS vs number of valid regions
# ============================================================

cor_meanMDS_nregions <- cor.test(
  sample_summary$mean_MDS,
  sample_summary$n_regions,
  method = "spearman",
  exact = FALSE
)

cor_meanMDS_nregions_table <- tibble(
  test = "Spearman correlation",
  variable_1 = "mean_MDS",
  variable_2 = "n_regions",
  rho = as.numeric(cor_meanMDS_nregions$estimate),
  p_value = cor_meanMDS_nregions$p.value
)

write_tsv(
  cor_meanMDS_nregions_table,
  file.path(sanity_dir, "sanity_cor_meanMDS_nregions.tsv")
)

p_mds_vs_nregions <- ggplot(
  sample_summary,
  aes(x = n_regions, y = mean_MDS, color = group)
) +
  geom_point(size = 2.5, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  scale_color_manual(values = group_colors) +
  labs(
    title = "Mean regional MDS vs number of valid regions",
    x = "Number of valid regions",
    y = "Mean regional MDS",
    color = "Group"
  )

ggsave(
  file.path(figure_dir, "sanity_mean_regional_mds_vs_n_regions.png"),
  p_mds_vs_nregions,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(figure_dir, "sanity_mean_regional_mds_vs_n_regions.pdf"),
  p_mds_vs_nregions,
  width = 7,
  height = 5,
  device = "pdf"
)

# ============================================================
# 9. PCA
# ============================================================

if (ncol(mds_matrix_complete) < 2) {
  warning("Fewer than 2 complete regions. PCA will use median-imputed filtered matrix.")
  pca_input <- median_impute_matrix(mds_matrix_filt)
} else {
  pca_input <- mds_matrix_complete
}

pca_res <- prcomp(pca_input, center = TRUE, scale. = TRUE)

percent_var <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 2)

sample_label_pca <- rownames(pca_input)
sample_label_pca <- gsub("^Cap", "", sample_label_pca)
sample_label_pca <- gsub("\\.dedup$", "", sample_label_pca)

pca_df <- data.frame(
  sample_id = rownames(pca_input),
  sample_label = sample_label_pca,
  group = group[match(rownames(pca_input), names(group))],
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2]
)

write_tsv(
  pca_df,
  file.path(table_dir, "regional_mds_pca_coordinates.tsv")
)

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_text_repel(
    aes(label = sample_label),
    size = 3,
    max.overlaps = Inf,
    show.legend = FALSE,
    family = base_font
  ) +
  stat_ellipse(
    aes(group = group),
    linewidth = 0.8,
    linetype = "dashed",
    show.legend = FALSE
  ) +
  scale_color_manual(values = group_colors) +
  labs(
    x = paste0("PC1 (", percent_var[1], "%)"),
    y = paste0("PC2 (", percent_var[2], "%)"),
    title = "PCA of regional MDS"
  )

ggsave(
  file.path(figure_dir, "regional_mds_pca.png"),
  plot = pca_plot,
  width = 7,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(figure_dir, "regional_mds_pca.pdf"),
  plot = pca_plot,
  width = 7,
  height = 7,
  device = "pdf"
)

# ============================================================
# 10. Jonckheere trend test per region
# Healthy -> Remission -> Relapse
# ============================================================

trend_results_regions <- lapply(colnames(mds_matrix), function(region) {
  x <- mds_matrix[[region]]
  keep <- !is.na(x)

  x_sub <- x[keep]
  g_sub <- as.numeric(group[rownames(mds_matrix)[keep]])

  p_inc <- tryCatch(
    jonckheere.test(x_sub, g_sub, alternative = "increasing")$p.value,
    error = function(e) NA_real_
  )

  p_dec <- tryCatch(
    jonckheere.test(x_sub, g_sub, alternative = "decreasing")$p.value,
    error = function(e) NA_real_
  )

  data.frame(
    region = region,
    p_increasing = p_inc,
    p_decreasing = p_dec
  )
}) %>%
  bind_rows()

trend_results_regions <- trend_results_regions %>%
  mutate(
    padj_increasing = p.adjust(p_increasing, method = "BH"),
    padj_decreasing = p.adjust(p_decreasing, method = "BH"),
    direction = ifelse(
      padj_increasing <= padj_decreasing,
      "Increasing",
      "Decreasing"
    ),
    best_padj = pmin(padj_increasing, padj_decreasing, na.rm = TRUE)
  ) %>%
  arrange(best_padj)

write_tsv(
  trend_results_regions,
  file.path(table_dir, "trend_results_regional_mds.tsv")
)

trend_direction_summary <- trend_results_regions %>%
  count(direction)

write_tsv(
  trend_direction_summary,
  file.path(sanity_dir, "sanity_trend_direction_summary.tsv")
)

region_group_means <- all_mds %>%
  group_by(region, group) %>%
  summarise(
    mean_MDS = mean(MDS, na.rm = TRUE),
    median_MDS = median(MDS, na.rm = TRUE),
    n_samples = sum(!is.na(MDS)),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = group,
    values_from = c(mean_MDS, median_MDS, n_samples)
  ) %>%
  mutate(
    mean_direction = case_when(
      mean_MDS_Healthy > mean_MDS_Remission &
        mean_MDS_Remission > mean_MDS_Relapse ~ "Decreasing",
      mean_MDS_Healthy < mean_MDS_Remission &
        mean_MDS_Remission < mean_MDS_Relapse ~ "Increasing",
      TRUE ~ "Non-monotonic"
    )
  )

write_tsv(
  region_group_means,
  file.path(sanity_dir, "sanity_region_group_mean_directions.tsv")
)

trend_results_regions_checked <- trend_results_regions %>%
  left_join(
    region_group_means %>%
      select(region, mean_direction),
    by = "region"
  ) %>%
  mutate(
    direction_consistent = direction == mean_direction
  )

write_tsv(
  trend_results_regions_checked,
  file.path(sanity_dir, "sanity_trend_results_regions_checked.tsv")
)

sig_region_summary <- trend_results_regions_checked %>%
  summarise(
    total_regions = n(),
    significant_best_padj_0_05 = sum(best_padj < 0.05, na.rm = TRUE),
    significant_increasing = sum(direction == "Increasing" & best_padj < 0.05, na.rm = TRUE),
    significant_decreasing = sum(direction == "Decreasing" & best_padj < 0.05, na.rm = TRUE),
    non_missing_best_padj = sum(!is.na(best_padj))
  )

write_tsv(
  sig_region_summary,
  file.path(sanity_dir, "sanity_significant_region_summary.tsv")
)

print(sig_region_summary)

# ============================================================
# 11. Heatmap of top regional MDS features
# ============================================================

top_n <- 50

top_regions <- trend_results_regions %>%
  filter(!is.na(best_padj)) %>%
  arrange(best_padj) %>%
  slice_head(n = top_n) %>%
  pull(region)

cat("\nNumber of selected regions for heatmap:\n")
print(length(top_regions))

heat_mat <- as.matrix(mds_matrix[, top_regions, drop = FALSE])

# Median imputation for heatmap visualization
for (j in seq_len(ncol(heat_mat))) {
  med_j <- median(heat_mat[, j], na.rm = TRUE)
  if (!is.finite(med_j)) med_j <- 0
  heat_mat[is.na(heat_mat[, j]), j] <- med_j
}

heat_mat_scaled <- scale(heat_mat)

bad_cols <- apply(heat_mat_scaled, 2, function(x) {
  all(is.na(x)) || any(!is.finite(x)) || sd(x, na.rm = TRUE) == 0
})

if (any(bad_cols)) {
  heat_mat_scaled <- heat_mat_scaled[, !bad_cols, drop = FALSE]
}

z_limit <- 2

heat_mat_scaled[heat_mat_scaled > z_limit] <- z_limit
heat_mat_scaled[heat_mat_scaled < -z_limit] <- -z_limit

annotation_row <- data.frame(
  Group = group[match(rownames(heat_mat_scaled), names(group))]
)

rownames(annotation_row) <- rownames(heat_mat_scaled)

ann_colors <- list(
  Group = group_colors
)

breaks <- seq(-z_limit, z_limit, length.out = 101)
heat_colors <- colorRampPalette(c("blue", "white", "red"))(100)

png(
  filename = file.path(figure_dir, "heatmap_top50_regional_MDS.png"),
  width = 3000,
  height = 2400,
  res = 300
)

pheatmap(
  heat_mat_scaled,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_row = annotation_row,
  annotation_colors = ann_colors,
  show_rownames = FALSE,
  show_colnames = TRUE,
  fontsize_col = 7,
  fontsize_row = 7,
  border_color = NA,
  color = heat_colors,
  breaks = breaks,
  main = "Heatmap of top regional MDS features"
)

dev.off()

pdf(
  file = file.path(figure_dir, "heatmap_top50_regional_MDS.pdf"),
  width = 10,
  height = 8,
  family = base_font
)

pheatmap(
  heat_mat_scaled,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_row = annotation_row,
  annotation_colors = ann_colors,
  show_rownames = FALSE,
  show_colnames = TRUE,
  fontsize_col = 7,
  fontsize_row = 7,
  border_color = NA,
  color = heat_colors,
  breaks = breaks,
  main = "Heatmap of top regional MDS features"
)

dev.off()

top_regions_table <- trend_results_regions %>%
  filter(region %in% top_regions) %>%
  arrange(best_padj)

write_tsv(
  top_regions_table,
  file.path(table_dir, "top50_regions_used_for_heatmap.tsv")
)

# ============================================================
# 12. t-SNE and UMAP
# ============================================================

X <- median_impute_matrix(mds_matrix_filt)
X_scaled <- scale(X)

sample_label <- rownames(X)
sample_label <- gsub("^Cap", "", sample_label)
sample_label <- gsub("\\.dedup$", "", sample_label)

# ----------------------------
# t-SNE
# ----------------------------

if (has_Rtsne) {
  set.seed(111)

  tsne_perplexity <- 15

  if (tsne_perplexity >= (nrow(X_scaled) - 1) / 3) {
    tsne_perplexity <- floor((nrow(X_scaled) - 1) / 3)
  }

  tsne_res <- Rtsne::Rtsne(
    X_scaled,
    dims = 2,
    perplexity = tsne_perplexity,
    pca = TRUE,
    max_iter = 1000,
    check_duplicates = FALSE,
    verbose = TRUE
  )

  tsne_df <- data.frame(
    sample_id = rownames(X),
    sample_label = sample_label,
    group = group[match(rownames(X), names(group))],
    tSNE1 = tsne_res$Y[, 1],
    tSNE2 = tsne_res$Y[, 2]
  )

  write_tsv(
    tsne_df,
    file.path(table_dir, "regional_mds_tsne_coordinates.tsv")
  )

  tsne_plot <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = group)) +
    geom_point(size = 3, alpha = 0.9) +
    geom_text_repel(
      aes(label = sample_label),
      size = 3,
      max.overlaps = Inf,
      show.legend = FALSE,
      family = base_font
    ) +
    stat_ellipse(
      aes(group = group),
      linewidth = 0.8,
      linetype = "dashed",
      show.legend = FALSE
    ) +
    scale_color_manual(values = group_colors) +
    labs(
      title = "t-SNE of regional MDS features",
      x = "t-SNE 1",
      y = "t-SNE 2"
    )

  ggsave(
    file.path(figure_dir, "regional_mds_tsne.png"),
    plot = tsne_plot,
    width = 8,
    height = 6,
    dpi = 300
  )

  ggsave(
    file.path(figure_dir, "regional_mds_tsne.pdf"),
    plot = tsne_plot,
    width = 8,
    height = 6,
    device = "pdf"
  )
} else {
  warning("Package Rtsne is not installed. Skipping t-SNE.")
}

# ----------------------------
# UMAP
# ----------------------------

if (has_uwot) {
  set.seed(111)

  umap_res <- uwot::umap(
    X_scaled,
    n_neighbors = 15,
    min_dist = 0.30,
    metric = "euclidean",
    n_components = 2,
    verbose = TRUE
  )

  umap_df <- data.frame(
    sample_id = rownames(X),
    sample_label = sample_label,
    group = group[match(rownames(X), names(group))],
    UMAP1 = umap_res[, 1],
    UMAP2 = umap_res[, 2]
  )

  write_tsv(
    umap_df,
    file.path(table_dir, "regional_mds_umap_coordinates.tsv")
  )

  umap_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = group)) +
    geom_point(size = 3, alpha = 0.9) +
    geom_text_repel(
      aes(label = sample_label),
      size = 3,
      max.overlaps = Inf,
      show.legend = FALSE,
      family = base_font
    ) +
    stat_ellipse(
      aes(group = group),
      linewidth = 0.8,
      linetype = "dashed",
      show.legend = FALSE
    ) +
    scale_color_manual(values = group_colors) +
    labs(
      title = "UMAP of regional MDS features",
      x = "UMAP 1",
      y = "UMAP 2"
    )

  ggsave(
    file.path(figure_dir, "regional_mds_umap.png"),
    plot = umap_plot,
    width = 8,
    height = 6,
    dpi = 300
  )

  ggsave(
    file.path(figure_dir, "regional_mds_umap.pdf"),
    plot = umap_plot,
    width = 8,
    height = 6,
    device = "pdf"
  )
} else {
  warning("Package uwot is not installed. Skipping UMAP.")
}

# ============================================================
# 13. Persian summary table
# ============================================================

regional_mds_summary_table <- sample_summary %>%
  group_by(group) %>%
  summarise(
    `تعداد نمونه` = n(),
    `میانگین MDS ناحیه‌ای` = round(mean(mean_MDS, na.rm = TRUE), 4),
    `میانه MDS ناحیه‌ای` = round(median(mean_MDS, na.rm = TRUE), 4),
    `میانگین تعداد نواحی معتبر` = round(mean(n_regions, na.rm = TRUE), 0),
    .groups = "drop"
  ) %>%
  mutate(
    group = recode(
      as.character(group),
      Healthy = "سالم",
      Remission = "خاموشی بیماری",
      Relapse = "عود بیماری"
    )
  ) %>%
  rename(`گروه` = group)

cor_row <- tibble::tibble(
  `گروه` = "همبستگی میانگین MDS با تعداد نواحی معتبر",
  `تعداد نمونه` = NA_real_,
  `میانگین MDS ناحیه‌ای` = NA_real_,
  `میانه MDS ناحیه‌ای` = NA_real_,
  `میانگین تعداد نواحی معتبر` = NA_real_,
  `توضیح` = paste0(
    "Spearman rho = ",
    round(unname(cor_meanMDS_nregions$estimate), 4),
    ", p-value = ",
    signif(cor_meanMDS_nregions$p.value, 4)
  )
)

regional_mds_summary_table <- regional_mds_summary_table %>%
  mutate(
    `توضیح` = "نمونه‌های باقی‌مانده پس از کنترل کیفیت"
  ) %>%
  bind_rows(cor_row)

write_tsv(
  regional_mds_summary_table,
  file.path(table_dir, "regional_mds_summary_table_persian.tsv")
)

# ============================================================
# 14. Save R objects
# ============================================================

saveRDS(
  mds_matrix,
  file.path(table_dir, "mds_matrix_numeric.rds")
)

saveRDS(
  list(
    mds_matrix = mds_matrix,
    mds_matrix_filt = mds_matrix_filt,
    mds_matrix_complete = mds_matrix_complete,
    all_mds = all_mds,
    sample_summary = sample_summary,
    trend_results_regions = trend_results_regions,
    pca_df = pca_df
  ),
  file.path(table_dir, "regional_mds_analysis_objects.rds")
)

# ============================================================
# 15. Final message
# ============================================================

cat("\n============================================================\n")
cat("Regional MDS analysis completed.\n")
cat("\nInput matrix:\n")
cat(matrix_file, "\n")

cat("\nTables saved in:\n")
cat(table_dir, "\n")

cat("\nFigures saved in:\n")
cat(figure_dir, "\n")

cat("\nMain outputs:\n")
cat("- sample_level_regional_mds_summary.tsv\n")
cat("- trend_results_regional_mds.tsv\n")
cat("- regional_mds_pca.png\n")
cat("- heatmap_top50_regional_MDS.png\n")
if (has_Rtsne) cat("- regional_mds_tsne.png\n")
if (has_uwot) cat("- regional_mds_umap.png\n")

cat("============================================================\n")
