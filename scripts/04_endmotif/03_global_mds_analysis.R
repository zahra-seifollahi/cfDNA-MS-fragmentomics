# ============================================================
# Global Motif Diversity Score (MDS) analysis
#
# This script performs:
#   1. Shannon entropy-based MDS calculation
#   2. Healthy-reference MDS z-score calculation
#   3. Group summary
#   4. Kruskal-Wallis test
#   5. Pairwise Wilcoxon tests with BH correction
#   6. MDS violin plots
#   7. MDS z-score violin plots
#
# Input:
#   results/intermediate/endmotif/global_endmotif_matrix.tsv
#   results/tables/endmotif/endmotif_sample_summary.tsv
#
# Output:
#   results/tables/endmotif/global_mds_scores.tsv
#   results/tables/endmotif/global_mds_summary_by_group.tsv
#   results/tables/endmotif/global_mds_kruskal_results.tsv
#   results/tables/endmotif/global_mds_pairwise_wilcoxon_results.tsv
#   results/figures/endmotif/global_mds_violin.png
#   results/figures/endmotif/global_mds_z_violin.png
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
})

# ============================================================
# 1. User settings
# ============================================================

motif_matrix_file <- "results/intermediate/endmotif/global_endmotif_matrix.tsv"
sample_summary_file <- "results/tables/endmotif/endmotif_sample_summary.tsv"

table_dir <- "results/tables/endmotif"
figure_dir <- "results/figures/endmotif"

dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

out_mds_file <- file.path(table_dir, "global_mds_scores.tsv")
out_summary_file <- file.path(table_dir, "global_mds_summary_by_group.tsv")
kruskal_out <- file.path(table_dir, "global_mds_kruskal_results.tsv")
pairwise_out <- file.path(table_dir, "global_mds_pairwise_wilcoxon_results.tsv")

plot_mds_out <- file.path(figure_dir, "global_mds_violin.png")
plot_mds_z_out <- file.path(figure_dir, "global_mds_z_violin.png")

group_levels <- c("Healthy", "Remission", "Relapse")

fill_cols <- c(
  "Healthy" = "darkseagreen3",
  "Remission" = "lightpink",
  "Relapse" = "lightblue"
)

point_cols <- c(
  "Healthy" = "darkgreen",
  "Remission" = "deeppink3",
  "Relapse" = "blue3"
)

# ============================================================
# 2. Helper functions
# ============================================================

calculate_shannon_entropy <- function(freq_vector) {
  freq_vector <- as.numeric(freq_vector)
  freq_vector <- freq_vector[!is.na(freq_vector)]
  freq_vector <- freq_vector[freq_vector > 0]

  -sum(freq_vector * log2(freq_vector))
}

p_to_label <- function(p) {
  if (is.na(p)) {
    return("NA")
  } else if (p < 0.001) {
    return("***")
  } else if (p < 0.01) {
    return("**")
  } else if (p < 0.05) {
    return("*")
  } else {
    return("ns")
  }
}

pairwise_to_table <- function(pairwise_result, variable_name) {
  as.data.frame(as.table(pairwise_result$p.value)) %>%
    rename(
      group2 = Var1,
      group1 = Var2,
      p_adj = Freq
    ) %>%
    filter(!is.na(p_adj)) %>%
    mutate(
      variable = variable_name,
      method = "Pairwise Wilcoxon rank-sum test",
      p_adjust_method = "BH",
      significance = vapply(p_adj, p_to_label, character(1))
    ) %>%
    select(variable, group1, group2, method, p_adjust_method, p_adj, significance)
}

add_pairwise_significance <- function(plot_object, df, y_col, pairwise_table) {
  y_max <- max(df[[y_col]], na.rm = TRUE)
  y_min <- min(df[[y_col]], na.rm = TRUE)
  y_range <- y_max - y_min

  if (y_range == 0) {
    y_range <- abs(y_max) * 0.1 + 1e-6
  }

  get_label <- function(g1, g2) {
    hit <- pairwise_table %>%
      filter(
        (group1 == g1 & group2 == g2) |
          (group1 == g2 & group2 == g1)
      )

    if (nrow(hit) == 0) {
      return("NA")
    }

    hit$significance[1]
  }

  lab_h_rem <- get_label("Healthy", "Remission")
  lab_h_rel <- get_label("Healthy", "Relapse")
  lab_rem_rel <- get_label("Remission", "Relapse")

  y1 <- y_max + 0.15 * y_range
  y2 <- y_max + 0.22 * y_range
  y3 <- y_max + 0.29 * y_range

  tick <- 0.025 * y_range
  text_offset <- 0.025 * y_range

  plot_object +
    annotate("segment", x = 1, xend = 2, y = y1, yend = y1, linewidth = 0.7) +
    annotate("segment", x = 1, xend = 1, y = y1, yend = y1 - tick, linewidth = 0.7) +
    annotate("segment", x = 2, xend = 2, y = y1, yend = y1 - tick, linewidth = 0.7) +
    annotate("text", x = 1.5, y = y1 + text_offset, label = lab_h_rem, size = 5) +

    annotate("segment", x = 2, xend = 3, y = y2, yend = y2, linewidth = 0.7) +
    annotate("segment", x = 2, xend = 2, y = y2, yend = y2 - tick, linewidth = 0.7) +
    annotate("segment", x = 3, xend = 3, y = y2, yend = y2 - tick, linewidth = 0.7) +
    annotate("text", x = 2.5, y = y2 + text_offset, label = lab_rem_rel, size = 5) +

    annotate("segment", x = 1, xend = 3, y = y3, yend = y3, linewidth = 0.7) +
    annotate("segment", x = 1, xend = 1, y = y3, yend = y3 - tick, linewidth = 0.7) +
    annotate("segment", x = 3, xend = 3, y = y3, yend = y3 - tick, linewidth = 0.7) +
    annotate("text", x = 2, y = y3 + text_offset, label = lab_h_rel, size = 5) +

    coord_cartesian(
      ylim = c(y_min - 0.10 * y_range, y_max + 0.35 * y_range),
      clip = "off"
    )
}

make_mds_violin_plot <- function(df, y_col, y_label, title_text, pairwise_table) {
  p <- ggplot(df, aes(x = group, y = .data[[y_col]], fill = group)) +
    geom_violin(
      trim = FALSE,
      color = "black",
      alpha = 0.65,
      width = 0.75,
      linewidth = 0.8
    ) +
    geom_boxplot(
      width = 0.22,
      outlier.shape = NA,
      fill = "white",
      color = "black",
      linewidth = 0.8
    ) +
    geom_jitter(
      aes(color = group),
      width = 0.08,
      size = 2,
      alpha = 0.85
    ) +
    scale_fill_manual(values = fill_cols) +
    scale_color_manual(values = point_cols) +
    labs(
      x = "Group",
      y = y_label,
      title = title_text
    ) +
    theme_classic(base_size = 16, base_family = "Times New Roman") +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5, size = 20),
      axis.text = element_text(color = "black", size = 13),
      axis.title = element_text(color = "black", size = 16),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.line = element_blank(),
      plot.margin = ggplot2::margin(10, 20, 10, 10)
    )

  add_pairwise_significance(
    plot_object = p,
    df = df,
    y_col = y_col,
    pairwise_table = pairwise_table
  )
}

# ============================================================
# 3. Check input files
# ============================================================

if (!file.exists(motif_matrix_file)) {
  stop("Motif matrix file not found: ", motif_matrix_file)
}

if (!file.exists(sample_summary_file)) {
  stop("Sample summary file not found: ", sample_summary_file)
}

# ============================================================
# 4. Load data
# ============================================================

motif_matrix_df <- read_tsv(motif_matrix_file, show_col_types = FALSE)

sample_summary <- read_tsv(sample_summary_file, show_col_types = FALSE) %>%
  mutate(
    group = factor(group, levels = group_levels)
  ) %>%
  filter(!is.na(group))

if (!"sample" %in% colnames(motif_matrix_df)) {
  stop("Motif matrix must contain a 'sample' column.")
}

required_summary_cols <- c("sample", "group")
missing_summary_cols <- setdiff(required_summary_cols, colnames(sample_summary))

if (length(missing_summary_cols) > 0) {
  stop(
    "Sample summary is missing required columns: ",
    paste(missing_summary_cols, collapse = ", ")
  )
}

motif_matrix_df <- motif_matrix_df %>%
  semi_join(sample_summary, by = "sample")

sample_summary <- sample_summary %>%
  semi_join(motif_matrix_df, by = "sample")

motif_only_df <- motif_matrix_df %>%
  select(-sample)

motif_matrix <- as.matrix(motif_only_df)
rownames(motif_matrix) <- motif_matrix_df$sample
storage.mode(motif_matrix) <- "numeric"

cat("\nLoaded motif matrix:\n")
cat("Samples:", nrow(motif_matrix), "\n")
cat("Motifs:", ncol(motif_matrix), "\n")

cat("\nSamples by group:\n")
print(sample_summary %>% count(group))

cat("\nRow sum summary:\n")
print(summary(rowSums(motif_matrix, na.rm = TRUE)))

# ============================================================
# 5. Calculate MDS
# ============================================================

mds_values <- apply(motif_matrix, 1, calculate_shannon_entropy)

mds_df <- tibble(
  sample = names(mds_values),
  MDS = as.numeric(mds_values)
) %>%
  left_join(
    sample_summary %>%
      select(sample, group, include_analysis, exclusion_reason),
    by = "sample"
  ) %>%
  mutate(
    group = factor(group, levels = group_levels)
  ) %>%
  arrange(group, sample)

# ============================================================
# 6. Calculate Healthy-reference MDS z-score
# ============================================================

healthy_mds <- mds_df %>%
  filter(group == "Healthy") %>%
  pull(MDS)

if (length(healthy_mds) < 2) {
  stop("Fewer than two Healthy samples found. Cannot calculate Healthy-reference z-score.")
}

healthy_mean <- mean(healthy_mds, na.rm = TRUE)
healthy_sd <- sd(healthy_mds, na.rm = TRUE)

if (is.na(healthy_sd) || healthy_sd == 0) {
  stop("Healthy MDS standard deviation is zero or NA. Cannot calculate z-score.")
}

mds_df <- mds_df %>%
  mutate(
    MDS_z = (MDS - healthy_mean) / healthy_sd
  )

# ============================================================
# 7. Summary by group
# ============================================================

mds_summary_by_group <- mds_df %>%
  group_by(group) %>%
  summarise(
    n = n(),
    MDS_mean = mean(MDS, na.rm = TRUE),
    MDS_sd = sd(MDS, na.rm = TRUE),
    MDS_median = median(MDS, na.rm = TRUE),
    MDS_min = min(MDS, na.rm = TRUE),
    MDS_max = max(MDS, na.rm = TRUE),
    MDS_z_mean = mean(MDS_z, na.rm = TRUE),
    MDS_z_sd = sd(MDS_z, na.rm = TRUE),
    MDS_z_median = median(MDS_z, na.rm = TRUE),
    MDS_z_min = min(MDS_z, na.rm = TRUE),
    MDS_z_max = max(MDS_z, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nMDS summary by group:\n")
print(mds_summary_by_group)

cat("\nHealthy reference values:\n")
cat("Healthy mean MDS:", healthy_mean, "\n")
cat("Healthy SD MDS:", healthy_sd, "\n")

# ============================================================
# 8. Kruskal-Wallis tests
# ============================================================

kruskal_mds <- kruskal.test(MDS ~ group, data = mds_df)
kruskal_mds_z <- kruskal.test(MDS_z ~ group, data = mds_df)

kruskal_results <- tibble(
  variable = c("MDS", "MDS_z"),
  method = "Kruskal-Wallis rank-sum test",
  statistic = c(
    as.numeric(kruskal_mds$statistic),
    as.numeric(kruskal_mds_z$statistic)
  ),
  df = c(
    as.numeric(kruskal_mds$parameter),
    as.numeric(kruskal_mds_z$parameter)
  ),
  p_value = c(
    kruskal_mds$p.value,
    kruskal_mds_z$p.value
  ),
  significance = vapply(
    c(kruskal_mds$p.value, kruskal_mds_z$p.value),
    p_to_label,
    character(1)
  )
)

cat("\nKruskal-Wallis results:\n")
print(kruskal_results)

# ============================================================
# 9. Pairwise Wilcoxon tests
# ============================================================

pairwise_mds <- pairwise.wilcox.test(
  x = mds_df$MDS,
  g = mds_df$group,
  p.adjust.method = "BH",
  exact = FALSE
)

pairwise_mds_z <- pairwise.wilcox.test(
  x = mds_df$MDS_z,
  g = mds_df$group,
  p.adjust.method = "BH",
  exact = FALSE
)

pairwise_results <- bind_rows(
  pairwise_to_table(pairwise_mds, "MDS"),
  pairwise_to_table(pairwise_mds_z, "MDS_z")
)

cat("\nPairwise Wilcoxon results:\n")
print(pairwise_results)

# ============================================================
# 10. Save tables
# ============================================================

write_tsv(mds_df, out_mds_file)
write_tsv(mds_summary_by_group, out_summary_file)
write_tsv(kruskal_results, kruskal_out)
write_tsv(pairwise_results, pairwise_out)

cat("\nSaved MDS scores to:\n")
cat(out_mds_file, "\n")

cat("\nSaved MDS summary to:\n")
cat(out_summary_file, "\n")

cat("\nSaved Kruskal-Wallis results to:\n")
cat(kruskal_out, "\n")

cat("\nSaved pairwise Wilcoxon results to:\n")
cat(pairwise_out, "\n")

# ============================================================
# 11. Plot MDS and MDS z-score
# ============================================================

pairwise_mds_table <- pairwise_results %>%
  filter(variable == "MDS")

pairwise_mds_z_table <- pairwise_results %>%
  filter(variable == "MDS_z")

p_mds <- make_mds_violin_plot(
  df = mds_df,
  y_col = "MDS",
  y_label = "Motif diversity score",
  title_text = "Motif Diversity Score",
  pairwise_table = pairwise_mds_table
)

p_mds_z <- make_mds_violin_plot(
  df = mds_df,
  y_col = "MDS_z",
  y_label = "MDS Z-score",
  title_text = "Motif Diversity Score Z-score",
  pairwise_table = pairwise_mds_z_table
)

ggsave(
  filename = plot_mds_out,
  plot = p_mds,
  width = 7,
  height = 5.5,
  dpi = 300
)

ggsave(
  filename = plot_mds_z_out,
  plot = p_mds_z,
  width = 7,
  height = 5.5,
  dpi = 300
)

cat("\nSaved MDS plot to:\n")
cat(plot_mds_out, "\n")

cat("\nSaved MDS z-score plot to:\n")
cat(plot_mds_z_out, "\n")

cat("\nDone.\n")
