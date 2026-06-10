# ============================================================
# Complete global end-motif analysis and plotting
#
# This script performs:
#   1. Ordered trend test: Healthy -> Remission -> Relapse
#   2. Top 10 increasing/decreasing trend motif plots
#   3. Kruskal-Wallis altered motif analysis
#   4. Heatmap of top 50 altered motifs
#   5. PCA of all motifs
#   6. PCA of significantly altered motifs
#   7. t-SNE and UMAP of all motifs
#   8. t-SNE and UMAP of top 50 altered motifs
#   9. Relapse vs Healthy weighted seqlogo
#   10. CpG-containing vs non-CpG motif burden
#
# Input:
#   results/intermediate/endmotif/global_endmotif_matrix.tsv
#   results/tables/endmotif/endmotif_sample_summary.tsv
#
# Output:
#   results/tables/endmotif/
#   results/figures/endmotif/
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(clinfun)
  library(pheatmap)
})

# Optional packages
has_Rtsne <- requireNamespace("Rtsne", quietly = TRUE)
has_uwot <- requireNamespace("uwot", quietly = TRUE)
has_ggseqlogo <- requireNamespace("ggseqlogo", quietly = TRUE)
has_patchwork <- requireNamespace("patchwork", quietly = TRUE)

# ============================================================
# 1. Settings
# ============================================================

motif_matrix_file <- "results/intermediate/endmotif/global_endmotif_matrix.tsv"
sample_summary_file <- "results/tables/endmotif/endmotif_sample_summary.tsv"

table_dir <- "results/tables/endmotif"
figure_dir <- "results/figures/endmotif"

dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

group_levels <- c("Healthy", "Remission", "Relapse")

group_colors <- c(
  "Healthy" = "darkgreen",
  "Remission" = "deeppink3",
  "Relapse" = "blue3"
)

fill_colors <- c(
  "Healthy" = "darkseagreen3",
  "Remission" = "lightpink",
  "Relapse" = "lightblue"
)

theme_endmotif <- function(base_size = 16) {
  theme_classic(base_size = base_size, base_family = "Times New Roman") +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, family = "Times New Roman"),
      axis.text = element_text(color = "black", family = "Times New Roman"),
      axis.title = element_text(color = "black", family = "Times New Roman"),
      legend.text = element_text(family = "Times New Roman"),
      legend.title = element_text(family = "Times New Roman"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.line = element_blank()
    )
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

# ============================================================
# 2. Load data
# ============================================================

if (!file.exists(motif_matrix_file)) {
  stop("Motif matrix file not found: ", motif_matrix_file)
}

if (!file.exists(sample_summary_file)) {
  stop("Sample summary file not found: ", sample_summary_file)
}

motif_matrix_df <- read_tsv(motif_matrix_file, show_col_types = FALSE)

sample_summary <- read_tsv(sample_summary_file, show_col_types = FALSE) %>%
  mutate(
    group = factor(group, levels = group_levels)
  ) %>%
  filter(!is.na(group))

motif_matrix_df <- motif_matrix_df %>%
  semi_join(sample_summary, by = "sample") %>%
  arrange(match(sample, sample_summary$sample))

motif_df <- motif_matrix_df %>%
  column_to_rownames("sample")

motif_df <- as.data.frame(lapply(motif_df, as.numeric))
rownames(motif_df) <- motif_matrix_df$sample

sample_summary <- sample_summary %>%
  filter(sample %in% rownames(motif_df)) %>%
  arrange(match(sample, rownames(motif_df)))

group <- sample_summary$group
names(group) <- sample_summary$sample

motif_matrix <- as.matrix(motif_df)
storage.mode(motif_matrix) <- "numeric"

cat("\nSamples used:\n")
print(sample_summary %>% count(group))

cat("\nMotif matrix dimensions:\n")
print(dim(motif_matrix))

cat("\nRow sum summary:\n")
print(summary(rowSums(motif_matrix)))

# ============================================================
# 3. Ordered trend test
# Healthy -> Remission -> Relapse
# ============================================================

safe_jt <- function(x, g, alternative) {
  out <- tryCatch(
    clinfun::jonckheere.test(x = x, g = g, alternative = alternative),
    error = function(e) NULL
  )

  if (is.null(out)) {
    return(NA_real_)
  } else {
    return(out$p.value)
  }
}

trend_results <- lapply(colnames(motif_df), function(m) {
  x <- motif_df[[m]]
  g_num <- as.integer(group)

  means <- tapply(x, group, mean, na.rm = TRUE)
  medians <- tapply(x, group, median, na.rm = TRUE)

  p_inc <- safe_jt(x, g_num, "increasing")
  p_dec <- safe_jt(x, g_num, "decreasing")

  data.frame(
    motif = m,
    Healthy = means["Healthy"],
    Remission = means["Remission"],
    Relapse = means["Relapse"],
    Healthy_median = medians["Healthy"],
    Remission_median = medians["Remission"],
    Relapse_median = medians["Relapse"],
    p_increasing = p_inc,
    p_decreasing = p_dec,
    stringsAsFactors = FALSE
  )
}) %>%
  bind_rows() %>%
  mutate(
    padj_increasing = p.adjust(p_increasing, method = "BH"),
    padj_decreasing = p.adjust(p_decreasing, method = "BH"),
    relapse_minus_healthy = Relapse - Healthy
  )

top10_increasing <- trend_results %>%
  filter(
    Healthy < Remission,
    Remission < Relapse,
    !is.na(padj_increasing)
  ) %>%
  arrange(padj_increasing) %>%
  slice_head(n = 10)

top10_decreasing <- trend_results %>%
  filter(
    Healthy > Remission,
    Remission > Relapse,
    !is.na(padj_decreasing)
  ) %>%
  arrange(padj_decreasing) %>%
  slice_head(n = 10)

top10_trend_table <- bind_rows(
  top10_increasing %>%
    transmute(
      motif,
      trend = "Increasing",
      Healthy_mean = Healthy,
      Remission_mean = Remission,
      Relapse_mean = Relapse,
      p_value = p_increasing,
      padj = padj_increasing
    ),
  top10_decreasing %>%
    transmute(
      motif,
      trend = "Decreasing",
      Healthy_mean = Healthy,
      Remission_mean = Remission,
      Relapse_mean = Relapse,
      p_value = p_decreasing,
      padj = padj_decreasing
    )
)

write_tsv(trend_results, file.path(table_dir, "endmotif_trend_test_results.tsv"))
write_tsv(top10_increasing, file.path(table_dir, "top10_increasing_trend_motifs.tsv"))
write_tsv(top10_decreasing, file.path(table_dir, "top10_decreasing_trend_motifs.tsv"))
write_tsv(top10_trend_table, file.path(table_dir, "top10_increasing_decreasing_trend_motifs_clean.tsv"))

cat("\nTop 10 increasing motifs:\n")
print(top10_increasing %>% select(motif, Healthy, Remission, Relapse, padj_increasing))

cat("\nTop 10 decreasing motifs:\n")
print(top10_decreasing %>% select(motif, Healthy, Remission, Relapse, padj_decreasing))

# ============================================================
# 4. Trend boxplot function: thesis-style grouped plot
# ============================================================

plot_trend_boxplot_checked <- function(motif_df,
                                       group,
                                       top_table,
                                       padj_col,
                                       title_text,
                                       direction = c("increasing", "decreasing")) {
  direction <- match.arg(direction)

  motifs <- top_table$motif

  if (length(motifs) == 0) {
    stop("No motifs available for plotting: ", title_text)
  }

  if (is.null(rownames(motif_df))) {
    stop("motif_df must have sample names as rownames.")
  }

  if (is.null(names(group))) {
    stop("group must be a named vector/factor.")
  }

  if (!all(rownames(motif_df) %in% names(group))) {
    missing_samples <- setdiff(rownames(motif_df), names(group))
    stop("Some samples are missing from names(group): ", paste(missing_samples, collapse = ", "))
  }

  if (!all(motifs %in% colnames(motif_df))) {
    missing_motifs <- setdiff(motifs, colnames(motif_df))
    stop("Some motifs are missing from motif_df: ", paste(missing_motifs, collapse = ", "))
  }

  plot_df <- motif_df %>%
    rownames_to_column("sample") %>%
    mutate(group = group[match(sample, names(group))]) %>%
    select(sample, group, all_of(motifs)) %>%
    pivot_longer(
      cols = all_of(motifs),
      names_to = "motif",
      values_to = "freq"
    ) %>%
    mutate(
      group = factor(group, levels = group_levels),
      motif = factor(motif, levels = motifs)
    )

  mean_check <- plot_df %>%
    group_by(motif, group) %>%
    summarise(plot_mean = mean(freq, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = group, values_from = plot_mean) %>%
    left_join(
      top_table %>%
        select(motif, Healthy, Remission, Relapse),
      by = "motif",
      suffix = c("_plot", "_table")
    )

  max_diff <- max(
    abs(mean_check$Healthy_plot - mean_check$Healthy_table),
    abs(mean_check$Remission_plot - mean_check$Remission_table),
    abs(mean_check$Relapse_plot - mean_check$Relapse_table),
    na.rm = TRUE
  )

  if (max_diff > 1e-10) {
    print(mean_check)
    stop("Plotted group means do not match trend table means.")
  }

  if (direction == "increasing") {
    direction_check <- mean_check %>%
      mutate(correct_direction = Healthy_plot < Remission_plot & Remission_plot < Relapse_plot)
  } else {
    direction_check <- mean_check %>%
      mutate(correct_direction = Healthy_plot > Remission_plot & Remission_plot > Relapse_plot)
  }

  if (!all(direction_check$correct_direction)) {
    print(direction_check)
    stop("Some motifs do not follow the expected trend direction.")
  }

  star_df <- top_table %>%
    mutate(
      motif = factor(motif, levels = motifs),
      label = p_to_star(.data[[padj_col]])
    ) %>%
    select(motif, label)

  global_top <- max(plot_df$freq, na.rm = TRUE) * 1.10
  global_bottom <- min(plot_df$freq, na.rm = TRUE)

  star_df <- star_df %>%
    mutate(y = global_top)

  p <- ggplot(plot_df, aes(x = motif, y = freq, color = group)) +
    geom_boxplot(
      aes(group = interaction(motif, group)),
      position = position_dodge(width = 0.72),
      width = 0.48,
      outlier.shape = NA,
      fill = "white",
      linewidth = 0.9
    ) +
    geom_jitter(
      aes(group = group),
      position = position_jitterdodge(jitter.width = 0.10, dodge.width = 0.72),
      size = 1.5,
      alpha = 0.7
    ) +
    geom_text(
      data = star_df,
      aes(x = motif, y = y, label = label),
      inherit.aes = FALSE,
      size = 5,
      family = "Times New Roman"
    ) +
    scale_color_manual(values = group_colors) +
    labs(
      x = NULL,
      y = "Fragment 5' end motif frequency",
      title = title_text,
      color = "*** p < 0.001\n** p < 0.01\n* p < 0.05"
    ) +
    coord_cartesian(
      ylim = c(global_bottom, global_top * 1.03),
      clip = "off"
    ) +
    theme_endmotif(base_size = 16) +
    theme(
      legend.position = c(1, 1.02),
      legend.justification = c(0, 1.2),
      legend.background = element_blank(),
      plot.margin = ggplot2::margin(10, 95, 10, 10)
    )

  return(p)
}

if (nrow(top10_increasing) > 0) {
  p_increasing <- plot_trend_boxplot_checked(
    motif_df = motif_df,
    group = group,
    top_table = top10_increasing,
    padj_col = "padj_increasing",
    title_text = "Top 10 increasing trend motifs",
    direction = "increasing"
  )

  ggsave(
    file.path(figure_dir, "top10_increasing_trend_motifs_grouped.png"),
    p_increasing,
    width = 10,
    height = 5.5,
    dpi = 300
  )
}

if (nrow(top10_decreasing) > 0) {
  p_decreasing <- plot_trend_boxplot_checked(
    motif_df = motif_df,
    group = group,
    top_table = top10_decreasing,
    padj_col = "padj_decreasing",
    title_text = "Top 10 decreasing trend motifs",
    direction = "decreasing"
  )

  ggsave(
    file.path(figure_dir, "top10_decreasing_trend_motifs_grouped.png"),
    p_decreasing,
    width = 10,
    height = 5.5,
    dpi = 300
  )
}

# ============================================================
# 5. Kruskal-Wallis altered motif analysis
# ============================================================

altered_results <- lapply(colnames(motif_df), function(m) {
  x <- motif_df[[m]]
  kw <- kruskal.test(x ~ group)

  group_means <- tapply(x, group, mean, na.rm = TRUE)
  group_medians <- tapply(x, group, median, na.rm = TRUE)

  data.frame(
    motif = m,
    Healthy_mean = group_means["Healthy"],
    Remission_mean = group_means["Remission"],
    Relapse_mean = group_means["Relapse"],
    Healthy_median = group_medians["Healthy"],
    Remission_median = group_medians["Remission"],
    Relapse_median = group_medians["Relapse"],
    p_overall = kw$p.value,
    stringsAsFactors = FALSE
  )
}) %>%
  bind_rows() %>%
  mutate(
    padj_overall = p.adjust(p_overall, method = "BH"),
    highest_group = apply(
      select(., Healthy_mean, Remission_mean, Relapse_mean),
      1,
      function(x) group_levels[which.max(x)]
    ),
    mean_range = apply(
      select(., Healthy_mean, Remission_mean, Relapse_mean),
      1,
      function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
    )
  ) %>%
  arrange(padj_overall)

write_tsv(altered_results, file.path(table_dir, "altered_motif_kruskal_results.tsv"))

top50_motifs <- altered_results %>%
  filter(!is.na(padj_overall)) %>%
  arrange(padj_overall) %>%
  slice_head(n = 50) %>%
  pull(motif)

sig_motifs <- altered_results %>%
  filter(padj_overall < 0.05) %>%
  pull(motif)

cat("\nNumber of significantly altered motifs, BH < 0.05:", length(sig_motifs), "\n")

# ============================================================
# 6. Heatmap of top 50 altered motifs
# ============================================================

if (length(top50_motifs) > 1) {
  heat_mat <- as.matrix(motif_df[, top50_motifs, drop = FALSE])

  annotation_row <- data.frame(Group = group)
  rownames(annotation_row) <- rownames(heat_mat)

  ann_colors <- list(
    Group = group_colors
  )

  png(
    filename = file.path(figure_dir, "heatmap_top50_altered_endmotifs.png"),
    width = 2200,
    height = 1700,
    res = 220
  )

  pheatmap(
    heat_mat,
    scale = "column",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_row = annotation_row,
    annotation_colors = ann_colors,
    show_rownames = FALSE,
    show_colnames = TRUE,
    fontsize_row = 7,
    fontsize_col = 8,
    fontfamily = "Times New Roman",
    color = colorRampPalette(c("blue", "white", "red"))(100),
    border_color = NA,
    main = "Heatmap of top 50 altered 5' end motifs"
  )

  dev.off()
}

# ============================================================
# 7. PCA helper function
# ============================================================

make_pca_plot <- function(X, group, title_text, output_file) {
  pca_res <- prcomp(X, center = TRUE, scale. = TRUE)
  percent_var <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 2)

  pca_df <- data.frame(
    sample = rownames(X),
    sample_label = sub("^Cap", "", rownames(X)),
    group = group,
    PC1 = pca_res$x[, 1],
    PC2 = pca_res$x[, 2]
  )

  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
    geom_point(size = 3, alpha = 0.9) +
    geom_text_repel(
      aes(label = sample_label),
      size = 3,
      max.overlaps = Inf,
      show.legend = FALSE,
      family = "Times New Roman"
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
      title = title_text
    ) +
    theme_endmotif(base_size = 16) +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black")
    )

  ggsave(output_file, p, width = 8, height = 6, dpi = 300)

  return(p)
}

pca_all_plot <- make_pca_plot(
  X = motif_matrix,
  group = group,
  title_text = "PCA of 5' end-motif frequencies",
  output_file = file.path(figure_dir, "pca_all_endmotifs.png")
)

if (length(sig_motifs) >= 2) {
  motif_matrix_sig <- motif_matrix[, sig_motifs, drop = FALSE]

  pca_sig_plot <- make_pca_plot(
    X = motif_matrix_sig,
    group = group,
    title_text = "PCA of significantly altered 5' end motifs",
    output_file = file.path(figure_dir, "pca_significantly_altered_endmotifs.png")
  )
} else {
  cat("\nSkipping PCA of significant motifs because fewer than 2 significant motifs were found.\n")
}

# ============================================================
# 8. t-SNE and UMAP helper functions
# ============================================================

make_tsne_plot <- function(X, group, title_text, output_file) {
  if (!has_Rtsne) {
    cat("\nPackage Rtsne is not installed. Skipping:", title_text, "\n")
    return(NULL)
  }

  X_scaled <- scale(X)
  tsne_perplexity <- min(15, floor((nrow(X_scaled) - 1) / 3))

  set.seed(111)
  tsne_res <- Rtsne::Rtsne(
    X_scaled,
    dims = 2,
    perplexity = tsne_perplexity,
    pca = TRUE,
    max_iter = 1000,
    check_duplicates = FALSE,
    verbose = FALSE
  )

  tsne_df <- data.frame(
    sample = rownames(X),
    sample_label = sub("^Cap", "", rownames(X)),
    group = group,
    tSNE1 = tsne_res$Y[, 1],
    tSNE2 = tsne_res$Y[, 2]
  )

  p <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = group)) +
    geom_point(size = 3, alpha = 0.9) +
    geom_text_repel(
      aes(label = sample_label),
      size = 3,
      max.overlaps = Inf,
      show.legend = FALSE,
      family = "Times New Roman"
    ) +
    stat_ellipse(
      aes(group = group),
      linewidth = 0.8,
      linetype = "dashed",
      show.legend = FALSE
    ) +
    scale_color_manual(values = group_colors) +
    labs(
      x = "t-SNE 1",
      y = "t-SNE 2",
      title = title_text
    ) +
    theme_endmotif(base_size = 16) +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black")
    )

  ggsave(output_file, p, width = 8, height = 6, dpi = 300)

  return(p)
}

make_umap_plot <- function(X, group, title_text, output_file) {
  if (!has_uwot) {
    cat("\nPackage uwot is not installed. Skipping:", title_text, "\n")
    return(NULL)
  }

  X_scaled <- scale(X)

  set.seed(111)
  umap_res <- uwot::umap(
    X_scaled,
    n_neighbors = 15,
    min_dist = 0.30,
    metric = "euclidean",
    n_components = 2,
    verbose = FALSE
  )

  umap_df <- data.frame(
    sample = rownames(X),
    sample_label = sub("^Cap", "", rownames(X)),
    group = group,
    UMAP1 = umap_res[, 1],
    UMAP2 = umap_res[, 2]
  )

  p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = group)) +
    geom_point(size = 3, alpha = 0.9) +
    geom_text_repel(
      aes(label = sample_label),
      size = 3,
      max.overlaps = Inf,
      show.legend = FALSE,
      family = "Times New Roman"
    ) +
    stat_ellipse(
      aes(group = group),
      linewidth = 0.8,
      linetype = "dashed",
      show.legend = FALSE
    ) +
    scale_color_manual(values = group_colors) +
    labs(
      x = "UMAP 1",
      y = "UMAP 2",
      title = title_text
    ) +
    theme_endmotif(base_size = 16) +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black")
    )

  ggsave(output_file, p, width = 8, height = 6, dpi = 300)

  return(p)
}

tsne_all_plot <- make_tsne_plot(
  X = motif_matrix,
  group = group,
  title_text = "t-SNE of 5' end-motif frequencies",
  output_file = file.path(figure_dir, "tsne_all_endmotifs.png")
)

umap_all_plot <- make_umap_plot(
  X = motif_matrix,
  group = group,
  title_text = "UMAP of 5' end-motif frequencies",
  output_file = file.path(figure_dir, "umap_all_endmotifs.png")
)

if (length(top50_motifs) >= 2) {
  X_top50 <- motif_matrix[, top50_motifs, drop = FALSE]

  tsne_top50_plot <- make_tsne_plot(
    X = X_top50,
    group = group,
    title_text = "t-SNE of top 50 altered 5' end motifs",
    output_file = file.path(figure_dir, "tsne_top50_altered_endmotifs.png")
  )

  umap_top50_plot <- make_umap_plot(
    X = X_top50,
    group = group,
    title_text = "UMAP of top 50 altered 5' end motifs",
    output_file = file.path(figure_dir, "umap_top50_altered_endmotifs.png")
  )
}

# ============================================================
# 9. Weighted seqlogo: Relapse vs Healthy
# ============================================================

if (has_ggseqlogo && has_patchwork) {
  library(ggseqlogo)
  library(patchwork)

  padj_cutoff <- 0.05
  log2fc_cutoff <- log2(1.25)
  top_n <- 30
  eps <- 1e-8

  relapse_healthy_df <- lapply(colnames(motif_df), function(m) {
    x_h <- motif_df[group == "Healthy", m]
    x_r <- motif_df[group == "Relapse", m]

    wt <- wilcox.test(x_r, x_h, exact = FALSE)

    mean_h <- mean(x_h, na.rm = TRUE)
    mean_r <- mean(x_r, na.rm = TRUE)

    data.frame(
      motif = m,
      Healthy_mean = mean_h,
      Relapse_mean = mean_r,
      log2FC = log2((mean_r + eps) / (mean_h + eps)),
      mean_difference = mean_r - mean_h,
      p_value = wt$p.value,
      stringsAsFactors = FALSE
    )
  }) %>%
    bind_rows() %>%
    mutate(
      padj = p.adjust(p_value, method = "BH"),
      neg_log10_padj = -log10(padj),
      direction = case_when(
        padj < padj_cutoff & log2FC >= log2fc_cutoff ~ "Higher in Relapse",
        padj < padj_cutoff & log2FC <= -log2fc_cutoff ~ "Higher in Healthy",
        TRUE ~ "Other"
      )
    )

  write_tsv(relapse_healthy_df, file.path(table_dir, "relapse_vs_healthy_endmotif_wilcoxon.tsv"))

  top_relapse <- relapse_healthy_df %>%
    filter(direction == "Higher in Relapse") %>%
    arrange(padj, desc(abs(log2FC))) %>%
    slice_head(n = top_n)

  top_healthy <- relapse_healthy_df %>%
    filter(direction == "Higher in Healthy") %>%
    arrange(padj, desc(abs(log2FC))) %>%
    slice_head(n = top_n)

  make_weighted_pwm <- function(motifs, weights) {
    motifs <- toupper(motifs)
    weights <- as.numeric(weights)
    weights[is.na(weights)] <- 0

    bases <- c("A", "C", "G", "T")
    k <- nchar(motifs[1])

    mat <- matrix(0, nrow = 4, ncol = k)
    rownames(mat) <- bases
    colnames(mat) <- seq_len(k)

    for (i in seq_along(motifs)) {
      chars <- strsplit(motifs[i], "")[[1]]
      for (pos in seq_len(k)) {
        mat[chars[pos], pos] <- mat[chars[pos], pos] + weights[i]
      }
    }

    mat <- sweep(mat, 2, colSums(mat), "/")
    mat[is.na(mat)] <- 0

    return(mat)
  }

  if (nrow(top_relapse) > 0 && nrow(top_healthy) > 0) {
    top_relapse <- top_relapse %>%
      mutate(logo_weight = abs(log2FC) * (-log10(padj)))

    top_healthy <- top_healthy %>%
      mutate(logo_weight = abs(log2FC) * (-log10(padj)))

    pwm_relapse <- make_weighted_pwm(top_relapse$motif, top_relapse$logo_weight)
    pwm_healthy <- make_weighted_pwm(top_healthy$motif, top_healthy$logo_weight)

    p_logo_healthy <- ggseqlogo::ggseqlogo(pwm_healthy, method = "custom") +
      labs(
        title = paste0("Top ", nrow(top_healthy), " motifs enriched in Healthy"),
        subtitle = "Relapse vs Healthy; ranked by adjusted p-value and fold-change",
        x = "Motif position",
        y = "Weighted base probability"
      ) +
      theme_endmotif(base_size = 16) +
      theme(panel.border = element_blank())

    p_logo_relapse <- ggseqlogo::ggseqlogo(pwm_relapse, method = "custom") +
      labs(
        title = paste0("Top ", nrow(top_relapse), " motifs enriched in Relapse"),
        subtitle = "Relapse vs Healthy; ranked by adjusted p-value and fold-change",
        x = "Motif position",
        y = "Weighted base probability"
      ) +
      theme_endmotif(base_size = 16) +
      theme(panel.border = element_blank())

    p_seqlogo_top <- p_logo_healthy / p_logo_relapse

    ggsave(
      file.path(figure_dir, "seqlogo_weighted_healthy_vs_relapse.png"),
      p_seqlogo_top,
      width = 8,
      height = 7,
      dpi = 300
    )
  } else {
    cat("\nSkipping seqlogo because enriched motifs were not found in both directions.\n")
  }

} else {
  cat("\nPackages ggseqlogo and/or patchwork are not installed. Skipping seqlogo.\n")
}

# ============================================================
# 10. CpG-containing vs non-CpG motif burden
# ============================================================

motif_category_df <- data.frame(
  motif = colnames(motif_df),
  motif_class = ifelse(
    grepl("CG", colnames(motif_df)),
    "CpG-containing motifs",
    "non-CpG motifs"
  ),
  stringsAsFactors = FALSE
)

write_tsv(motif_category_df, file.path(table_dir, "motif_cpg_category.tsv"))

cpg_cols <- motif_category_df$motif[motif_category_df$motif_class == "CpG-containing motifs"]
non_cpg_cols <- motif_category_df$motif[motif_category_df$motif_class == "non-CpG motifs"]

cpg_burden_df <- data.frame(
  sample = rownames(motif_df),
  group = group,
  CpG_containing_motifs = rowSums(motif_df[, cpg_cols, drop = FALSE], na.rm = TRUE),
  non_CpG_motifs = rowSums(motif_df[, non_cpg_cols, drop = FALSE], na.rm = TRUE)
)

cpg_burden_long <- cpg_burden_df %>%
  pivot_longer(
    cols = c(CpG_containing_motifs, non_CpG_motifs),
    names_to = "motif_class",
    values_to = "burden"
  ) %>%
  mutate(
    group = factor(group, levels = group_levels),
    motif_class = recode(
      motif_class,
      CpG_containing_motifs = "CpG-containing motifs",
      non_CpG_motifs = "non-CpG motifs"
    )
  )

cpg_summary <- cpg_burden_long %>%
  group_by(group, motif_class) %>%
  summarise(
    n = n(),
    mean = mean(burden, na.rm = TRUE),
    sd = sd(burden, na.rm = TRUE),
    median = median(burden, na.rm = TRUE),
    min = min(burden, na.rm = TRUE),
    max = max(burden, na.rm = TRUE),
    .groups = "drop"
  )

cpg_tests <- cpg_burden_long %>%
  filter(group %in% c("Healthy", "Relapse")) %>%
  group_by(motif_class) %>%
  summarise(
    p_value = wilcox.test(burden ~ group, exact = FALSE)$p.value,
    Healthy_median = median(burden[group == "Healthy"], na.rm = TRUE),
    Relapse_median = median(burden[group == "Relapse"], na.rm = TRUE),
    median_difference = Relapse_median - Healthy_median,
    .groups = "drop"
  ) %>%
  mutate(
    padj = p.adjust(p_value, method = "BH"),
    label = p_to_star(padj)
  )

write_tsv(cpg_burden_df, file.path(table_dir, "cpg_non_cpg_burden_by_sample.tsv"))
write_tsv(cpg_summary, file.path(table_dir, "cpg_non_cpg_burden_summary.tsv"))
write_tsv(cpg_tests, file.path(table_dir, "cpg_non_cpg_tests_healthy_vs_relapse.tsv"))

sig_df <- cpg_burden_long %>%
  group_by(motif_class) %>%
  summarise(
    y_max = max(burden, na.rm = TRUE),
    y_min = min(burden, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    y_range = ifelse(y_max - y_min == 0, abs(y_max) * 0.1 + 1e-6, y_max - y_min),
    y_sig = y_max + 0.28 * y_range,
    y_text = y_max + 0.34 * y_range,
    tick = 0.04 * y_range
  ) %>%
  left_join(cpg_tests %>% select(motif_class, label, padj), by = "motif_class")

p_cpg_burden <- ggplot(
  cpg_burden_long,
  aes(x = group, y = burden, fill = group)
) +
  geom_violin(
    trim = FALSE,
    alpha = 0.65,
    color = "black",
    linewidth = 0.7
  ) +
  geom_boxplot(
    width = 0.18,
    outlier.shape = NA,
    fill = "white",
    color = "black",
    linewidth = 0.7
  ) +
  geom_jitter(
    aes(color = group),
    width = 0.08,
    size = 2,
    alpha = 0.85
  ) +
  geom_segment(
    data = sig_df,
    aes(x = 1, xend = 2, y = y_sig, yend = y_sig),
    inherit.aes = FALSE,
    linewidth = 0.7
  ) +
  geom_segment(
    data = sig_df,
    aes(x = 1, xend = 1, y = y_sig, yend = y_sig - tick),
    inherit.aes = FALSE,
    linewidth = 0.7
  ) +
  geom_segment(
    data = sig_df,
    aes(x = 2, xend = 2, y = y_sig, yend = y_sig - tick),
    inherit.aes = FALSE,
    linewidth = 0.7
  ) +
  geom_text(
    data = sig_df,
    aes(x = 1.5, y = y_text, label = label),
    inherit.aes = FALSE,
    size = 5,
    family = "Times New Roman"
  ) +
  facet_wrap(~ motif_class, scales = "free_y") +
  scale_fill_manual(values = fill_colors) +
  scale_color_manual(values = group_colors) +
  labs(
    x = "Group",
    y = "Sum of end-motif frequencies",
    title = "CpG-containing and non-CpG end-motif burden"
  ) +
  theme_endmotif(base_size = 16) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 14, family = "Times New Roman")
  )

ggsave(
  file.path(figure_dir, "cpg_non_cpg_burden_plot.png"),
  p_cpg_burden,
  width = 8,
  height = 5,
  dpi = 300
)

cat("\nAll complete global end-motif analysis outputs saved to:\n")
cat(table_dir, "\n")
cat(figure_dir, "\n")

cat("\nDone.\n")
