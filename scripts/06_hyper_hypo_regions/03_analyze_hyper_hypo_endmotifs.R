# ============================================================
# Hyper vs hypo methylated-region end-motif analysis
#
# Purpose:
#   Compare cfDNA end-motif profiles between hypermethylated
#   and hypomethylated regions within QC-passed Healthy samples.
#
# Design:
#   Healthy samples only.
#   Paired comparison:
#     each sample has one hyper file and one hypo file.
#
# Expected input files:
#   Cap01_hyper_endmotif.tsv
#   Cap01_hypo_endmotif.tsv
#   ...
#
# Expected FinaleToolkit columns:
#   contig/start/stop/name/count + 256 motif columns
#   or seq/start/stop/name/count + 256 motif columns
#
# Usage:
#   Rscript scripts/06_hyper_hypo_regions/03_analyze_hyper_hypo_endmotifs.R /path/to/interval_end_motifs
#
# Default input:
#   results/intermediate/hyper_hypo_regions/interval_end_motifs
#
# Outputs:
#   results/tables/hyper_hypo_regions/
#   results/figures/hyper_hypo_regions/
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(pheatmap)
  library(tibble)
})

graphics.off()

# ============================================================
# 1. Paths and settings
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

motif_dir <- ifelse(
  length(args) >= 1,
  args[1],
  "results/intermediate/hyper_hypo_regions/interval_end_motifs"
)

qc_metadata_file <- "results/tables/qc/final_sample_metadata.csv"

table_dir <- "results/tables/hyper_hypo_regions"
figure_dir <- "results/figures/hyper_hypo_regions"
sanity_dir <- file.path(table_dir, "sanity_checks")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sanity_dir, recursive = TRUE, showWarnings = FALSE)

expected_regions_per_class <- 52
expected_motif_files_after_qc <- 46

padj_cutoff <- 0.05
log2fc_cutoff <- 0.25
eps <- 1e-10

plot_family <- "serif"

class_levels <- c("hyper", "hypo")

class_colors <- c(
  hyper = "lightskyblue",
  hypo = "violetred3"
)

class_light_colors <- c(
  hyper = "lightblue1",
  hypo = "lavenderblush1"
)

volcano_colors <- c(
  "Hypo-enriched" = "steelblue4",
  "Hyper-enriched" = "maroon",
  "Significant, small effect" = "thistle",
  "Not significant" = "grey55"
)

theme_set(
  theme_classic(base_size = 15, base_family = plot_family) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, family = plot_family),
      axis.title = element_text(color = "black", family = plot_family),
      axis.text = element_text(color = "black", family = plot_family),
      legend.title = element_blank(),
      legend.text = element_text(family = plot_family)
    )
)

# ============================================================
# 2. Helper functions
# ============================================================

p_to_label <- function(p) {
  case_when(
    is.na(p) ~ "NA",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

calculate_mds <- function(x) {
  x <- as.numeric(x)
  x <- x[!is.na(x)]
  x <- x[x > 0]
  
  if (length(x) == 0) {
    return(NA_real_)
  }
  
  x <- x / sum(x)
  
  -sum(x * log2(x)) / log2(256)
}

safe_weighted_mean <- function(x, w) {
  x <- as.numeric(x)
  w <- as.numeric(w)
  
  keep <- !is.na(x) & !is.na(w) & w > 0
  
  if (sum(keep) == 0) {
    return(NA_real_)
  }
  
  weighted.mean(x[keep], w = w[keep])
}

safe_paired_wilcox <- function(x, y) {
  tryCatch(
    wilcox.test(x, y, paired = TRUE, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
}

safe_paired_ttest <- function(x, y) {
  tryCatch(
    t.test(x, y, paired = TRUE)$p.value,
    error = function(e) NA_real_
  )
}

make_pwm_from_motif_profile <- function(motif_freq) {
  motif_names <- names(motif_freq)
  
  motif_freq <- as.numeric(motif_freq)
  names(motif_freq) <- motif_names
  
  motif_freq[is.na(motif_freq)] <- 0
  
  valid_motifs <- !is.na(names(motif_freq)) &
    str_detect(names(motif_freq), "^[ACGT]{4}$")
  
  motif_freq <- motif_freq[valid_motifs]
  
  if (length(motif_freq) == 0) {
    stop("No valid 4-mer motifs found in motif profile.")
  }
  
  if (sum(motif_freq, na.rm = TRUE) <= 0) {
    stop("Motif profile has zero total frequency.")
  }
  
  motif_freq <- motif_freq / sum(motif_freq, na.rm = TRUE)
  
  bases <- c("A", "C", "G", "T")
  
  pwm <- matrix(
    0,
    nrow = 4,
    ncol = 4,
    dimnames = list(bases, as.character(1:4))
  )
  
  for (motif in names(motif_freq)) {
    chars <- strsplit(motif, "")[[1]]
    
    if (length(chars) != 4) {
      next
    }
    
    for (pos in 1:4) {
      if (chars[pos] %in% bases) {
        pwm[chars[pos], pos] <- pwm[chars[pos], pos] + motif_freq[[motif]]
      }
    }
  }
  
  col_sums <- colSums(pwm)
  col_sums[col_sums == 0] <- 1
  
  pwm <- sweep(pwm, 2, col_sums, "/")
  
  if (all(pwm == 0)) {
    stop("PWM is all zero after motif parsing.")
  }
  
  pwm
}

# ============================================================
# 3. Read QC metadata
# ============================================================

qc_healthy_samples <- NULL

if (file.exists(qc_metadata_file)) {
  qc_metadata <- read_csv(qc_metadata_file, show_col_types = FALSE) %>%
    mutate(sample = as.character(sample))
  
  if (all(c("sample", "group", "include_analysis") %in% colnames(qc_metadata))) {
    qc_healthy_samples <- qc_metadata %>%
      filter(group == "Healthy", include_analysis == "yes") %>%
      pull(sample)
    
    cat("\nQC metadata found.\n")
    cat("QC-passed Healthy samples:", length(qc_healthy_samples), "\n")
    print(qc_healthy_samples)
  } else {
    warning("QC metadata exists but does not contain sample, group, include_analysis columns.")
  }
} else {
  warning("QC metadata not found. All hyper/hypo endmotif files will be used.")
}

# ============================================================
# 4. Read motif files
# ============================================================

if (!dir.exists(motif_dir)) {
  stop("Motif directory does not exist: ", motif_dir)
}

motif_files <- list.files(
  motif_dir,
  pattern = "_endmotif\\.tsv$",
  full.names = TRUE
)

cat("\nNumber of motif files found before QC filtering:", length(motif_files), "\n")

if (length(motif_files) == 0) {
  stop("No *_endmotif.tsv files found in: ", motif_dir)
}

file_info <- tibble(
  file = motif_files,
  file_name = basename(motif_files),
  sample = str_extract(file_name, "Cap\\d+"),
  methylation_class = case_when(
    str_detect(file_name, "_hyper_") ~ "hyper",
    str_detect(file_name, "_hypo_") ~ "hypo",
    TRUE ~ NA_character_
  )
)

bad_file_info <- file_info %>%
  filter(is.na(sample) | is.na(methylation_class))

if (nrow(bad_file_info) > 0) {
  print(bad_file_info)
  stop("Some files do not have clear sample/class names.")
}

if (!is.null(qc_healthy_samples)) {
  file_info <- file_info %>%
    filter(sample %in% qc_healthy_samples)
}

cat("\nNumber of motif files after QC Healthy filtering:", nrow(file_info), "\n")
print(file_info %>% count(methylation_class))

if (nrow(file_info) != expected_motif_files_after_qc) {
  warning(
    "Expected ", expected_motif_files_after_qc,
    " files for 23 QC-passed Healthy samples x 2 classes, but found ",
    nrow(file_info), ". This may be okay if QC sample count changed."
  )
}

sample_file_check <- file_info %>%
  count(sample, methylation_class) %>%
  pivot_wider(
    names_from = methylation_class,
    values_from = n,
    values_fill = 0
  ) %>%
  arrange(sample)

print(sample_file_check, n = 100)

bad_sample_files <- sample_file_check %>%
  filter(hyper != 1 | hypo != 1)

if (nrow(bad_sample_files) > 0) {
  print(bad_sample_files)
  stop("Some samples do not have exactly one hyper and one hypo file.")
}

write_tsv(
  file_info,
  file.path(sanity_dir, "hyper_hypo_endmotif_file_info.tsv")
)

write_tsv(
  sample_file_check,
  file.path(sanity_dir, "hyper_hypo_sample_file_check.tsv")
)

# ============================================================
# 5. Load motif tables
# ============================================================

motif_df <- file_info$file %>%
  map_dfr(function(f) {
    file_name <- basename(f)
    
    sample <- str_extract(file_name, "Cap\\d+")
    
    methylation_class <- case_when(
      str_detect(file_name, "_hyper_") ~ "hyper",
      str_detect(file_name, "_hypo_") ~ "hypo",
      TRUE ~ NA_character_
    )
    
    tmp <- read_tsv(f, show_col_types = FALSE, progress = FALSE)
    
    if (!"contig" %in% colnames(tmp) && "seq" %in% colnames(tmp)) {
      tmp <- tmp %>%
        rename(contig = seq)
    }
    
    if (!"stop" %in% colnames(tmp) && "end" %in% colnames(tmp)) {
      tmp <- tmp %>%
        rename(stop = end)
    }
    
    required_cols <- c("contig", "start", "stop", "name", "count")
    missing_cols <- setdiff(required_cols, colnames(tmp))
    
    if (length(missing_cols) > 0) {
      stop(
        "Missing required columns in ",
        file_name,
        ": ",
        paste(missing_cols, collapse = ", ")
      )
    }
    
    tmp %>%
      mutate(
        file_name = file_name,
        sample = sample,
        methylation_class = methylation_class,
        beta = as.numeric(name),
        region_id = paste(contig, start, stop, sep = ":")
      )
  })

motif_df <- motif_df %>%
  mutate(
    methylation_class = factor(methylation_class, levels = class_levels),
    count = as.numeric(count)
  )

# ============================================================
# 6. Sanity checks
# ============================================================

non_motif_cols <- c(
  "contig", "start", "stop", "name", "count",
  "file_name", "sample", "methylation_class",
  "beta", "region_id"
)

motif_cols <- setdiff(colnames(motif_df), non_motif_cols)
motif_cols <- motif_cols[str_detect(motif_cols, "^[ACGT]{4}$")]

cat("\nNumber of motif columns:", length(motif_cols), "\n")

if (length(motif_cols) != 256) {
  warning("Expected 256 4-mer motif columns, found ", length(motif_cols))
}

bad_motif_names <- motif_cols[!str_detect(motif_cols, "^[ACGT]{4}$")]

if (length(bad_motif_names) > 0) {
  print(bad_motif_names)
  stop("Some motif columns are not valid 4-mers.")
}

motif_df <- motif_df %>%
  mutate(
    across(all_of(motif_cols), as.numeric)
  )

basic_summary <- motif_df %>%
  summarise(
    n_rows = n(),
    n_samples = n_distinct(sample),
    n_unique_regions = n_distinct(region_id),
    min_beta = min(beta, na.rm = TRUE),
    max_beta = max(beta, na.rm = TRUE),
    min_count = min(count, na.rm = TRUE),
    max_count = max(count, na.rm = TRUE),
    missing_beta = sum(is.na(beta)),
    missing_count = sum(is.na(count))
  )

cat("\nBasic summary:\n")
print(basic_summary)

row_balance <- motif_df %>%
  count(sample, methylation_class) %>%
  arrange(sample, methylation_class)

cat("\nRows per sample/class:\n")
print(row_balance, n = 100)

bad_row_balance <- row_balance %>%
  filter(n != expected_regions_per_class)

if (nrow(bad_row_balance) > 0) {
  print(bad_row_balance)
  warning("Some sample/class groups do not have exactly ", expected_regions_per_class, " regions.")
}

motif_df <- motif_df %>%
  rowwise() %>%
  mutate(
    motif_sum = sum(c_across(all_of(motif_cols)), na.rm = TRUE)
  ) %>%
  ungroup()

motif_sum_summary <- motif_df %>%
  summarise(
    min_motif_sum = min(motif_sum, na.rm = TRUE),
    q1_motif_sum = quantile(motif_sum, 0.25, na.rm = TRUE),
    median_motif_sum = median(motif_sum, na.rm = TRUE),
    mean_motif_sum = mean(motif_sum, na.rm = TRUE),
    q3_motif_sum = quantile(motif_sum, 0.75, na.rm = TRUE),
    max_motif_sum = max(motif_sum, na.rm = TRUE)
  )

beta_summary <- motif_df %>%
  group_by(methylation_class) %>%
  summarise(
    n = n(),
    min_beta = min(beta, na.rm = TRUE),
    median_beta = median(beta, na.rm = TRUE),
    mean_beta = mean(beta, na.rm = TRUE),
    max_beta = max(beta, na.rm = TRUE),
    .groups = "drop"
  )

count_summary <- motif_df %>%
  group_by(sample, methylation_class) %>%
  summarise(
    total_count = sum(count, na.rm = TRUE),
    mean_count = mean(count, na.rm = TRUE),
    median_count = median(count, na.rm = TRUE),
    .groups = "drop"
  )

write_tsv(basic_summary, file.path(sanity_dir, "hyper_hypo_basic_summary.tsv"))
write_tsv(row_balance, file.path(sanity_dir, "hyper_hypo_row_balance.tsv"))
write_tsv(motif_sum_summary, file.path(sanity_dir, "hyper_hypo_motif_sum_summary.tsv"))
write_tsv(beta_summary, file.path(sanity_dir, "hyper_hypo_beta_summary.tsv"))
write_tsv(count_summary, file.path(sanity_dir, "hyper_hypo_count_summary_by_sample.tsv"))

# ============================================================
# 7. Sanity plots
# ============================================================

p_motif_sum <- ggplot(motif_df, aes(x = motif_sum)) +
  geom_histogram(bins = 40, fill = "blue3", color = "white") +
  labs(
    x = "Sum of 256 motif frequencies per region",
    y = "Number of regions",
    title = "Sanity check: motif frequencies should sum near 1"
  )

ggsave(
  file.path(figure_dir, "sanity_motif_frequency_sum.png"),
  p_motif_sum,
  width = 7,
  height = 5,
  dpi = 300
)

p_beta <- ggplot(
  motif_df,
  aes(x = methylation_class, y = beta, fill = methylation_class)
) +
  geom_boxplot(width = 0.50, outlier.shape = NA, alpha = 0.9, color = "black") +
  geom_jitter(width = 0.10, alpha = 0.35, size = 0.35, color = "grey50") +
  scale_fill_manual(values = class_colors) +
  labs(
    x = "Methylation class",
    y = "Beta value",
    title = "Selected methylation-defined regions"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(figure_dir, "sanity_selected_region_beta_values.png"),
  p_beta,
  width = 6,
  height = 5,
  dpi = 300
)

p_count <- count_summary %>%
  ggplot(aes(x = methylation_class, y = total_count, group = sample)) +
  geom_line(color = "grey60", alpha = 0.5) +
  geom_point(aes(fill = methylation_class), shape = 21, size = 2.5, color = "black") +
  scale_fill_manual(values = class_colors) +
  labs(
    x = "Methylation class",
    y = "Total end-motif count",
    title = "End-motif fragment representation"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(figure_dir, "sanity_total_endmotif_count.png"),
  p_count,
  width = 6,
  height = 5,
  dpi = 300
)

# ============================================================
# 8. Sample-level weighted motif summary
# ============================================================

sample_motif_summary <- motif_df %>%
  group_by(sample, methylation_class) %>%
  summarise(
    across(
      all_of(motif_cols),
      ~ safe_weighted_mean(.x, w = count)
    ),
    total_count = sum(count, na.rm = TRUE),
    mean_beta = mean(beta, na.rm = TRUE),
    n_regions = n(),
    .groups = "drop"
  ) %>%
  mutate(
    methylation_class = factor(methylation_class, levels = class_levels)
  )

write_tsv(
  sample_motif_summary,
  file.path(table_dir, "sample_level_weighted_endmotif_profiles.tsv")
)

write_tsv(
  sample_motif_summary %>%
    select(sample, methylation_class, total_count, mean_beta, n_regions),
  file.path(table_dir, "sample_level_weighted_endmotif_summary.tsv")
)

# ============================================================
# 9. MDS analysis
# ============================================================

sample_mds <- sample_motif_summary %>%
  rowwise() %>%
  mutate(
    MDS = calculate_mds(c_across(all_of(motif_cols)))
  ) %>%
  ungroup() %>%
  select(sample, methylation_class, n_regions, total_count, mean_beta, MDS)

mds_wide <- sample_mds %>%
  select(sample, methylation_class, total_count, MDS) %>%
  pivot_wider(
    names_from = methylation_class,
    values_from = c(total_count, MDS)
  ) %>%
  mutate(
    delta_MDS_hypo_minus_hyper = MDS_hypo - MDS_hyper,
    count_ratio_hypo_to_hyper = total_count_hypo / total_count_hyper
  )

mds_wilcox_p <- safe_paired_wilcox(mds_wide$MDS_hypo, mds_wide$MDS_hyper)
mds_ttest_p <- safe_paired_ttest(mds_wide$MDS_hypo, mds_wide$MDS_hyper)

mds_wilcox_stat <- tryCatch(
  unname(wilcox.test(mds_wide$MDS_hypo, mds_wide$MDS_hyper, paired = TRUE, exact = FALSE)$statistic),
  error = function(e) NA_real_
)

mds_ttest_stat <- tryCatch(
  unname(t.test(mds_wide$MDS_hypo, mds_wide$MDS_hyper, paired = TRUE)$statistic),
  error = function(e) NA_real_
)

mds_tests <- tibble(
  test = c("paired_wilcoxon", "paired_t_test"),
  statistic = c(mds_wilcox_stat, mds_ttest_stat),
  p_value = c(mds_wilcox_p, mds_ttest_p),
  mean_hyper = mean(mds_wide$MDS_hyper, na.rm = TRUE),
  mean_hypo = mean(mds_wide$MDS_hypo, na.rm = TRUE),
  median_hyper = median(mds_wide$MDS_hyper, na.rm = TRUE),
  median_hypo = median(mds_wide$MDS_hypo, na.rm = TRUE),
  mean_delta_hypo_minus_hyper = mean(mds_wide$delta_MDS_hypo_minus_hyper, na.rm = TRUE),
  median_delta_hypo_minus_hyper = median(mds_wide$delta_MDS_hypo_minus_hyper, na.rm = TRUE)
)

write_tsv(sample_mds, file.path(table_dir, "hyper_hypo_sample_mds.tsv"))
write_tsv(mds_wide, file.path(table_dir, "hyper_hypo_paired_mds_wide.tsv"))
write_tsv(mds_tests, file.path(table_dir, "hyper_hypo_mds_paired_tests.tsv"))

sample_mds <- sample_mds %>%
  mutate(methylation_class = factor(methylation_class, levels = class_levels))

mds_sig_label <- p_to_label(mds_wilcox_p)

y_max <- max(sample_mds$MDS, na.rm = TRUE)
y_min <- min(sample_mds$MDS, na.rm = TRUE)
y_range <- y_max - y_min

if (!is.finite(y_range) || y_range == 0) {
  y_range <- abs(y_max) * 0.1 + 1e-6
}

y_sig <- y_max + 0.20 * y_range
tick <- 0.025 * y_range
text_offset <- 0.025 * y_range

p_mds <- ggplot(sample_mds, aes(x = methylation_class, y = MDS)) +
  geom_violin(
    aes(fill = methylation_class),
    trim = FALSE,
    color = "black",
    alpha = 0.6
  ) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    fill = "white",
    color = "black"
  ) +
  geom_line(
    aes(group = sample),
    color = "grey60",
    alpha = 0.45
  ) +
  geom_point(
    aes(fill = methylation_class),
    shape = 21,
    size = 2.3,
    color = "black"
  ) +
  scale_fill_manual(values = class_colors) +
  annotate("segment", x = 1, xend = 2, y = y_sig, yend = y_sig, linewidth = 0.7) +
  annotate("segment", x = 1, xend = 1, y = y_sig, yend = y_sig - tick, linewidth = 0.7) +
  annotate("segment", x = 2, xend = 2, y = y_sig, yend = y_sig - tick, linewidth = 0.7) +
  annotate(
    "text",
    x = 1.5,
    y = y_sig + text_offset,
    label = mds_sig_label,
    size = 5,
    family = plot_family
  ) +
  coord_cartesian(
    ylim = c(y_min - 0.10 * y_range, y_max + 0.25 * y_range),
    clip = "off"
  ) +
  labs(
    x = "Methylation class",
    y = "Motif diversity score",
    title = "End-motif diversity: hyper vs hypo"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(figure_dir, "hyper_hypo_mds_paired_violin.png"),
  p_mds,
  width = 6,
  height = 5.5,
  dpi = 300
)

# ============================================================
# 10. Motif-level paired tests
# ============================================================

motif_long <- sample_motif_summary %>%
  select(sample, methylation_class, all_of(motif_cols)) %>%
  pivot_longer(
    cols = all_of(motif_cols),
    names_to = "motif",
    values_to = "frequency"
  )

motif_wide_by_sample <- motif_long %>%
  pivot_wider(
    id_cols = c(sample, motif),
    names_from = methylation_class,
    values_from = frequency
  ) %>%
  mutate(
    diff_hypo_minus_hyper = hypo - hyper,
    log2FC_hypo_vs_hyper = log2((hypo + eps) / (hyper + eps))
  )

motif_tests <- motif_wide_by_sample %>%
  group_by(motif) %>%
  summarise(
    n_samples = sum(!is.na(hyper) & !is.na(hypo)),
    mean_hyper = mean(hyper, na.rm = TRUE),
    mean_hypo = mean(hypo, na.rm = TRUE),
    median_hyper = median(hyper, na.rm = TRUE),
    median_hypo = median(hypo, na.rm = TRUE),
    mean_diff_hypo_minus_hyper = mean(diff_hypo_minus_hyper, na.rm = TRUE),
    median_diff_hypo_minus_hyper = median(diff_hypo_minus_hyper, na.rm = TRUE),
    mean_log2FC_hypo_vs_hyper = mean(log2FC_hypo_vs_hyper, na.rm = TRUE),
    p_value = safe_paired_wilcox(hypo, hyper),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    neg_log10_padj = -log10(p_adj_BH + 1e-300),
    direction = case_when(
      p_adj_BH < padj_cutoff & mean_log2FC_hypo_vs_hyper > log2fc_cutoff ~ "Hypo-enriched",
      p_adj_BH < padj_cutoff & mean_log2FC_hypo_vs_hyper < -log2fc_cutoff ~ "Hyper-enriched",
      p_adj_BH < padj_cutoff ~ "Significant, small effect",
      TRUE ~ "Not significant"
    )
  ) %>%
  arrange(p_adj_BH, p_value)

motif_direction_summary <- motif_tests %>%
  count(direction)

write_tsv(motif_long, file.path(table_dir, "hyper_hypo_motif_long.tsv"))
write_tsv(motif_wide_by_sample, file.path(table_dir, "hyper_hypo_motif_wide_by_sample.tsv"))
write_tsv(motif_tests, file.path(table_dir, "hyper_hypo_motif_paired_tests.tsv"))
write_tsv(motif_direction_summary, file.path(table_dir, "hyper_hypo_motif_direction_summary.tsv"))

# ============================================================
# 11. Volcano plot
# ============================================================

p_volcano <- ggplot(
  motif_tests,
  aes(x = mean_log2FC_hypo_vs_hyper, y = neg_log10_padj, color = direction)
) +
  geom_point(size = 2, alpha = 0.90) +
  geom_vline(
    xintercept = c(-log2fc_cutoff, log2fc_cutoff),
    linetype = "dashed",
    color = "grey35",
    linewidth = 0.4
  ) +
  geom_hline(
    yintercept = -log10(padj_cutoff),
    linetype = "dashed",
    color = "grey35",
    linewidth = 0.4
  ) +
  scale_color_manual(values = volcano_colors) +
  coord_cartesian(xlim = c(-2, 2)) +
  labs(
    x = "Mean log2FC: hypo vs hyper",
    y = "-log10 adjusted p-value",
    color = "Motif category",
    title = "End-motif differences between hypo and hyper regions"
  )

ggsave(
  file.path(figure_dir, "hyper_hypo_motif_volcano.png"),
  p_volcano,
  width = 8,
  height = 6,
  dpi = 300
)

# ============================================================
# 12. Top motif plots
# ============================================================

top10_significant <- motif_tests %>%
  filter(p_adj_BH < padj_cutoff, abs(mean_log2FC_hypo_vs_hyper) > log2fc_cutoff) %>%
  arrange(p_adj_BH, desc(abs(mean_log2FC_hypo_vs_hyper))) %>%
  slice_head(n = 10) %>%
  pull(motif)

if (length(top10_significant) < 10) {
  top10_significant <- motif_tests %>%
    filter(p_adj_BH < padj_cutoff) %>%
    arrange(p_adj_BH) %>%
    slice_head(n = 10) %>%
    pull(motif)
}

top10_hypo <- motif_tests %>%
  filter(p_adj_BH < padj_cutoff) %>%
  arrange(desc(mean_log2FC_hypo_vs_hyper)) %>%
  slice_head(n = 10) %>%
  pull(motif)

top10_hyper <- motif_tests %>%
  filter(p_adj_BH < padj_cutoff) %>%
  arrange(mean_log2FC_hypo_vs_hyper) %>%
  slice_head(n = 10) %>%
  pull(motif)

write_tsv(
  motif_tests %>%
    arrange(p_adj_BH) %>%
    slice_head(n = 30),
  file.path(table_dir, "top30_most_significant_motifs.tsv")
)

write_tsv(
  motif_tests %>%
    arrange(desc(mean_log2FC_hypo_vs_hyper)) %>%
    slice_head(n = 30),
  file.path(table_dir, "top30_hypo_enriched_motifs.tsv")
)

write_tsv(
  motif_tests %>%
    arrange(mean_log2FC_hypo_vs_hyper) %>%
    slice_head(n = 30),
  file.path(table_dir, "top30_hyper_enriched_motifs.tsv")
)

plot_top_motif_set <- function(motifs, plot_title, output_prefix) {
  if (length(motifs) == 0) {
    warning("No motifs found for ", plot_title)
    return(NULL)
  }
  
  plot_df <- motif_long %>%
    filter(motif %in% motifs) %>%
    group_by(methylation_class, motif) %>%
    summarise(
      mean_frequency = mean(frequency, na.rm = TRUE),
      se_frequency = sd(frequency, na.rm = TRUE) / sqrt(sum(!is.na(frequency))),
      .groups = "drop"
    ) %>%
    mutate(
      motif = factor(motif, levels = rev(motifs))
    )
  
  p <- ggplot(
    plot_df,
    aes(x = motif, y = mean_frequency, fill = methylation_class)
  ) +
    geom_col(
      position = position_dodge(width = 0.8),
      width = 0.7,
      color = "black",
      linewidth = 0.25
    ) +
    geom_errorbar(
      aes(
        ymin = mean_frequency - se_frequency,
        ymax = mean_frequency + se_frequency
      ),
      position = position_dodge(width = 0.8),
      width = 0.2
    ) +
    coord_flip() +
    scale_fill_manual(values = class_colors) +
    labs(
      x = "4-mer motif",
      y = "Mean motif frequency",
      fill = "Methylation class",
      title = plot_title
    )
  
  ggsave(
    file.path(figure_dir, paste0(output_prefix, ".png")),
    p,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  return(p)
}

plot_top_motif_set(
  top10_hypo,
  "Top 10 hypo-enriched end motifs",
  "top10_hypo_enriched_end_motifs"
)

plot_top_motif_set(
  top10_hyper,
  "Top 10 hyper-enriched end motifs",
  "top10_hyper_enriched_end_motifs"
)

if (length(top10_significant) > 0) {
  top10_paired_df <- motif_wide_by_sample %>%
    filter(motif %in% top10_significant) %>%
    pivot_longer(
      cols = c(hyper, hypo),
      names_to = "methylation_class",
      values_to = "frequency"
    ) %>%
    mutate(
      methylation_class = factor(methylation_class, levels = class_levels),
      motif = factor(motif, levels = top10_significant)
    )
  
  p_top10_paired <- ggplot(
    top10_paired_df,
    aes(x = methylation_class, y = frequency, group = sample)
  ) +
    geom_line(color = "grey55", alpha = 0.45) +
    geom_point(aes(fill = methylation_class), shape = 21, size = 2, color = "black", alpha = 0.9) +
    facet_wrap(~ motif, scales = "free_y", ncol = 5) +
    scale_fill_manual(values = class_colors) +
    labs(
      x = "Methylation class",
      y = "Motif frequency",
      fill = "Class",
      title = "Paired sample-level changes in top altered motifs"
    )
  
  ggsave(
    file.path(figure_dir, "top10_altered_motifs_paired_samples.png"),
    p_top10_paired,
    width = 11,
    height = 6,
    dpi = 300
  )
}

# ============================================================
# 13. Heatmaps
# ============================================================

top_heatmap_motifs <- motif_tests %>%
  filter(p_adj_BH < padj_cutoff) %>%
  arrange(p_adj_BH, desc(abs(mean_log2FC_hypo_vs_hyper))) %>%
  slice_head(n = 50) %>%
  pull(motif)

if (length(top_heatmap_motifs) >= 2) {
  heatmap_df <- sample_motif_summary %>%
    select(sample, methylation_class, all_of(top_heatmap_motifs)) %>%
    unite("sample_class", sample, methylation_class, remove = FALSE) %>%
    column_to_rownames("sample_class") %>%
    select(all_of(top_heatmap_motifs))
  
  heatmap_mat <- as.matrix(heatmap_df)
  heatmap_mat_scaled <- scale(heatmap_mat)
  
  annotation_df <- sample_motif_summary %>%
    select(sample, methylation_class) %>%
    unite("sample_class", sample, methylation_class, remove = FALSE) %>%
    column_to_rownames("sample_class") %>%
    select(methylation_class)
  
  annotation_colors <- list(
    methylation_class = class_colors
  )
  
  png(
    filename = file.path(figure_dir, "heatmap_top50_differential_end_motifs.png"),
    width = 3000,
    height = 2400,
    res = 300
  )
  
  pheatmap(
    heatmap_mat_scaled,
    annotation_row = annotation_df,
    annotation_colors = annotation_colors,
    color = colorRampPalette(c("lightblue1", "white", "lightpink"))(100),
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    fontsize = 10,
    fontsize_row = 7,
    fontsize_col = 9,
    border_color = NA,
    main = "Top 50 differential end motifs"
  )
  
  dev.off()
  
  diff_heatmap_df <- motif_wide_by_sample %>%
    filter(motif %in% top_heatmap_motifs) %>%
    select(sample, motif, diff_hypo_minus_hyper) %>%
    pivot_wider(
      names_from = motif,
      values_from = diff_hypo_minus_hyper
    ) %>%
    column_to_rownames("sample")
  
  diff_mat <- as.matrix(diff_heatmap_df)
  
  png(
    filename = file.path(figure_dir, "heatmap_paired_motif_differences_hypo_minus_hyper.png"),
    width = 3000,
    height = 2200,
    res = 300
  )
  
  pheatmap(
    diff_mat,
    color = colorRampPalette(c("lightblue1", "white", "lightpink"))(100),
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    fontsize = 9,
    border_color = NA,
    main = "Paired motif differences: hypo minus hyper"
  )
  
  dev.off()
} else {
  warning("Fewer than 2 significant motifs. Skipping heatmaps.")
}

# ============================================================
# 14. PCA of motif profiles
# ============================================================

pca_input <- sample_motif_summary %>%
  select(sample, methylation_class, all_of(motif_cols))

pca_mat <- pca_input %>%
  select(all_of(motif_cols)) %>%
  as.matrix()

pca_res <- prcomp(pca_mat, center = TRUE, scale. = TRUE)

pca_df <- as_tibble(pca_res$x[, 1:4]) %>%
  bind_cols(pca_input %>% select(sample, methylation_class))

var_explained <- summary(pca_res)$importance[2, 1:2] * 100

write_tsv(
  pca_df,
  file.path(table_dir, "hyper_hypo_motif_pca_coordinates.tsv")
)

p_pca <- ggplot(
  pca_df,
  aes(x = PC1, y = PC2, fill = methylation_class)
) +
  geom_point(shape = 21, size = 3.2, color = "black", alpha = 0.9) +
  stat_ellipse(aes(color = methylation_class), linewidth = 0.8, alpha = 0.8) +
  scale_fill_manual(values = class_colors) +
  scale_color_manual(values = class_colors) +
  labs(
    x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "%)"),
    fill = "Methylation class",
    color = "Methylation class",
    title = "PCA of 4-mer end-motif profiles"
  )

ggsave(
  file.path(figure_dir, "hyper_hypo_motif_pca.png"),
  p_pca,
  width = 7,
  height = 6,
  dpi = 300
)

# ============================================================
# 15. Base probability / PWM analysis
# ============================================================

class_profiles <- sample_motif_summary %>%
  group_by(methylation_class) %>%
  summarise(
    across(
      all_of(motif_cols),
      ~ weighted.mean(.x, w = total_count, na.rm = TRUE)
    ),
    .groups = "drop"
  )

hyper_profile <- class_profiles %>%
  filter(methylation_class == "hyper") %>%
  select(all_of(motif_cols)) %>%
  unlist()

hypo_profile <- class_profiles %>%
  filter(methylation_class == "hypo") %>%
  select(all_of(motif_cols)) %>%
  unlist()

names(hyper_profile) <- motif_cols
names(hypo_profile) <- motif_cols

hyper_pwm <- make_pwm_from_motif_profile(hyper_profile)
hypo_pwm <- make_pwm_from_motif_profile(hypo_profile)
diff_pwm <- hypo_pwm - hyper_pwm

write_tsv(
  as.data.frame(hyper_pwm) %>% rownames_to_column("base"),
  file.path(table_dir, "hyper_pwm.tsv")
)

write_tsv(
  as.data.frame(hypo_pwm) %>% rownames_to_column("base"),
  file.path(table_dir, "hypo_pwm.tsv")
)

write_tsv(
  as.data.frame(diff_pwm) %>% rownames_to_column("base"),
  file.path(table_dir, "hypo_minus_hyper_pwm.tsv")
)

diff_pwm_long <- as.data.frame(diff_pwm) %>%
  rownames_to_column("base") %>%
  pivot_longer(
    cols = -base,
    names_to = "position",
    values_to = "difference"
  ) %>%
  mutate(
    position = factor(position, levels = as.character(1:4)),
    base = factor(base, levels = c("A", "C", "G", "T"))
  )

p_pwm_diff <- ggplot(
  diff_pwm_long,
  aes(x = position, y = base, fill = difference)
) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(difference, 3)), size = 3.5, family = plot_family) +
  scale_fill_gradient2(
    low = "lightblue1",
    mid = "white",
    high = "lightpink",
    midpoint = 0
  ) +
  labs(
    x = "Position in 4-mer",
    y = "Base",
    fill = "Hypo - hyper",
    title = "Difference in base probability by position"
  )

ggsave(
  file.path(figure_dir, "base_probability_difference_hypo_minus_hyper.png"),
  p_pwm_diff,
  width = 6,
  height = 4.5,
  dpi = 300
)

pwm_long <- bind_rows(
  as.data.frame(hyper_pwm) %>%
    rownames_to_column("base") %>%
    pivot_longer(cols = -base, names_to = "position", values_to = "probability") %>%
    mutate(methylation_class = "hyper"),
  as.data.frame(hypo_pwm) %>%
    rownames_to_column("base") %>%
    pivot_longer(cols = -base, names_to = "position", values_to = "probability") %>%
    mutate(methylation_class = "hypo")
) %>%
  mutate(
    methylation_class = factor(methylation_class, levels = class_levels),
    position = factor(position, levels = as.character(1:4)),
    base = factor(base, levels = c("A", "C", "G", "T"))
  )

p_pwm_prob <- ggplot(
  pwm_long,
  aes(x = position, y = probability, fill = base)
) +
  geom_col(position = "stack", color = "black", linewidth = 0.2) +
  facet_wrap(~ methylation_class) +
  labs(
    x = "Position in 4-mer",
    y = "Base probability",
    fill = "Base",
    title = "5' end-motif base composition by methylation class"
  )

ggsave(
  file.path(figure_dir, "base_probability_hyper_vs_hypo.png"),
  p_pwm_prob,
  width = 8,
  height = 4.5,
  dpi = 300
)

# ============================================================
# 16. Motif SD / variability summary
# ============================================================

motif_sd_summary <- motif_long %>%
  group_by(motif, methylation_class) %>%
  summarise(
    n_samples = sum(!is.na(frequency)),
    mean_frequency = mean(frequency, na.rm = TRUE),
    sd_frequency = sd(frequency, na.rm = TRUE),
    se_frequency = sd_frequency / sqrt(n_samples),
    cv_frequency = sd_frequency / mean_frequency,
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = methylation_class,
    values_from = c(
      n_samples,
      mean_frequency,
      sd_frequency,
      se_frequency,
      cv_frequency
    )
  ) %>%
  mutate(
    mean_diff_hypo_minus_hyper = mean_frequency_hypo - mean_frequency_hyper,
    sd_diff_hypo_minus_hyper = sd_frequency_hypo - sd_frequency_hyper,
    abs_sd_diff = abs(sd_diff_hypo_minus_hyper),
    max_sd = pmax(sd_frequency_hyper, sd_frequency_hypo, na.rm = TRUE)
  ) %>%
  arrange(desc(max_sd))

motif_tests_with_sd <- motif_tests %>%
  left_join(motif_sd_summary, by = "motif") %>%
  arrange(p_adj_BH, p_value)

write_tsv(motif_sd_summary, file.path(table_dir, "motif_sd_summary_hyper_hypo.tsv"))
write_tsv(motif_tests_with_sd, file.path(table_dir, "motif_tests_with_sd_hyper_hypo.tsv"))

p_sd_scatter <- ggplot(
  motif_sd_summary,
  aes(x = sd_frequency_hyper, y = sd_frequency_hypo)
) +
  geom_point(alpha = 0.7, size = 2, color = "blue3") +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "grey40"
  ) +
  labs(
    x = "SD of motif frequency in hyper regions",
    y = "SD of motif frequency in hypo regions",
    title = "Motif-level variability across samples"
  )

ggsave(
  file.path(figure_dir, "motif_sd_hyper_vs_hypo_scatter.png"),
  p_sd_scatter,
  width = 6,
  height = 5,
  dpi = 300
)

# ============================================================
# 17. CpG-containing vs non-CpG motif burden
# ============================================================

motif_category_df <- tibble(
  motif = motif_cols,
  motif_class = ifelse(
    grepl("CG", motif_cols),
    "CpG-containing motifs",
    "non-CpG motifs"
  )
)

cpg_cols <- motif_category_df$motif[
  motif_category_df$motif_class == "CpG-containing motifs"
]

non_cpg_cols <- motif_category_df$motif[
  motif_category_df$motif_class == "non-CpG motifs"
]

cpg_burden_df <- sample_motif_summary %>%
  mutate(
    CpG_containing = rowSums(across(all_of(cpg_cols)), na.rm = TRUE),
    non_CpG = rowSums(across(all_of(non_cpg_cols)), na.rm = TRUE),
    total_burden_check = CpG_containing + non_CpG
  ) %>%
  select(
    sample,
    methylation_class,
    CpG_containing,
    non_CpG,
    total_burden_check
  )

cpg_burden_long <- cpg_burden_df %>%
  pivot_longer(
    cols = c(CpG_containing, non_CpG),
    names_to = "motif_class",
    values_to = "burden"
  ) %>%
  mutate(
    methylation_class = factor(methylation_class, levels = class_levels),
    motif_class = recode(
      motif_class,
      CpG_containing = "CpG-containing motifs",
      non_CpG = "non-CpG motifs"
    ),
    motif_class = factor(
      motif_class,
      levels = c("CpG-containing motifs", "non-CpG motifs")
    )
  )

cpg_summary <- cpg_burden_long %>%
  group_by(methylation_class, motif_class) %>%
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
  pivot_wider(
    id_cols = c(sample, motif_class),
    names_from = methylation_class,
    values_from = burden
  ) %>%
  group_by(motif_class) %>%
  summarise(
    n_samples = sum(!is.na(hyper) & !is.na(hypo)),
    mean_hyper = mean(hyper, na.rm = TRUE),
    mean_hypo = mean(hypo, na.rm = TRUE),
    median_hyper = median(hyper, na.rm = TRUE),
    median_hypo = median(hypo, na.rm = TRUE),
    mean_diff_hypo_minus_hyper = mean(hypo - hyper, na.rm = TRUE),
    p_value = safe_paired_wilcox(hypo, hyper),
    .groups = "drop"
  ) %>%
  mutate(
    padj = p.adjust(p_value, method = "BH"),
    label = p_to_label(padj)
  )

write_tsv(motif_category_df, file.path(table_dir, "motif_cpg_category.tsv"))
write_tsv(cpg_burden_df, file.path(table_dir, "cpg_non_cpg_burden_by_sample.tsv"))
write_tsv(cpg_summary, file.path(table_dir, "cpg_non_cpg_burden_summary.tsv"))
write_tsv(cpg_tests, file.path(table_dir, "cpg_non_cpg_burden_paired_tests.tsv"))

sig_df <- cpg_burden_long %>%
  group_by(motif_class) %>%
  summarise(
    y_max = max(burden, na.rm = TRUE),
    y_min = min(burden, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    y_range = ifelse(
      y_max - y_min == 0,
      abs(y_max) * 0.1 + 1e-6,
      y_max - y_min
    ),
    y_sig = y_max + 0.28 * y_range,
    y_text = y_max + 0.34 * y_range,
    tick = 0.04 * y_range
  ) %>%
  left_join(
    cpg_tests %>% select(motif_class, label, padj),
    by = "motif_class"
  )

p_cpg_burden <- ggplot(
  cpg_burden_long,
  aes(x = methylation_class, y = burden, fill = methylation_class)
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
  geom_line(
    aes(group = sample),
    color = "grey60",
    alpha = 0.35
  ) +
  geom_point(
    aes(color = methylation_class),
    position = position_jitter(width = 0.04, height = 0),
    size = 2,
    alpha = 0.85
  ) +
  facet_wrap(~ motif_class, scales = "free_y") +
  scale_fill_manual(values = class_colors) +
  scale_color_manual(values = class_colors) +
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
    family = plot_family
  ) +
  labs(
    x = "Methylation class",
    y = "Sum of end-motif frequencies",
    title = "CpG-containing and non-CpG end-motif burden"
  ) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 13)
  )

ggsave(
  file.path(figure_dir, "cpg_non_cpg_burden_hyper_vs_hypo.png"),
  p_cpg_burden,
  width = 9,
  height = 5.5,
  dpi = 300,
  bg = "white"
)

# ============================================================
# 18. Save R objects
# ============================================================

saveRDS(
  list(
    motif_df = motif_df,
    sample_motif_summary = sample_motif_summary,
    sample_mds = sample_mds,
    mds_wide = mds_wide,
    motif_long = motif_long,
    motif_wide_by_sample = motif_wide_by_sample,
    motif_tests = motif_tests,
    motif_tests_with_sd = motif_tests_with_sd,
    hyper_pwm = hyper_pwm,
    hypo_pwm = hypo_pwm,
    diff_pwm = diff_pwm,
    cpg_burden_df = cpg_burden_df,
    cpg_tests = cpg_tests
  ),
  file.path(table_dir, "hyper_hypo_endmotif_analysis_objects.rds")
)

# ============================================================
# 19. Final message
# ============================================================

cat("\n============================================================\n")
cat("Hyper/hypo end-motif analysis completed.\n")

cat("\nInput motif directory:\n")
cat(motif_dir, "\n")

cat("\nSamples used:\n")
print(sample_file_check)

cat("\nBasic summary:\n")
print(basic_summary)

cat("\nMDS paired tests:\n")
print(mds_tests)

cat("\nMotif direction summary:\n")
print(motif_direction_summary)

cat("\nCpG/non-CpG burden tests:\n")
print(cpg_tests)

cat("\nTables saved in:\n")
cat(table_dir, "\n")

cat("\nFigures saved in:\n")
cat(figure_dir, "\n")

cat("============================================================\n")