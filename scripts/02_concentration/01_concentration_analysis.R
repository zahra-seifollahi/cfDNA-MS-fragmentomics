# ============================================================
# cfDNA concentration / quantity statistical analysis
#
# Purpose:
#   Analyze Qubit cfDNA quantity/concentration data.
#
# Main analysis:
#   DNA Quantity (ng)
#
# Additional recorded metric:
#   Concentration (ng/ul)
#
# Important:
#   This script uses all available Qubit samples.
#   It does not apply BAM-derived QC filtering.
#
# Usage:
#   Rscript scripts/02_concentration/01_concentration_analysis.R /path/to/cfdna_quant.xlsx
#
# If no input path is provided, default is:
#   data/example/cfdna_quant.xlsx
#
# Output:
#   results/tables/concentration/
#   results/figures/concentration/
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(writexl)
  library(tibble)
})

graphics.off()

# ============================================================
# 1. Settings
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

file_path <- ifelse(
  length(args) >= 1,
  args[1],
  "data/example/cfdna_quant.xlsx"
)

table_dir <- "results/tables/concentration"
figure_dir <- "results/figures/concentration"

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c("Healthy", "Remission", "Relapse")

group_cols <- c(
  Healthy = "darkgreen",
  Remission = "deeppink3",
  Relapse = "blue3"
)

group_fills <- c(
  Healthy = "darkseagreen3",
  Remission = "lightpink",
  Relapse = "lightblue"
)

ms_cols <- c(
  Healthy = "cornflowerblue",
  MS = "tomato"
)

ms_fills <- c(
  Healthy = "lightblue",
  MS = "mistyrose"
)

plot_family <- "serif"

base_theme <- theme_classic(base_size = 14, base_family = plot_family) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, family = plot_family),
    axis.title = element_text(color = "black", family = plot_family),
    axis.text = element_text(color = "black", family = plot_family),
    legend.position = "none"
  )

main_metric <- "dna_quantity_ng"
main_metric_label <- "DNA Quantity (ng)"

# ============================================================
# 2. Helper functions
# ============================================================

clean_group <- function(group) {
  case_when(
    str_detect(group, regex("healthy", ignore_case = TRUE)) ~ "Healthy",
    str_detect(group, regex("remission", ignore_case = TRUE)) ~ "Remission",
    str_detect(group, regex("relapse", ignore_case = TRUE)) ~ "Relapse",
    TRUE ~ as.character(group)
  )
}

p_to_star <- function(p) {
  case_when(
    is.na(p) ~ "NA",
    p <= 0.0001 ~ "****",
    p <= 0.001 ~ "***",
    p <= 0.01 ~ "**",
    p <= 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

format_p <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p < 0.001, "p < 0.001", paste0("p = ", signif(p, 2))))
}

safe_shapiro <- function(x) {
  x <- x[!is.na(x)]

  if (length(x) < 3) {
    return(NA_real_)
  }

  if (length(unique(x)) < 3) {
    return(NA_real_)
  }

  tryCatch(
    shapiro.test(x)$p.value,
    error = function(e) NA_real_
  )
}

pairwise_to_table <- function(pairwise_result, variable_name) {
  as.data.frame(as.table(pairwise_result$p.value)) %>%
    rename(
      group2 = Var1,
      group1 = Var2,
      p_adjusted = Freq
    ) %>%
    filter(!is.na(p_adjusted)) %>%
    mutate(
      variable = variable_name,
      test = "Pairwise Wilcoxon rank-sum test",
      p_value = p_adjusted,
      p_adjust_method = "BH",
      significance = p_to_star(p_adjusted)
    ) %>%
    select(
      variable,
      test,
      group1,
      group2,
      p_value,
      p_adjusted,
      p_adjust_method,
      significance
    )
}

get_pw <- function(pw, row, col) {
  value <- tryCatch(
    pw$p.value[row, col],
    error = function(e) NA_real_
  )

  as.numeric(value)
}

make_summary <- function(data, group_col, value_col) {
  data %>%
    group_by(.data[[group_col]]) %>%
    summarise(
      n = n(),
      mean = mean(.data[[value_col]], na.rm = TRUE),
      median = median(.data[[value_col]], na.rm = TRUE),
      sd = sd(.data[[value_col]], na.rm = TRUE),
      min = min(.data[[value_col]], na.rm = TRUE),
      max = max(.data[[value_col]], na.rm = TRUE),
      IQR = IQR(.data[[value_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(group = .data[[group_col]])
}

# ============================================================
# 3. Load data
# ============================================================

if (!file.exists(file_path)) {
  stop("Input concentration file not found: ", file_path)
}

df_raw <- read_excel(file_path)

required_cols <- c("Sample ID", "DNA Quantity (ng)", "Group", "Concentration (ng/ul)")
missing_cols <- setdiff(required_cols, colnames(df_raw))

if (length(missing_cols) > 0) {
  stop(
    "Input file is missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

df <- df_raw %>%
  transmute(
    sample_id = as.character(`Sample ID`),
    dna_quantity_ng = readr::parse_number(as.character(`DNA Quantity (ng)`)),
    concentration_ng_ul = readr::parse_number(as.character(`Concentration (ng/ul)`)),
    group = clean_group(as.character(Group))
  ) %>%
  mutate(
    group = factor(group, levels = group_levels),
    cfdna_value = .data[[main_metric]]
  ) %>%
  filter(!is.na(sample_id), !is.na(group), !is.na(cfdna_value))

cat("\nSamples loaded:\n")
print(df %>% count(group))

write_tsv(
  df,
  file.path(table_dir, "cfdna_concentration_clean_data.tsv")
)

# ============================================================
# 4. Summary statistics
# ============================================================

group_summary_quantity <- make_summary(df, "group", "dna_quantity_ng") %>%
  mutate(metric = "DNA Quantity (ng)") %>%
  select(metric, everything())

group_summary_concentration <- make_summary(df, "group", "concentration_ng_ul") %>%
  mutate(metric = "Concentration (ng/ul)") %>%
  select(metric, everything())

group_summary <- bind_rows(
  group_summary_quantity,
  group_summary_concentration
)

write_tsv(
  group_summary,
  file.path(table_dir, "cfdna_group_summary.tsv")
)

# ============================================================
# 5. Normality test per group
# ============================================================

normality_results <- df %>%
  group_by(group) %>%
  summarise(
    variable = main_metric,
    n = n(),
    shapiro_p_value = safe_shapiro(cfdna_value),
    normality = ifelse(shapiro_p_value < 0.05, "Non-normal", "Not rejected"),
    .groups = "drop"
  )

write_tsv(
  normality_results,
  file.path(table_dir, "cfdna_shapiro_normality_results.tsv")
)

# ============================================================
# 6. Three-group comparison
# ============================================================

kruskal_base <- kruskal.test(cfdna_value ~ group, data = df)

kruskal_result <- tibble(
  variable = main_metric,
  comparison = "Healthy vs Remission vs Relapse",
  test = "Kruskal-Wallis rank-sum test",
  n = nrow(df),
  statistic = as.numeric(kruskal_base$statistic),
  df = as.numeric(kruskal_base$parameter),
  p_value = kruskal_base$p.value,
  significance = p_to_star(p_value)
)

pairwise_base <- pairwise.wilcox.test(
  x = df$cfdna_value,
  g = df$group,
  p.adjust.method = "BH",
  exact = FALSE
)

pairwise_three_groups <- pairwise_to_table(
  pairwise_base,
  variable_name = main_metric
)

write_tsv(
  kruskal_result,
  file.path(table_dir, "cfdna_kruskal_three_group_result.tsv")
)

write_tsv(
  pairwise_three_groups,
  file.path(table_dir, "cfdna_pairwise_wilcoxon_three_groups.tsv")
)

# ============================================================
# 7. Healthy vs MS combined
# ============================================================

df_ms <- df %>%
  mutate(
    disease_group = ifelse(group == "Healthy", "Healthy", "MS"),
    disease_group = factor(disease_group, levels = c("Healthy", "MS"))
  )

ms_summary <- make_summary(df_ms, "disease_group", "cfdna_value") %>%
  mutate(metric = main_metric_label) %>%
  select(metric, everything())

healthy_vs_ms_base <- wilcox.test(
  cfdna_value ~ disease_group,
  data = df_ms,
  exact = FALSE
)

healthy_vs_ms <- tibble(
  variable = main_metric,
  comparison = "Healthy vs MS",
  test = "Wilcoxon rank-sum test",
  group1 = "Healthy",
  group2 = "MS",
  n1 = sum(df_ms$disease_group == "Healthy"),
  n2 = sum(df_ms$disease_group == "MS"),
  statistic = as.numeric(healthy_vs_ms_base$statistic),
  p_value = healthy_vs_ms_base$p.value,
  significance = p_to_star(p_value)
)

write_tsv(
  ms_summary,
  file.path(table_dir, "cfdna_healthy_vs_ms_summary.tsv")
)

write_tsv(
  healthy_vs_ms,
  file.path(table_dir, "cfdna_healthy_vs_ms_wilcoxon.tsv")
)

# ============================================================
# 8. Healthy vs Relapse only
# ============================================================

df_healthy_relapse <- df %>%
  filter(group %in% c("Healthy", "Relapse")) %>%
  mutate(
    group = factor(group, levels = c("Healthy", "Relapse"))
  )

healthy_relapse_summary <- make_summary(df_healthy_relapse, "group", "cfdna_value") %>%
  mutate(metric = main_metric_label) %>%
  select(metric, everything())

healthy_vs_relapse_base <- wilcox.test(
  cfdna_value ~ group,
  data = df_healthy_relapse,
  exact = FALSE
)

healthy_vs_relapse <- tibble(
  variable = main_metric,
  comparison = "Healthy vs Relapse",
  test = "Wilcoxon rank-sum test",
  group1 = "Healthy",
  group2 = "Relapse",
  n1 = sum(df_healthy_relapse$group == "Healthy"),
  n2 = sum(df_healthy_relapse$group == "Relapse"),
  statistic = as.numeric(healthy_vs_relapse_base$statistic),
  p_value = healthy_vs_relapse_base$p.value,
  significance = p_to_star(p_value)
)

write_tsv(
  healthy_relapse_summary,
  file.path(table_dir, "cfdna_healthy_vs_relapse_summary.tsv")
)

write_tsv(
  healthy_vs_relapse,
  file.path(table_dir, "cfdna_healthy_vs_relapse_wilcoxon.tsv")
)

# ============================================================
# 9. Combined report-ready tables
# ============================================================

descriptive_statistics_with_p_values <- bind_rows(
  group_summary_quantity %>%
    mutate(
      analysis = "Three-group comparison",
      comparison = "Healthy vs Remission vs Relapse",
      test = "Kruskal-Wallis rank-sum test",
      p_value = kruskal_result$p_value,
      p_adjusted = NA_real_
    ),
  ms_summary %>%
    mutate(
      analysis = "Two-group comparison",
      comparison = "Healthy vs MS",
      test = "Wilcoxon rank-sum test",
      p_value = healthy_vs_ms$p_value,
      p_adjusted = NA_real_
    ),
  healthy_relapse_summary %>%
    mutate(
      analysis = "Two-group comparison",
      comparison = "Healthy vs Relapse",
      test = "Wilcoxon rank-sum test",
      p_value = healthy_vs_relapse$p_value,
      p_adjusted = NA_real_
    )
) %>%
  select(
    analysis,
    comparison,
    metric,
    group,
    n,
    mean,
    median,
    sd,
    min,
    max,
    IQR,
    test,
    p_value,
    p_adjusted
  )

p_value_summary_table <- bind_rows(
  kruskal_result %>%
    transmute(
      analysis = "Three-group comparison",
      comparison,
      test,
      p_value,
      p_adjusted = NA_real_,
      significance
    ),
  pairwise_three_groups %>%
    transmute(
      analysis = "Pairwise three-group comparison",
      comparison = paste(group1, "vs", group2),
      test,
      p_value,
      p_adjusted,
      significance
    ),
  healthy_vs_ms %>%
    transmute(
      analysis = "Two-group comparison",
      comparison,
      test,
      p_value,
      p_adjusted = NA_real_,
      significance
    ),
  healthy_vs_relapse %>%
    transmute(
      analysis = "Two-group comparison",
      comparison,
      test,
      p_value,
      p_adjusted = NA_real_,
      significance
    )
)

write_tsv(
  descriptive_statistics_with_p_values,
  file.path(table_dir, "cfdna_descriptive_statistics_with_p_values.tsv")
)

write_tsv(
  p_value_summary_table,
  file.path(table_dir, "cfdna_p_value_summary.tsv")
)

results_list <- list(
  "Clean_data" = df,
  "Group_summary" = group_summary,
  "Normality_test" = normality_results,
  "Kruskal_3_groups" = kruskal_result,
  "Pairwise_Wilcoxon_3_groups" = pairwise_three_groups,
  "Healthy_vs_MS_summary" = ms_summary,
  "Healthy_vs_MS" = healthy_vs_ms,
  "Healthy_vs_Relapse_summary" = healthy_relapse_summary,
  "Healthy_vs_Relapse" = healthy_vs_relapse,
  "Descriptive_p_values" = descriptive_statistics_with_p_values,
  "P_value_summary" = p_value_summary_table
)

write_xlsx(
  results_list,
  file.path(table_dir, "cfdna_concentration_statistical_results.xlsx")
)

# ============================================================
# 10. Plots
# ============================================================

p_three <- ggplot(
  df,
  aes(x = group, y = cfdna_value, fill = group)
) +
  geom_violin(trim = FALSE, alpha = 0.55, color = NA) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.85, color = "black") +
  geom_jitter(aes(color = group), width = 0.12, size = 2, alpha = 0.85) +
  scale_fill_manual(values = group_fills) +
  scale_color_manual(values = group_cols) +
  labs(
    title = "cfDNA quantity across Healthy, Remission, and Relapse groups",
    x = "",
    y = main_metric_label
  ) +
  base_theme

ggsave(
  filename = file.path(figure_dir, "cfdna_three_group_comparison.png"),
  plot = p_three,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  filename = file.path(figure_dir, "cfdna_three_group_comparison.pdf"),
  plot = p_three,
  width = 7,
  height = 5,
  device = "pdf"
)

p_ms <- ggplot(
  df_ms,
  aes(x = disease_group, y = cfdna_value, fill = disease_group)
) +
  geom_violin(trim = FALSE, alpha = 0.55, color = NA) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.85, color = "black") +
  geom_jitter(aes(color = disease_group), width = 0.12, size = 2, alpha = 0.85) +
  scale_fill_manual(values = ms_fills) +
  scale_color_manual(values = ms_cols) +
  labs(
    title = "cfDNA quantity: Healthy vs MS",
    x = "",
    y = main_metric_label
  ) +
  base_theme

ggsave(
  filename = file.path(figure_dir, "cfdna_healthy_vs_ms_comparison.png"),
  plot = p_ms,
  width = 6,
  height = 5,
  dpi = 300
)

ggsave(
  filename = file.path(figure_dir, "cfdna_healthy_vs_ms_comparison.pdf"),
  plot = p_ms,
  width = 6,
  height = 5,
  device = "pdf"
)

p_healthy_relapse <- ggplot(
  df_healthy_relapse,
  aes(x = group, y = cfdna_value, fill = group)
) +
  geom_violin(trim = FALSE, alpha = 0.55, color = NA) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.85, color = "black") +
  geom_jitter(aes(color = group), width = 0.12, size = 2, alpha = 0.85) +
  scale_fill_manual(values = group_fills[c("Healthy", "Relapse")]) +
  scale_color_manual(values = group_cols[c("Healthy", "Relapse")]) +
  labs(
    title = "cfDNA quantity: Healthy vs Relapse",
    x = "",
    y = main_metric_label
  ) +
  base_theme

ggsave(
  filename = file.path(figure_dir, "cfdna_healthy_vs_relapse_comparison.png"),
  plot = p_healthy_relapse,
  width = 6,
  height = 5,
  dpi = 300
)

ggsave(
  filename = file.path(figure_dir, "cfdna_healthy_vs_relapse_comparison.pdf"),
  plot = p_healthy_relapse,
  width = 6,
  height = 5,
  device = "pdf"
)

# ============================================================
# 11. Final report
# ============================================================

cat("\ncfDNA concentration / quantity analysis completed.\n")

cat("\nInput file:\n")
cat(file_path, "\n")

cat("\nSamples used:\n")
print(df %>% count(group))

cat("\nThree-group Kruskal-Wallis result:\n")
print(kruskal_result)

cat("\nPairwise Wilcoxon results:\n")
print(pairwise_three_groups)

cat("\nHealthy vs MS:\n")
print(healthy_vs_ms)

cat("\nHealthy vs Relapse:\n")
print(healthy_vs_relapse)

cat("\nOutputs saved to:\n")
cat(table_dir, "\n")
cat(figure_dir, "\n")

cat("\nDone.\n")
