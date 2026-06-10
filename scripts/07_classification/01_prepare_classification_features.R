# ============================================================
# Prepare classification feature sets
#
# Feature sets:
#   1. Global end-motif frequencies
#   2. Regional MDS
#   3. Bioanalyzer 10-bp fragment-length features
#
# Important:
#   - Classification uses QC-passed samples only.
#   - Length features are taken from Bioanalyzer, not BAM-derived length.
#   - Hyper/hypo features are not used because they are Healthy-only.
#
# Usage:
#   Rscript scripts/07_classification/01_prepare_classification_features.R
#
# Optional Bioanalyzer input:
#   Rscript scripts/07_classification/01_prepare_classification_features.R /path/to/bioanalyzer_bins_50_300_intervals.xlsx
#
# Outputs:
#   results/intermediate/classification/feature_sets/
#   results/tables/classification/
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(readxl)
  library(stringr)
  library(tibble)
})

# ============================================================
# 1. Paths
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

qc_file <- "results/tables/qc/final_sample_metadata.csv"

global_motif_file <- "results/intermediate/endmotif/global_endmotif_matrix.tsv"

regional_mds_file <- "results/tables/regional_mds/regional_MDS_matrix_samples_by_regions.tsv"

bioanalyzer_clean_file <- "results/tables/fragment_length/bioanalyzer_clean_interval_data.tsv"

bioanalyzer_input_file <- ifelse(
  length(args) >= 1,
  args[1],
  NA_character_
)

feature_dir <- "results/intermediate/classification/feature_sets"
table_dir <- "results/tables/classification"

dir.create(feature_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 2. Helper functions
# ============================================================

check_file <- function(path) {
  if (!file.exists(path)) {
    stop("Required file not found: ", path)
  }
}

clean_sample_id <- function(x) {
  y <- x %>%
    as.character() %>%
    str_remove("\\.dedup$") %>%
    str_remove("\\.frag$") %>%
    str_remove("\\.frag\\.gz$") %>%
    str_remove("\\.tsv$") %>%
    str_trim()

  num <- suppressWarnings(as.integer(str_extract(y, "\\d+")))

  y_clean <- ifelse(
    !is.na(num),
    paste0("Cap", str_pad(num, width = 2, pad = "0")),
    y
  )

  y_clean
}

sample_to_number <- function(x) {
  suppressWarnings(as.integer(str_extract(as.character(x), "\\d+")))
}

write_feature_set <- function(df, feature_set_name) {

  out_file <- file.path(
    feature_dir,
    paste0(feature_set_name, ".tsv")
  )

  write_tsv(df, out_file)

  summary_df <- tibble(
    feature_set = feature_set_name,
    output_file = out_file,
    n_samples = nrow(df),
    n_columns_total = ncol(df),
    n_features = ncol(df) - 2,
    n_missing_values = sum(is.na(df)),
    groups = paste(
      names(table(df$group)),
      as.integer(table(df$group)),
      sep = "=",
      collapse = "; "
    )
  )

  summary_df
}

# ============================================================
# 3. Check required files
# ============================================================

check_file(qc_file)
check_file(global_motif_file)
check_file(regional_mds_file)

# ============================================================
# 4. Load QC metadata
# ============================================================

metadata <- read_csv(qc_file, show_col_types = FALSE) %>%
  mutate(
    sample = clean_sample_id(sample),
    sample_num = sample_to_number(sample),
    group = as.character(group),
    include_analysis = as.character(include_analysis)
  ) %>%
  filter(include_analysis == "yes") %>%
  mutate(
    group = factor(group, levels = c("Healthy", "Remission", "Relapse"))
  ) %>%
  arrange(group, sample_num) %>%
  select(sample, sample_num, group, include_analysis, exclusion_reason)

cat("\nQC-passed samples used for classification:\n")
print(metadata %>% count(group))

write_tsv(
  metadata,
  file.path(table_dir, "classification_qc_passed_samples.tsv")
)

# ============================================================
# 5. Feature set 1: Global end-motif frequencies
# ============================================================

global_motif <- read_tsv(global_motif_file, show_col_types = FALSE) %>%
  rename(sample = 1) %>%
  mutate(sample = clean_sample_id(sample))

motif_cols <- setdiff(colnames(global_motif), "sample")
motif_cols <- motif_cols[str_detect(motif_cols, "^[ACGT]{4}$")]

if (length(motif_cols) != 256) {
  warning("Expected 256 motif columns, found: ", length(motif_cols))
}

global_motif <- global_motif %>%
  select(sample, all_of(motif_cols)) %>%
  mutate(across(all_of(motif_cols), as.numeric))

global_endmotif_256 <- metadata %>%
  select(sample, group) %>%
  inner_join(global_motif, by = "sample") %>%
  arrange(group, sample)

cat("\nFeature set: global_endmotif_256\n")
print(global_endmotif_256 %>% count(group))
cat("Features:", ncol(global_endmotif_256) - 2, "\n")
cat("Missing values:", sum(is.na(global_endmotif_256)), "\n")

# ============================================================
# 6. Feature set 2: Regional MDS
# ============================================================

regional_mds <- read_tsv(regional_mds_file, show_col_types = FALSE)

if ("sample_id" %in% colnames(regional_mds)) {
  regional_mds <- regional_mds %>%
    rename(sample = sample_id)
}

regional_mds <- regional_mds %>%
  mutate(sample = clean_sample_id(sample))

regional_feature_cols <- setdiff(colnames(regional_mds), "sample")

regional_mds <- regional_mds %>%
  mutate(across(all_of(regional_feature_cols), as.numeric))

regional_mds <- regional_mds %>%
  rename_with(
    .fn = ~ paste0("regional_MDS_", make.names(.x)),
    .cols = all_of(regional_feature_cols)
  )

regional_feature_cols <- setdiff(colnames(regional_mds), "sample")

regional_mds_167 <- metadata %>%
  select(sample, group) %>%
  inner_join(regional_mds, by = "sample") %>%
  arrange(group, sample)

cat("\nFeature set: regional_mds_167\n")
print(regional_mds_167 %>% count(group))
cat("Features:", ncol(regional_mds_167) - 2, "\n")
cat("Missing values:", sum(is.na(regional_mds_167)), "\n")

# ============================================================
# 7. Feature set 3: Bioanalyzer 10-bp length features
# ============================================================

# This section first tries to use the cleaned Bioanalyzer table generated
# by the fragment-length script. If that file does not exist, it reads the
# original Bioanalyzer Excel file provided as a command-line argument.

if (file.exists(bioanalyzer_clean_file)) {

  cat("\nReading Bioanalyzer clean table:\n")
  cat(bioanalyzer_clean_file, "\n")

  bio_clean <- read_tsv(bioanalyzer_clean_file, show_col_types = FALSE)

  # Expected columns from the cleaned script:
  # sample_id, group, range, range_start, range_end, range_mid, percent
  if (!all(c("sample_id", "range_start", "range_end", "percent") %in% colnames(bio_clean))) {
    stop(
      "Bioanalyzer clean file exists but does not contain required columns: ",
      "sample_id, range_start, range_end, percent"
    )
  }

  bio_clean <- bio_clean %>%
    mutate(
      sample = clean_sample_id(sample_id),
      range_start = as.numeric(range_start),
      range_end = as.numeric(range_end),
      percent = as.numeric(percent)
    )

} else {

  if (is.na(bioanalyzer_input_file)) {
    stop(
      "Bioanalyzer clean table not found:\n",
      bioanalyzer_clean_file,
      "\n\nProvide the original Bioanalyzer Excel file as an argument:\n",
      "Rscript scripts/07_classification/01_prepare_classification_features.R /path/to/bioanalyzer_bins_50_300_intervals.xlsx"
    )
  }

  check_file(bioanalyzer_input_file)

  cat("\nReading Bioanalyzer Excel file:\n")
  cat(bioanalyzer_input_file, "\n")

  bio_raw <- read_excel(bioanalyzer_input_file)

  required_bio_cols <- c("Sample ID", "Range", "%")

  missing_bio_cols <- setdiff(required_bio_cols, colnames(bio_raw))

  if (length(missing_bio_cols) > 0) {
    stop(
      "Bioanalyzer Excel file is missing required columns: ",
      paste(missing_bio_cols, collapse = ", ")
    )
  }

  range_parts <- str_match(
    as.character(bio_raw$Range),
    "^\\s*(\\d+(?:\\.\\d+)?)\\s*-\\s*(\\d+(?:\\.\\d+)?)\\s*$"
  )

  bio_clean <- bio_raw %>%
    transmute(
      sample = clean_sample_id(`Sample ID`),
      range = as.character(`Range`),
      range_start = as.numeric(range_parts[, 2]),
      range_end = as.numeric(range_parts[, 3]),
      percent = readr::parse_number(as.character(`%`))
    ) %>%
    filter(
      !is.na(sample),
      !is.na(range_start),
      !is.na(range_end),
      !is.na(percent)
    )
}

# Keep only 50-300 bp Bioanalyzer intervals and collapse to 10-bp bins.
# If original data are 5-bp bins, each two adjacent bins become one 10-bp feature.

bio_10bp_long <- bio_clean %>%
  filter(
    range_start >= 50,
    range_end <= 300
  ) %>%
  mutate(
    bin10_start = floor((range_start - 50) / 10) * 10 + 50,
    bin10_end = bin10_start + 10,
    bin10_mid = (bin10_start + bin10_end) / 2,
    bin10 = paste0("length_", bin10_start, "_", bin10_end, "bp")
  ) %>%
  group_by(sample, bin10, bin10_start, bin10_end, bin10_mid) %>%
  summarise(
    percent_10bp = mean(percent, na.rm = TRUE),
    n_original_bins = n(),
    .groups = "drop"
  ) %>%
  arrange(sample, bin10_start)

write_tsv(
  bio_10bp_long,
  file.path(table_dir, "bioanalyzer_10bp_length_long.tsv")
)

bioanalyzer_10bp_matrix <- bio_10bp_long %>%
  select(sample, bin10, percent_10bp) %>%
  pivot_wider(
    names_from = bin10,
    values_from = percent_10bp
  )

bio_length_cols <- setdiff(colnames(bioanalyzer_10bp_matrix), "sample")

bioanalyzer_10bp_matrix <- bioanalyzer_10bp_matrix %>%
  mutate(across(all_of(bio_length_cols), as.numeric))

bioanalyzer_10bp_length <- metadata %>%
  select(sample, group) %>%
  inner_join(bioanalyzer_10bp_matrix, by = "sample") %>%
  arrange(group, sample)

cat("\nFeature set: bioanalyzer_10bp_length\n")
print(bioanalyzer_10bp_length %>% count(group))
cat("Features:", ncol(bioanalyzer_10bp_length) - 2, "\n")
cat("Missing values:", sum(is.na(bioanalyzer_10bp_length)), "\n")

# ============================================================
# 8. Check sample consistency
# ============================================================

sample_check <- bind_rows(
  global_endmotif_256 %>%
    select(sample, group) %>%
    mutate(feature_set = "global_endmotif_256"),
  regional_mds_167 %>%
    select(sample, group) %>%
    mutate(feature_set = "regional_mds_167"),
  bioanalyzer_10bp_length %>%
    select(sample, group) %>%
    mutate(feature_set = "bioanalyzer_10bp_length")
) %>%
  count(feature_set, group, name = "n_samples") %>%
  arrange(feature_set, group)

write_tsv(
  sample_check,
  file.path(table_dir, "classification_feature_set_sample_check.tsv")
)

cat("\nSample consistency check:\n")
print(sample_check)

# ============================================================
# 9. Save feature sets
# ============================================================

feature_summaries <- bind_rows(
  write_feature_set(global_endmotif_256, "global_endmotif_256"),
  write_feature_set(regional_mds_167, "regional_mds_167"),
  write_feature_set(bioanalyzer_10bp_length, "bioanalyzer_10bp_length")
)

write_tsv(
  feature_summaries,
  file.path(table_dir, "classification_feature_set_summary.tsv")
)

classification_labels <- metadata %>%
  select(sample, group) %>%
  arrange(group, sample)

write_tsv(
  classification_labels,
  file.path(table_dir, "classification_sample_labels.tsv")
)

# ============================================================
# 10. Final report
# ============================================================

cat("\n============================================================\n")
cat("Classification feature preparation completed.\n")

cat("\nFeature sets saved in:\n")
cat(feature_dir, "\n")

cat("\nSummary table saved to:\n")
cat(file.path(table_dir, "classification_feature_set_summary.tsv"), "\n")

cat("\nFeature set summary:\n")
print(feature_summaries)

cat("============================================================\n")
