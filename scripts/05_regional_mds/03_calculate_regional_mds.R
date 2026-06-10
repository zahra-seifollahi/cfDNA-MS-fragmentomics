# ============================================================
# Calculate regional MDS matrix from FinaleToolkit
# interval-end-motifs output files
#
# Purpose:
#   Calculate regional motif diversity score (MDS) from
#   4-mer end-motif frequencies for each sample and region.
#
# Input:
#   Directory containing *.interval_endmotifs.tsv files
#
# Usage:
#   Rscript scripts/05_regional_mds/03_calculate_regional_mds.R /path/to/interval_endmotifs
#
# If no input path is provided, default is:
#   results/intermediate/regional_mds/interval_endmotifs
#
# Output:
#   results/tables/regional_mds/
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
})

graphics.off()

# ============================================================
# 1. Settings
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

endmotif_dir <- ifelse(
  length(args) >= 1,
  args[1],
  "results/intermediate/regional_mds/interval_endmotifs"
)

out_dir <- "results/tables/regional_mds"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

k <- 4
min_count <- 30
sum_tolerance <- 1e-4

qc_metadata_file <- "results/tables/qc/final_sample_metadata.csv"

# ============================================================
# 2. Helper functions
# ============================================================

clean_sample_id <- function(file_path) {
  basename(file_path) %>%
    str_remove("\\.interval_endmotifs\\.tsv$") %>%
    str_remove("\\.tsv$") %>%
    str_remove("\\.dedup$")
}

calculate_mds <- function(x, k = 4) {
  x <- as.numeric(x)
  x[is.na(x)] <- 0

  total <- sum(x, na.rm = TRUE)

  if (!is.finite(total) || total <= 0) {
    return(NA_real_)
  }

  p <- x / total
  p <- p[p > 0]

  shannon_entropy <- -sum(p * log(p))
  max_entropy <- log(4^k)

  shannon_entropy / max_entropy
}

process_one_file <- function(file_path, k = 4, min_count = 30) {
  sample_id <- clean_sample_id(file_path)
  file_name <- basename(file_path)

  message("Processing: ", sample_id)

  df <- read_tsv(
    file_path,
    show_col_types = FALSE,
    progress = FALSE
  )

  # FinaleToolkit may use contig instead of seq
  if (!"seq" %in% colnames(df) && "contig" %in% colnames(df)) {
    df <- df %>%
      rename(seq = contig)
  }

  # Some versions may use end/stop; standardize to stop
  if (!"stop" %in% colnames(df) && "end" %in% colnames(df)) {
    df <- df %>%
      rename(stop = end)
  }

  required_region_cols <- c("seq", "start", "stop", "name", "count")
  missing_region_cols <- setdiff(required_region_cols, colnames(df))

  if (length(missing_region_cols) > 0) {
    stop(
      "Missing required columns in ",
      basename(file_path),
      ": ",
      paste(missing_region_cols, collapse = ", ")
    )
  }

  motif_cols <- colnames(df)[str_detect(colnames(df), "^[ACGT]{4}$")]

  if (length(motif_cols) == 0) {
    stop("No 4-mer motif columns found in file: ", basename(file_path))
  }

  if (length(motif_cols) != 4^k) {
    warning(
      "Expected ", 4^k, " motif columns but found ",
      length(motif_cols), " in file: ", basename(file_path)
    )
  }

  df <- df %>%
    mutate(
      across(all_of(motif_cols), as.numeric),
      count = as.numeric(count),
      start = as.integer(start),
      stop = as.integer(stop)
    )

  motif_sum <- rowSums(
    df[, motif_cols, drop = FALSE],
    na.rm = TRUE
  )

  mds_values <- apply(
    df[, motif_cols, drop = FALSE],
    1,
    calculate_mds,
    k = k
  )

  mds_df <- df %>%
    mutate(
      sample_id = sample_id,
      file_name = file_name,
      region_id = ifelse(
        !is.na(name) & name != "",
        as.character(name),
        paste(seq, start, stop, sep = "_")
      ),
      region_coord = paste(seq, start, stop, sep = "_"),
      motif_frequency_sum = motif_sum,
      motif_sum_close_to_1 = abs(motif_frequency_sum - 1) <= sum_tolerance,
      low_count = count < min_count,
      MDS_raw = mds_values,
      MDS = ifelse(low_count, NA_real_, MDS_raw)
    ) %>%
    select(
      sample_id,
      file_name,
      region_id,
      region_coord,
      seq,
      start,
      stop,
      name,
      count,
      motif_frequency_sum,
      motif_sum_close_to_1,
      low_count,
      MDS_raw,
      MDS
    )

  mds_df
}

# ============================================================
# 3. Load QC metadata if available
# ============================================================

qc_metadata <- NULL
included_samples <- NULL

if (file.exists(qc_metadata_file)) {
  qc_metadata <- read_csv(qc_metadata_file, show_col_types = FALSE) %>%
    mutate(sample = as.character(sample))

  if (all(c("sample", "include_analysis") %in% colnames(qc_metadata))) {
    included_samples <- qc_metadata %>%
      filter(include_analysis == "yes") %>%
      pull(sample)

    cat("\nQC metadata found.\n")
    cat("QC-passed samples:", length(included_samples), "\n")
  } else {
    warning("QC metadata file exists but does not contain sample/include_analysis columns.")
  }
} else {
  warning("QC metadata not found. All interval-endmotif files will be used.")
}

# ============================================================
# 4. Process all files
# ============================================================

files <- list.files(
  endmotif_dir,
  pattern = "\\.interval_endmotifs\\.tsv$",
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No .interval_endmotifs.tsv files found in: ", endmotif_dir)
}

file_sample_ids <- vapply(files, clean_sample_id, character(1))

if (!is.null(included_samples)) {
  files <- files[file_sample_ids %in% included_samples]
  file_sample_ids <- file_sample_ids[file_sample_ids %in% included_samples]
}

if (length(files) == 0) {
  stop("No interval-endmotif files remain after QC filtering.")
}

cat("\nNumber of interval-endmotif files used:\n")
print(length(files))

all_mds_long <- map_dfr(
  files,
  process_one_file,
  k = k,
  min_count = min_count
)

# ============================================================
# 5. Sanity checks
# ============================================================

cat("\nNumber of samples:\n")
print(n_distinct(all_mds_long$sample_id))

cat("\nNumber of regions:\n")
print(n_distinct(all_mds_long$region_id))

cat("\nMDS summary after count filtering:\n")
print(summary(all_mds_long$MDS))

cat("\nRaw MDS summary before count filtering:\n")
print(summary(all_mds_long$MDS_raw))

cat("\nNumber of missing MDS values after count filtering:\n")
print(sum(is.na(all_mds_long$MDS)))

cat("\nRows with count < min_count:\n")
print(sum(all_mds_long$low_count, na.rm = TRUE))

cat("\nMotif frequency sum summary:\n")
print(summary(all_mds_long$motif_frequency_sum))

cat("\nRows where motif frequency sum is not close to 1:\n")
print(sum(!all_mds_long$motif_sum_close_to_1, na.rm = TRUE))

region_count_per_sample <- all_mds_long %>%
  count(sample_id, name = "n_regions") %>%
  summarise(
    min_regions = min(n_regions),
    max_regions = max(n_regions),
    mean_regions = mean(n_regions)
  )

cat("\nRegion count per sample:\n")
print(region_count_per_sample)

valid_regions_per_sample <- all_mds_long %>%
  group_by(sample_id) %>%
  summarise(
    n_regions_total = n(),
    n_regions_with_MDS = sum(!is.na(MDS)),
    n_low_count = sum(low_count, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nNon-missing MDS regions per sample summary:\n")
print(
  valid_regions_per_sample %>%
    summarise(
      min_regions_with_MDS = min(n_regions_with_MDS),
      max_regions_with_MDS = max(n_regions_with_MDS),
      mean_regions_with_MDS = mean(n_regions_with_MDS),
      min_low_count = min(n_low_count),
      max_low_count = max(n_low_count),
      mean_low_count = mean(n_low_count)
    )
)

# ============================================================
# 6. Save sanity check tables
# ============================================================

motif_sum_sanity <- all_mds_long %>%
  select(
    sample_id,
    file_name,
    region_id,
    region_coord,
    count,
    motif_frequency_sum,
    motif_sum_close_to_1,
    low_count,
    MDS_raw,
    MDS
  )

motif_sum_suspicious <- motif_sum_sanity %>%
  filter(
    !motif_sum_close_to_1 |
      is.na(motif_frequency_sum) |
      motif_frequency_sum <= 0
  )

write_tsv(
  motif_sum_sanity,
  file.path(out_dir, "regional_mds_motif_frequency_sum_sanity_check.tsv")
)

write_tsv(
  motif_sum_suspicious,
  file.path(out_dir, "regional_mds_motif_frequency_sum_suspicious_rows.tsv")
)

write_tsv(
  valid_regions_per_sample,
  file.path(out_dir, "regional_mds_valid_regions_per_sample.tsv")
)

# ============================================================
# 7. Save long-format MDS table
# ============================================================

write_tsv(
  all_mds_long,
  file.path(out_dir, "regional_MDS_long.tsv")
)

# ============================================================
# 8. Make sample × region MDS matrix
# ============================================================

mds_matrix <- all_mds_long %>%
  select(sample_id, region_id, MDS) %>%
  distinct(sample_id, region_id, .keep_all = TRUE) %>%
  pivot_wider(
    names_from = region_id,
    values_from = MDS
  ) %>%
  arrange(sample_id)

write_tsv(
  mds_matrix,
  file.path(out_dir, "regional_MDS_matrix_samples_by_regions.tsv")
)

write_csv(
  mds_matrix,
  file.path(out_dir, "regional_MDS_matrix_samples_by_regions.csv")
)

mds_matrix_numeric <- as.data.frame(mds_matrix)
rownames(mds_matrix_numeric) <- mds_matrix_numeric$sample_id
mds_matrix_numeric$sample_id <- NULL

# ============================================================
# 9. Matrix sanity checks
# ============================================================

cat("\nMDS matrix dimensions:\n")
print(dim(mds_matrix))

cat("\nFirst columns of MDS matrix:\n")
print(colnames(mds_matrix)[1:min(5, ncol(mds_matrix))])

cat("\nFirst sample IDs in MDS matrix:\n")
print(head(mds_matrix$sample_id))

cat("\nNumeric MDS matrix dimensions:\n")
print(dim(mds_matrix_numeric))

cat("\nTotal NA values in numeric MDS matrix:\n")
print(sum(is.na(mds_matrix_numeric)))

cat("\nNA values per region summary:\n")
print(summary(colSums(is.na(mds_matrix_numeric))))

cat("\nNA values per sample summary:\n")
print(summary(rowSums(is.na(mds_matrix_numeric))))

region_na_summary <- data.frame(
  region_id = colnames(mds_matrix_numeric),
  n_missing = colSums(is.na(mds_matrix_numeric)),
  n_non_missing = colSums(!is.na(mds_matrix_numeric)),
  missing_fraction = colMeans(is.na(mds_matrix_numeric))
) %>%
  arrange(desc(n_missing))

sample_na_summary <- data.frame(
  sample_id = rownames(mds_matrix_numeric),
  n_missing = rowSums(is.na(mds_matrix_numeric)),
  n_non_missing = rowSums(!is.na(mds_matrix_numeric)),
  missing_fraction = rowMeans(is.na(mds_matrix_numeric))
) %>%
  arrange(desc(n_missing))

write_tsv(
  region_na_summary,
  file.path(out_dir, "regional_MDS_region_NA_summary.tsv")
)

write_tsv(
  sample_na_summary,
  file.path(out_dir, "regional_MDS_sample_NA_summary.tsv")
)

# ============================================================
# 10. Save region annotation and sample ID check
# ============================================================

region_annotation <- all_mds_long %>%
  select(region_id, region_coord, seq, start, stop, name) %>%
  distinct() %>%
  arrange(seq, start, stop)

sample_id_check <- all_mds_long %>%
  select(sample_id, file_name) %>%
  distinct() %>%
  arrange(sample_id)

write_tsv(
  region_annotation,
  file.path(out_dir, "regional_MDS_region_annotation.tsv")
)

write_tsv(
  sample_id_check,
  file.path(out_dir, "regional_MDS_sample_id_check.tsv")
)

# ============================================================
# 11. Final message
# ============================================================

cat("\n============================================================\n")
cat("Regional MDS calculation finished.\n")
cat("Input interval-endmotif directory:\n")
cat(endmotif_dir, "\n")

cat("\nOutputs saved in:\n")
cat(out_dir, "\n")

cat("\nMain matrix file:\n")
cat(file.path(out_dir, "regional_MDS_matrix_samples_by_regions.tsv"), "\n")

cat("\nMain long-format file:\n")
cat(file.path(out_dir, "regional_MDS_long.tsv"), "\n")

cat("\nSample ID check file:\n")
cat(file.path(out_dir, "regional_MDS_sample_id_check.tsv"), "\n")

cat("============================================================\n")
