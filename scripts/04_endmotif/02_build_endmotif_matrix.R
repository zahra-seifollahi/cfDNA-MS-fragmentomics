# ============================================================
# Build global end-motif matrix
#
# Purpose:
#   Read FinaleToolkit end-motif frequency files and build a
#   sample-by-motif matrix for downstream cfDNA fragmentomics.
#
# Input:
#   FinaleToolkit end-motif output files:
#   Cap01.dedup.endmotif.tsv
#   Cap02.dedup.endmotif.tsv
#   ...
#
# Output:
#   results/intermediate/endmotif/global_endmotif_matrix.tsv
#   results/tables/endmotif/endmotif_sample_summary.tsv
#
# Usage:
#   Rscript scripts/04_endmotif/02_build_endmotif_matrix.R /path/to/endmotif/files
#
# If no path is provided, the script uses:
#   results/intermediate/endmotif
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(stringr)
  library(tibble)
})

# ============================================================
# 1. User settings
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

input_dir <- ifelse(
  length(args) >= 1,
  args[1],
  "results/intermediate/endmotif"
)

matrix_out <- "results/intermediate/endmotif/global_endmotif_matrix.tsv"
sample_summary_out <- "results/tables/endmotif/endmotif_sample_summary.tsv"
qc_metadata_file <- "results/tables/qc/final_sample_metadata.csv"

dir.create(dirname(matrix_out), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(sample_summary_out), showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 2. Helper functions
# ============================================================

clean_sample_name <- function(file_path) {
  basename(file_path) %>%
    str_remove("\\.dedup\\.endmotif\\.tsv$") %>%
    str_remove("\\.endmotif\\.tsv$") %>%
    str_remove("\\.tsv$")
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

read_endmotif_file <- function(file_path) {
  df <- read.delim(
    file_path,
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE
  )

  if (ncol(df) < 2) {
    stop("End-motif file has fewer than 2 columns: ", file_path)
  }

  df <- df[, 1:2]
  colnames(df) <- c("motif", "freq")

  df %>%
    mutate(
      motif = as.character(motif),
      freq = as.numeric(freq)
    ) %>%
    filter(!is.na(motif), !is.na(freq))
}

# ============================================================
# 3. Find end-motif files
# ============================================================

files <- list.files(
  input_dir,
  pattern = "\\.endmotif\\.tsv$",
  full.names = TRUE
)

if (length(files) == 0) {
  stop(
    "No .endmotif.tsv files found in: ",
    input_dir,
    "\nProvide the folder path as an argument, for example:\n",
    "Rscript scripts/04_endmotif/02_build_endmotif_matrix.R /path/to/endmotif/files"
  )
}

files <- sort(files)

cat("\nInput directory:\n")
cat(input_dir, "\n\n")

cat("Number of end-motif files found:", length(files), "\n\n")

# ============================================================
# 4. Read files
# ============================================================

endmotif_list <- map(files, read_endmotif_file)
sample_names <- map_chr(files, clean_sample_name)

names(endmotif_list) <- sample_names

# ============================================================
# 5. Check motif consistency across samples
# ============================================================

reference_motifs <- endmotif_list[[1]]$motif

motif_consistency <- map_lgl(endmotif_list, function(df) {
  identical(df$motif, reference_motifs)
})

if (!all(motif_consistency)) {
  problematic_samples <- names(motif_consistency)[!motif_consistency]

  stop(
    "Motif order or motif set is not identical across all files. Problematic samples: ",
    paste(problematic_samples, collapse = ", "),
    "\nCheck FinaleToolkit outputs before building the matrix."
  )
}

# ============================================================
# 6. Build sample-by-motif matrix
# ============================================================

motif_matrix <- do.call(
  rbind,
  lapply(endmotif_list, function(df) df$freq)
)

colnames(motif_matrix) <- reference_motifs
rownames(motif_matrix) <- sample_names

motif_matrix_df <- as.data.frame(motif_matrix) %>%
  rownames_to_column("sample")

# ============================================================
# 7. Create sample summary
# ============================================================

sample_summary <- tibble(
  sample = sample_names,
  sample_id = as.integer(str_extract(sample_names, "\\d+")),
  group = assign_group_from_sample(sample_names),
  n_motifs = ncol(motif_matrix),
  motif_frequency_sum = rowSums(motif_matrix, na.rm = TRUE)
) %>%
  mutate(
    group = factor(group, levels = c("Healthy", "Remission", "Relapse", "Unknown"))
  ) %>%
  arrange(sample_id)

# ============================================================
# 8. Apply QC sample filtering if metadata exists
# ============================================================

if (file.exists(qc_metadata_file)) {
  cat("QC metadata found:\n")
  cat(qc_metadata_file, "\n\n")

  qc_metadata <- read_csv(qc_metadata_file, show_col_types = FALSE) %>%
    mutate(sample = as.character(sample))

  required_cols <- c("sample", "include_analysis")
  missing_cols <- setdiff(required_cols, colnames(qc_metadata))

  if (length(missing_cols) > 0) {
    stop(
      "QC metadata is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  sample_summary <- sample_summary %>%
    left_join(
      qc_metadata %>%
        select(sample, include_analysis, exclusion_reason),
      by = "sample"
    )

  missing_qc <- sample_summary %>%
    filter(is.na(include_analysis))

  if (nrow(missing_qc) > 0) {
    warning(
      "Some end-motif samples were not found in QC metadata: ",
      paste(missing_qc$sample, collapse = ", ")
    )
  }

  keep_samples <- sample_summary %>%
    filter(include_analysis == "yes") %>%
    pull(sample)

  motif_matrix_df <- motif_matrix_df %>%
    filter(sample %in% keep_samples)

  sample_summary <- sample_summary %>%
    filter(sample %in% keep_samples)

  cat("Applied QC filtering.\n")
  cat("Samples kept after QC:", length(keep_samples), "\n\n")

} else {
  cat("No QC metadata file found. All end-motif samples will be kept.\n\n")

  sample_summary <- sample_summary %>%
    mutate(
      include_analysis = "yes",
      exclusion_reason = NA_character_
    )
}

# ============================================================
# 9. Final checks
# ============================================================

cat("Final sample counts by group:\n")
print(sample_summary %>% count(group))

cat("\nMotif matrix dimensions:\n")
cat("Samples:", nrow(motif_matrix_df), "\n")
cat("Motifs:", ncol(motif_matrix_df) - 1, "\n\n")

cat("Range of motif frequencies:\n")
print(range(as.matrix(motif_matrix_df[, -1]), na.rm = TRUE))

cat("\nSummary of row sums:\n")
print(
  motif_matrix_df %>%
    select(-sample) %>%
    as.matrix() %>%
    rowSums(na.rm = TRUE) %>%
    summary()
)

# ============================================================
# 10. Save outputs
# ============================================================

write_tsv(motif_matrix_df, matrix_out)
write_tsv(sample_summary, sample_summary_out)

cat("\nSaved motif matrix to:\n")
cat(matrix_out, "\n")

cat("\nSaved sample summary to:\n")
cat(sample_summary_out, "\n")

cat("\nDone.\n")
