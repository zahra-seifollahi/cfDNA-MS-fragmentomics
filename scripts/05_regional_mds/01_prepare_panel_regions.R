# ============================================================
# Prepare 167 marker regions for regional MDS analysis
#
# Purpose:
#   Read the 167-marker Excel file
#   Extract marker coordinates
#   Save BED file for FinaleToolkit interval-end-motifs
#
# Input:
#   12-07-2024_167_markers.xlsx
#
# Required columns:
#   marker_ID
#   CHR
#   START(hg19_pos)
#   END(hg19_pos)
#   hyper_hypo
#
# Usage:
#   Rscript scripts/05_regional_mds/01_prepare_panel_regions.R /path/to/12-07-2024_167_markers.xlsx
#
# Output:
#   metadata/167_original_markers.bed
#   metadata/ms_panel_regions_unique.bed
#   results/tables/regional_mds/167_markers_clean.tsv
#   results/tables/regional_mds/167_markers_summary.tsv
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(readr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

marker_file <- ifelse(
  length(args) >= 1,
  args[1],
  "metadata/12-07-2024_167_markers.xlsx"
)

bed_out_original <- "metadata/167_original_markers.bed"
bed_out_regional_mds <- "metadata/ms_panel_regions_unique.bed"

table_dir <- "results/tables/regional_mds"

dir.create("metadata", recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(marker_file)) {
  stop("Marker file not found: ", marker_file)
}

markers <- read_excel(marker_file)

required_cols <- c(
  "marker_ID",
  "CHR",
  "START(hg19_pos)",
  "END(hg19_pos)",
  "hyper_hypo"
)

missing_cols <- setdiff(required_cols, colnames(markers))

if (length(missing_cols) > 0) {
  stop(
    "Marker file is missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

markers_clean <- markers %>%
  transmute(
    marker_ID = as.character(marker_ID),
    chr = as.character(CHR),
    start = as.integer(round(as.numeric(`START(hg19_pos)`))),
    end = as.integer(round(as.numeric(`END(hg19_pos)`))),
    hyper_hypo = as.character(hyper_hypo),
    marker_key = paste(chr, start, end, sep = "_")
  ) %>%
  filter(
    !is.na(marker_ID),
    !is.na(chr),
    !is.na(start),
    !is.na(end),
    end > start
  ) %>%
  arrange(chr, start, end)

duplicated_coordinates <- markers_clean %>%
  group_by(chr, start, end) %>%
  summarise(
    n_markers = n_distinct(marker_ID),
    marker_IDs = paste(sort(unique(marker_ID)), collapse = ";"),
    hyper_hypo_values = paste(sort(unique(hyper_hypo)), collapse = ";"),
    .groups = "drop"
  ) %>%
  filter(n_markers > 1) %>%
  arrange(desc(n_markers))

markers_unique <- markers_clean %>%
  group_by(chr, start, end) %>%
  summarise(
    marker_ID = paste(sort(unique(marker_ID)), collapse = ";"),
    hyper_hypo = paste(sort(unique(hyper_hypo)), collapse = ";"),
    .groups = "drop"
  ) %>%
  arrange(chr, start, end)

marker_summary <- data.frame(
  metric = c(
    "original_rows",
    "valid_marker_rows",
    "unique_marker_IDs",
    "unique_coordinates",
    "zero_or_negative_width_removed",
    "duplicated_coordinates"
  ),
  value = c(
    nrow(markers),
    nrow(markers_clean),
    n_distinct(markers_clean$marker_ID),
    nrow(markers_unique),
    nrow(markers) - nrow(markers_clean),
    nrow(duplicated_coordinates)
  )
)

write_tsv(
  markers_clean,
  file.path(table_dir, "167_markers_clean.tsv")
)

write_tsv(
  markers_unique,
  file.path(table_dir, "167_markers_unique_regions.tsv")
)

write_tsv(
  duplicated_coordinates,
  file.path(table_dir, "167_markers_duplicated_coordinates.tsv")
)

write_tsv(
  marker_summary,
  file.path(table_dir, "167_markers_summary.tsv")
)

# BED for original marker rows
write.table(
  markers_clean %>%
    select(chr, start, end, marker_ID),
  file = bed_out_original,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# BED used by the rest of the regional MDS workflow
# This keeps one row per unique coordinate.
write.table(
  markers_unique %>%
    select(chr, start, end, marker_ID),
  file = bed_out_regional_mds,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

cat("\nRegional marker BED preparation completed.\n")

cat("\nInput file:\n")
cat(marker_file, "\n")

cat("\nOriginal rows:\n")
cat(nrow(markers), "\n")

cat("\nValid marker rows:\n")
cat(nrow(markers_clean), "\n")

cat("\nUnique marker IDs:\n")
cat(n_distinct(markers_clean$marker_ID), "\n")

cat("\nUnique BED regions:\n")
cat(nrow(markers_unique), "\n")

cat("\nDuplicated coordinate groups:\n")
cat(nrow(duplicated_coordinates), "\n")

cat("\nBED saved to:\n")
cat(bed_out_regional_mds, "\n")

cat("\nTables saved to:\n")
cat(table_dir, "\n")

cat("\nDone.\n")
