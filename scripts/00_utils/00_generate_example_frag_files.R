# ============================================================
# Generate synthetic example fragment files for GitHub
#
# Purpose:
#   Create 10 fictional cfDNA fragment files in .frag.gz format
#   plus a metadata table.
#
# Output:
#   data/example/fragments/CapXX_example.frag.gz
#   data/example/example_sample_metadata.csv
#
# Notes:
#   - Synthetic data only
#   - No real patient/sample information
#   - Fragment length maximum = 500 bp
#   - Columns: chr, start, end, fragment_length, sample
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

set.seed(123)

# ============================================================
# 1. Output paths
# ============================================================

example_dir <- "data/example"
frag_dir <- file.path(example_dir, "fragments")

dir.create(example_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(frag_dir, recursive = TRUE, showWarnings = FALSE)

metadata_file <- file.path(example_dir, "example_sample_metadata.csv")

# ============================================================
# 2. Synthetic sample metadata
# ============================================================

sample_metadata <- tibble(
  sample = sprintf("Cap%02d_example", 1:10),
  group = c(
    "Healthy", "Healthy", "Healthy", "Healthy",
    "Remission", "Remission", "Remission",
    "Relapse", "Relapse", "Relapse"
  )
) %>%
  mutate(
    frag_file = file.path(frag_dir, paste0(sample, ".frag.gz")),
    synthetic = TRUE
  )

# ============================================================
# 3. Chromosome sizes
# Approximate hg19 chromosome lengths
# ============================================================

chrom_sizes <- tibble(
  chr = c(paste0("chr", 1:22), "chrX", "chrY"),
  size = c(
    249250621, 243199373, 198022430, 191154276, 180915260,
    171115067, 159138663, 146364022, 141213431, 135534747,
    135006516, 133851895, 115169878, 107349540, 102531392,
    90354753, 81195210, 78077248, 59128983, 63025520,
    48129895, 51304566, 155270560, 59373566
  )
)

# ============================================================
# 4. Fragment generator
# ============================================================

generate_fragments <- function(sample_id,
                               group_name,
                               n_fragments = 1000,
                               min_len = 50,
                               max_len = 500) {
  
  # Mild group-specific length tendency only for realistic example data
  # All fragments are still capped at 500 bp.
  if (group_name == "Healthy") {
    frag_len <- round(rnorm(n_fragments, mean = 170, sd = 45))
  } else if (group_name == "Remission") {
    frag_len <- round(rnorm(n_fragments, mean = 160, sd = 55))
  } else if (group_name == "Relapse") {
    frag_len <- round(rnorm(n_fragments, mean = 145, sd = 65))
  } else {
    frag_len <- round(rnorm(n_fragments, mean = 160, sd = 55))
  }
  
  frag_len <- pmax(frag_len, min_len)
  frag_len <- pmin(frag_len, max_len)
  
  chr_sample <- sample(
    chrom_sizes$chr,
    size = n_fragments,
    replace = TRUE,
    prob = chrom_sizes$size / sum(chrom_sizes$size)
  )
  
  chr_size <- chrom_sizes$size[match(chr_sample, chrom_sizes$chr)]
  
  start_pos <- mapply(
    function(cs, fl) {
      sample.int(cs - fl, size = 1)
    },
    cs = chr_size,
    fl = frag_len
  )
  
  end_pos <- start_pos + frag_len
  
  frag_df <- tibble(
    chr = chr_sample,
    start = as.integer(start_pos),
    end = as.integer(end_pos),
    fragment_length = as.integer(frag_len),
    sample = sample_id
  ) %>%
    arrange(chr, start, end)
  
  return(frag_df)
}

# ============================================================
# 5. Write .frag.gz files
# ============================================================

for (i in seq_len(nrow(sample_metadata))) {
  
  sample_id <- sample_metadata$sample[i]
  group_name <- sample_metadata$group[i]
  output_file <- sample_metadata$frag_file[i]
  
  frag_df <- generate_fragments(
    sample_id = sample_id,
    group_name = group_name,
    n_fragments = 1000,
    min_len = 50,
    max_len = 500
  )
  
  write_tsv(
    frag_df,
    gzfile(output_file)
  )
  
  cat("Written:", output_file, "\n")
}

# ============================================================
# 6. Write metadata
# ============================================================

sample_metadata <- sample_metadata %>%
  mutate(
    n_fragments = 1000,
    file_format = "frag.gz",
    genome_build = "synthetic_hg19_like",
    note = "Synthetic example data for GitHub testing only"
  )

write_csv(sample_metadata, metadata_file)

cat("\nMetadata written:\n")
cat(metadata_file, "\n")

cat("\nGroup summary:\n")
print(sample_metadata %>% count(group))

cat("\nDone. Synthetic example fragment files created in:\n")
cat(frag_dir, "\n")
