# FinaleToolkit notes

This project used FinaleToolkit for cfDNA fragmentomic feature extraction, including global and region-based end-motif analyses.

The input files were compressed and indexed fragment files generated from deduplicated BAM files.

Because the sequencing data were targeted rather than whole-genome, some genomic regions or contigs may be absent from individual fragment files. During the original analysis, FinaleToolkit behavior was adjusted to skip missing contigs instead of stopping with a `ValueError`.

For reproducibility, this repository does not include a modified conda environment or edited package files. Users should either:

1. use a FinaleToolkit version that safely skips missing contigs, or
2. pre-filter analysis regions to contigs present in each fragment file, or
3. apply an equivalent local patch and document the exact change.

The original adjustment was conceptually equivalent to returning/skipping when `tbx.fetch(contig, start, stop)` raises a `ValueError`.
