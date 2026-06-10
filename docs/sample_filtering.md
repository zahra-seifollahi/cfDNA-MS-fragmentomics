# Sample filtering and QC

Sample filtering was performed before downstream cfDNA fragmentomic analyses.

The main exclusion criterion was evidence of fragment degradation based on BAM-derived fragment-length profiles. Fragment length distributions were generated from fragment files and compared with Bioanalyzer profiles. Samples with excessive short-fragment enrichment were considered unreliable for end-motif analysis because severe degradation can distort cfDNA fragmentomic signatures.

The length-based QC script calculates sample-level metrics from FinaleToolkit fragment-length output files, including:

- weighted median fragment length below 300 bp
- fraction of fragments shorter than 130 bp
- fraction of fragments in the 120–220 bp mononucleosomal window
- short-to-mononucleosomal fragment ratio

Samples were removed if they met at least one of the following degradation-based criteria:

1. weighted median fragment length below 300 bp < 125 bp
2. fraction of fragments in the 120–220 bp window < 0.50
3. fraction of fragments shorter than 130 bp > 0.50

cfDNA concentration was also checked using Qubit measurements. The concentration values were within the expected range and were therefore used as a QC confirmation rather than as the main exclusion criterion.

The final included sample list should be used consistently in all downstream analyses, including global end-motif analysis, regional motif diversity score analysis, hyper/hypomethylated region analyses, fragment length analysis, and classification.
