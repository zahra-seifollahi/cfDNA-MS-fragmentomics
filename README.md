# cfDNA-MS-fragmentomics

This repository contains the analysis workflow for a cell-free DNA (cfDNA) fragmentomics study in multiple sclerosis (MS). The project analyzes plasma cfDNA features across three groups: Healthy, Remission, and Relapse.

The workflow includes sample quality control, cfDNA concentration analysis, fragment length analysis, global end-motif analysis, regional motif diversity score analysis, hyper/hypomethylated region analysis, and machine-learning classification.

Raw sequencing files, private metadata, BED files, intermediate outputs, and generated results are not included in this public repository.

---

## Project overview

The main goal of this project is to investigate whether cfDNA fragmentomic features can distinguish between healthy individuals and MS-related clinical groups.

The analyzed feature categories include:

1. **cfDNA concentration**

   * Qubit-based cfDNA concentration comparison across groups.

2. **Fragment length features**

   * Fragment length analysis from BAM-derived fragment files.
   * Bioanalyzer-derived fragment length distribution.
   * Bioanalyzer 10-bp length bins used for classification.

3. **Global end-motif frequencies**

   * 4-mer end-motif frequencies extracted using FinaleToolkit.
   * 256 global end-motif features per sample.
   * Motif diversity score calculated from end-motif distributions.

4. **Regional motif diversity**

   * End-motif profiles calculated in selected genomic regions.
   * Regional motif diversity score matrix generated for downstream analysis and classification.

5. **Hyper/hypomethylated region analysis**

   * End-motif analysis repeated in selected hypermethylated and hypomethylated regions.

6. **Classification**

   * Feature sets:

     * global end-motif frequencies
     * regional motif diversity scores
     * Bioanalyzer 10-bp length bins
   * Classifiers:

     * Support Vector Machine
     * Random Forest
     * Elastic Net
   * Evaluation:

     * nested cross-validation
     * multiple binary classification tasks
     * classifier comparison
     * consensus feature ranking
     * Random Forest using consensus-selected features

---

## Repository structure

```text
cfDNA-MS-fragmentomics/
├── docs/
│   ├── finaletoolkit_notes.md
│   └── sample_filtering.md
│
├── metadata/
│   └── .gitkeep
│
├── scripts/
│   ├── 01_qc/
│   │   ├── 01_make_fragment_files_from_bam.sh
│   │   └── 02_length_based_sample_qc.R
│   │
│   ├── 02_concentration/
│   │   └── 01_concentration_analysis.R
│   │
│   ├── 03_fragment_length/
│   │   ├── 01_bam_fragment_length_analysis.R
│   │   └── 02_bioanalyzer_fragment_length_analysis.R
│   │
│   ├── 04_endmotif/
│   │   ├── 01_run_global_endmotifs.sh
│   │   ├── 02_build_endmotif_matrix.R
│   │   ├── 03_global_mds_analysis.R
│   │   └── 04_global_endmotif_complete_analysis.R
│   │
│   ├── 05_regional_mds/
│   │   ├── 01_prepare_panel_regions.R
│   │   ├── 02_run_regional_endmotifs.sh
│   │   ├── 03_calculate_regional_mds.R
│   │   └── 04_analyze_regional_mds.R
│   │
│   ├── 06_hyper_hypo_regions/
│   │   ├── 01_select_hyper_hypo_regions.sh
│   │   ├── 02_run_hyper_hypo_endmotifs.sh
│   │   └── 03_analyze_hyper_hypo_endmotifs.R
│   │
│   └── 07_classification/
│       ├── 01_prepare_classification_features.R
│       ├── 02_run_svm_nested_cv.R
│       ├── 03_run_rf_nested_cv.R
│       ├── 04_run_elastic_net_nested_cv.R
│       ├── 05_compare_classifiers.R
│       ├── 06_create_consensus_feature_lists.R
│       └── 07_run_consensus_top50_rf.R
│
├── .gitignore
└── README.md
```

---

## Input data

This repository does not include raw or private data.

Expected input files include:

* BAM files or paired-end alignment files for fragment generation
* fragment files generated from BAM files
* Bioanalyzer interval data
* cfDNA concentration table
* sample metadata
* genomic region BED files
* FinaleToolkit end-motif output files
* reference genome file for FinaleToolkit, such as a `.2bit` file

BED files used for regional analysis and hyper/hypomethylated region analysis should be provided locally by the user and should not be committed to the public repository.

---

## Software requirements

The workflow uses:

* R
* Bash
* FinaleToolkit
* command-line tools for file handling and genomic interval processing

Main R packages used across scripts include:

```r
readr
readxl
dplyr
tidyr
stringr
ggplot2
pROC
caret
randomForest
glmnet
limma
clinfun
scales
forcats
```

FinaleToolkit is used for extracting cfDNA end-motif features and regional motif profiles.

---

## Workflow

### 1. Quality control

```bash
bash scripts/01_qc/01_make_fragment_files_from_bam.sh
Rscript scripts/01_qc/02_length_based_sample_qc.R
```

This step generates fragment files and identifies QC-passed samples for downstream analysis.

---

### 2. cfDNA concentration analysis

```bash
Rscript scripts/02_concentration/01_concentration_analysis.R
```

This step compares cfDNA concentration across Healthy, Remission, and Relapse groups.

---

### 3. Fragment length analysis

```bash
Rscript scripts/03_fragment_length/01_bam_fragment_length_analysis.R
Rscript scripts/03_fragment_length/02_bioanalyzer_fragment_length_analysis.R
```

The Bioanalyzer-derived 10-bp length bins are used later as one of the classification feature sets.

---

### 4. Global end-motif analysis

```bash
bash scripts/04_endmotif/01_run_global_endmotifs.sh
Rscript scripts/04_endmotif/02_build_endmotif_matrix.R
Rscript scripts/04_endmotif/03_global_mds_analysis.R
Rscript scripts/04_endmotif/04_global_endmotif_complete_analysis.R
```

This step builds the 256-feature global end-motif matrix and performs global motif diversity and differential motif analyses.

---

### 5. Regional motif diversity analysis

```bash
Rscript scripts/05_regional_mds/01_prepare_panel_regions.R
bash scripts/05_regional_mds/02_run_regional_endmotifs.sh
Rscript scripts/05_regional_mds/03_calculate_regional_mds.R
Rscript scripts/05_regional_mds/04_analyze_regional_mds.R
```

This step calculates motif diversity scores across selected genomic regions.

---

### 6. Hyper/hypomethylated region analysis

```bash
bash scripts/06_hyper_hypo_regions/01_select_hyper_hypo_regions.sh
bash scripts/06_hyper_hypo_regions/02_run_hyper_hypo_endmotifs.sh
Rscript scripts/06_hyper_hypo_regions/03_analyze_hyper_hypo_endmotifs.R
```

This step evaluates end-motif patterns in selected hypermethylated and hypomethylated regions.

---

### 7. Classification analysis

First prepare the standardized classification feature sets:

```bash
Rscript scripts/07_classification/01_prepare_classification_features.R
```

This creates three standardized feature matrices:

```text
global_endmotif_256
regional_mds_167
bioanalyzer_10bp_length
```

Run SVM:

```bash
Rscript scripts/07_classification/02_run_svm_nested_cv.R global_endmotif_256
Rscript scripts/07_classification/02_run_svm_nested_cv.R regional_mds_167
Rscript scripts/07_classification/02_run_svm_nested_cv.R bioanalyzer_10bp_length
```

Run Random Forest:

```bash
Rscript scripts/07_classification/03_run_rf_nested_cv.R global_endmotif_256
Rscript scripts/07_classification/03_run_rf_nested_cv.R regional_mds_167
Rscript scripts/07_classification/03_run_rf_nested_cv.R bioanalyzer_10bp_length
```

Run Elastic Net:

```bash
Rscript scripts/07_classification/04_run_elastic_net_nested_cv.R global_endmotif_256
Rscript scripts/07_classification/04_run_elastic_net_nested_cv.R regional_mds_167
Rscript scripts/07_classification/04_run_elastic_net_nested_cv.R bioanalyzer_10bp_length
```

Compare classifiers:

```bash
Rscript scripts/07_classification/05_compare_classifiers.R
```

Create consensus feature lists:

```bash
Rscript scripts/07_classification/06_create_consensus_feature_lists.R
```

Run Random Forest using consensus-selected Top50 features from all feature sets:

```bash
Rscript scripts/07_classification/07_run_consensus_top50_rf.R
```

---

## Classification tasks

The classification scripts evaluate the following binary tasks:

```text
Healthy vs MS
Healthy vs Remission
Healthy vs Relapse
NonRelapse vs Relapse
Remission vs Relapse
```

For SVM and Random Forest, multiple feature-selection methods are evaluated:

```text
RFE
RF importance
Limma FDR<0.05
Top50 Limma
Top50 trend
```

Elastic Net performs embedded feature selection through regularization.

---

## Cross-validation design

The classification workflow uses nested cross-validation:

* outer 5-fold cross-validation for unbiased test evaluation
* inner 3-fold cross-validation for model tuning, feature selection, and threshold selection

Data imputation, filtering, feature selection, model tuning, and threshold selection are performed within cross-validation folds to reduce data leakage.

---

## Output organization

Generated results are written under:

```text
results/tables/
results/figures/
results/intermediate/
```

These folders are excluded from version control by default.

---

## Notes on public repository use

This repository is intended to share analysis code and workflow structure, not raw data or private metadata.

Before running the scripts, users should adapt:

* input file locations
* sample metadata paths
* BED file paths
* reference genome paths
* FinaleToolkit reference `.2bit` path
* output directory permissions

Private paths, raw sequencing data, generated results, and BED files should not be committed.
## R environment

The R package environment is documented in `renv.lock`. To restore the environment, open R in the project directory and run:

```r
install.packages("renv")
renv::restore()
```

---

## Author

Zahra Seifollahi
MSc thesis project
School of Biotechnology, University of Tehran


