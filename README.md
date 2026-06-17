Last login: Wed Jun 17 15:30:56 on ttys000
zahraseyfollahi@Zahras-MacBook ~ % nano downloads/README.md
zahraseyfollahi@Zahras-MacBook ~ % /Users/zahraseyfollahi/Downloads/README_updated.md 
zsh: permission denied: /Users/zahraseyfollahi/Downloads/README_updated.md
zahraseyfollahi@Zahras-MacBook ~ % nano /Users/zahraseyfollahi/Downloads/README_updated.md






























































  UW PICO 5.09                        File: /Users/zahraseyfollahi/Downloads/README_updated.md                          

# cfDNA-MS-fragmentomics

This repository contains the analysis workflow for a cell-free DNA (cfDNA) fragmentomics study in multiple sclerosis (M$

The workflow includes sample quality control, cfDNA concentration analysis, fragment length analysis, global end-motif $

Raw sequencing files, private metadata, original BED files, intermediate outputs, and generated results are **not inclu$

---

## Project overview

The main goal of this project is to evaluate whether cfDNA fragmentomic features can distinguish healthy individuals fr$

The analyzed feature categories include:

1. **cfDNA concentration**
   - Qubit-based cfDNA concentration comparison across groups.

2. **Fragment length features**
   - Fragment length distributions from BAM-derived fragment files.
   - Bioanalyzer-derived fragment length distributions.
   - Bioanalyzer 10-bp length bins used as one classification feature set.

3. **Global end-motif frequencies**
   - 4-mer end-motif frequencies extracted using FinaleToolkit.
   - 256 global end-motif features per sample.
   - Motif diversity score calculated from global end-motif distributions.

4. **Regional motif diversity**
   - End-motif profiles calculated across selected genomic regions.
   - Regional motif diversity score matrix generated for downstream analysis and classification.

5. **Hyper/hypomethylated region analysis**
   - End-motif analysis repeated in selected hypermethylated and hypomethylated regions.
   - Comparison of methylation-associated end-motif patterns.

6. **Machine-learning classification**
   - Feature sets:
     - global end-motif frequencies
     - regional motif diversity scores
     - Bioanalyzer 10-bp length bins
   - Classifiers:
     - Support Vector Machine
     - Random Forest
     - Elastic Net
   - Evaluation:
     - nested cross-validation
     - multiple binary classification tasks
     - classifier comparison
     - consensus feature ranking
     - Random Forest using consensus-selected features

---

## Repository structure

```text
cfDNA-MS-fragmentomics/
├── data/
│   └── example/
│       ├── README.md

