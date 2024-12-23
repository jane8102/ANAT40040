# [24234227] ANAT40040 Assignment 2: Gene Expression Analysis and Interpretation (2024)

### Overview
This repository contains code to perform a comprehensive analysis of RNA sequencing (RNA-seq) data to identify differentially expressed (DE) genes and evaluate their prognostic significance in breast cancer (BC). The analysis includes data extraction, preprocessing, differential expression analysis, pathway enrichment analyses, and survival modelling using Lasso Regularised Cox Regression. The key outputs include pathway enrichment plots, Principal Component Analysis (PCA) plots, heatmaps, Kaplan-Meier survival curves, and hazard ratios.

### Code Components
1. **Data Extraction**
- Extracts RNA-seq data, clinical patient, and copy number aberration (CNA) data from a .tar.gz archive.
- Verifies extracted files and cleans column names for further analysis.


2. **Data Preprocessing**
- Matches and subsets RNA-seq, clinical patient, and CNA data based on shared patient IDs.
- Creates a metadata table including ERBB2 amplification status and cancer stage.

3. **Differential Expression Analysis**
- Uses the DESeq2 package to identify DE genes between ERBB2-amplified and ERBB2-not amplified samples.
- Filters significant DE genes with absolute _log2_ fold change values and adjusted _p_-value < 0.05.
- Outputs counts of upregulated and downregulated genes.

4. **Pathway Enrichment Analyses**
- Utilises clusterProfiler, ReactomePA, and pathview packages for Gene Ontology (GO), KEGG, and Reactome pathway enrichment analyses.
- Generates dot plots and tree plots for visualising top pathways for over- and under-expressed genes.

5. **Data Visualisation**
- Plots histograms and bar plots to analyse the distributions of clinical patient data.
- Produces PCA plots to explore variance in expression profiles by ERBB2 amplification status and cancer stage.
- Creates a heatmap of the top 20 DE genes ranked by adjusted _p_-value to visualise the expression patterns.

6. **Survival Analysis Using Lasso Regularised Cox Regression**
- Implements Lasso Regularised Cox Regression using glmnet package on the top 200 DE genes to identify prognostic genes.
- Performs 5-fold cross-validation to determine the optimal penalty parameter λ.
- Plots the Lasso path, cross-validation plot.
- Outputs the significant DE genes selected at the optimal λ.
- Plots Kaplan-Meier survival curves to visualise the overal survival (OS) probabilities for high- and low-expression groups of significant genes (genes with _p_-values < 0.05 are considered statistically significant).
- Calculates hazard ratios (HR) for significant genes using Cox proportional hazards models and displays them as forest plots for easy interpretation.

### Dependencies
Ensure the following R packages are installed:
- Core Analysis: DESeq2, glmnet, survival
- Pathway Enrichment: clusterProfiler, ReactomePA, pathview, org.Hs.eg.db
- Data Visualisation: ggplot2, ggpubr, pheatmap, survminer
- Utilities: BiocManager, gridExtra

### How to Run?
1. **Prepare Data**
- Insert the .tar.gz file in the directory specified in the script.
- Ensure file paths in the script are correctly set.

2. **Install Required Packages**
- Run the provided package installation checks to ensure all dependencies are installed.

3. **Execute Analysis**
- Run the codes sequentially to extract, preprocess, and analyse data.
- Outputs will include visualisations (plots), DE gene lists, and pathway enrichment results.

4. **Visualise Outputs**
- Observe the histograms and bar plots to analyse the distribution of clinical patient data.
- Examine pathway enrichment dot plots and tree plots.
- Review Kaplan-Meier survival curves and hazard ratio forest plots.

### Contact
For questions or issues, contact: Chek, Zi Yan Jane \
Email: zi.chek@ucdconnect.ie

