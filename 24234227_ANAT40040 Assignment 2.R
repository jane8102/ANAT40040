#------------------------------------------
# title: "24234227_ANAT40040 Assignment 2"
# author: "Chek, Zi Yan Jane"
#------------------------------------------
  
########## Data Extraction
  
# Define the path to the tar.gz file.
tar_file  = "C:/Users/chekj/Documents/UCD/Academic/MSc AI in Medicine & Medical Research/Autumn/ANAT40040 Biological Principles & Cellular Organisation/Assignments/Gene Expression Analysis and Interpretation/brca_tcga_pan_can_atlas_2018.tar.gz"


# Define the destination folder where files will be extracted.
output_dir = "C:/Users/chekj/Documents/UCD/Academic/MSc AI in Medicine & Medical Research/Autumn/ANAT40040 Biological Principles & Cellular Organisation/Assignments/Gene Expression Analysis and Interpretation"


# Check if the output directory exists, create it if necessary.
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE) # Create directory if not exists
  cat("Output directory created at:", output_dir, "\n")
} else {
  cat("Output directory already exists:", output_dir, "\n")
}


# Use untar() to extract the tar.gz achieve.
untar(tarfile = tar_file, exdir = output_dir)
cat("Extraction completed successfully.\n")


# Display the files extracted for verification.
extracted_files = list.files(output_dir, recursive = TRUE) # List all files extracted
print(extracted_files) # Print files that have been extracted


# Define the folder path relative to the current working directory.
folder_path = paste(getwd(),"brca_tcga_pan_can_atlas_2018", sep = "/" )
cat("The folder path is:", folder_path, "\n") # Print folder path to verify


# Define file paths for:
# - RNA-seq file: "data_mrna_seq_v2_rsem.txt"
# - Patient Data: "data_clinical_patient.txt"
# - Copy Number Aberrations Data (CNA): "data_cna.txt"

# Define file paths
rnaSeq_path = paste(folder_path,"data_mrna_seq_v2_rsem.txt", sep = "/")
clinPatient_path = paste(folder_path, "data_clinical_patient.txt", sep = "/")
cna_path = paste(folder_path, "data_cna.txt", sep = "/")


# Read files using read.delim().

# Read rna-seq data
rnaSeq_data = read.delim(rnaSeq_path) 

# Read patient data
clinPatient_data = read.delim(clinPatient_path) 
clinPatient_data = clinPatient_data[5:dim(clinPatient_data)[1],] # Skip 5 rows of column descriptions
rownames(clinPatient_data) = NULL # Reset row indices

# Read cna data
cna_data = read.delim(cna_path) 


# Define a function to replace '.' with '-' in column names.
clean_colNames = function(df) {
  colnames(df) = gsub("\\.", "-", colnames(df))  # Replace '.' with '-'
  return(df)
}


# Apply clean_colNames function to RNA-seq and CNA data.
rnaSeq_data = clean_colNames(rnaSeq_data)
cna_data = clean_colNames(cna_data)


# Display each dataframe.

# Display rna-seq data
cat("RNA-seq data:\n") 
print(rnaSeq_data)

# Display patient data
cat("Clinical data:\n")
print(clinPatient_data)

# Display cna data
cat("CNA data:\n")
print(cna_data)


########## Data Preprocessing

# Extract Patient IDs.

# Extract patient IDs from rna-seq data
rnaSeq_ids = colnames(rnaSeq_data)[3:ncol(rnaSeq_data)] # Patient IDs are in the column headers starting from columns 3 to end
rnaSeq_cIds = gsub("-[0-9]+$", "", rnaSeq_ids) # Remove any suffix after the last dash (e.g., -01)
rnaSeq_mapping = setNames(rnaSeq_ids, rnaSeq_cIds) # Create a mapping of cleaned IDs (xxxx_cIds) to original IDs (xxxx_oIDs)

# Extract patient IDs from cna_data
cna_ids = colnames(cna_data)[3:ncol(cna_data)] # patient IDs are in the column headers starting from columns 3 to end
cna_cIds = gsub("-[0-9]+$", "", cna_ids) # Remove any suffix after the last dash (e.g., -01)
cna_mapping = setNames(cna_ids, cna_cIds) # Create a mapping of cleaned IDs (xxxx_cIds) to original IDs

# Extract patient IDs from clinPatient_data
clinPatient_ids = clinPatient_data[, 1] # the patient IDs are stored in column 1

# Check the number of IDs
cat("Number of RNA-seq Patient IDs:", length(rnaSeq_cIds), "\n")
cat("Number of CNA Patient IDs:", length(cna_cIds), "\n")
cat("Number of Patient Data IDs:", length(clinPatient_ids), "\n")


# Match RNA-seq patient IDs with the CNA IDs and patient data IDs.

# Find common patient IDs
common_ids = Reduce(intersect, list(rnaSeq_cIds, cna_cIds, clinPatient_ids))

cat("First 6 common IDs: ", head(common_ids), "\n")
cat("Number of matching patient IDs across all datasets:", length(common_ids), "\n")


# Subset the data based on common IDs.

# Get original column names that match the cleaned common IDs
rnaSeq_oIds = rnaSeq_mapping[common_ids]
cna_oIds = cna_mapping[common_ids]

# Subset rna-seq data (columns with matching original IDs)
rnaSeq_mData = rnaSeq_data[, c(1, 2, which(colnames(rnaSeq_data) %in% rnaSeq_oIds))]

# Subset cna data (columns with matching original IDs)
cna_mData = cna_data[, c(1, 2, which(colnames(cna_data) %in% cna_oIds))]

# Subset patient data (rows with matching original IDs)
clinPatient_mData = clinPatient_data[clinPatient_data[, 1] %in% common_ids, ]

cat("Dimensions of RNA-seq subset:", dim(rnaSeq_mData), "\n")
cat("Dimensions of CNA subset:", dim(cna_mData), "\n")
cat("Dimensions of Patient Data subset:", dim(clinPatient_mData), "\n")


########## Data Analysis

# Descriptive statistics of the patient data.

# Set up the plotting area to display 3 plots in a single row
par(mfrow = c(1,3))

# Plot a histogram for age
col_age = which(colnames(clinPatient_mData)=='Diagnosis.Age')

hist(as.numeric(clinPatient_mData[, col_age]), # Plot histogram
     main = "Histogram of Age at Diagnosis", 
     xlab = "Age at Diagnosis",
     breaks = 50,
     col = "darkturquoise")
mtext("A", side = 3, line = 1, adj = 0, cex = 1.5)

# Plot a bar plot for stage
col_stage = which(colnames(clinPatient_mData) == "Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code")

# Remove subcategories A, B, C for simplicity
stage = clinPatient_mData[, col_stage]
stage = gsub("I[^VI]", "I", stage)
stage = gsub("V[^\\s]", "V", stage)

bar_positions = barplot(table(as.factor(stage)),
                        main = "Barplot of Stage",
                        ylim = c(0, max(table(as.factor(stage))) + 50), col = "salmon") # Plot barplot
mtext("B", side = 3, line = 1, adj = 0, cex = 1.5)
text(x = bar_positions, y = table(as.factor(stage)), label = table(as.factor(stage)), pos = 3, cex = 0.8, col = "black")

# Histogram of (uncensored) survival times 
col_osStatus = which(colnames(clinPatient_mData) == "Overall.Survival.Status")
uncensored = which(clinPatient_mData[, col_osStatus] == "1:DECEASED")

col_osMonths = which(colnames(clinPatient_mData) == "Overall.Survival..Months.")

hist(as.numeric(clinPatient_mData[uncensored, col_osMonths]), # Plot histogram
     main = "Histogram of (Uncensored) Overall Survival (OS) Months", 
     xlab = "Survival (Months)",
     breaks = 50,
     col = "aquamarine2")
mtext("C", side = 3, line = 1, adj = 0, cex = 1.5)

# Reset to default plotting parameters to avoid affecting future plots
par(mfrow = c(1,1))


# Create a metadata using the CNA level of ERBB2+, where:
# - CNA level of ERBB2+ greater than 0 means amplified
# - CNA level of ERBB2+ less than or equal to 0 means not amplified

# Convert rna-seq data to assay matrix
# Columns 1 and 2 are gene names
assay = round(as.matrix(rnaSeq_mData[, -c(1,2)]))
rownames(assay) = rnaSeq_mData[,1] # Set gene names as rownames

# Extract ERBB2 CNA levels
erbb2_row = which(cna_mData[, 1] == "ERBB2") # Find the row corresponding to ERBB2 in cna data
erbb2_cna = as.numeric(cna_mData[erbb2_row, -c(1, 2)]) # Extract cna values and skip first two columns

pat_ids = clinPatient_mData[, 1]  # Get patient IDs from patient data
metadata = matrix(0, ncol(assay), 2)  # Initialise metadata matrix

# For loop through the columns of the assay to create metadata
for (i in 1:ncol(assay)) {
  pat_barcode = colnames(assay)[i]
  pat_barcode = substr(pat_barcode, 1, 12) # Keep first 12 characters of the barcode
  pat_barcode = gsub("-[0-9]+$", "", pat_barcode) # Remove any suffix after the last dash (e.g., -01)
  
  # Match the barcode to the patient IDs
  idx = which(pat_barcode == pat_ids)
  
  # Extract the ERBB2 cna value for this column
  cna_status = erbb2_cna[i]
  
  # Assign "amplified" or "not amplified" based on CNA level
  metadata[i, 1] = ifelse(cna_status > 0, "Amplified", "Not Amplified")
  
  # Add stage factor from patient data
  metadata[i, 2] = ifelse(length(idx) > 0, as.character(stage[idx]), NA)
}


# Replace NAs in metadata with default values
metadata[is.na(metadata)] = "NA"

# Set column names for metadata
colnames(metadata) = c("ERBB2_Amplification", "Stage")

# Convert metadata to data frame
metadata = as.data.frame(metadata)

print(head(metadata))
cat("Dimension of Metadata: ", dim(metadata))


# Install Deseq2 using Bioconductor.

# Check if the BiocManager package is installed
if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager") # Install BiocManager if not already installed
}

# Check if the DESeq2 package is installed
if (!require("DESeq2", quietly = TRUE)){
  BiocManager::install("DESeq2")  # Install DESeq2 if not already installed 
}

library(DESeq2) # Load DESeq2 library


# Create an DESeq object and normalise the data.

assay[is.na(assay)] = 0  # Impute NAs with zeros
assay[assay < 0] = 0 # Replace negative values with 0


# Filter out genes with low expression
smallestGroupSize = 3
keep = rowSums(assay >= 10) >= smallestGroupSize # Keep genes expressed in at least 'smallestGroupSize' samples
assay = assay[keep,]

# Relevel the factor to make "Not Amplified" the reference level
metadata$ERBB2_Amplification = factor(metadata$ERBB2_Amplification, level = c("Not Amplified", "Amplified"))

# Create DESeq2 dataset
dds =  DESeqDataSetFromMatrix(
  countData = assay,             # Assay matrix with integer counts
  colData = metadata,            # Metadata with amplification information
  design = ~ ERBB2_Amplification # Design formula for differential expression analysis
)


dds = DESeq(dds)  # Run the DESeq pipeline


### Perform Differential Expression Analysis

# Extract the results after running the DESeq pipeline and obtain differentially expressed (DE) genes ranked by fold change.

resultsNames(dds) # Lists the coefficients
cat("\n")

# Obtain the results of the differential expression analysis
res = results(dds)


# Calculate total DE genes with adjusted p-value<0.05
total_pvalue_DE = nrow(res[res$padj < 0.05, ])

# Calculate upregulated genes (positive log2 fold change and adjusted p-value<0.05)
upregulated_genes = nrow(res[res$log2FoldChange > 0 & res$padj < 0.05, ])

# Calculate downregulated genes (negative log2 fold change and adjusted p-value<0.05)
downregulated_genes = nrow(res[res$log2FoldChange < 0 & res$padj < 0.05, ])

# Print results
cat("============================================================================ \n")
cat("Differentially Expressed Genes (DEGs) Summary \n")
cat("============================================================================ \n")
cat("Total Differentially Expressed (DE) Genes (adjusted p-value < 0.05):", total_pvalue_DE, "\n")
cat("Upregulated Genes (log2FoldChange > 0 & adjusted p-value < 0.05):", upregulated_genes, "\n")
cat("Downregulated Genes (log2FoldChange < 0 & adjusted p-value < 0.05):", downregulated_genes, "\n")

# Rank genes by log2FoldChange (fold change)
# Sorts the genes by the magnitude of their fold change (regardless of direction)
# Order by absolute log2FoldChange to get the most differentially expressed genes
rank_foldChange = res[order(abs(res$log2FoldChange), decreasing = TRUE), ] 
cat("============================================================================ \n")
cat("Top 10 Differentially Expressed (DE) Genes Ranked by Absolute Fold Change \n")
cat("============================================================================ \n")
head(rank_foldChange, 10) # Print top 10 most differentially expressed genes by fold change
cat("\n")


### Perform Pathway Enrichment Analyses

# Perform Gene Ontology (GO) Enrichment Analysis on over- and under-expressed genes.

# Check if the clusterProfiler package is installed
if (!requireNamespace("clusterProfiler", quietly = TRUE))
  BiocManager::install("clusterProfiler") # Install clusterProfiler if not already installed

# Check if the org.Hs.eg.db package is installed
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
  BiocManager::install("org.Hs.eg.db") # Install org.Hs.eg.db if not already installed

# Check if the enrichplot package is installed
if (!requireNamespace("enrichplot", quietly = TRUE))
  install.packages("enrichplot") # Install enrichplot if not already installed 

# Check if the pathview package is installed
if (!requireNamespace("pathview", quietly = TRUE))
  BiocManager::install("pathview") # Install pathview if not already installed

# Check if the ReactomePA package is installed
if (!requireNamespace("ReactomePA", quietly = TRUE))
  BiocManager::install("ReactomePA") # Install ReactomePA if not already installed

# Load the required libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(pathview)
library(ReactomePA)


# Check if ggplot2 package is installed
if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages("ggplot2") # Install ggplot2 if not already installed

# Check if ggpubr package is installed
if (!requireNamespace("ggpubr", quietly = TRUE)) {
  install.packages("ggpubr") # Install ggpubr if not already installed
}

# Load the libraries
library(ggplot2)
library(ggpubr)


# Obtain the subset of differentially expressed genes that are statistically significant
# Filtering genes with an adjusted p-value (padj) < 0.05
res_signf = res[res$padj < 0.05, ]

# Separate the significant genes into over-expressed and under-expressed categories
# Over-expressed: Genes with positive log2 fold change
# Under-expressed: Genes with negative log2 fold change
over_DE = rownames(res_signf[res_signf$log2FoldChange > 0, ])
under_DE = rownames(res_signf[res_signf$log2FoldChange < 0, ])


# Perform GO enrichment analysis for under-expressed genes
go_results_under = enrichGO(
  gene          = under_DE,      # List of under-expressed gene symbols
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",      # Input gene IDs are gene symbols
  ont           = "BP",          # Ontology type is biological process
  pAdjustMethod = "BH",          # Adjust p-values using Benjamini-Hochberg
  pvalueCutoff  = 0.05,          # p-value cutoff for significance
  qvalueCutoff  = 0.05           # q-value cutoff for significance
)

# Perform GO enrichment analysis for over-expressed genes
go_results_over = enrichGO(
  gene          = over_DE,        # List of over-expressed gene symbols
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",       # Input gene IDs are gene symbols
  ont           = "BP",           # Ontology type is biological process
  pAdjustMethod = "BH",           # Adjust p-values using Benjamini-Hochberg
  pvalueCutoff  = 0.05,           # p-value cutoff for significance
  qvalueCutoff  = 0.05            # q-value cutoff for significance
)

# Print the top results for GO enrichment of under-expressed genes
print(head(go_results_under))

# Print the top results for GO enrichment of over-expressed genes
print(head(go_results_over))

# Plot dotplot for the top 10 GO terms for under-expressed genes
dotplotGO_under = dotplot(go_results_under, showCategory = 10) + 
  ggtitle("Top 10 Gene Ontology (GO) Enrichment for Underexpressed Genes") +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5), # Title font size and style
    axis.title.x = element_text(size = 12),  # X-axis title font size
    axis.title.y = element_text(size = 12),  # Y-axis title font size
    axis.text.x = element_text(size = 9),    # X-axis text size
    axis.text.y = element_text(size = 9),    # Y-axis text size
    legend.title = element_text(size = 12),  # Legend title font size
    legend.text = element_text(size = 9)     # Legend text font size
  )

# Plot dotplot for the top 10 GO terms for over-expressed genes
dotplotGO_over = dotplot(go_results_over, showCategory = 10) +
  ggtitle("Top 10 Gene Ontology (GO) Enrichment for Overexpressed Genes") +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5), # Title font size and style
    axis.title.x = element_text(size = 12),  # X-axis title font size
    axis.title.y = element_text(size = 12),  # Y-axis title font size
    axis.text.x = element_text(size = 9),    # X-axis text size
    axis.text.y = element_text(size = 9),    # Y-axis text size
    legend.title = element_text(size = 12),  # Legend title font size
    legend.text = element_text(size = 9)     # Legend text font size
  )

# Arrange the plots side by side
ggarrange(dotplotGO_under, dotplotGO_over, ncol = 2,
          labels = c("A","B"))


# Perform KEGG Pathway Enrichment Analysis on over- and under-expressed genes.

# Map under-expressed genes (under_DE) from SYMBOL to ENTREZID
gene_entrez_under = bitr(
  under_DE,                  # List of under-expressed gene symbols
  fromType = "SYMBOL",       # Input type is gene symbols
  toType   = "ENTREZID",     # Outout type is Entrez IDs
  OrgDb    = org.Hs.eg.db
)

# Map over-expressed genes (over_DE) from SYMBOL to ENTREZID
gene_entrez_over = bitr(
  over_DE,                   # List of over-expressed gene symbols
  fromType = "SYMBOL",       # Input type is gene symbols
  toType   = "ENTREZID",     # Outout type is Entrez IDs
  OrgDb    = org.Hs.eg.db
)


# Perform KEGG pathway enrichment analysis for under-expressed genes
kegg_results_under =  enrichKEGG(
  gene          = gene_entrez_under[,2], # Entrez IDs of under-expressed genes
  organism      = "human",   
  pAdjustMethod = "BH",                  # Adjust p-values using Benjamini-Hochberg
  pvalueCutoff  = 0.05,                  # p-value cutoff for significance
  qvalueCutoff  = 0.05                   # q-value cutoff for significance
)

# Perform KEGG pathway enrichment analysis for over-expressed genes
kegg_results_over =  enrichKEGG(
  gene          = gene_entrez_over[,2], # Entrez IDs of over-expressed genes
  organism      = "human",              
  pAdjustMethod = "BH",                 # Adjust p-values using Benjamini-Hochberg
  pvalueCutoff  = 0.05,                 # p-value cutoff for significance
  qvalueCutoff  = 0.05                  # q-value cutoff for significance
)

# Print the top results for KEGG pathway enrichment of under-expressed genes
print(head(kegg_results_under))

# Print the top results for KEGG pathway enrichment of over-expressed genes
print(head(kegg_results_over))

# Plot dotplot for the top 10 KEGG pathway for under-expressed genes
dotplotKEGG_under = dotplot(kegg_results_under, showCategory = 10) + 
  ggtitle("Top 10 KEGG Pathway Enrichment for Underexpressed Genes") +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5), # Title font size and style
    axis.title.x = element_text(size = 10),   # X-axis title font size
    axis.title.y = element_text(size = 10),   # Y-axis title font size
    axis.text.x = element_text(size = 8),     # X-axis text size
    axis.text.y = element_text(size = 8),     # Y-axis text size
    legend.title = element_text(size = 10),   # Legend title font size
    legend.text = element_text(size = 8)      # Legend text font size
  )

# Plot dotplot for the top 10 KEGG pathway for over-expressed genes
dotplotKEGG_over = dotplot(kegg_results_over, showCategory = 10) + 
  ggtitle("Top 10 KEGG Pathway Enrichment for Overexpressed Genes") +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5), # Title font size and style
    axis.title.x = element_text(size = 10),   # X-axis title font size
    axis.title.y = element_text(size = 10),   # Y-axis title font size
    axis.text.x = element_text(size = 8),     # X-axis text size
    axis.text.y = element_text(size = 8),     # Y-axis text size
    legend.title = element_text(size = 10),   # Legend title font size
    legend.text = element_text(size = 8)      # Legend text font size
  )

# Arrange the plots side by side
ggarrange(dotplotKEGG_under, dotplotKEGG_over, ncol = 2,
          labels = c("A","B"))

# Perform Reactome Pathway Enrichment Analysis on over- and under-expressed genes.

# Perform Reactome pathway enrichment analysis for under-expressed genes
reactome_results_under =  enrichPathway(
  gene          = gene_entrez_under[ ,2], # Entrez IDs of over-expressed genes
  organism      = "human",   
  pAdjustMethod = "BH",                   # Adjust p-values using Benjamini-Hochberg
  pvalueCutoff  = 0.05,                   # p-value cutoff for significance
  qvalueCutoff  = 0.05,                   # q-value cutoff for significance
)

# Perform Reactome pathway enrichment analysis for over-expressed genes
reactome_results_over =  enrichPathway(
  gene          = gene_entrez_over[ ,2], # Entrez IDs of over-expressed genes
  organism      = "human",   
  pAdjustMethod = "BH",                  # Adjust p-values using Benjamini-Hochberg
  pvalueCutoff  = 0.05,                  # p-value cutoff for significance
  qvalueCutoff  = 0.05,                  # q-value cutoff for significance
)

# Print the top results for Reactome pathway enrichment of under-expressed genes
print(head(reactome_results_under))

# Print the top results for Reactome pathway enrichment of over-expressed genes
print(head(reactome_results_over))

# Plot dotplot for the top 10 Reactome pathway for under-expressed genes
dotplotReactome_under = dotplot(reactome_results_under, showCategory = 10) + 
  ggtitle("Top 10 Reactome Pathway Enrichment for Underexpressed Genes")+
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5), # Title font size and style
    axis.title.x = element_text(size = 10),   # X-axis title font size
    axis.title.y = element_text(size = 10),   # Y-axis title font size
    axis.text.x = element_text(size = 8),     # X-axis text size
    axis.text.y = element_text(size = 8),     # Y-axis text size
    legend.title = element_text(size = 10),   # Legend title font size
    legend.text = element_text(size = 8)      # Legend text font size
  )

# Plot dotplot for the top 10 Reactome pathway for over-expressed genes
dotplotReactome_over = dotplot(reactome_results_over, showCategory = 10) + 
  ggtitle("Top 10 Reactome Pathway Enrichment for Overexpressed Genes") +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5), # Title font size and style
    axis.title.x = element_text(size = 10),   # X-axis title font size
    axis.title.y = element_text(size = 10),   # Y-axis title font size
    axis.text.x = element_text(size = 8),     # X-axis text size
    axis.text.y = element_text(size = 8),     # Y-axis text size
    legend.title = element_text(size = 10),   # Legend title font size
    legend.text = element_text(size = 8)      # Legend text font size
  )

# Arrange the plots side by side
ggarrange(dotplotReactome_under, dotplotReactome_over, ncol = 2,
          labels = c("A","B"))

# It is observed that there is a significant enrichment of differentially expressed genes in immunity-related pathways across all the datasets.
# We can visualise these findings effectively using a tree plot, which clusters related pathways and provides an overview of their relationships.

# Gene Ontology (GO) Pathway Enrichment Treeplots

# Calculate pairwise similarity between enriched GO terms
go_results_under_pw = pairwise_termsim(go_results_under) # For under-expressed genes
go_results_over_pw = pairwise_termsim(go_results_over) # For over-expressed genes

# Plot treeplot for GO
treeplotGO_under = treeplot(go_results_under_pw) + # For under-expressed genes
  ggtitle("Treeplot of Gene Ontology (GO) Pathway Enrichment for Underexpressed Genes") +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5), # Title font size and style
    legend.title = element_text(size = 8), # Legend title font size
    legend.text = element_text(size = 8)   # Legend text font size
  )

treeplotGO_over = treeplot(go_results_over_pw) + # For over-expressed genes
  ggtitle("Treeplot of Gene Ontology (GO) Pathway Enrichment for Overexpressed Genes") +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5), # Title font size and style
    legend.title = element_text(size = 8), # Legend title font size
    legend.text = element_text(size = 8),  # Legend text font size
  )

# Arrange the plots
ggarrange(treeplotGO_under, treeplotGO_over, ncol = 1,
          labels = c("A","B"))


# KEGG Pathway Enrichment Treeplots

# Calculate pairwise similarity
kegg_results_under_pw = pairwise_termsim(kegg_results_under) # For under-expressed genes
kegg_results_over_pw = pairwise_termsim(kegg_results_over) # For over-expressed genes

# Plot treeplots for KEGG
treeplotKEGG_under = treeplot(kegg_results_under_pw) + # For under-expressed genes
  ggtitle("Treeplot of KEGG Pathway Enrichment for Underexpressed Genes") +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5), # Title font size and style
    legend.title = element_text(size = 10), # Legend title font size
    legend.text = element_text(size = 10)   # Legend text font size
  )

treeplotKEGG_over = treeplot(kegg_results_over_pw) + # For over-expressed genes
  ggtitle("Treeplot of KEGG Pathway Enrichment for Overexpressed Genes") +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5), # Title font size and style
    legend.title = element_text(size = 10), # Legend title font size
    legend.text = element_text(size = 10),  # Legend text font size
  )

# Arrange the plots
ggarrange(treeplotKEGG_under, treeplotKEGG_over, ncol = 1,
          labels = c("A","B"))


# Reactome Pathway Enrichment Treeplots

# Calculate pairwise similarity
reactome_results_under_pw = pairwise_termsim(reactome_results_under) # For under-expressed genes
reactome_results_over_pw = pairwise_termsim(reactome_results_over) # For over-expressed genes

# Plot treeplots for Reactome
treeplotReactome_under = treeplot(reactome_results_under_pw) + # For under-expressed genes
  ggtitle("Treeplot of Reactome Pathway Enrichment for Underexpressed Genes") +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5), # Title font size and style
    legend.title = element_text(size = 10), # Legend title font size
    legend.text = element_text(size = 10)   # Legend text font size
  )

treeplotReactome_over = treeplot(reactome_results_over_pw) + # For over-expressed genes
  ggtitle("Treeplot of Reactome Pathway Enrichment for Overexpressed Genes") +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5), # Title font size and style
    legend.title = element_text(size = 10), # Legend title font size
    legend.text = element_text(size = 10),  # Legend text font size
  )

# Arrange the plots
ggarrange(treeplotReactome_under, treeplotReactome_over, ncol = 1,
          labels = c("A","B"))


### Data Visualisation

# Obtain Variance Stabilized Transformed (VST) expression values for visualisations.

# Perform variance stabilising transformation
vsd = vst(dds)

# Extract transformed expression values as a matrix of variance-stabilised expression values
vst_values = assay(vsd)


# Plot Principal Component Analysis (PCA) Plots with the VST values obtained.

# PCA Plot grouped by ERBB2_Amplification
pca1 = plotPCA(vsd, intgroup = c("ERBB2_Amplification")) +
  ggtitle("Principal Component Analysis (PCA) Plot: \nAmplified vs. Not Amplified (ERBB2 Status)") +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5), # Adjust title size and alignment
  )

# PCA Plot grouped by "Stage"
pca2 = plotPCA(vsd, intgroup = c("Stage")) +
  ggtitle("Principal Component Analysis (PCA) Plot: \nStage of Cancer Progression") +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),  # Adjust title size and alignment
  )

# Arrange the plots side by side
ggarrange(pca1, pca2, ncol = 2,
          labels = c("A","B"))


# Plot heatmap with the VST values obtained.

# Check if the pheatmap package is installed
if (!requireNamespace("pheatmap", quietly = TRUE)){
  install.packages("pheatmap") # Install pheatmap if not already installed
}

library(pheatmap) # Load pheatmap library


# Subset the dataset on differentially expressed genes ranked by adjusted p-value
top_DE = order(res$padj) # Order by adjusted p-value
vsd_DE_visual = assay(vsd)[top_DE[1:20], ] # For visualisation, extract the top 20 DE genes


# Create a column annotation based on the Amplified/Not Amplified status
annotation_col = data.frame(Amplification = metadata$ERBB2_Amplification)
rownames(annotation_col) = colnames(vsd)  # Match rownames with column names of the data

# Define custom colors for Amplified and Not Amplified
annotation_colors = list(
  Amplification = c(Not.Amplified = "salmon", Amplified = "darkturquoise")
)

# Plot heatmap for the top 20 differentially expressed genes ranked by P-value
pheatmap(
  vsd_DE_visual,
  cluster_rows = TRUE, # Cluster genes based on expression patterns
  cluster_cols = TRUE, # Cluster samples based on similarity
  scale = 'row', # Scale rows for better visualization
  show_colnames = FALSE, # Show sample names
  show_rownames = FALSE, # Show gene names
  annotation_col = annotation_col, # Add column annotations for Amplification status
  main = "Heatmap of Top 20 DE Genes Ranked by Adjusted P-value",
  fontsize = 8
)


### Generates an overall survival model using the glmnet package with the VST values of the DE genes.

# Check if glmnet package is installed
if (!requireNamespace("glmnet", quietly = TRUE)) {
  install.packages("glmnet") # Install glmnet if not already installed
}

# Load the glmnet package
library(glmnet)


# Prepare the overall survival model data

# Extract overall survival status and overall survival months from matched patient data
survival_data = clinPatient_mData[, c("Overall.Survival.Status", "Overall.Survival..Months.")]
survival_data

# Check the structure of Overall.Survival..Months. variable
str(survival_data$Overall.Survival..Months.) # Character

# Convert Overall.Survival..Months. to numeric
survival_data$Overall.Survival..Months. = as.numeric(survival_data$Overall.Survival..Months.)

# Check for non-positive values for Overall.Survival..Months.
table(survival_data$Overall.Survival..Months. <= 0) # 1055 positive values, 13 non-positive values

# Get index numbers of rows where Overall.Survival..Months. > 0
positive_indices = which(survival_data$Overall.Survival..Months. > 0)

# Remove rows with non-positive event times
survival_data = survival_data[survival_data$Overall.Survival..Months. > 0, ]

# Check the structure and summary of Overall.Survival..Months.
str(survival_data$Overall.Survival..Months.) # Numeric
summary(survival_data$Overall.Survival..Months.)

# Check the structure of Overall.Survival.Status variable
str(survival_data$Overall.Survival.Status) # Character
unique(survival_data$Overall.Survival.Status) # 0:LIVING, 1:DECEASED
table(survival_data$Overall.Survival.Status) # 907 zeros, 148 ones

# Recode Overall.Survival.Status variable
survival_data$Overall.Survival.Status = ifelse(survival_data$Overall.Survival.Status == "0:LIVING", 0, 1)

# Check the structure and frequency (0s and 1s) of Overall.Survival.Status
str(survival_data$Overall.Survival.Status) # Numeric
table(survival_data$Overall.Survival.Status) # 907 zeros, 148 ones

# Check got NA values
sum(is.na(survival_data)) # No NA values

# Subset vst values for the top 200 DE genes
vsd_DE_pred = assay(vsd)[top_DE[1:200], ] # For predictive modelling, consider top 200 DE genes
vst_subset = t(vsd_DE_pred)  # Transpose vst values so rows = samples, columns = genes
vst_subset = vst_subset[positive_indices, ]
as.data.frame(vst_subset)

# Combine survival data with expression data
survivalModel_data = cbind(survival_data, vst_subset)
survivalModel_data


# Prepare the matrix of expression data (x) and the survival object (y) for glmnet.

# Check if survival package is installed
if (!requireNamespace("survival", quietly = TRUE)) {
  install.packages("survival") # Install survival if not already installed
}

# Load the survival package
library(survival)

# Design matrix (expression data only)
x_train = as.matrix(survivalModel_data[, -c(1, 2)]) # Exclude Overall.Survival.Status and Overall.Survival..Months.

# Check the dimension of x_train
dim(x_train)

# Survival response variable
y_train = Surv(time = survivalModel_data$Overall.Survival..Months., event = survivalModel_data$Overall.Survival.Status)

# Check the dimension of y_train
dim(y_train)


# Use glmnet to fit a Cox proportional hazards model with Lasso regularisation.

# Check if glmnet package is installed
if (!requireNamespace("glmnet", quietly = TRUE)) {
  install.packages("glmnet") # Install glmnet if not already installed
}

# Check if plotmo package is installed
if (!requireNamespace("plotmo", quietly = TRUE)) {
  install.packages("plotmo") # Install plotmo if not already installed
}

# Load the packages
library(glmnet)
library(plotmo)


# Fit Lasso Regularised Cox model
cox_model = glmnet(
  x = x_train, 
  y = y_train, 
  family = "cox", # Specify Cox proportional hazards model
  alpha = 1 # Alpha = 1 for LASSO regularisation
)


# Perform cross-validation to find the best value of lambda (regularisation parameter).

set.seed(4321)  # Set seed for reproducibility

# Perform cross-validation to select the best lambda
cv_fit = cv.glmnet(
  x = x_train, 
  y = y_train, 
  family = "cox", 
  alpha = 1, # LASSO regularisation
  nfolds = 5 # 5-fold
)


# Plot LASSO path and cross-validation plots for visualisation

# Set up the plotting area to display 2 plots in a single row
par(mfrow = c(1, 2))

# Plot the Lasso path
plot_glmnet(cox_model, xvar = "lambda", label = F, main = "LASSO Path for Cox Model with Top 200 DE Genes Ranked by Adjusted P-value")

# Highlight lambda.min and lambda.1se
abline(v = log(cv_fit$lambda.min), col = "orange", lty = 2)
abline(v = log(cv_fit$lambda.1se), col = "blue", lty = 2)

# Plot cross-validation results
plot(cv_fit)
mtext("Cross-Validation Results for LASSO Cox Model", side = 3, line = 2, cex = 1.2)

par(mfrow = c(1,1))


# Choose the best lambda.
# Identify genes with non-zero coefficients at the optimal lambda.
# Selected genes: Genes with non-zero coefficients are predictive of overall survival.
# Significance: Genes selected by LASSO can provide biological insights into survival.

# Coefficients at lambda.min (Minimum Cross-Validation Error)
coef_min = coef(cox_model, s = cv_fit$lambda.min)
sigfGenes_min = rownames(coef_min)[as.numeric(coef_min) != 0]
cat("Number of significant genes at lambda.min,", cv_fit$lambda.min,":", length(sigfGenes_min), "\n")
cat("Significant genes at lambda.min,", cv_fit$lambda.min, ":", sigfGenes_min, "\n")

# Coefficients at lambda.1se (One Standard Error Rule)
coef_1se = coef(cox_model, s = cv_fit$lambda.1se)
sigfGenes_1se = rownames(coef_1se)[as.numeric(coef_1se) != 0]
cat("Number of significant genes at lambda.1se,", cv_fit$lambda.1se, ":", length(sigfGenes_1se), "\n")
cat("Significant genes at lambda.1se,", cv_fit$lambda.1se, ":", sigfGenes_1se, "\n")


# Compare the cross-validation (CV) errors at the lambda values.

# Cross-validation error at lambda.min
cv_error_min = min(cv_fit$cvm)

# Cross-validation error at lambda.1se
cv_error_1se = cv_fit$cvm[cv_fit$lambda == cv_fit$lambda.1se]

cat("CV error at lambda.min: ", cv_error_min, "\n")
cat("CV error at lambda.1se: ", cv_error_1se, "\n")


# Based on the numbers of significant genes and the CV errors calculated for the possible lambdas, the best lambda seems to be lambda.min.

# This lambda strikes a balance:
# - Selects 22 significant genes.
# - CV error (12.80188) which is the lowest.
# - A model with 22 genes is interpretable and biologically meaningful, while maintaining predictive accuracy.

# Hence, we use lambda.min
coef_min = coef(cox_model, s = cv_fit$lambda.min)
sigfGenes_min = rownames(coef_min)[as.numeric(coef_min) != 0]
cat("Number of significant genes at lambda.min,", cv_fit$lambda.min,":", length(sigfGenes_min), "\n")
cat("Significant genes at lambda.min,", cv_fit$lambda.min, ":", sigfGenes_min, "\n")


# Plot Kaplan Meier Curves and Hazard Ratios

# Check if survminer package is installed
if (!requireNamespace("survminer", quietly = TRUE)) {
  install.packages("survminer") # Install survminer if not already installed
}

# Load the survminer package
library(survminer)


# Create a list to store the plots
kmplot_list = list()

# Loop through each significant gene
for (gene in sigfGenes_min){
  
  # Add gene expression to survival data
  survival_data$gene_expression = vst_subset[, gene]
  
  # Divide samples into high and low expression groups (based on median)
  survival_data$group = ifelse(survival_data$gene_expression > median(survival_data$gene_expression), "High", "Low")
  
  # Create a Kaplan-Meier survival object
  km_fit = survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ group, data = survival_data)
  
  # Plot Kaplan-Meier curve
  km_plot = ggsurvplot(
    km_fit,
    data = survival_data,
    pval = TRUE,                          # Add p-value for log-rank test
    risk.table = TRUE,                    # Add risk table
    legend.labs = c("High Expression", "Low Expression"),
    title = paste("Kaplan-Meier Curve for", gene),
    xlab = "Time (Months)",
    ylab = "Survival Probability"
  )
  
  km_plot$plot = km_plot$plot +
    theme(plot.title = element_text(hjust = 0.5))
  
  print(km_plot)
  
  # Store the plot in the list
  kmplot_list[[gene]] = km_plot
}

# Genes with index numbers 3, 6, 9, 11, 13, 16, 18 are statistically significant (<0.05).
# The rest are not significant.

# Index numbers of genes to include
indices_sigfGenes = c(3, 6, 9, 11, 13, 16, 18)

# Include genes based on index numbers
sigfGenes_min_filtered = sigfGenes_min[indices_sigfGenes]

# Create a list of only the plot components from the ggsurvplot objects
kmplot_filtered = lapply(sigfGenes_min_filtered, function(gene) {
  kmplot_list[[gene]]$plot  # Extract the plot component
})

# Arrange the plots in a 3x3 grid
grid.arrange(grobs = kmplot_list_filtered, ncol = 3, nrow = 3)


# Plot Hazard Ratios as forest plots.

# Create an empty list to store forest plots
forestplot_list = list()

# Loop through each significant gene
for (gene in sigfGenes_min_filtered) {
  
  # Add gene expression to survival data
  survival_data$gene_expression = vst_subset[, gene]
  
  # Fit a Cox proportional hazards model
  cox_model_gene = coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ gene_expression, data = survival_data)
  
  # Generate the forest plot
  forest_plot = ggforest(
    cox_model_gene,
    data = survival_data,
    main = paste("Hazard Ratio for", gene)
  )
  
  # Store the forest plot in the list
  forestplot_list[[gene]] = forest_plot
}

# Arrange the forest plots in a 3x3 grid
grid.arrange(grobs = forestplot_list, ncol = 3, nrow = 3)



