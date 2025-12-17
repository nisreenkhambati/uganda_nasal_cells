# Upload rnaseq counts data to R and subset for protein coding genes.

# This R script uploads the gene count data following STAR alignment
# of all fastq files and subsets for protein coding genes.
# It creates a gene annotation file converting ensembl_gene_id_version to external_gene_name
# It also changes sample names to make more readable.

#### Load libraries ####
library(biomaRt)

#### Load data  ####
# Import counts data for 75 fastq samples post STAR allignment
counts_data <-read.delim("data/count_table_final.csv", header = TRUE, sep = "\t", dec=".")
colnames(counts_data) # see column names are sample names
# 62754 genes

# Rename samples 
column_names <- colnames(counts_data)
for (old_name in column_names) {
  new_name <- sub(".*_(S\\d+)_.*", "\\1", old_name)
  colnames(counts_data)[colnames(counts_data) == old_name] <- new_name
}

# get an annotated gene list 
ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  version = 110   # Ensembl release used in Dec 2023 at time of analysis
)

genes <- getBM(
  attributes = c(
    'ensembl_gene_id_version',
    'external_gene_name',
    'gene_biotype'
  ),
  filters = "ensembl_gene_id_version",
  values = counts_data$ID,
  mart = ensembl
)
# 62754

#### Filter for protein coding #####
protein <- c("protein_coding") 
genes <- genes[genes$gene_biotype %in% protein ,]
# 20070

counts_data_protein <- counts_data[counts_data$ID %in% genes$ensembl_gene_id_version,]
# 20070

rnaseq_counts_protein_coding_final <- counts_data_protein
genes_annotation_final <- genes






