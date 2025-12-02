## **Uganda_nasal_cells**

This repository contains data and code used in the analysis of the RNA sequencing nasal and blood Ugandan data 

### **Code**

1.	*qc_pca_github.R*: processes the nasal and blood counts data and meta files and performs quality control, heatmaps and PCA.

2.	*DEG_visualisation_github.R*: performs differential gene expression between TB cases and controls for nasal and blood samples and creates volcano plots, boxplots with jitter, heatmaps based on DEG, scatter plots and GSEA enrichment plots

3.	*enrichment_cibersortx_github.R*: performs gsea and ORA enrichment analysis for nasal and blood samples and estimates immune cell composition using cibersortx

4.	*roc_models_github.R*: trains six machine learning models using nasal and blood DEG and calculates cross-validated AUC for predicting TB

5.	*TB_signatureprofiler_github.R*: evaluates the performance of previously published diagnostic TB blood signatures in the nasal and blood datasets. 

### **Data** 

*full_metadata.csv*: metadata for nasal and blood samples 

*rnaseq_counts_protein_coding_final.txt*: raw RNAseq counts of protein coding genes for nasal and blood samples

*genes_annotation_final.txt*: ensemble gene id to gene symbol conversion

*nasalcounts.csv*: raw counts data for nasal samples that have had mitochondrial, haemoglobin and ribosomal genes, as well as low count genes filtered out, sample outliers excluded and genes id converted to symbol (32 nasal samples)

*nasalcoldata.csv*: meta data to match nasalcounts.csv

*bloodcounts.csv*: raw counts data for blood samples that have had mitochondrial, haemoglobin and ribosomal genes, as well as low count genes filtered out, and genes id converted to symbol (40 blood samples)

*bloodcoldata.csv*: meta data to match bloodcounts.csv:

*em_nasal.csv*: table of normalised nasal counts following DESeq2

*de_nasal.csv*: table of nasal differential gene expression results using DESeq2

*sample_info_nasal.csv*: sample info for nasal samples 

*em_blood.csv*: table of normalised blood counts following DESeq2

*de_blood.csv*: table of blood differential genesexpression results using DESeq2

*sample_info_blood.csv*: sample info for blood samples 

*h.all.v2023.2.Hs.symbols.gmt.txt*: Hallmark gene sets (MSigDB)

*c2.cp.kegg_legacy.v2023.2.Hs.symbols.gmt.txt*: Canonical Pathways gene sets derived from the KEGG pathway database (MSigDB)

*c2.cp.reactome.v2023.2.Hs.symbols.gmt.txt*: Canonical Pathways gene sets derived from the Reactome pathway database (MSigDB)

*CIBERSORTx_Job2_Adjusted.txt*: CIBERSORTx results for nasal samples 

*CIBERSORTx_Job3_Adjusted.txt*: CIBERSORTx results for blood samples 

*sig_genes_nasal.rds*: vector of 44 nasal DEG

*sig_genes_blood rds*: vector of 238 blood DEG



