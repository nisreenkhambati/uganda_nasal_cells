#enrichment and cibersort #

# This R script
# does fgsea  and ORA enrichment analysis for nasal and blood
# does cibersort x for estimating immune cell compositions

#### Load libraries ####

library(tidyverse)
library(DESeq2)
library(org.Hs.eg.db)
library(fgsea)
library(reactome.db)
library(clusterProfiler)
library(ggplot2)
library(tidyverse)
library(ggrepel)
library(BiocManager)
library(reshape2)
library(circlize)
library(ggrepel)
library(ComplexHeatmap)
library(GenomicFeatures)
library(ggarchery)
library(txdbmaker)

#### GSEA ENRICHMENT ANALYSIS ####

#### NASAL DATA ####

# read in your DEG dataset
res <- read.csv("data/de_nasal.csv", header=TRUE)
head(res)
colnames(res)[1] <- "symbol"
head(res)
res <- as_tibble(res)
head(res)
res2 <- res %>%
  dplyr::select(symbol, stat) %>%
  na.omit() %>%
  distinct()
head(res2)
ranks <- deframe(res2)
head(ranks, 20)
class(ranks)

#### Hallmark gene sets ####
pathways.hallmark <- gmtPathways("data/h.all.v2023.2.Hs.symbols.gmt.txt")
head(pathways.hallmark)

# run the fgsea function
set.seed(42)
fgseaRes <- fgsea(pathways = pathways.hallmark,
                  stats    = ranks,
                  minSize  = 15,
                  maxSize  = 500)

fgseaRes
# arrange ascending pvalue (smallest first)
fgseaRes.p<- fgseaRes %>%
  dplyr::select(-leadingEdge, -ES) %>%
  arrange(padj)

fgseaRes.p.sub <- fgseaRes.p[1:20,]

#### KEGG pathways ####
pathways.kegg <- gmtPathways("data/c2.cp.kegg_legacy.v2023.2.Hs.symbols.gmt.txt")
head(pathways.kegg)

# run the fgsea function
set.seed(42)
fgseaRes <- fgsea(pathways = pathways.kegg,
                  stats    = ranks,
                  minSize  = 15,
                  maxSize  = 500)

fgseaRes

# arrange ascending pvalue (smallest first)
fgseaRes.p<- fgseaRes %>%
  dplyr::select(-leadingEdge, -ES) %>%
  arrange(padj)

fgseaRes.p.sub <- fgseaRes.p[1:20,]


#### Reactome database ####

pathways.reactome <- gmtPathways("data/c2.cp.reactome.v2023.2.Hs.symbols.gmt.txt")
head(pathways.reactome)

# run the fgsea function
set.seed(42)
fgseaRes <- fgsea(pathways = pathways.reactome,
                  stats    = ranks,
                  minSize  = 15,
                  maxSize  = 500)
fgseaRes

# arrange ascending pvalue (smallest first)
fgseaRes.p<- fgseaRes %>%
  dplyr::select(-leadingEdge, -ES) %>%
  arrange(padj)

fgseaRes.p.sub <- fgseaRes.p[1:20,]

#### BLOOD DATA ####

# read in your DEG dataset
res <- read.csv("data/de_blood.csv", header=TRUE)
colnames(res)[1] <- "symbol"
res <- as_tibble(res)
head(res)
res2 <- res %>%
  dplyr::select(symbol, stat) %>%
  na.omit() %>%
  distinct()
head(res2)
ranks <- deframe(res2)
head(ranks, 20)
class(ranks)

#### Hallmark gene sets ####
pathways.hallmark <- gmtPathways("data/h.all.v2023.2.Hs.symbols.gmt.txt")
head(pathways.hallmark)


# run the fgsea function
set.seed(42)
fgseaRes <- fgsea(pathways = pathways.hallmark,
                  stats    = ranks,
                  minSize  = 15,
                  maxSize  = 500)

fgseaRes

# arrrange ascending pvalue (smallest first)
fgseaRes.p<- fgseaRes %>%
  dplyr::select(-leadingEdge, -ES) %>%
  arrange(padj)

fgseaRes.p.sub <- fgseaRes.p[1:20,]

#### KEGG pathways ####
pathways.kegg <- gmtPathways("data/c2.cp.kegg_legacy.v2023.2.Hs.symbols.gmt.txt")
head(pathways.kegg)

# run the fgsea function
set.seed(42)
fgseaRes <- fgsea(pathways = pathways.kegg,
                  stats    = ranks,
                  minSize  = 15,
                  maxSize  = 500)

fgseaRes

# arrrange ascending pvalue (smallest first)
fgseaRes.p<- fgseaRes %>%
  dplyr::select(-leadingEdge, -ES) %>%
  arrange(padj)

fgseaRes.p.sub <- fgseaRes.p[1:20,]

#### Reactome database ####

pathways.reactome <- gmtPathways("data/c2.cp.reactome.v2023.2.Hs.symbols.gmt.txt")
head(pathways.reactome)

# run the fgsea function
set.seed(42)
fgseaRes <- fgsea(pathways = pathways.reactome,
                  stats    = ranks,
                  eps = 0.0,
                  minSize  = 15,
                  maxSize  = 500)


fgseaRes

# arrrange ascending pvalue (smallest first)
fgseaRes.p<- fgseaRes %>%
  dplyr::select(-leadingEdge, -ES) %>%
  arrange(padj)

fgseaRes.p.sub <- fgseaRes.p[1:20,]


#### ORA ANALYSIS ####

#### NASAL DATA ####

# normalised gene expression table
em <- read.csv("data/em_nasal.csv", header=TRUE, row.names = 1)
# differential expression table
de <- read.csv("data/de_nasal.csv", header=TRUE, row.names = 1)
de <- de[,-c(1,3,4)]
colnames(de) <- c("log2fold", "p", "p.adj")
#sample info table
ss <- read.csv("data/sample_info_nasal.csv", header =TRUE, row.names = 1)
ss <- ss[, c(2,3)]
colnames(ss) <- c("sample", "sample_group")
ss$sample <- row.names(ss)


#### Manipulate files ####
master <- merge(em, de, by.x=0, by.y=0)
row.names(master) = master[,"Row.names"]
names(master)[1] <- "gene_name"
master <- na.omit(master)
em <- na.omit(em)
sorted_order <- order(master[,"p.adj"], decreasing=FALSE)
master <- master[sorted_order,]
master$mean_counts <- rowMeans(master[,ss$sample])
master$mlog10p.adj <- -log10(master$p.adj)
master$sig <- as.factor(master$p.adj < 0.05 & abs(master$log2fold) > 1.0)
scale(em)
scale(t(em))
t(scale(t(em)))
data.frame(t(scale(t(em))))
em_scaled <- data.frame(t(scale(t(em))))
em_scaled <-  na.omit(em_scaled)
master_sig <- subset(master, p.adj < 0.05 & abs(log2fold) >1.0)
sig_genes <- row.names(master_sig)
class(sig_genes)
em_sig <- em[sig_genes,]
em_scaled_sig <- em_scaled[sig_genes,]


# Create a vector of gene names
master_sig_up <- subset(master, p.adj< 0.05 & log2fold >1)
master_sig_down <- subset(master, p.adj< 0.05 & log2fold < -1)
sig_genes_up <- row.names(master_sig_up)
sig_genes_down <- row.names(master_sig_down)


#### PATHWAY ANALYSIS ####
# padjusted value as 0.1.

do_pathway <- function(organism_db, genes)
{
  # libraries
  library(clusterProfiler)

  # converts from ensembl Symbols to Entrez
  sig_genes_entrez <- bitr(genes,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = organism_db)

  # gets the enrichment
  pathway_data <- enrichGO(gene = sig_genes_entrez$ENTREZID,
      OrgDb = organism_db,
      readable = T,
      ont = "BP",
      pvalueCutoff = 0.1,
      qvalueCutoff = 0.10)

  # add this check right here
  if (is.null(pathway_data) || nrow(as.data.frame(pathway_data)) == 0) {
    message("No enriched terms found at this cutoff.")
    return(NULL)
  }

  # bar and dotplot
  ggp.barplot <- barplot(pathway_data, showCategory=10)
  ggp.dotplot <- dotplot(pathway_data, showCategory=10, font.size = 10)
  ggp.goplot <-  goplot(pathway_data, showCategory = 10)
  ggp.cnetplot <- cnetplot(pathway_data,
                          categorySize= "geneNum",
                          color_category='firebrick',
                          color_gene='steelblue',  cex_label_gene = 0.7,
                          cex_gene = 0.3)


  # extract the enrichment table from go enrich
  gene_sets <- pathway_data$geneID
  description <- pathway_data$Description
  p.adj <- pathway_data$p.adjust
  ora_results <- data.frame(cbind(gene_sets, description, p.adj))

  # Create the list to store the results
  pathway_results <- list("tables" = list(), "plots" = list())
  pathway_results$tables$pathway_data <- pathway_data
  pathway_results$tables$ora_results <- ora_results
  pathway_results$plots$ggp.barplot <- ggp.barplot
  pathway_results$plots$ggp.dotplot <- ggp.dotplot
  pathway_results$plots$ggp.goplot <- ggp.goplot
  pathway_results$plots$ggp.cnetplot <- ggp.cnetplot

  # return the results
  return(pathway_results)
}


#pathway for upregulated DEG
pathway_results_up <- do_pathway(org.Hs.eg.db, sig_genes_up)
# view results
pathway_results_up$plots$ggp.dotplot
ora_results_up <- pathway_results_up$tables$ora_results
pathway_data_up <- pathway_results_up$tables$pathway_data

#pathway for downregulated DEG
pathway_results_down <- do_pathway(org.Hs.eg.db, sig_genes_down)
# view results
pathway_results_down$plots$ggp.dotplot
ora_results_down <- pathway_results_down$tables$ora_results
pathway_data_down <- pathway_results_down$tables$pathway_data

nasal_up_dotplot <- pathway_results_up$plots$ggp.dotplot
nasal_down_dotplot <-pathway_results_down$plots$ggp.dotplot


#### BLOOD DATA ####

#### Load datasets ####

# normalised gene expression table
em <- read.csv("data/em_blood.csv", header=TRUE, row.names = 1)
# differential expression table
de_blood <- read.csv("data/de_blood.csv", header=TRUE, row.names = 1)
de_blood <- de_blood[,-c(1,3,4)]
colnames(de_blood) <- c("log2fold", "p", "p.adj")
#sample info table
ss <- read.csv("data/sample_info_blood.csv", header =TRUE, row.names = 1)
ss <- ss[, c(2,3)]
colnames(ss) <- c("sample", "sample_group")
ss$sample <- row.names(ss)


#### Manipulate files ####
master <- merge(em, de_blood, by.x=0, by.y=0)
row.names(master) <- master[,"Row.names"]
names(master)[1] <- "gene_name"
nrow(master[master$p.adj == "NA",])
master <- na.omit(master)
em <- na.omit(em)
sorted_order <- order(master[,"p.adj"], decreasing=FALSE)
master <- master[sorted_order,]
master$mean_counts <- rowMeans(master[,ss$sample])
master$mlog10p.adj <- -log10(master$p.adj)
master$sig <- as.factor(master$p.adj < 0.05 & abs(master$log2fold) > 1.0)

scale(em)
scale(t(em))
t(scale(t(em)))
data.frame(t(scale(t(em))))
em_scaled <- data.frame(t(scale(t(em))))
em_scaled <-  na.omit(em_scaled)
master_sig <- subset(master, p.adj < 0.05 & abs(log2fold) >1.0)
sig_genes_blood <- row.names(master_sig)
class(sig_genes_blood)
em_sig <- em[sig_genes_blood,]
em_scaled_sig <- em_scaled[sig_genes_blood,]

# Create a vector of gene names
master_sig_up <- subset(master, p.adj< 0.05 & log2fold >1)
master_sig_down <- subset(master, p.adj< 0.05 & log2fold < -1)
sig_genes_up <- row.names(master_sig_up)
sig_genes_down <- row.names(master_sig_down)

#### PATHWAY ANALYSIS -  ####

#cut off is 0.05
do_pathway <- function(organism_db, genes)
{
  # libraries
  library(clusterProfiler)

  # converts from ensembl Symbols to Entrez
  sig_genes_entrez <- bitr(genes,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = organism_db)

  # gets the enrichment
  pathway_data <- enrichGO(gene = sig_genes_entrez$ENTREZID,
      OrgDb = organism_db,
      readable = T,
      ont = "BP",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.10)

  # add this check right here
  if (is.null(pathway_data) || nrow(as.data.frame(pathway_data)) == 0) {
    message("No enriched terms found at this cutoff.")
    return(NULL)
  }

  # bar and dotplot
  ggp.barplot <- barplot(pathway_data, showCategory=10)
  ggp.dotplot <- dotplot(pathway_data, showCategory=10, font.size = 10)
  ggp.goplot <- goplot(pathway_data, showCategory = 10)
  ggp.cnetplot <- cnetplot(pathway_data,
                          categorySize= "geneNum",
                          color_category='firebrick',
                          color_gene='steelblue',  cex_label_gene = 0.7)

  # extract the enrichment table from go enrich
  gene_sets <- pathway_data$geneID
  description <- pathway_data$Description
  p.adj <- pathway_data$p.adjust
  ora_results <- data.frame(cbind(gene_sets, description, p.adj))

  # Create the list to store the results
  pathway_results <- list("tables" = list(), "plots" = list())
  pathway_results$tables$pathway_data <- pathway_data
  pathway_results$tables$ora_results <- ora_results
  pathway_results$plots$ggp.barplot <- ggp.barplot
  pathway_results$plots$ggp.dotplot <- ggp.dotplot
  pathway_results$plots$ggp.goplot <- ggp.goplot
  pathway_results$plots$ggp.cnetplot <- ggp.cnetplot

  # return the results
  return(pathway_results)
}


#pathway for upregulated DEG
pathway_results_up <- do_pathway(org.Hs.eg.db, sig_genes_up)
# view results
pathway_results_up$plots$ggp.dotplot
ora_results_up <- pathway_results_up$tables$ora_results
pathway_data_up <- pathway_results_up$tables$pathway_data

#pathway for downregulated DEG
pathway_results_down <- do_pathway(org.Hs.eg.db, sig_genes_down)
# view results - start here
pathway_results_down$plots$ggp.dotplot
ora_results_down <- pathway_results_down$tables$ora_results
pathway_data_down <- pathway_results_down$tables$pathway_data

blood_up_dotplot<- pathway_results_up$plots$ggp.dotplot
blood_down_dotplot<-pathway_results_down$plots$ggp.dotplot


#### panelled fig ####

# Dotplot #

nasal_up_dotplot
nasal_down_dotplot
blood_up_dotplot
blood_down_dotplot

cowplot::plot_grid(nasal_up_dotplot,
                   blood_up_dotplot,
                   nasal_down_dotplot,
                   blood_down_dotplot,
                   ncol=2,
                   labels=LETTERS[1:4])



title_row <- cowplot::plot_grid(
  cowplot::ggdraw() + cowplot::draw_label("Nasal", fontface = "bold", size = 14),
  cowplot::ggdraw() + cowplot::draw_label("Blood", fontface = "bold", size = 14),
  ncol = 2)

plot_grid_with_labels <- cowplot::plot_grid(nasal_up_dotplot,
                                            blood_up_dotplot,
                                            nasal_down_dotplot,
                                            blood_down_dotplot,
                                            ncol=2,
                                            labels=LETTERS[1:4])


final_plot_1 <- cowplot::plot_grid(title_row, plot_grid_with_labels,
  ncol = 1, rel_heights = c(0.05, 1))

print(final_plot_1)




#### CIBERSORTX: Immune cell composition  ####

# NASAL #

#### Preparing data  ####
counts_data <- read.csv("data/nasalcounts.csv", header=TRUE, row.names = 1)
coldata <- read.csv("data/nasalcoldata.csv", header=TRUE, row.names = 1)
str(coldata)
coldata <- coldata %>%
  dplyr::mutate(
    status = factor(status, levels = c("control", "case")),
    sex = factor(sex, levels = c("male", "female"))
  )
str(coldata)


# gene length

# need to download GTF annotation file from
# https://www.gencodegenes.org/human/release_44.html)
# GTF file is comprehensive gene annotation on the primary assembly (chromosomes
# and scaffolds - PRI Regions) sequence regions (Release 44 -GRCh38.p14)

txdb <- txdbmaker::makeTxDbFromGFF(
    "gencode.v44.primary_assembly.annotation.gtf", format="gtf")
exons.list.per.gene <- exonsBy(txdb,by="gene")
head(exons.list.per.gene )
exonic.gene.sizes <- as.data.frame(sum(width(GenomicRanges::reduce(exons.list.per.gene))))
head(exonic.gene.sizes)
names(exonic.gene.sizes) <- "gene_length"
exonic.gene.sizes.df <- data.frame(
    ensembl_gene_id_version = rownames(exonic.gene.sizes),
    exonic.gene.sizes,
    row.names = NULL)
head(exonic.gene.sizes.df)
genes <- read.delim("data/genes_annotation_final.txt",
    header = TRUE, sep = "\t", dec=".")
genes_test <- merge(genes, exonic.gene.sizes.df,
    by = "ensembl_gene_id_version", all.x = TRUE)
any(is.na(genes_test$gene_length)) #FALSE


# conversion to TPM #
gene_lengths_kb <- genes_test$gene_length[match(row.names(counts_data),
    genes_test$external_gene_name)] / 1000
any(is.na(gene_lengths_kb))
identical(rownames(counts_data),
    genes_test$external_gene_name[match(rownames(counts_data),
        genes_test$external_gene_name)]) #TRUE

# Calculate RPK (Reads Per Kilobase)
rpk <- sweep(counts_data, 1, gene_lengths_kb, FUN = "/")
scaling_factor <- colSums(rpk) / 1e6
tpm <- sweep(rpk, 2, scaling_factor, FUN = "/")
counts_data_tpm <- data.frame(Genes = rownames(counts_data), tpm)
counts_data_tpm  <- na.omit(counts_data_tpm)


#### Run cibersort ####
# Run this count_data_tpm matrix using the "Impute Cell Fractions" analysis
# module on the CIBERSORTX website
# Run with B-mode batch correction enabled, disable quantile normalization true,
# run in relative mode, and set permutations 500 (recommended is >=100)
# Source GEP file used for batch correction: LM22.update-gene-symbols.txt
# The output includes a p-value for the global deconvolution of each sample. A
# p-value threshold <0.05 is recommended.

#### Analyse NASAL CIBERSORT ####
cibersort <- read.delim("data/CIBERSORTx_Job2_Adjusted.txt",
    header = TRUE, sep = "\t", dec = ".")
names(cibersort) <- gsub(x = names(cibersort), pattern = "\\.", replacement = "_")
row.names(cibersort) <- cibersort$Mixture
cibersort <- cibersort[,c(2:23)]
row_sums <- rowSums(cibersort)
print(row_sums)

# Combining columns of LM22 cell sub-types
cibersort$B_cells <- cibersort$B_cells_naive + cibersort$B_cells_memory
cibersort$T_cells_CD4 <- cibersort$T_cells_CD4_naive +
    cibersort$T_cells_CD4_memory_resting +
    cibersort$T_cells_CD4_memory_activated +
    cibersort$T_cells_follicular_helper +
    cibersort$T_cells_regulatory__Tregs_
cibersort$T_cells_gamma_delta <- cibersort$T_cells_gamma_delta
cibersort$NK_cells <- cibersort$NK_cells_resting + cibersort$NK_cells_activated
cibersort$Mono_Macrophages <- cibersort$Monocytes +
    cibersort$Macrophages_M0 +
    cibersort$Macrophages_M1 +
    cibersort$Macrophages_M2
cibersort$Dendritic_cells <- cibersort$Dendritic_cells_resting +
    cibersort$Dendritic_cells_activated
cibersort$Mast_cells <- cibersort$Mast_cells_resting +
    cibersort$Mast_cells_activated
cibersort_subset <- cibersort %>%
    dplyr::select(B_cells,
        T_cells_CD8,
        T_cells_CD4,
        T_cells_gamma_delta,
        Mono_Macrophages,
        Neutrophils,
        NK_cells,
        Dendritic_cells,
        Eosinophils,
        Mast_cells,
        Plasma_cells)
row_sums <- rowSums(cibersort_subset)
print(row_sums)

# Merging with Metadata
cibersort_metadata <- merge(coldata, cibersort_subset, by='row.names')
row.names(cibersort_metadata) <- cibersort_metadata$Row.names
cibersort_metadata$Row.names <- NULL
cibersort_metadata <- cibersort_metadata %>%
  dplyr::select (status,
      B_cells,
      T_cells_CD8,
      T_cells_CD4,
      T_cells_gamma_delta,
      Mono_Macrophages,
      Neutrophils,
      NK_cells,
      Dendritic_cells,
      Eosinophils,
      Mast_cells,
      Plasma_cells)


# summarise the data
summary_of_cells <- summary(cibersort_metadata)
summary_of_cells
cells_control <- cibersort_metadata[cibersort_metadata$status == "control",]
cells_case <- cibersort_metadata[cibersort_metadata$status == "case",]
summary_of_cells_control <- summary(cells_control)
summary_of_cells_case <- summary(cells_case)
summary_of_cells_control
summary_of_cells_case


str(cibersort_metadata)
# filter the cibersort_metadata dataframe to keep only the cell types that are
# prevalent i.e., appear in at least 50% of the samples.
check_prevalence <- function(cell_data, threshold = 0.50) {
  mean(cell_data != 0) >= threshold
}

cell_types <- colnames(cibersort_metadata)[2:ncol(cibersort_metadata)]
prevalent_cells <- sapply(cibersort_metadata[, cell_types], check_prevalence)
cell_types_to_keep <- names(prevalent_cells[prevalent_cells])
test <- cibersort_metadata[, c("status", cell_types_to_keep)]
cibersort_metadata <- test



cibersort_metadata_long <- cibersort_metadata %>%
    rownames_to_column(var = "Sample")
cibersort_metadata_long <- cibersort_metadata_long %>%
  pivot_longer(
    cols = -c(Sample, status),
    names_to = "CellType",
    values_to = "Proportion")


cibersort_metadata_long$cell_type2[cibersort_metadata_long$CellType=="B_cells"] <- "B cells"
cibersort_metadata_long$cell_type2[cibersort_metadata_long$CellType=="Plasma_cells"] <- "Plasma cells"
cibersort_metadata_long$cell_type2[cibersort_metadata_long$CellType=="T_cells_CD8"] <- "CD8+ T cells"
cibersort_metadata_long$cell_type2[cibersort_metadata_long$CellType=="T_cells_CD4"] <- "CD4+ T cells"
cibersort_metadata_long$cell_type2[cibersort_metadata_long$CellType=="NK_cells"] <- "NK cells"
cibersort_metadata_long$cell_type2[cibersort_metadata_long$CellType=="Mono_Macrophages"] <- "Monocytes/macrophages"
cibersort_metadata_long$cell_type2[cibersort_metadata_long$CellType=="Dendritic_cells"] <- "Dendritic cells"
cibersort_metadata_long$cell_type2[cibersort_metadata_long$CellType=="Mast_cells"] <- "Mast cells"
cibersort_metadata_long$cell_type2[cibersort_metadata_long$CellType=="Eosinophils"] <- "Eosinophils"
cibersort_metadata_long$cell_type2[cibersort_metadata_long$CellType=="Neutrophils"] <- "Neutrophils"


# check for normality

# controls
shapiro.test(cells_control$B_cells)
shapiro.test(cells_control$T_cells_CD8)
shapiro.test(cells_control$T_cells_CD4)
shapiro.test(cells_control$Mono_Macrophages)
shapiro.test(cells_control$Neutrophils)
shapiro.test(cells_control$NK_cells)
shapiro.test(cells_control$Dendritic_cells)
shapiro.test(cells_control$Eosinophils)
shapiro.test(cells_control$Mast_cells)
shapiro.test(cells_control$Plasma_cells)

shapiro.test(cells_case$B_cells)
shapiro.test(cells_case$T_cells_CD8)
shapiro.test(cells_case$T_cells_CD4)
shapiro.test(cells_case$Mono_Macrophages)
shapiro.test(cells_case$Neutrophils)
shapiro.test(cells_case$NK_cells)
shapiro.test(cells_case$Dendritic_cells)
shapiro.test(cells_case$Eosinophils)
shapiro.test(cells_case$Mast_cells)
shapiro.test(cells_case$Plasma_cells)


# statistical tests

wilcox_test_cd8 <- wilcox.test(cells_control$T_cells_CD8, cells_case$T_cells_CD8)
CD8_p <- wilcox_test_cd8$p.value

wilcox_test_den <- wilcox.test(cells_control$Dendritic_cells, cells_case$Dendritic_cells)
den_p <- wilcox_test_den$p.value

wilcox_test_eos <- wilcox.test(cells_control$Eosinophils, cells_case$Eosinophils)
eos_p <- wilcox_test_eos$p.value

wilcox_test_mast <- wilcox.test(cells_control$Mast_cells, cells_case$Mast_cells)
mast_p <- wilcox_test_mast$p.value

wilcox_test_plasma <- wilcox.test(cells_control$Plasma_cells, cells_case$Plasma_cells)
plasma_p<- wilcox_test_plasma$p.value

ttest_B_cells <- t.test(cells_control$B_cells, cells_case$B_cells)
B_cells_p <- ttest_B_cells$p.value

ttest_cd4 <- t.test(cells_control$T_cells_CD4, cells_case$T_cells_CD4)
CD4_p <- ttest_cd4$p.value


ttest_mono <- t.test(cells_control$Mono_Macrophages, cells_case$Mono_Macrophages)
macro_p <- ttest_mono$p.value

t_test_neutrophils <- t.test(cells_control$Neutrophils, cells_case$Neutrophils)
Neutrophils_p <- t_test_neutrophils$p.value

t_test_nk <- t.test(cells_control$NK_cells, cells_case$NK_cells)
nk_p <-t_test_nk$p.value


all_p_values <- c(B_cells_p, Neutrophils_p, CD4_p, mast_p, macro_p, CD8_p,
    eos_p, plasma_p, den_p, nk_p)
adjusted_p_values <- p.adjust(all_p_values, method = "BH")
adjusted_p_values_results <- data.frame(
  Test = c("B_cells", "Neutrophils", "T_cells_CD4", "Mast_cells",
      "Mono_Macrophages", "T_cells_CD8", "Eosinophils", "Plasma_cells",
      "Dendritic_cells", "NK_cells"),
  Adjusted_p_value = adjusted_p_values
)
adjusted_p_values_results

# fig
nasal_imps <- ggplot(cibersort_metadata_long,
  aes(x = cell_type2, y = Proportion, fill=status)) +  ylim(0, 0.99) +
  geom_boxplot() +
  labs(x = "Immune Cell Type", y = "Imputed immune cell proportions") +
  ggtitle("Nasal samples") + # No significant differences in cell populations
  annotate("text", x=1, y=0.96, size=3, label= "p = 0.25") +
  annotate("text", x=2, y=0.96, size=3, label= "p = 0.47") +
  annotate("text", x=3, y=0.96, size=3, label= "p = 0.75") +
  annotate("text", x=4, y=0.96, size=3, label= "p = 0.94") +
  annotate("text", x=5, y=0.96, size=3, label= "p = 0.66") +
  annotate("text", x=6, y=0.96, size=3, label= "p = 0.66") +
  annotate("text", x=7, y=0.96, size=3, label= "p = 0.66") +
  annotate("text", x=8, y=0.96, size=3, label= "p = 0.66") +
  annotate("text", x=9, y=0.96, size=3, label= "p = 0.66") +
  annotate("text", x=10, y=0.96, size=3, label= "p = 0.66") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12))




# BLOOD #
counts_data_blood <- read.csv("data/bloodcounts.csv",
    header=TRUE, row.names = 1)
coldata_blood <- read.csv("data/bloodcoldata.csv", header=TRUE, row.names = 1)
str(coldata_blood)
coldata_blood <- coldata_blood %>%
  dplyr::mutate(
    status = factor(status, levels = c("control", "case")),
    sex = factor(sex, levels = c("male", "female"))
  )
str(coldata_blood)

# conversion to TPM
gene_lengths_kb <- genes_test$gene_length[match(row.names(counts_data_blood),
    genes_test$external_gene_name)] / 1000
any(is.na(gene_lengths_kb))
identical(rownames(counts_data_blood),
    genes_test$external_gene_name[match(rownames(counts_data_blood),
        genes_test$external_gene_name)])

# Calculate RPK (Reads Per Kilobase)
rpk <- sweep(counts_data_blood, 1, gene_lengths_kb, FUN = "/")
scaling_factor <- colSums(rpk) / 1e6
tpm <- sweep(rpk, 2, scaling_factor, FUN = "/")
counts_data_blood_tpm <- data.frame(Genes = rownames(counts_data_blood), tpm)
counts_data_blood_tpm  <- na.omit(counts_data_blood_tpm)


#### Run cibersort ####
# Run this count_data_blood_tpm matrix using the "Impute Cell Fractions"
# analysis module on the CIBERSORTx website
# Run with B-mode batch correction enabled, disable quantile normalization true,
# run in relative mode, and set permutations 500 (recommended is >=100)
# Source GEP file used for batch correction: LM22.update-gene-symbols.txt
# The output includes a p-value for the global deconvolution of each sample. A
# p-value threshold <0.05 is recommended.

#### Analyse BLOOD CIBERSORT ####
cibersort <- read.delim("data/CIBERSORTx_Job3_Adjusted.txt",
    header = TRUE, sep = "\t", dec = ".")
names(cibersort) <- gsub(x = names(cibersort), pattern = "\\.", replacement = "_")
row.names(cibersort) <- cibersort$Mixture
cibersort <- cibersort[,c(2:23)]
row_sums <- rowSums(cibersort)
print(row_sums)


# Combining columns of LM22 cell types
cibersort$B_cells <- cibersort$B_cells_naive + cibersort$B_cells_memory
cibersort$T_cells_CD4 <- cibersort$T_cells_CD4_naive +
    cibersort$T_cells_CD4_memory_resting +
    cibersort$T_cells_CD4_memory_activated +
    cibersort$T_cells_follicular_helper +
    cibersort$T_cells_regulatory__Tregs_
cibersort$T_cells_gamma_delta <- cibersort$T_cells_gamma_delta
cibersort$NK_cells <- cibersort$NK_cells_resting + cibersort$NK_cells_activated
cibersort$Mono_Macrophages <- cibersort$Monocytes +
    cibersort$Macrophages_M0 +
    cibersort$Macrophages_M1 +
    cibersort$Macrophages_M2
cibersort$Dendritic_cells <- cibersort$Dendritic_cells_resting +
    cibersort$Dendritic_cells_activated
cibersort$Mast_cells <- cibersort$Mast_cells_resting +
    cibersort$Mast_cells_activated

cibersort_subset <- cibersort %>%
    dplyr::select(B_cells,
        T_cells_CD8,
        T_cells_CD4,
        T_cells_gamma_delta,
        Mono_Macrophages,
        Neutrophils,
        NK_cells,
        Dendritic_cells,
        Eosinophils,
        Mast_cells,
        Plasma_cells)
row_sums <- rowSums(cibersort_subset)
print(row_sums)

cibersort_metadata <- merge(coldata_blood, cibersort_subset, by='row.names')
row.names(cibersort_metadata) <- cibersort_metadata$Row.names
cibersort_metadata$Row.names <- NULL
cibersort_metadata <- cibersort_metadata %>%
  dplyr::select (status,
      B_cells,
      T_cells_CD8,
      T_cells_CD4,
      T_cells_gamma_delta,
      Mono_Macrophages,
      Neutrophils,
      NK_cells,
      Dendritic_cells,
      Eosinophils,
      Mast_cells,
      Plasma_cells)

# summarise the blood data
summary_of_cells <- summary(cibersort_metadata)
summary_of_cells
cells_control <- cibersort_metadata[cibersort_metadata$status == "control",]
cells_case <- cibersort_metadata[cibersort_metadata$status == "case",]
summary_of_cells_control <- summary(cells_control)
summary_of_cells_case <- summary(cells_case)
summary_of_cells_control
summary_of_cells_case


str(cibersort_metadata)
# filter the cibersort_metadata dataframe to keep only the cell types that are
# prevalent i.e., appear in at least 50% of the samples.
check_prevalence <- function(cell_data, threshold = 0.50) {
  mean(cell_data != 0) >= threshold
}

cell_types <- colnames(cibersort_metadata)[2:ncol(cibersort_metadata)]
prevalent_cells <- sapply(cibersort_metadata[, cell_types], check_prevalence)
cell_types_to_keep <- names(prevalent_cells[prevalent_cells])
test <- cibersort_metadata[, c("status", cell_types_to_keep)]
cibersort_metadata <- test


cibersort_metadata_long <- cibersort_metadata %>%
    rownames_to_column(var = "Sample")
cibersort_metadata_long <- cibersort_metadata_long %>%
  pivot_longer(
    cols = -c(Sample, status),
    names_to = "CellType",
    values_to = "Proportion"
  )

cibersort_metadata_long$cell_type2[cibersort_metadata_long$CellType=="B_cells"] <- "B cells"
cibersort_metadata_long$cell_type2[cibersort_metadata_long$CellType=="T_cells_CD8"] <- "CD8+ T cells"
cibersort_metadata_long$cell_type2[cibersort_metadata_long$CellType=="T_cells_CD4"] <- "CD4+ T cells"
cibersort_metadata_long$cell_type2[cibersort_metadata_long$CellType=="NK_cells"] <- "NK cells"
cibersort_metadata_long$cell_type2[cibersort_metadata_long$CellType=="Mono_Macrophages"] <- "Monocytes/macrophages"
cibersort_metadata_long$cell_type2[cibersort_metadata_long$CellType=="Dendritic_cells"] <- "Dendritic cells"
cibersort_metadata_long$cell_type2[cibersort_metadata_long$CellType=="Mast_cells"] <- "Mast cells"
cibersort_metadata_long$cell_type2[cibersort_metadata_long$CellType=="Neutrophils"] <- "Neutrophils"

cell_types_to_keep


# check for normality
#controls
shapiro.test(cells_control$B_cells)
shapiro.test(cells_control$T_cells_CD8)
shapiro.test(cells_control$T_cells_CD4)
shapiro.test(cells_control$Mono_Macrophages)
shapiro.test(cells_control$Neutrophils)
shapiro.test(cells_control$NK_cells)
shapiro.test(cells_control$Dendritic_cells)
shapiro.test(cells_control$Mast_cells)


#cases
shapiro.test(cells_case$B_cells)
shapiro.test(cells_case$T_cells_CD8)
shapiro.test(cells_case$T_cells_CD4)
shapiro.test(cells_case$Mono_Macrophages)
shapiro.test(cells_case$Neutrophils)
shapiro.test(cells_case$NK_cells)
shapiro.test(cells_case$Dendritic_cells)
shapiro.test(cells_case$Mast_cells)

#statistical tests
wilcox_test_cd4 <- wilcox.test(cells_control$T_cells_CD4, cells_case$T_cells_CD4)
CD4_p <- wilcox_test_cd4$p.value

wilcox_test_macro <- wilcox.test(cells_control$Mono_Macrophages, cells_case$Mono_Macrophages)
macro_p <- wilcox_test_macro $p.value

wilcox_test_den <- wilcox.test(cells_control$Dendritic_cells, cells_case$Dendritic_cells)
den_p <- wilcox_test_den$p.value

wilcox_test_mast <- wilcox.test(cells_control$Mast_cells, cells_case$Mast_cells)
mast_p <- wilcox_test_mast$p.value


ttest_B_cells <- t.test(cells_control$B_cells, cells_case$B_cells)
B_cells_p <- ttest_B_cells$p.value


ttest_cd8 <- t.test(cells_control$T_cells_CD8, cells_case$T_cells_CD8)
CD8_p <- ttest_cd8$p.value


t_test_neutrophils <- t.test(cells_control$Neutrophils, cells_case$Neutrophils)
Neutrophils_p <- t_test_neutrophils$p.value


t_test_nk <- t.test(cells_control$NK_cells, cells_case$NK_cells)
nk_p <-t_test_nk$p.value

all_p_values <- c(B_cells_p, Neutrophils_p, CD4_p, mast_p, macro_p, CD8_p, den_p, nk_p)
adjusted_p_values <- p.adjust(all_p_values, method = "BH")
adjusted_p_values_results <- data.frame(
  Test = c("B_cells", "Neutrophils", "T_cells_CD4", "Mast_cells",
      "Mono_Macrophages", "T_cells_CD8", "Dendritic_cells", "NK_cells"),
  Adjusted_p_value = adjusted_p_values)

adjusted_p_values_results


#fig
blood_imps <- ggplot(cibersort_metadata_long,
  aes(x = cell_type2, y = Proportion, fill=status)) +  ylim(0, 0.99) +
  geom_boxplot() +
  labs(x = "Immune Cell Type", y = "Imputed immune cell proportions") +
  ggtitle("Blood samples") + # No significant differences in cell populations
  annotate("text", x=1, y=0.96, size=3, label= "p = 0.96") +
  annotate("text", x=2, y=0.96, size=3, label= "p = 0.25") +
  annotate("text", x=3, y=0.96, size=3, label= "p = 0.22") + # CD8+ T cells:
  annotate("text", x=4, y=0.96, size=3, label= "p = 0.97") +
  annotate("text", x=5, y=0.96, size=3, label= "p = 0.68") +
  annotate("text", x=6, y=0.96, size=3, label= "p = 0.97") +
  annotate("text", x=7, y=0.96, size=3, label= "p = 0.04") + #
  annotate("text", x=8, y=0.96, size=3, label= "p = 0.04") + #
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12))


# paneled figs

nasal_imps
blood_imps
cowplot::plot_grid(nasal_imps, blood_imps, ncol=1, labels=LETTERS[1:2])

nasal_imps <- nasal_imps + labs(fill = "TB status")
blood_imps <- blood_imps + labs(fill = "TB status")

legend_theme <- theme(
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 12)
)

nasal_imps2 <- nasal_imps + legend_theme + theme(legend.position = "none")
blood_imps2 <- blood_imps + legend_theme + theme(legend.position = "none")

legend <- cowplot::get_legend(
  blood_imps2 + theme(legend.position = "right")
)

combined_imputation <- cowplot::plot_grid(
  cowplot::plot_grid(nasal_imps2, blood_imps2, ncol=1, labels=LETTERS[1:2]),
  legend,
  ncol = 2,
  rel_widths = c(1.5, 0.2)
)
