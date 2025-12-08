# QC and PCA
# This R script processes the counts data and metafiles and does QC, heatmaps and pca.

#### Load libraries ####
library(ggplot2)
library(DESeq2)
library(tidyverse)
library(ggrepel)
library(ComplexHeatmap)
library(grid)
library(RColorBrewer)
library(pheatmap)
library(sva)
library(org.Hs.eg.db)
library(BiocManager)
library(reshape2)
library(circlize)
library(ggpubr)


#### Preparing nasal data  ####
coldata<- read.csv("data/full_metadata.csv", header = TRUE, row.names = 4, sep = ",")
coldata <- coldata[,-1]
coldata <- coldata[coldata$type == "nasal",]
row.names(coldata)
row_order <- paste0("S", 1:35)
coldata <- coldata[row_order,]

counts_data <-read.delim("data/rnaseq_counts_protein_coding_final.txt", header = TRUE, sep = "\t", dec=".")
row.names(counts_data) <- counts_data$ID #make ID rownames
counts_data <- counts_data[, -1]
colnames(counts_data)
counts_data <- counts_data[colnames(counts_data) %in% rownames(coldata)]
colnames(counts_data)
counts_data <- counts_data[, rownames(coldata)]
colnames(counts_data)
all(rownames(coldata) %in% colnames(counts_data))
all(rownames(coldata) == colnames(counts_data))
all(colnames(counts_data) %in% rownames(coldata))
all(colnames(counts_data) == rownames(coldata))


# Remove non-useful genes
genes <-read.delim("data/genes_annotation_final.txt", header = TRUE, sep = "\t", dec=".")
genes[(duplicated(genes$ensembl_gene_id_version) == TRUE),]
all(rownames(counts_data) %in% genes$ensembl_gene_id_version)
counts_data <- merge(counts_data, genes[1:2], by.x = "row.names", by.y= "ensembl_gene_id_version", all= FALSE)
row.names(counts_data) <- counts_data$Row.names
counts_data <- counts_data[, -1]

counts_data<- counts_data[!grepl("HB", counts_data$external_gene_name),]
counts_data<- counts_data[!grepl("^MT", counts_data$external_gene_name),]
counts_data<- counts_data[!grepl("RPL", counts_data$external_gene_name),]
counts_data<- counts_data[!grepl("RPS", counts_data$external_gene_name),]
counts_data<- counts_data[!grepl("RNA", counts_data$external_gene_name),]
counts_data <-  counts_data[,-36]
colnames(counts_data)


# Filter for genes with a median >= 0
counts_data = subset(counts_data,apply(counts_data, 1, mean) >= 1)


str(coldata) #status and sex are character
coldata <- coldata %>%
  dplyr::mutate(
    status = factor(status, levels = c("control", "case")),
    sex = factor(sex, levels = c("male", "female"))
  )


# Add gene names to your counts data
genes <-read.delim("data/genes_annotation_final.txt", header = TRUE, sep = "\t", dec=".")
genes[(duplicated(genes$ensembl_gene_id_version) == TRUE),]
all(rownames(counts_data) %in% genes$ensembl_gene_id_version)
genes <- genes[genes$ensembl_gene_id_version %in% rownames(counts_data) ,]
all(rownames(counts_data) %in% genes$ensembl_gene_id_version)
test <- merge(counts_data, genes[1:2], by.x = "row.names", by.y= "ensembl_gene_id_version", all= FALSE)
row.names(test) <- test$Row.names
test <- test[, -1]
missing_names<- genes[genes$external_gene_name == "", ]
missing_names2 <- genes[is.na(genes$external_gene_name),]
duplicate_names<- test[duplicated(test$external_gene_name),] #160
duplicate_names$external_gene_name
non_duplicates_idx <- which(duplicated(test$external_gene_name) == FALSE)
test <- test[non_duplicates_idx, ]
duplicate_names<- test[duplicated(test$external_gene_name),] # 0
duplicate_names$external_gene_name
row.names(test) = test[,"external_gene_name"]

missing_names<- test[test$external_gene_name == "", ] # 1
missing_names2 <- test[is.na(test$external_gene_name),] #-
non_empty_idx <- which(test$external_gene_name != "")
test <- test[non_empty_idx,]
counts_data <- test
counts_data <- counts_data[,-36]

# Convert to matrix
counts_data = as.matrix(counts_data)

all(rownames(coldata) %in% colnames(counts_data))
all(rownames(coldata) == colnames(counts_data))
all(colnames(counts_data) %in% rownames(coldata))
all(colnames(counts_data) == rownames(coldata))


####  Construct a DESeqDataSet ####
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = coldata,
                              design = ~ status)
dds
mcols(dds)

colSums(counts_data)
counts_table <- as.data.frame(colSums(counts_data))
colnames(counts_table) <- ("total counts")
summary(counts_table)

# Make euclidean heatmap for identifying outliers
vsd <- vst(dds, blind = TRUE)
vsd_mat <- assay(vsd)
sampleDists <- dist(t(vsd_mat), method = "euclidean")
dm <- as.matrix(sampleDists)
rownames(dm) <- colnames(vsd_mat)
colnames(dm) <- colnames(vsd_mat)

df <- data.frame(status = colData(dds)$status, row.names = colnames(vsd_mat))

ha <- HeatmapAnnotation(
  df = df,
  col = list(status = c("control" = "salmon",
                        "case" = "turquoise")),
  simple_anno_size = unit(0.7, "cm"),
  show_legend = TRUE,
  annotation_legend_param = list(
    status = list(
      title = "TB status",
      title_gp = gpar(fontsize = 16, fontface = "bold"),
      labels_gp = gpar(fontsize = 12)
    )
  ),
  annotation_name_gp = gpar(fontsize = 0)
)


pheatmap_cols <- colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(255)
hc <- hclust(sampleDists, method = "complete")
ord <- hc$order


euclid_ht_no <- Heatmap(
  dm,
  name = "Distance",
  col = pheatmap_cols,
  column_title = "Sample-to-sample distance heatmap",
  column_title_gp = gpar(fontsize = 20, fontface = "bold"),
  show_row_names = TRUE,
  show_column_names = TRUE,
  top_annotation = ha,
  show_column_dend = FALSE,
  show_row_dend = TRUE,
  column_dend_height = unit(1, "cm"),
  column_names_side = "top",
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  row_order = ord,
  column_order = ord,
  heatmap_legend_param = list(
    title = "",
    title_gp = gpar(fontsize = 16, fontface = "bold"),
    labels_gp = gpar(fontsize = 14)
  )
)

euclid_ht_no


#### Preparing blood data  ####

coldata<- read.csv("data/full_metadata.csv", header = TRUE, row.names = 4, sep = ",")
coldata <- coldata[,-1]
coldata <- coldata[coldata$type == "blood",] #subset blood only
row.names(coldata)
row_order <- paste0("S", 36:75)
coldata <- coldata[row_order,]

counts_data <-read.delim("data/rnaseq_counts_protein_coding_final.txt", header = TRUE, sep = "\t", dec=".")
row.names(counts_data) <- counts_data$ID
counts_data <- counts_data[, -1]
colnames(counts_data)
counts_data <- counts_data[colnames(counts_data) %in% rownames(coldata)]
colnames(counts_data)
counts_data <- counts_data[, rownames(coldata)]
colnames(counts_data)

all(rownames(coldata) %in% colnames(counts_data))
all(rownames(coldata) == colnames(counts_data))
all(colnames(counts_data) %in% rownames(coldata))
all(colnames(counts_data) == rownames(coldata))

#### Remove non-useful genes  ####
genes <-read.delim("data/genes_annotation_final.txt", header = TRUE, sep = "\t", dec=".")
genes[(duplicated(genes$ensembl_gene_id_version) == TRUE),]
all(rownames(counts_data) %in% genes$ensembl_gene_id_version)
counts_data <- merge(counts_data, genes[1:2], by.x = "row.names", by.y= "ensembl_gene_id_version", all= FALSE)
row.names(counts_data) <- counts_data$Row.names
counts_data <- counts_data[, -1]

counts_data<- counts_data[!grepl("HB", counts_data$external_gene_name),]
counts_data<- counts_data[!grepl("^MT", counts_data$external_gene_name),]
counts_data<- counts_data[!grepl("RPL", counts_data$external_gene_name),]
counts_data<- counts_data[!grepl("RPS", counts_data$external_gene_name),]
counts_data<- counts_data[!grepl("RNA", counts_data$external_gene_name),]
counts_data <-  counts_data[,-41]
colnames(counts_data)


# Filter for genes with a median >= 0
counts_data = subset(counts_data,apply(counts_data, 1, mean) >= 1)
# 15946

str(coldata)
coldata <- coldata %>%
  dplyr::mutate(
    status = factor(status, levels = c("control", "case")),
    sex = factor(sex, levels = c("male", "female"))
  )

#### Add gene names to your counts data ####
# Add gene names to your counts data
genes <-read.delim("data/genes_annotation_final.txt", header = TRUE, sep = "\t", dec=".")
genes[(duplicated(genes$ensembl_gene_id_version) == TRUE),]
all(rownames(counts_data) %in% genes$ensembl_gene_id_version)
genes <- genes[genes$ensembl_gene_id_version %in% rownames(counts_data) ,]
all(rownames(counts_data) %in% genes$ensembl_gene_id_version)
test <- merge(counts_data, genes[1:2], by.x = "row.names", by.y= "ensembl_gene_id_version", all= FALSE)
row.names(test) <- test$Row.names
test <- test[, -1]
missing_names<- genes[genes$external_gene_name == "", ]
missing_names2 <- genes[is.na(genes$external_gene_name),]
duplicate_names<- test[duplicated(test$external_gene_name),] #160
duplicate_names$external_gene_name
non_duplicates_idx <- which(duplicated(test$external_gene_name) == FALSE)
test <- test[non_duplicates_idx, ]
duplicate_names<- test[duplicated(test$external_gene_name),] # 0
duplicate_names$external_gene_name
row.names(test) = test[,"external_gene_name"]

missing_names<- test[test$external_gene_name == "", ] # 1
missing_names2 <- test[is.na(test$external_gene_name),] #-
non_empty_idx <- which(test$external_gene_name != "")
test <- test[non_empty_idx,]
counts_data <- test
counts_data <- counts_data[,-41]

# Convert to matrix
counts_data = as.matrix(counts_data)

all(rownames(coldata) %in% colnames(counts_data))
all(rownames(coldata) == colnames(counts_data))
all(colnames(counts_data) %in% rownames(coldata))
all(colnames(counts_data) == rownames(coldata))


####  Construct a DESeqDataSet ####
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = coldata,
                              design = ~ status)
dds
mcols(dds)


colSums(counts_data)
counts_table <- as.data.frame(colSums(counts_data))
colnames(counts_table) <- ("total counts")
summary(counts_table)

# euclidean heatmap
vsd <- vst(dds, blind = TRUE)
vsd_mat <- assay(vsd)
sampleDists <- dist(t(vsd_mat), method = "euclidean")
dm <- as.matrix(sampleDists)
rownames(dm) <- colnames(vsd_mat)
colnames(dm) <- colnames(vsd_mat)


df <- data.frame(status = colData(dds)$status, row.names = colnames(vsd_mat))

ha <- HeatmapAnnotation(
  df = df,
  col = list(status = c("control" = "salmon",
                        "case" = "turquoise")),
  simple_anno_size = unit(0.7, "cm"),
  show_legend = TRUE,
  annotation_legend_param = list(
    status = list(
      title = "TB status",
      title_gp = gpar(fontsize = 16, fontface = "bold"),
      labels_gp = gpar(fontsize = 12)
    )
  ),
  annotation_name_gp = gpar(fontsize = 0)
)

pheatmap_cols <- colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(255)
hc <- hclust(sampleDists, method = "complete")
ord <- hc$order

euclid_ht_no <- Heatmap(
  dm,
  name = "Distance",
  col = pheatmap_cols,
  column_title = "Sample-to-sample distance heatmap",
  column_title_gp = gpar(fontsize = 20, fontface = "bold"),
  show_row_names = TRUE,
  show_column_names = TRUE,
  top_annotation = ha,
  show_column_dend = FALSE,
  show_row_dend = TRUE,
  column_dend_height = unit(1, "cm"),
  column_names_side = "top",
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  row_order = ord,
  column_order = ord,
  heatmap_legend_param = list(
    title = "",
    title_gp = gpar(fontsize = 16, fontface = "bold"),
    labels_gp = gpar(fontsize = 14)
  )
)

euclid_ht_no


#### PCA ####

#### NASAL PCA ####

#### Preparing data  ####

#Bring in your counts data and coldata
# coldata: 32 obs 19 variables
# counts : 16889 obs, 32 columns
#have been wrangled with following steps - removed low count genes, non useful genes, removed 3 outlier samples, have given them gene names
counts_data <- read.csv("data/nasalcounts.csv", header=TRUE, row.names = 1)
coldata <- read.csv("data/nasalcoldata.csv", header=TRUE, row.names = 1)
coldata <- coldata %>%
  dplyr::mutate(
    status = factor(status, levels = c("control", "case")),
    sex = factor(sex, levels = c("male", "female"))
  )
str(coldata)
dim(coldata)
dim(counts_data)

#  Now Construct a DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = coldata,
                              design = ~ status)

dds
mcols(dds)

# Make a VSD object using VST transformation
vsd <- vst(dds)
ntop <- 500
rv <- rowVars(assay(vsd))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
mat <- t( assay(vsd)[select, ] )
pca<-prcomp(mat)


#### biological variables ####
pc1_pc2_noeclip = function(colour_groups, e_data)
{
  ntop <- 500
  rv <- rowVars(assay(e_data))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  mat <- t( assay(e_data)[select, ] )
  pca<-prcomp(mat)
  pca_coordinates <-  as.data.frame(pca$x)

  vars = apply(pca$x, 2, var)
  prop_x = round(vars["PC1"] / sum(vars) * 100, 0)
  prop_y = round(vars["PC2"] / sum(vars) * 100, 0)
  x_axis_label = paste("PC1 (", prop_x, "%)", sep = "")
  y_axis_label = paste("PC2 (", prop_y, "%)", sep = "")

  ggp = ggplot(pca_coordinates, aes(x=PC1, y=PC2, colour = colour_groups)) +
    geom_point() +
    labs(title = "", x = x_axis_label , y= y_axis_label, colour = NULL) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.text = element_text(size = 16),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          plot.margin  = margin(2, 2, 2, 2)
    )

  return(ggp)
}


#### PCA function PCA3/4 ####
pc3_pc4_noeclip = function(colour_groups, e_data)
{
  ntop <- 500
  rv <- rowVars(assay(e_data))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  mat <- t( assay(e_data)[select, ] )
  pca<-prcomp(mat)
  pca_coordinates <-  as.data.frame(pca$x)

  vars = apply(pca$x, 2, var)
  prop_x = round(vars["PC3"] / sum(vars) * 100, 0)
  prop_y = round(vars["PC4"] / sum(vars) * 100, 0)
  x_axis_label = paste("PC3 (", prop_x, "%)", sep = "")
  y_axis_label = paste("PC4 (", prop_y, "%)", sep = "")

  ggp = ggplot(pca_coordinates, aes(x=PC3, y=PC4, colour = colour_groups)) +
    geom_point() +
    labs(title = "", x = x_axis_label , y= y_axis_label, colour = NULL) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.text = element_text(size = 16),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          plot.margin  = margin(2, 2, 2, 2)
    )

  return(ggp)
}


##### tb #####
n_pca_tb1_2 <- pc1_pc2_noeclip(coldata$status, vsd)
n_pca_tb3_4 <- pc3_pc4_noeclip(coldata$status, vsd)

##### sex #####
n_pca_gender1_2 <- pc1_pc2_noeclip(coldata$sex, vsd)
n_pca_gender3_4 <- pc3_pc4_noeclip(coldata$sex, vsd)

##### Hiv stage  #####
n_pca_stage1_2 <- pc1_pc2_noeclip(coldata$HIV_stage, vsd)
n_pca_stage3_4 <- pc3_pc4_noeclip(coldata$HIV_stage, vsd)


make_pc1_pc2_cont = function(colour_groups, e_data)
{
  ntop <- 500
  rv <- rowVars(assay(e_data))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  mat <- t( assay(e_data)[select, ] )
  pca<-prcomp(mat)
  pca_coordinates <-  as.data.frame(pca$x)

  vars = apply(pca$x, 2, var)
  prop_x = round(vars["PC1"] / sum(vars) * 100, 0)
  prop_y = round(vars["PC2"] / sum(vars) * 100, 0)
  x_axis_label = paste("PC1 (", prop_x, "%)", sep = "")
  y_axis_label = paste("PC2 (", prop_y, "%)", sep = "")

  ggp = ggplot(pca_coordinates, aes(x=PC1, y=PC2, colour = colour_groups)) +
    geom_point(size = 2) +
    scale_colour_gradientn(
      colours = c("#1b7837", "#92c5de", "#762a83")
    ) +
    labs(title = "", x = x_axis_label , y= y_axis_label, colour = NULL) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.text = element_text(size = 16),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          plot.margin  = margin(2, 2, 2, 2)
    )

  return(ggp)
}

make_pc3_pc4_cont = function(colour_groups, e_data)
{
  ntop <- 500
  rv <- rowVars(assay(e_data))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  mat <- t( assay(e_data)[select, ] )
  pca<-prcomp(mat)
  pca_coordinates <-  as.data.frame(pca$x)

  vars = apply(pca$x, 2, var)
  prop_x = round(vars["PC3"] / sum(vars) * 100, 0)
  prop_y = round(vars["PC4"] / sum(vars) * 100, 0)
  x_axis_label = paste("PC3 (", prop_x, "%)", sep = "")
  y_axis_label = paste("PC4 (", prop_y, "%)", sep = "")

  ggp = ggplot(pca_coordinates, aes(x=PC3, y=PC4, colour = colour_groups)) +
    geom_point(size = 2) +
    scale_colour_gradientn(
      colours = c("#1b7837", "#92c5de", "#762a83")
    ) +
    labs(title = "", x = x_axis_label , y= y_axis_label, colour = NULL ) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.text = element_text(size = 16),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          plot.margin  = margin(2, 2, 2, 2)
    )

  return(ggp)
}

##### CD4 count  #####

n_pca_cd4cont1_2 <- make_pc1_pc2_cont(coldata$CD4_count, vsd)
n_pca_cd4cont3_4 <- make_pc3_pc4_cont(coldata$CD4_count, vsd)

##### age  #####

n_pca_agecont1_2 <- make_pc1_pc2_cont(coldata$age, vsd)
n_pca_agecont3_4 <- make_pc3_pc4_cont(coldata$age, vsd)

# sex has a strong visual relationship on PCA

# pca scores
pc_scores <- pca$x
coldata_subset <- coldata[, c(3, 7), drop = FALSE]
pc_scores <- merge(pc_scores, coldata_subset, by = 0, all.x = TRUE)
pc_scores <- pc_scores %>%
  as_tibble() %>%
  rename("sample" = "Row.names")

pc_scores <- pc_scores %>%
  relocate(c(status, sex), .after = sample)
pc_scores

pc_scores_control = pc_scores[pc_scores$status == "control",]
pc_scores_case = pc_scores[pc_scores$status == "case",]

run_pca_tests2 <- function(pc_scores, group_var) {
  pc_scores[[group_var]] <- as.factor(pc_scores[[group_var]])
  group1 <- pc_scores[pc_scores[[group_var]] == levels(pc_scores[[group_var]])[1],]
  group2 <- pc_scores[pc_scores[[group_var]] == levels(pc_scores[[group_var]])[2],]
  t_test_results <- c()
  wilcox_test_results <- c()
  for (i in 1:4) {
    pc_col <- paste0("PC", i)
    shapiro_group1 <- shapiro.test(group1[[pc_col]])$p.value
    shapiro_group2 <- shapiro.test(group2[[pc_col]])$p.value
    if (shapiro_group1 > 0.05 && shapiro_group2 > 0.05) {
      t_test_results[pc_col] <- t.test(group1[[pc_col]], group2[[pc_col]])$p.value
      wilcox_test_results[pc_col] <- NA  # NA because t-test is appropriate
    } else {
      wilcox_test_results[pc_col] <- wilcox.test(group1[[pc_col]], group2[[pc_col]])$p.value
      t_test_results[pc_col] <- NA  # NA because Wilcoxon test is appropriate
    }
  }

  list(
    t_test_p_values = t_test_results,
    wilcox_test_p_values = wilcox_test_results
  )
}

# Run the tests for TB status (main factor of interest) and sex (strong visual relationship on P3 and P4)
pca_tb_results <- run_pca_tests2(pc_scores, "status")
pca_tb_results
pca_sex_results <- run_pca_tests2(pc_scores, "sex")
pca_sex_results

# Combine the results into a table
p_values_table <- data.frame(
  PC = paste0("PC", 1:4),
  TB_t_test = unlist(pca_tb_results$t_test_p_values),
  TB_wilcox_test = unlist(pca_tb_results$wilcox_test_p_values),
  Sex_t_test = unlist(pca_sex_results$t_test_p_values),
  Sex_wilcox_test = unlist(pca_sex_results$wilcox_test_p_values)
)

# Print the table
print(p_values_table)


#### PCA for panelled figs ####

#### TB status ####
make_pc3_pc4 = function(colour_groups, e_data)
{
  ntop <- 500
  rv <- rowVars(assay(e_data))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  mat <- t( assay(e_data)[select, ] )
  pca<-prcomp(mat)
  pca_coordinates <-  as.data.frame(pca$x)
  vars = apply(pca$x, 2, var)
  prop_x = round(vars["PC3"] / sum(vars) * 100, 0)
  prop_y = round(vars["PC4"] / sum(vars) * 100, 0)
  x_axis_label = paste("PC3 (", prop_x, "%)", sep = "")
  y_axis_label = paste("PC4 (", prop_y, "%)", sep = "")

  ggp = ggplot(pca_coordinates, aes(x=PC3, y=PC4, colour = colour_groups)) +
    geom_point() +
    labs(title = "Nasal", x = x_axis_label , y= y_axis_label, colour = NULL ) +
    stat_ellipse() +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.text = element_text(size = 12)
    )

  return(ggp)
}

n_pca_p3and4 <- make_pc3_pc4(coldata$status, vsd)
n_bp_pc4 <- ggboxplot(pc_scores,x="status",y="PC4", color="status",
                      palette=c("#00BFC4","#F8766D"),
                      order=c("case","control"),
                      add="jitter",xlab=F,outlier.shape=NA)

n_bp_pc4<-n_bp_pc4+stat_compare_means(method = "t.test") +
  theme(legend.position = "none", text = element_text(size =12))

n_bp_pc4



#### BLOOD PCA ####

#Bring in your counts data and coldata
# coldata_blood: 40 obs 19 variables
# counts_blood : 15787 obs, 40 columns
#have been wrangled with following steps - removed low count genes, non useful genes, removed 3 outlier samples, have given them gene names
counts_data_blood <- read.csv("data/bloodcounts.csv", header=TRUE, row.names = 1)
coldata_blood <- read.csv("data/bloodcoldata.csv", header=TRUE, row.names = 1)
coldata_blood <- coldata_blood %>%
  dplyr::mutate(
    status = factor(status, levels = c("control", "case")),
    sex = factor(sex, levels = c("male", "female"))
  )
str(coldata_blood)
dim(coldata_blood)
dim(counts_data_blood)

dds <- DESeqDataSetFromMatrix(countData = counts_data_blood,
                              colData = coldata_blood,
                              design = ~ status)

dds
mcols(dds)

# Make a VSD object using VST transformation
vsd <- vst(dds)#
rv <- rowVars(assay(vsd))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
mat <- t( assay(vsd)[select, ] )
pca<-prcomp(mat)




#### biological variables- blood  ####
##### tb status #####
b_pca_tb1_2 <- pc1_pc2_noeclip(coldata_blood$status, vsd)
b_pca_tb3_4 <- pc3_pc4_noeclip(coldata_blood$status, vsd)


##### sex #####
b_pca_gender1_2 <- pc1_pc2_noeclip(coldata_blood$sex, vsd)
b_pca_gender3_4 <- pc3_pc4_noeclip(coldata_blood$sex, vsd)


##### Hiv stage  #####
b_pca_stage1_2 <- pc1_pc2_noeclip(coldata_blood$HIV_stage, vsd)
b_pca_stage3_4 <- pc3_pc4_noeclip(coldata_blood$HIV_stage, vsd)


##### CD4 count  #####

b_pca_cd4cont1_2 <- make_pc1_pc2_cont(coldata_blood$CD4_count, vsd)
b_pca_cd4cont3_4 <- make_pc3_pc4_cont(coldata_blood$CD4_count, vsd)

##### AGE  #####

b_pca_agecont1_2 <- make_pc1_pc2_cont(coldata_blood$age, vsd)
b_pca_agecont3_4 <- make_pc3_pc4_cont(coldata_blood$age, vsd)

pc_scores <- pca$x
coldata_blood_subset <- coldata_blood[, c(3, 7), drop = FALSE]
pc_scores <- merge(pc_scores, coldata_blood_subset, by = 0, all.x = TRUE)
pc_scores <- pc_scores %>%
  as_tibble() %>%
  rename("sample" = "Row.names")

pc_scores <- pc_scores %>%
  relocate(c(status, sex), .after = sample)
pc_scores

pc_scores_control = pc_scores[pc_scores$status == "control",]
pc_scores_case = pc_scores[pc_scores$status == "case",]

# Run the tests for TB status (main factor of interest) and sex (strong relationship on P3 and P4)
pca_tb_results <- run_pca_tests2(pc_scores, "status")
pca_tb_results
pca_sex_results <- run_pca_tests2(pc_scores, "sex")
pca_sex_results

# Combine the results into a table
p_values_table <- data.frame(
  PC = paste0("PC", 1:4),
  TB_t_test = unlist(pca_tb_results$t_test_p_values),
  TB_wilcox_test = unlist(pca_tb_results$wilcox_test_p_values),
  Sex_t_test = unlist(pca_sex_results$t_test_p_values),
  Sex_wilcox_test = unlist(pca_sex_results$wilcox_test_p_values)
)

# Print the table
print(p_values_table)


#### PCA for panelled figs - blood ####

make_pc1_pc2 = function(colour_groups, e_data)
{
  ntop <- 500
  rv <- rowVars(assay(e_data))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  mat <- t( assay(e_data)[select, ] )
  pca<-prcomp(mat)
  pca_coordinates <-  as.data.frame(pca$x)

  vars = apply(pca$x, 2, var)
  prop_x = round(vars["PC1"] / sum(vars) * 100, 0)
  prop_y = round(vars["PC2"] / sum(vars) * 100, 0)
  x_axis_label = paste("PC1 (", prop_x, "%)", sep = "")
  y_axis_label = paste("PC2 (", prop_y, "%)", sep = "")

  ggp = ggplot(pca_coordinates, aes(x=PC1, y=PC2, colour = colour_groups)) +
    geom_point() +
    labs(title = "Blood", x = x_axis_label , y= y_axis_label, colour = NULL) +
    stat_ellipse() +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.text = element_text(size = 12)
    )

  return(ggp)
}


b_pca_pc1and2 <- make_pc1_pc2(coldata_blood$status, vsd)
b_bp_pc1 <- ggboxplot(pc_scores,x="status",y="PC1", color="status",
                      palette=c("#00BFC4","#F8766D"),
                      order=c("case","control"),
                      add="jitter",xlab=F,outlier.shape=NA)

b_bp_pc1 <-b_bp_pc1+stat_compare_means(method = "t.test", label.y = max(pc_scores$PC1) + 2) +
  theme(legend.position = "none", text = element_text(size =12))


#### PANELLED FIGURES  ####

# TB status
n_pca_p3and4
n_bp_pc4
b_pca_pc1and2
b_bp_pc1


Fig1 <- cowplot::plot_grid(n_pca_p3and4,
                           b_pca_pc1and2,
                           n_bp_pc4,
                           b_bp_pc1,
                           ncol=2, labels=LETTERS[1:2:3:4])

Fig1


# PCA of ALL NON TB variables

# nasal
n_pca_gender1_2
n_pca_stage1_2
n_pca_cd4cont1_2
n_pca_agecont1_2

n_pca_gender3_4
n_pca_stage3_4
n_pca_cd4cont3_4
n_pca_agecont3_4

stack_pair_no_legend <- function(p_top, p_bottom, align = "hv") {
  p_top_nl <- p_top + theme(legend.position = "none")
  p_bottom_nl <- p_bottom + theme(legend.position = "none")
  pair <- cowplot::plot_grid(p_top_nl, p_bottom_nl, ncol = 1, align = align)
  return(pair)
}

pair_gender  <- stack_pair_no_legend(n_pca_gender1_2,  n_pca_gender3_4)
pair_stage  <- stack_pair_no_legend(n_pca_stage1_2,  n_pca_stage3_4)
pair_cd4cont  <- stack_pair_no_legend(n_pca_cd4cont1_2,  n_pca_cd4cont3_4)
pair_age  <- stack_pair_no_legend(n_pca_agecont1_2,  n_pca_agecont3_4)


clinical <- cowplot::plot_grid(
  pair_gender, pair_stage, pair_cd4cont, pair_age,
  ncol = 4, align = "hv"
)


title_row <- ggdraw() + draw_label("PCA for nasal variables", size = 20, x = 0.5, hjust = 0.5)


sub1 <- ggdraw() + draw_label("Sex",                fontface = "plain", size = 14, x = 0.5, hjust = 0.5)
sub2 <- ggdraw() + draw_label("HIV stage",                fontface = "plain", size = 14, x = 0.5, hjust = 0.5)
sub3 <- ggdraw() + draw_label("CD4 count",                fontface = "plain", size = 14, x = 0.5, hjust = 0.5)
sub4 <- ggdraw() + draw_label("Age",         fontface = "plain", size = 14, x = 0.5, hjust = 0.5)


subtitle_row <- cowplot::plot_grid(sub1, sub2, sub3, sub4, ncol = 4)

clinical_with_title <- cowplot::plot_grid(
  title_row,
  subtitle_row,
  clinical,
  ncol = 1,
  rel_heights = c(0.08, 0.05, 1)
)


clinical_with_title
dev.off()


#blood uncorrected

b_pca_gender1_2
b_pca_stage1_2
b_pca_cd4cont1_2
b_pca_agecont1_2

b_pca_gender3_4
b_pca_stage3_4
b_pca_cd4cont3_4
b_pca_agecont3_4

pair_gender  <- stack_pair_no_legend(b_pca_gender1_2,  b_pca_gender3_4)
pair_stage  <- stack_pair_no_legend(b_pca_stage1_2,  b_pca_stage3_4)
pair_cd4cont  <- stack_pair_no_legend(b_pca_cd4cont1_2,  b_pca_cd4cont3_4)
pair_age  <- stack_pair_no_legend(b_pca_agecont1_2,  b_pca_agecont3_4)

clinical <- cowplot::plot_grid(
  pair_gender, pair_stage, pair_cd4cont, pair_age,
  ncol = 4, align = "hv"
)

title_row <- ggdraw() + draw_label("PCA for blood variables", size = 20, x = 0.5, hjust = 0.5)


sub1 <- ggdraw() + draw_label("Sex",                fontface = "plain", size = 14, x = 0.5, hjust = 0.5)
sub2 <- ggdraw() + draw_label("HIV stage",                fontface = "plain", size = 14, x = 0.5, hjust = 0.5)
sub3 <- ggdraw() + draw_label("CD4 count",                fontface = "plain", size = 14, x = 0.5, hjust = 0.5)
sub4 <- ggdraw() + draw_label("Age",         fontface = "plain", size = 14, x = 0.5, hjust = 0.5)



subtitle_row <- cowplot::plot_grid(sub1, sub2, sub3, sub4, ncol = 4)

clinical_with_title_blood <- cowplot::plot_grid(
  title_row,
  subtitle_row,
  clinical,
  ncol = 1,
  rel_heights = c(0.08, 0.05, 1)
)

clinical_with_title_blood
dev.off()


#### NASAL PCA post sex correction ####

coldata
counts_data
str(coldata)
batch = factor(coldata$sex)
sample_group = factor(coldata$status)
str(coldata)
counts_corrected = ComBat_seq(as.matrix(counts_data), batch=batch, group=sample_group)


dds <- DESeqDataSetFromMatrix(countData = counts_corrected,
                              colData = coldata,
                              design = ~ status)
dds
mcols(dds)
vsd <- vst(dds)
ntop <- 500
rv <- rowVars(assay(vsd))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
mat <- t( assay(vsd)[select, ] )
pca<-prcomp(mat)

# tb
pc1_pc2_noeclip(coldata$status, vsd)
pc3_pc4_noeclip(coldata$status, vsd)
# sex
pc1_pc2_noeclip(coldata$sex, vsd)
pc3_pc4_noeclip(coldata$sex, vsd)

# pca scores
pc_scores <- pca$x
coldata_subset <- coldata[, c(3, 7), drop = FALSE]
pc_scores <- merge(pc_scores, coldata_subset, by = 0, all.x = TRUE)
pc_scores <- pc_scores %>%
  as_tibble() %>%
  rename("sample" = "Row.names")
pc_scores <- pc_scores %>%
  relocate(c(status, sex), .after = sample)
pc_scores
pc_scores_control = pc_scores[pc_scores$status == "control",]
pc_scores_case = pc_scores[pc_scores$status == "case",]
pca_tb_results <- run_pca_tests2(pc_scores, "status")
pca_tb_results
pca_sex_results <- run_pca_tests2(pc_scores, "sex")
pca_sex_results
p_values_table <- data.frame(
  PC = paste0("PC", 1:4),
  TB_t_test = unlist(pca_tb_results$t_test_p_values),
  TB_wilcox_test = unlist(pca_tb_results$wilcox_test_p_values),
  Sex_t_test = unlist(pca_sex_results$t_test_p_values),
  Sex_wilcox_test = unlist(pca_sex_results$wilcox_test_p_values)
)

# Print the table
print(p_values_table)

#### PCA for panelled figs - nasal ####

n_pca_p3and4_corr <- make_pc3_pc4(coldata$status, vsd)
n_bp_pc3_corr <- ggboxplot(pc_scores,x="status",y="PC3", color="status",
                           palette=c("#00BFC4","#F8766D"),
                           order=c("case","control"),
                           add="jitter",xlab=F,outlier.shape=NA)

n_bp_pc3_corr<-n_bp_pc3_corr+stat_compare_means(method = "t.test", label.y = max(pc_scores$PC3) + 2.5) +
  theme(legend.position = "none", text = element_text(size =12))

n_bp_pc3_corr

#### Blood PCA post sex correction ####

coldata_blood
counts_data_blood
str(coldata_blood)
batch = factor(coldata_blood$sex)
sample_group = factor(coldata_blood$status)
str(coldata_blood)
counts_corrected_blood = ComBat_seq(as.matrix(counts_data_blood), batch=batch, group=sample_group)
dds <- DESeqDataSetFromMatrix(countData = counts_corrected_blood,
                              colData = coldata_blood,
                              design = ~ status)
dds
mcols(dds)
vsd <- vst(dds)
ntop <- 500
rv <- rowVars(assay(vsd))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
mat <- t( assay(vsd)[select, ] )
pca<-prcomp(mat)

# tb
pc1_pc2_noeclip(coldata_blood$status, vsd)
pc3_pc4_noeclip(coldata_blood$status, vsd)
# sex
pc1_pc2_noeclip(coldata_blood$sex, vsd)
pc3_pc4_noeclip(coldata_blood$sex, vsd)


# pca scores
pc_scores <- pca$x
coldata_blood_subset <- coldata_blood[, c(3, 7), drop = FALSE]
pc_scores <- merge(pc_scores, coldata_blood_subset, by = 0, all.x = TRUE)
pc_scores <- pc_scores %>%
  as_tibble() %>%
  rename("sample" = "Row.names")
pc_scores <- pc_scores %>%
  relocate(c(status, sex), .after = sample)
pc_scores

pc_scores_control = pc_scores[pc_scores$status == "control",]
pc_scores_case = pc_scores[pc_scores$status == "case",]
pca_tb_results <- run_pca_tests2(pc_scores, "status")
pca_tb_results
pca_sex_results <- run_pca_tests2(pc_scores, "sex")
pca_sex_results

p_values_table <- data.frame(
  PC = paste0("PC", 1:4),
  TB_t_test = unlist(pca_tb_results$t_test_p_values),
  TB_wilcox_test = unlist(pca_tb_results$wilcox_test_p_values),
  Sex_t_test = unlist(pca_sex_results$t_test_p_values),
  Sex_wilcox_test = unlist(pca_sex_results$wilcox_test_p_values)
)

print(p_values_table)


#### PCA for panelled figs - blood ####
b_pca_pc1and2_corr <- make_pc1_pc2(coldata_blood$status, vsd)
b_bp_pc1_corr <- ggboxplot(pc_scores,x="status",y="PC1", color="status",
                           palette=c("#00BFC4","#F8766D"),
                           order=c("case","control"),
                           add="jitter",xlab=F,outlier.shape=NA)

b_bp_pc1_corr<-b_bp_pc1_corr+stat_compare_means(method = "t.test", label.y = max(pc_scores$PC1) + 2.5) +
  theme(legend.position = "none", text = element_text(size =12))


#### PANELLED FIGURE ####

n_pca_p3and4_corr
n_bp_pc3_corr
b_pca_pc1and2_corr
b_bp_pc1_corr


Fig2 <- cowplot::plot_grid(n_pca_p3and4_corr,
                           b_pca_pc1and2_corr,
                           n_bp_pc3_corr,
                           b_bp_pc1_corr, ncol=2, labels=LETTERS[1:2:3:4])

