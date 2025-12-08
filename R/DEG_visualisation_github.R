# DEG and visualization

# This R script processes and wrangles counts data and coldata for nasal and blood samples
# Runs DESEQ2 to get DEG for nasal and blood samples
# Creates volcano plots, boxplots with jitter, heatmaps based on DEG, scatter plots for genes and GSEA enrichment plots


#### NASAL ####

#### Load libraries ####
library(ggplot2)
library(DESeq2)
library(tidyverse)
library(ggrepel)
library(ComplexHeatmap)
library(sva)
library(reshape2)
library(circlize)
library(cowplot)
library(amap)
library(org.Hs.eg.db)
library(fgsea)


#### Preparing data  ####
counts_data <- read.csv("data/nasalcounts.csv", header=TRUE, row.names = 1)
coldata <- read.csv("data/nasalcoldata.csv", header=TRUE, row.names = 1)
coldata <- coldata %>%
  dplyr::mutate(
    status = factor(status, levels = c("control", "case")),
    sex = factor(sex, levels = c("male", "female"))
  )
str(coldata)

counts_data = as.matrix(counts_data)
all(rownames(coldata) %in% colnames(counts_data))
all(rownames(coldata) == colnames(counts_data))
all(colnames(counts_data) %in% rownames(coldata))
all(colnames(counts_data) == rownames(coldata))


#sex correction
batch = factor(coldata$sex)
sample_group = factor(coldata$status)
counts_corrected = ComBat_seq(as.matrix(counts_data), batch=batch, group=sample_group)


#### DESeq2  ####
dds <- DESeqDataSetFromMatrix(countData = counts_corrected,
                              colData = coldata,
                              design = ~ status)
dds
mcols(dds)
dds$status <- relevel(dds$status, ref = "control")
levels(dds$status)
dds <- DESeq(dds)
res05 <- results(dds, alpha=0.05)
resOrdered_05 <- res05[order(res05$pvalue),]
de <- as.data.frame(resOrdered_05)
signif_de <- subset(de, padj < 0.05 & abs(log2FoldChange) >1.0)
sig_genes <- row.names(signif_de)
sig_genes
sig_genes_nasal <- sig_genes
table_counts_normalized <- counts(dds, normalized=TRUE)

#### Parse files ####
# normalised gene expression table
em = read.csv("data/em_nasal.csv", header=TRUE, row.names = 1)

# differential expression table
de = read.csv("data/de_nasal.csv", header=TRUE, row.names = 1)
res_nasal <- de
de <- de[,-c(1,3,4)]
colnames(de) <- c("log2fold", "p", "p.adj")

#sample info table
ss = read.csv("data/sample_info_nasal.csv", header =TRUE, row.names = 1)
ss <- ss[, c(2,3)]
colnames(ss) <- c("sample", "sample_group")
ss$sample <- row.names(ss)
str(ss)

master <- merge(em, de, by.x=0, by.y=0)
row.names(master) = master[,"Row.names"]
names(master)[1] = "gene_name"
master = na.omit(master)
em = na.omit(em)
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

#### Volcano plot ####
master_sig_up = subset(master, p.adj< 0.05 & log2fold >1)
master_sig_down = subset(master, p.adj< 0.05 & log2fold < -1)
de_sig_up_top5 <- head(subset(de, p.adj < 0.05 & log2fold > 1), 5)
de_sig_down_top5 <-head(subset(de, p.adj< 0.05 & log2fold < -1), 5)
master_sig_up_top5 <- master[master$gene_name %in% rownames(de_sig_up_top5), ]
master_sig_down_top5 <- master[master$gene_name %in% rownames(de_sig_down_top5), ]

nasalvp <- ggplot(master, aes(x=log2fold, y=mlog10p.adj)) +
  geom_point(aes(colour = "a")) +
  geom_point(data = master_sig_up, aes(colour = "b")) +
  geom_point(data = master_sig_down, aes(colour = "c")) +
  labs(title = "", x= "Log2 fold change", y= "-Log10 p") +
  theme_classic() +
  geom_vline(xintercept= -1, linetype="dashed", colour="grey", size=0.5) +
  geom_vline(xintercept= 1, linetype="dashed", colour="grey", size=0.5) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="grey", size=0.5) +
  geom_text_repel(data=master_sig_up_top5, aes(label=gene_name, colour= "b"), show.legend = FALSE) +
  geom_text_repel(data=master_sig_down_top5, aes(label=gene_name, colour= "c"), show.legend = FALSE) +
  scale_colour_manual(values = c("black", "red", "blue"), labels=c("No change", "up","down"), name="")


####  Boxplot with jitter ####

master10.s <- de[1:10,] #selection
candidate_genes <- row.names(master10.s)
em[candidate_genes, ]
gene_data <- em[candidate_genes, ]
gene_data <-  data.frame(t(gene_data))
gene_data$sample_group = ss$sample_group
gene_data.m = melt(gene_data, id.vars = "sample_group")
sum(gene_data.m$value == 0, na.rm = TRUE)
which(gene_data.m$value == 0)
table(gene_data.m$variable[gene_data.m$value == 0])
gene_data.m$value_adj <- gene_data.m$value + 1

nasalbp_dot <- ggplot(gene_data.m,
                      aes(x = variable, y = value_adj, fill = sample_group)) +
  geom_boxplot(width = 0.7, outlier.shape = NA) +
  geom_point(aes(group = sample_group),
             position = position_jitterdodge(jitter.width = 0.15,
                                             jitter.height = 0,
                                             dodge.width = 0.75),
             color = "black", size = 1, alpha = 0.8,
             show.legend = FALSE) +
  scale_fill_manual(values = c("#00BFC4", "#F8766D"),
                    labels = c("cases", "controls"), name = "") +
  scale_y_log10() +
  xlab("Genes") +
  ylab("Normalized counts (log10 scale)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#### Nasal Heatmaps ####

mat <- as.matrix(em_scaled_sig) # needs a numeric matrix
df <- data.frame(status = ss$sample_group)
col_fun <- colorRamp2(c(-2, 0, 2), c("#0000FF", "white", "#FF0000"))
col_fun(seq(-3, 3)) # Not sure about the importance

ha <- HeatmapAnnotation(
  df = df,
  col = list(status = c("case" = "#00BFC4", "control" = "#F8766D")),
  simple_anno_size = unit(0.7, "cm"),
  show_legend = TRUE,
  annotation_legend_param = list(
    status = list(title = "TB status", # changes legend name to TB status and
                  title_gp = gpar(fontsize = 16, fontface = "bold"),
                  labels_gp = gpar(fontsize = 14))
  ),
  annotation_name_gp = gpar(fontsize = 0)
)


nasal_heatmap_DEG <- Heatmap(mat,
                             name="Scaled value",
                             col=col_fun,
                             column_title = "Heatmap of nasal DEG",
                             column_title_gp = gpar(fontsize = 16, fontface = "bold"),
                             show_row_names = F, #no gene names
                             show_column_names = T,
                             top_annotation = ha,
                             show_column_dend = F, show_row_dend = F, column_dend_height = unit(2, "cm"),
                             column_names_side = "top",
                             row_names_gp = gpar(fontsize=6),
                             column_names_gp = gpar(fontsize=12),
                             heatmap_legend_param = list(
                               title = "Scaled value",
                               title_gp = gpar(fontsize = 16, fontface = "bold"),
                               labels_gp = gpar(fontsize = 14)
                             )
)


#### BLOOD ####
counts_data_blood <- read.csv("data/bloodcounts.csv", header=TRUE, row.names = 1)
coldata_blood <- read.csv("data/bloodcoldata.csv", header=TRUE, row.names = 1)
coldata_blood <- coldata_blood %>%
  dplyr::mutate(
    status = factor(status, levels = c("control", "case")),
    sex = factor(sex, levels = c("male", "female"))
  )
str(coldata_blood)

counts_data_blood = as.matrix(counts_data_blood)
all(rownames(coldata_blood) %in% colnames(counts_data_blood))
all(rownames(coldata_blood) == colnames(counts_data_blood))
all(colnames(counts_data_blood) %in% rownames(coldata_blood))
all(colnames(counts_data_blood) == rownames(coldata_blood))

#sex correction
batch = factor(coldata_blood$sex)
sample_group = factor(coldata_blood$status)
counts_corrected_blood = ComBat_seq(as.matrix(counts_data_blood), batch=batch, group=sample_group)


#### DESeq2  ####
dds <- DESeqDataSetFromMatrix(countData = counts_corrected_blood,
                              colData = coldata_blood,
                              design = ~ status)
dds$status <- relevel(dds$status, ref = "control")
levels(dds$status)
dds <- DESeq(dds)
res05 <- results(dds, alpha=0.05)
resOrdered_05 <- res05[order(res05$pvalue),]
de_blood <- as.data.frame(resOrdered_05)
signif_de_blood <- subset(de_blood, padj < 0.05 & abs(log2FoldChange) >1.0)
sig_genes_blood <- row.names(signif_de_blood)
table_counts_normalized_blood <- counts(dds, normalized=TRUE)



#### Parse files ####

# normalised gene expression table
em = read.csv("data/em_blood.csv", header=TRUE, row.names = 1)

# differential expression table
de = read.csv("data/de_blood.csv", header=TRUE, row.names = 1)
res_blood <- de_blood
de_blood <- de_blood[,-c(1,3,4)]
colnames(de_blood) <- c("log2fold", "p", "p.adj")

#sample info table
ss = read.csv("data/sample_info_blood.csv", header =TRUE, row.names = 1)
ss <- ss[, c(2,3)]
colnames(ss) <- c("sample", "sample_group")
ss$sample <- row.names(ss)

master <- merge(em, de_blood, by.x=0, by.y=0)
row.names(master) = master[,"Row.names"]
names(master)[1] = "gene_name"
nrow(master[master$p.adj == "NA",])
master = na.omit(master)
em = na.omit(em)
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

#### Volcano plot ####
master_sig_up = subset(master, p.adj< 0.05 & log2fold >1)
master_sig_down = subset(master, p.adj< 0.05 & log2fold < -1)
de_sig_up_top5 <- head(subset(de_blood, p.adj < 0.05 & log2fold > 1), 5)
de_sig_down_top5 <-head(subset(de_blood, p.adj< 0.05 & log2fold < -1), 5)
master_sig_up_top5 <- master[master$gene_name %in% rownames(de_sig_up_top5), ]
master_sig_down_top5 <- master[master$gene_name %in% rownames(de_sig_down_top5), ]

bloodvp <- ggplot(master, aes(x=log2fold, y=mlog10p.adj)) +
  geom_point(aes(colour = "a")) +
  geom_point(data = master_sig_up, aes(colour = "b")) +
  geom_point(data = master_sig_down, aes(colour = "c")) +
  labs(title = "", x= "Log2 fold change", y= "-Log10 p") +
  theme_classic() +
  geom_vline(xintercept= -1, linetype="dashed", colour="grey", size=0.5) +
  geom_vline(xintercept= 1, linetype="dashed", colour="grey", size=0.5) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="grey", size=0.5) +
  geom_text_repel(data=master_sig_up_top5, aes(label=gene_name, colour= "b"), show.legend = FALSE) +
  geom_text_repel(data=master_sig_down_top5, aes(label=gene_name, colour= "c"), show.legend = FALSE) +
  scale_colour_manual(values = c("black", "red", "blue"), labels=c("No change", "up","down"), name="")



####  Boxplot with jitter ####

master10.s <- de_blood[1:10,] #selection
candidate_genes <- row.names(master10.s)
em[candidate_genes, ]
gene_data <- em[candidate_genes, ]
gene_data <-  data.frame(t(gene_data))
gene_data$sample_group = ss$sample_group
gene_data.m = melt(gene_data, id.vars = "sample_group")
sum(gene_data.m$value == 0, na.rm = TRUE)
which(gene_data.m$value == 0)
table(gene_data.m$variable[gene_data.m$value == 0])

bloodbp_dot <- ggplot(gene_data.m,
                      aes(x=variable,y=value, fill=sample_group)) +
  geom_boxplot(width = 0.7, outlier.shape = NA) +
  geom_point(aes(group = sample_group),
             position = position_jitterdodge(jitter.width = 0.15,
                                             jitter.height = 0,
                                             dodge.width = 0.75),
             color = "black", size = 1, alpha = 0.8,
             show.legend = FALSE) +
  scale_fill_manual(values = c("#00BFC4", "#F8766D"),
                    labels=c("cases", "controls"), name ="") +
  scale_y_log10() +
  xlab("Genes") +
  ylab("Normalized counts (log10 scale)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#### Blood Heatmaps ####


mat <- as.matrix(em_scaled_sig)
df <- data.frame(status = ss$sample_group)
col_fun <- colorRamp2(c(-2, 0, 2), c("#0000FF", "white", "#FF0000"))
col_fun(seq(-3, 3)) # Not sure about the importance

ha <- HeatmapAnnotation(
  df = df,
  col = list(status = c("case" = "#00BFC4", "control" = "#F8766D")),
  simple_anno_size = unit(0.7, "cm"),
  show_legend = TRUE,
  annotation_legend_param = list(
    status = list(title = "TB status",
                  title_gp = gpar(fontsize = 16, fontface = "bold"),
                  labels_gp = gpar(fontsize = 14))
  ),
  annotation_name_gp = gpar(fontsize = 0)  # Hides the "status" label
)

blood_heatmap_DEG <- Heatmap(mat,
                             name="Scaled value",
                             col=col_fun,
                             column_title = "Heatmap of blood DEG",
                             column_title_gp = gpar(fontsize = 16, fontface = "bold"),
                             show_row_names = F, #no gene names
                             show_column_names = T,
                             top_annotation = ha,
                             show_column_dend = F, show_row_dend = F, column_dend_height = unit(2, "cm"),
                             column_names_side = "top",
                             row_names_gp = gpar(fontsize=6),
                             column_names_gp = gpar(fontsize=12),
                             heatmap_legend_param = list(
                               title = "Scaled value",
                               title_gp = gpar(fontsize = 16, fontface = "bold"),
                               labels_gp = gpar(fontsize = 14)
                             )
)


#### Comparing blood and nose ####

#### Panelled vp and bp  ####

nasalbp_dot
nasalvp
bloodbp_dot
bloodvp

VP_BP_dot <- cowplot::plot_grid(nasalvp, bloodvp, nasalbp_dot, bloodbp_dot, ncol=2, labels=LETTERS[1:2:3:4])
VP_BP_dot
title_row <- cowplot::plot_grid(
  cowplot::ggdraw() + cowplot::draw_label("Nasal", fontface = "bold", size = 16),
  cowplot::ggdraw() + cowplot::draw_label("Blood", fontface = "bold", size = 16),
  ncol = 2
)

VP_BP_titles_dot <- cowplot::plot_grid(
  title_row, VP_BP_dot,
  ncol = 1, rel_heights = c(0.05, 1)
)
VP_BP_titles_dot


#### Scatter plot ####

# here include all genes under evaluation and use colours to code for nasal deg, blood deg, and both deg.
# wrangle data
de$Gene <- rownames(de)
de_blood$Gene <- rownames(de_blood)
merged_data <- merge(de, de_blood, by = "Gene")

merged_data$label <- ifelse((merged_data$p.adj.x < 0.05 & abs(merged_data$log2fold.x) >1.0)  & (merged_data$p.adj.y < 0.05 & abs(merged_data$log2fold.y) >1.0),
                            as.character(merged_data$Gene), NA)

merged_data$label_blood_alone <- ifelse(
  (merged_data$p.adj.x >= 0.05 | abs(merged_data$log2fold.x) <= 1.0) &
    (merged_data$p.adj.y < 0.05 & abs(merged_data$log2fold.y) > 1.0),
  as.character(merged_data$Gene),
  NA
)

merged_data$label_nasal_alone <- ifelse(
  (merged_data$p.adj.y >= 0.05 | abs(merged_data$log2fold.y) <= 1.0) &
    (merged_data$p.adj.x < 0.05 & abs(merged_data$log2fold.x) > 1.0),
  as.character(merged_data$Gene),
  NA
)

merged_data <- merged_data %>%
  mutate(category = case_when(
    !is.na(label_blood_alone) ~ "Blood DEG Only",
    !is.na(label_nasal_alone) ~ "Nasal DEG Only",
    !is.na(label) ~ "DEG in both nose and blood",
    TRUE ~ "Not Significant"
  )) %>%
  mutate(category = factor(category,
                           levels = c("Not Significant",
                                      "Nasal DEG Only",
                                      "Blood DEG Only",
                                      "DEG in both nose and blood")))


agreement_plot <- ggplot(merged_data, aes(x = log2fold.x, y = log2fold.y, color = category)) +
  geom_point(
    data = subset(merged_data, category == "Not Significant"),
    size = 2.5, alpha = 0.6
  ) +
  geom_point(
    data = subset(merged_data, category == "Nasal DEG Only"),
    size = 2.8
  ) +
  geom_point(
    data = subset(merged_data, category == "Blood DEG Only"),
    size = 2.8
  ) +
  geom_point(
    data = subset(merged_data, category == "DEG in both nose and blood"),
    size = 3
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, color = "black", linetype = "solid") +
  geom_vline(xintercept = 0, color = "black", linetype = "solid") +
  geom_hline(yintercept = c(-1, 1), color = "black", linetype = "dotted") +
  geom_vline(xintercept = c(-1, 1), color = "black", linetype = "dotted") +
  geom_text_repel(
    data = subset(merged_data, category == "DEG in both nose and blood"),
    aes(label = Gene),
    show.legend = FALSE,
    max.overlaps = 30,
    size = 4
  ) +
  coord_cartesian(xlim = c(-5, 5), ylim = c(-3, 3)) +
  scale_color_manual(
    name = "",
    values = c(
      "Blood DEG Only" = "red",
      "Nasal DEG Only" = "blue",
      "DEG in both nose and blood" = "black",
      "Not Significant" = "grey70"
    )
  ) +
  theme_minimal() +
  labs(
    x = "log2 FC (Nasal)",
    y = "log2 FC (Blood)",
    title = "Gene Expression in Blood and Nasal Samples"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )


agreement_plot


#### Enrichment plots ####

#Nasal

#prepare ranked blood list

res_blood$symbol <- rownames(res_blood)
res_blood <- res_blood %>% relocate(symbol, .before = baseMean)
res <- res_blood
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

# prepare nasal deg gene set
nasal_genes <-sig_genes_nasal
head(names(ranks))
present_nasal <- intersect(nasal_genes, names(ranks))
missing_nasal <- setdiff(nasal_genes, names(ranks))

cat(length(present_nasal), "of", length(nasal_genes), "nasal genes found in 'ranks'.\n")
if(length(missing_nasal) > 0){
  cat("Missing genes (not found in ranks):", paste(missing_nasal, collapse = ", "), "\n")
}
nasal_pathway <- list(Nasal_44 = present_nasal)
set.seed(42)
fgsea_nasal <- fgsea::fgsea(pathways = nasal_pathway,
                            stats    = ranks,
                            minSize  = 15,
                            maxSize  = 500)


fgsea_nasal

enrich_nasal <- fgsea::plotEnrichment(nasal_pathway[[1]], ranks) +
  labs(title = "Enrichment plot: Nasal DEG",
       subtitle = paste0("NES = ", round(fgsea_nasal$NES, 3),
                         "   padj = ", signif(fgsea_nasal$padj, 3)))
enrich_nasal

#Blood

#prepare ranked nasal list

res_nasal$symbol <- rownames(res_nasal)
res_nasal<- res_nasal%>% relocate(symbol, .before = baseMean)
res <- res_nasal
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

# prep blood deg gene set
blood_genes <- sig_genes_blood
head(names(ranks))
present_blood <- intersect(blood_genes, names(ranks))
missing_blood <- setdiff(blood_genes, names(ranks))

cat(length(present_blood), "of", length(blood_genes), "blood genes found in 'ranks'.\n")
if(length(missing_blood) > 0){
  cat("Missing genes (not found in ranks):", paste(missing_blood, collapse = ", "), "\n")
}

blood_pathway <- list(Blood_238 = present_blood)
set.seed(42)
fgsea_blood <- fgsea::fgsea(pathways = blood_pathway,
                            stats    = ranks,
                            minSize  = 15,
                            maxSize  = 500)

fgsea_blood
enrich_blood  <- fgsea::plotEnrichment(blood_pathway[[1]], ranks) +
  labs(title = "Enrichment plot: Blood DEG",
       subtitle = paste0("NES = ", round(fgsea_blood$NES, 3),
                         "   padj = ", signif(fgsea_blood$padj, 3)))
enrich_blood

enrichmentplot <- cowplot::plot_grid(
  enrich_nasal,
  enrich_blood,
  ncol = 2,
  labels = LETTERS[1:2]
)



#### Figures @ end ####
nasalbp_dot
nasalvp
bloodbp_dot
bloodvp
nasal_heatmap_DEG
blood_heatmap_DEG
VP_BP_titles_dot
agreement_plot
enrichmentplot
