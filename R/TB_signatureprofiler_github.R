# TB signature profiler

# This R script evaluates the performance of previously published diagnostic TB blood signatures in our nasal  and blood datasets.
# Signature performance scores were calculated using gene set variation analysis (GSVA) and single-sample GSEA (ssGSEA)


#### Load libraries ####
library(ggplot2)
library(DESeq2)
library(tidyverse)
library(TBSignatureProfiler)
library(SummarizedExperiment)
library(pROC)
library(HGNChelper)


# update boxplot function
signatureBoxplot <- function (inputData, annotationData, signatureColNames,
    annotationColName, name = "Signatures", scale = FALSE, violinPlot = FALSE,
    includePoints = TRUE, notch = FALSE, rotateLabels = FALSE, nrow = NULL,
    ncol = NULL, fill_colors = NULL) {

  if (methods::is(inputData, "SummarizedExperiment")) {
    if (any(duplicated(signatureColNames))) {
      stop("Duplicate signature column name is not supported.")
    }
    if (!all(signatureColNames %in% colnames(SummarizedExperiment::colData(inputData)))) {
      stop("Signature column name not found in inputData.")
    }
    if (!all(annotationColName %in% colnames(SummarizedExperiment::colData(inputData)))) {
      stop("Annotation column name not found in inputData.")
    }
    annotationData <- data.frame(
        SummarizedExperiment::colData(inputData)[,annotationColName, drop = FALSE])
    inputData <- data.frame(
        SummarizedExperiment::colData(inputData)[,signatureColNames, drop = FALSE])
  }
  else {
    if (ncol(annotationData) != 1) {
      stop("annotationData must have only one column.")
    }
    annotationColName <- colnames(annotationData)
  }
  if (length(annotationColName) != 1) {
    stop("You must specify a single annotation column name to color boxplots by.")
  }
  if (!is.factor(annotationData[, 1])) {
    annotationData[, 1] <- as.factor(annotationData[, 1])
  }
  n <- length(levels(annotationData[, 1]))
  if (n > 9) {
    stop("Too many levels in the annotation data. The boxplot can contain a maximum of 9 levels")
  }
  if (nrow(annotationData) == nrow(inputData)) {
    if (!all(rownames(annotationData) == rownames(inputData))) {
      stop("Annotation data and signature data does not match.")
    }
  }
  else if (nrow(annotationData) == ncol(inputData)) {
    if (!all(rownames(annotationData) == colnames(inputData))) {
      stop("Annotation data and signature data does not match.")
    }
    inputData <- t(inputData)
  }
  else {
    stop("Annotation data and signature data does not match.")
  }
  pathwaydata <- t(inputData)
  if (scale) {
    pathwaydata <- t(scale(t(pathwaydata)))
  }
  boxplotdf <- data.frame(t(pathwaydata), Group = annotationData[,
                                                                 1])
  boxplotdfm <- reshape2::melt(boxplotdf, value.name = "Score",
                               variable.name = "Signature", id.vars = "Group")
  theplot <- ggplot2::ggplot(boxplotdfm, ggplot2::aes(Group, Score)) +
      ggplot2::facet_wrap(~Signature, scales = "free", nrow = nrow, ncol = ncol)
  if (violinPlot) {
    theplot <- theplot +
        ggplot2::geom_violin(ggplot2::aes(fill = Group)) +
        ggplot2::theme_classic()
  }
  else {
    theplot <- theplot +
        ggplot2::geom_boxplot(outlier.shape = NA,
            ggplot2::aes(fill = Group), notch = notch) +
        ggplot2::theme_classic()
  }
  if (includePoints) {
    theplot <- theplot +
        ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.1))
  }
  if (rotateLabels) {
    theplot <- theplot +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
  }
  # Check and apply custom colors
  if (!is.null(fill_colors)) {
    if (!all(levels(annotationData[, 1]) %in% names(fill_colors))) {
      stop("fill_colors must be a named vector with names matching the annotation levels.")
    }
    theplot <- theplot + ggplot2::scale_fill_manual(values = fill_colors)
  } else {
    if (n < 3) n <- 3
    fill_colors <- RColorBrewer::brewer.pal(n, "Set1")
    theplot <- theplot + ggplot2::scale_fill_manual(values = fill_colors)
  }

  return(theplot + ggplot2::ggtitle(name))
}

# Define custom colors
custom_colors <- c("case" = "#00BFC4", "control" = "#F8766D")


# Nasal #
counts_data <- read.csv("data/nasalcounts.csv", header=TRUE, row.names = 1)
coldata <- read.csv("data/nasalcoldata.csv", header=TRUE, row.names = 1)
coldata <- coldata[, c(2,3)]
colnames(coldata) <- c("sample", "sample_group")
coldata$sample <- row.names(coldata)
coldata <- coldata %>%
  dplyr::mutate(sample_group = factor(sample_group,
      levels = c("control", "case")))
str(coldata)

all(rownames(coldata) %in% colnames(counts_data))
all(rownames(coldata) == colnames(counts_data))
all(colnames(counts_data) %in% rownames(coldata))
all(colnames(counts_data) == rownames(coldata))
counts <- counts_data


#### Make summarised experiment ####
hivtb_data<- SummarizedExperiment(assays = list(counts = as.matrix(counts)),
    colData = coldata)
hivtb_data <- mkAssay(hivtb_data, log = TRUE, counts_to_CPM = TRUE)

#### Subset signatures for Diagnosis only ####
# which signatures in blood aim to diagnose TB, i.e not risk, failure, treatment response
sigAnnotData[sigAnnotData$disease %in% c("Disease", "OD", "HIV"), ]
diag_subset <- sigAnnotData[sigAnnotData$disease %in% c("Disease", "OD", "HIV"), ]
names_diag <- diag_subset$names

#### Run the TBSignatureProfiler: ssGSEA method ####
out <- capture.output(ssgsea_result <- runTBsigProfiler(input = hivtb_data,
                                                        useAssay = "log_counts_cpm",
                                                        signatures = TBsignatures,
                                                        algorithm = "ssGSEA",
                                                        combineSigAndAlgorithm = TRUE,
                                                        parallel.sz = 1))

# Remove any signatures that were not scored in RNA seq dataset
names_diag <- names_diag[!(names_diag %in% c("Chendi_HIV_2", "Hoang_OD_3"))]
TBsignatures_diag <- TBsignatures[names_diag]


# AUC table
set.seed(0)
table_AUC_diag_ssgsea59 <- tableAUC(ssgsea_result,
                                    annotationColName = "sample_group",
                                    signatureColNames = names(TBsignatures_diag),
                                    num.boot = 100,
                                    pb.show = FALSE,
                                    output = "data.frame")


# Which stat signif signatures with AUC >0.7
table_AUC_diag_ssgsea59 %>%  filter(AUC >= 0.70 & P.value <0.05)
sigs_ssgsea_seq <- table_AUC_diag_ssgsea59 %>%  filter(AUC >= 0.70 & P.value <0.05)


# Boxplot of stat signif with signatures AUC >0.7
sigs <- c("Maertzdorf_OD_100", "Sloot_HIV_2", "Lee_4")
set.seed(0)
nasal_ssGSEA_seq <- signatureBoxplot(
  inputData = ssgsea_result,
  name = "Significant Diagnostic Signatures in Nasal Dataset, ssGSEA",
  signatureColNames = sigs,
  annotationColName = "sample_group",
  rotateLabels = FALSE,
  fill_colors = custom_colors # Pass the custom colors
)

#### Run the TBSignatureProfiler: GSVA method ####
out <- capture.output(gsva_result <- runTBsigProfiler(input = hivtb_data,
                                                      useAssay = "log_counts_cpm",
                                                      signatures = TBsignatures,
                                                      algorithm = "GSVA",
                                                      combineSigAndAlgorithm = TRUE,
                                                      parallel.sz = 1))


# Remove any signatures that were not scored in RNA seq dataset
names_diag <- names_diag[!(names_diag %in% c("Chendi_HIV_2", "Hoang_OD_3"))]
TBsignatures_diag <- TBsignatures[names_diag]

# AUC table
set.seed(0)
table_AUC_diag_gsva59 <- tableAUC(gsva_result,
                                  annotationColName = "sample_group",
                                  signatureColNames = names(TBsignatures_diag),
                                  num.boot = 100,
                                  pb.show = FALSE,
                                  output = "data.frame")

# Which stat signif signatures AUC >0.7
table_AUC_diag_gsva59 %>%  filter(AUC >= 0.7 & P.value <0.05)
sigs_gsva_seq<- table_AUC_diag_gsva59 %>%  filter(AUC >= 0.7 & P.value <0.05)

# Boxplot of stat signif signatures AUC >0.7
sigs <- c("Natarajan_7", "Sloot_HIV_2", "Jacobsen_3", "Gjoen_10" )
set.seed(0)
nasal_gsva_seq <- signatureBoxplot(inputData = gsva_result,
    name = "Significant Diagnostic Signatures in Nasal Dataset, GSVA",
    signatureColNames = sigs,
    annotationColName = "sample_group",
    rotateLabels = FALSE,
    fill_colors = custom_colors # Pass the custom colors
)




# figures from nasal
nasal_gsva_seq
nasal_ssGSEA_seq

#### BLOOD  ####
counts_data_blood <- read.csv("data/bloodcounts.csv", header=TRUE, row.names = 1)
coldata_blood <- read.csv("data/bloodcoldata.csv", header=TRUE, row.names = 1)

coldata_blood <- coldata_blood[, c(2,3)]
colnames(coldata_blood) <- c("sample", "sample_group")
coldata_blood$sample <- row.names(coldata_blood)

coldata_blood <- coldata_blood %>%
  dplyr::mutate(sample_group = factor(sample_group,
      levels = c("control", "case")))
str(coldata_blood)

# final check
all(rownames(coldata_blood) %in% colnames(counts_data_blood))
all(rownames(coldata_blood) == colnames(counts_data_blood))
all(colnames(counts_data_blood) %in% rownames(coldata_blood))
all(colnames(counts_data_blood) == rownames(coldata_blood))

counts_blood <- counts_data_blood


#### Make summarised experiment ####

hivtb_data<- SummarizedExperiment(assays = list(counts = as.matrix(counts_blood)),
    colData = coldata_blood)
hivtb_data <- mkAssay(hivtb_data, log = TRUE, counts_to_CPM = TRUE)


#### Subset signatures for Diagnosis only ####
sigAnnotData[sigAnnotData$disease %in% c("Disease", "OD", "HIV"), ]
diag_subset <- sigAnnotData[sigAnnotData$disease %in% c("Disease", "OD", "HIV"), ]
names_diag <- diag_subset$names

#### Run the TBSignatureProfiler: ssGSEA method ####

out <- capture.output(ssgsea_result <- runTBsigProfiler(input = hivtb_data,
                                                        useAssay = "log_counts_cpm",
                                                        signatures = TBsignatures,
                                                        algorithm = "ssGSEA",
                                                        combineSigAndAlgorithm = TRUE,
                                                        parallel.sz = 1))

# Remove any signatures that were not scored in RNA seq dataset
names_diag <- names_diag[!(names_diag %in% c("Chendi_HIV_2", "Hoang_OD_3"))]
TBsignatures_diag <- TBsignatures[names_diag] # excludes signatures that were not scored
# 59 diagnostic signatures

# Get an AUC table
set.seed(0)
table_AUC_diag_ssgsea59_blood <- tableAUC(ssgsea_result,
                                          annotationColName = "sample_group",
                                          signatureColNames = names(TBsignatures_diag),
                                          num.boot = 100,
                                          pb.show = FALSE,
                                          output = "data.frame")

# Which stat signif signatures AUC >0.7
table_AUC_diag_ssgsea59_blood %>%  filter(AUC >= 0.70 & P.value <0.05)
sigs_ssgsea_seq_blood <- table_AUC_diag_ssgsea59_blood %>%  filter(AUC >= 0.70 & P.value <0.05)


#### Boxplot of stat signif signatures AUC >0.7
#  only show the top 10 for the plot
sigs_ssgsea_seq_10_blood <- table_AUC_diag_ssgsea59_blood %>%
  arrange(desc(AUC)) %>%
  filter(AUC > 0.7) %>%
  slice_head(n = 10) %>%
  dplyr::select(Signature)

custom_colors <- c("case" = "#00BFC4", "control" = "#F8766D")
set.seed(0)
blood_ssGSEA_seq <- signatureBoxplot(inputData = ssgsea_result,
    name = "Significant Diagnostic Signatures in Blood Dataset*, ssGSEA",
    signatureColNames = sigs_ssgsea_seq_10_blood$Signature,
    annotationColName = "sample_group",
    rotateLabels = FALSE,
    fill_colors = custom_colors)


#### Run the TBSignatureProfiler: GSVA method ####
out <- capture.output(gsva_result <- runTBsigProfiler(input = hivtb_data,
                                                      useAssay = "log_counts_cpm", #log counts cpm
                                                      signatures = TBsignatures,
                                                      algorithm = "GSVA",
                                                      combineSigAndAlgorithm = TRUE,
                                                      parallel.sz = 1))


# Remove any signatures that were not scored in RNA seq dataset
names_diag <- names_diag[!(names_diag %in% c("Chendi_HIV_2", "Hoang_OD_3"))]
TBsignatures_diag <- TBsignatures[names_diag]

# Get an AUC table
set.seed(0)
table_AUC_diag_gsva59_blood <- tableAUC(gsva_result,
                                        annotationColName = "sample_group",
                                        signatureColNames = names(TBsignatures_diag),
                                        num.boot = 100,
                                        pb.show = FALSE,
                                        output = "data.frame")


#### Which stat signif signatures AUC >0.7
table_AUC_diag_gsva59_blood%>%  filter(AUC >= 0.70 & P.value <0.05)
sigs_gsva_seq_blood<- table_AUC_diag_gsva59_blood%>%  filter(AUC >= 0.70 & P.value <0.05)



# Boxplot of stat signif signatures AUC >0.7
# only show the top 10 for the plot
sigs_gsva_seq_10_blood <- table_AUC_diag_gsva59_blood %>%
  arrange(desc(AUC)) %>%
  filter(AUC > 0.7) %>%
  slice_head(n = 10) %>%
  dplyr::select(Signature)

set.seed(0)
blood_gsva_seq <- signatureBoxplot(inputData = gsva_result,
    name = "Significant Diagnostic Signatures in Blood Dataset*, GSVA",
    signatureColNames = sigs_gsva_seq_10_blood$Signature,
    annotationColName = "sample_group",
    rotateLabels = FALSE,
    fill_colors = custom_colors)

blood_gsva_seq
blood_ssGSEA_seq

#### overlap blood and nose  ####

#what overlap between the sig nasal and sig blood for seq ssGSEA
sigs_ssgsea_seq_blood
sigs_ssgsea_seq
intersect(sigs_ssgsea_seq$Signature, sigs_ssgsea_seq_blood$Signature)

#what overlap between the sig nasal and sig blood for seq GSVA
sigs_gsva_seq_blood
sigs_gsva_seq
intersect(sigs_gsva_seq$Signature, sigs_gsva_seq_blood$Signature)

#### Final figures ####

nasal_gsva_seq
nasal_ssGSEA_seq
blood_gsva_seq
blood_ssGSEA_seq


legend_theme <- theme(
  legend.title = element_text(size = 16),  # Adjust legend title size
  legend.text = element_text(size = 14)    # Adjust legend text size
)

nasal_ssGSEA_seq2 <- nasal_ssGSEA_seq +
    legend_theme +
    theme(legend.position = "none")
blood_ssGSEA_seq2 <- blood_ssGSEA_seq +
    legend_theme +
    theme(legend.position = "none")
nasal_gsva_seq2 <- nasal_gsva_seq +
    legend_theme +
    theme(legend.position = "none")
blood_gsva_seq2 <- blood_gsva_seq +
    legend_theme +
    theme(legend.position = "none")


legend <- cowplot::get_legend(
  nasal_ssGSEA_seq2 + theme(legend.position = "right")
)


tbsp_nasal_blood <- cowplot::plot_grid(
  cowplot::plot_grid(nasal_ssGSEA_seq2,
      blood_ssGSEA_seq2,
      nasal_gsva_seq2,
      blood_gsva_seq2,
      ncol = 2,
      labels = LETTERS[1:4]),
  legend,
  ncol = 2,
  rel_widths = c(1.5, 0.2)
)

print(tbsp_nasal_blood)
