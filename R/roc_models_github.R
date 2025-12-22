# ROC and models

# This R script uses nasal and blood DEGs for training six machine learning
# models and calculates cross-validated AUC for predicting TB based on all 
# 44 nasal and 238 blood DEGs
# Compares nasal_44 and blood_238 missed TB cases for best performing model
# Derives nasal parsimonious gene signature based on variable importance (NASAL_4)
# and refits models

####   Libraries   ####
library(doParallel)
library(caret)
library(devEMF)
library(ggplot2)
library(MLeval)
library(tidyverse)
library(ggplot2)
library(DESeq2)
library(sva)
library(patchwork)
library(cowplot)
library(ranger)
library(glmnet)
library(kernlab)
library(pls)


#### Theme   ####
theme_SL2 <- function () {
  theme_bw() %+replace%
    theme(panel.grid = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      plot.background = element_blank(),
      legend.background = element_rect(fill="transparent", colour=NA),
      legend.key = element_rect(fill="transparent", colour=NA),
      legend.text=element_text(size=12, family="Arial", face="bold"),
      legend.title=element_blank())
}


#### NASAL ####

#### Modelling dataframe ####
counts_data <- read.csv("data/nasalcounts.csv", header=TRUE, row.names = 1)
coldata <- read.csv("data/nasalcoldata.csv", header=TRUE, row.names = 1)
coldata <- coldata %>%
  dplyr::mutate(
    status = factor(status, levels = c("control", "case")),
    sex = factor(sex, levels = c("male", "female")))
str(coldata)

# Convert to matrix
counts_data <- as.matrix(counts_data)
all(rownames(coldata) %in% colnames(counts_data))
all(rownames(coldata) == colnames(counts_data))
all(colnames(counts_data) %in% rownames(coldata))
all(colnames(counts_data) == rownames(coldata))

batch <- factor(coldata$sex)
sample_group <- factor(coldata$status)
counts_corrected <- ComBat_seq(as.matrix(counts_data),
    batch=batch,
    group=sample_group)


# DESeq2
dds <- DESeqDataSetFromMatrix(countData = counts_corrected,
                              colData = coldata,
                              design = ~ status)

levels(dds$status)
dds
mcols(dds)
vsd <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd)
dim(vst_mat)
head(vst_mat[, 1:5])
rownames(vst_mat)
colnames(vst_mat)
all(colnames(vst_mat) == rownames(coldata))

sig_genes <- readRDS("data/sig_genes_nasal.rds")
sig_genes_present <- intersect(sig_genes, rownames(vst_mat))
if(length(sig_genes_present) < length(sig_genes)) {
  warning("Some sig_genes not found in VST = matrix")
}

feats <- t(vst_mat[sig_genes_present, , drop = FALSE]) %>% as.data.frame()
head(feats)
feats$status <- factor(coldata[rownames(coldata), "status"])
colnames(feats)
head(feats)
stopifnot(!anyNA(feats$status))
levels(feats$status)
table(feats$status)
feats_nasal<- feats

#### Run models ALL DEG ####

#### RANGER ####
set.seed(42)
fit_control <- trainControl(method = "repeatedcv",
                            number = 10,
                            repeats = 10, 
                            summaryFunction = twoClassSummary,
                            savePredictions = TRUE,
                            classProbs = TRUE,
                            verboseIter = TRUE,
                            search = 'random')

set.seed(42)
r_fit <- caret::train(status~.,
                      data = feats,
                      method = 'ranger',
                      metric = 'ROC',
                      tuneLength = 15,
                      trControl = fit_control,
                      importance = 'permutation')


fm_model_r <- evalm(r_fit, gnames = 'random forest', plots = "r", fsize = 11)
ggr <- fm_model_r$roc + theme_SL2() + theme(legend.position = "bottom")
dev.off()
fm_model_r

r_fit_nasal44 <- r_fit

#### ELASTIC NET ####
set.seed(42)
glmnet_fit <- train(
  status ~ .,
  data = feats,
  method = "glmnet",
  metric = "ROC",
  trControl = fit_control,
  preProcess = c("center","scale"))


fm_model_e <- evalm(glmnet_fit, gnames = 'elastic net', plots = "r", fsize = 11)
gge <- fm_model_e$roc + theme_SL2() + theme(legend.position = "bottom")
dev.off()
fm_model_e

#### SVM radial ####

set.seed(42)
svmr_fit <- train(status ~ .,
  data = feats,
  method = "svmRadial",
  metric = "ROC",
  trControl = fit_control,
  tuneLength = 15,
  preProcess = c("center","scale")
)


fm_model_svmr <- evalm(svmr_fit, gnames = 'SVM (radial)', plots = "r", fsize = 11)
ggsvmr <- fm_model_svmr$roc + theme_SL2() + theme(legend.position = "bottom")
dev.off()
fm_model_svmr

#### SVM linear ####
set.seed(42)
svm_fit <- train(status ~ .,
  data = feats,
  method = "svmLinear",
  metric = "ROC",
  trControl = fit_control,
  tuneLength = 15,
  preProcess = c("center","scale")
)


fm_model_svm <- evalm(svm_fit, gnames='SVM (linear)', plots="r", fsize=11)
ggsvm <- fm_model_svm$roc + theme_SL2() + theme(legend.position = "bottom")
dev.off()
fm_model_svm

#### KNN ####

set.seed(42)
knn_fit <- train(status ~ .,
  data = feats,
  method = "knn",
  metric = "ROC",
  trControl = fit_control,
  tuneLength=15,
  preProcess = c("center","scale"))


fm_model_knn <- evalm(knn_fit, gnames='kNN', plots="r", fsize=11)
ggk <- fm_model_knn$roc + theme_SL2() + theme(legend.position = "bottom")
dev.off()
fm_model_knn

#### PLS####
set.seed(42)
pls_fit <- train(status ~ .,
  data = feats,
  method   = "pls",                 # PLS-DA via caret
  metric   = "ROC",
  trControl = fit_control,          # your 10x10 repeated CV
  preProcess = c("center","scale"))


fm_model_pls <- evalm(pls_fit, gnames='PLS', plots="r", fsize=11)
ggp <- fm_model_pls$roc + theme_SL2() + theme(legend.position = "bottom")
dev.off()
fm_model_pls

#all models
fm_model_r$roc
fm_model_e$roc
fm_model_svm$roc
fm_model_svmr$roc
fm_model_knn$roc
fm_model_pls$roc


all_plots <- list(ggr, gge, ggk, ggsvm, ggsvmr, ggp)


all_plots <- lapply(all_plots, function(p) {
  p + theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6)
  )
})

all_plots <- wrap_plots(all_plots, ncol = 3, nrow = 2)



####  VARIABLE IMPORTANCE  ####

plot_importance <- function(fit, model_name, top = 10) {
  imp <- varImp(fit)
  ggplot(imp, top = top) +
    theme_SL2() +
    ggtitle(paste("Top", top, "Genes by Importance in", model_name, "Model (Nasal)")) +
    xlab("Genes") +
    ylab("Variable Importance") +
    theme(axis.text.x = element_text(size = 16, angle = 45, hjust = 1))
}


# generate importance plots
nasal_rf    <- plot_importance(r_fit,       "Random Forest")
nasal_knn   <- plot_importance(knn_fit,     "kNN")
nasal_svmr  <- plot_importance(svmr_fit,    "SVM (Radial)")
nasal_svml  <- plot_importance(svm_fit,     "SVM (Linear)")
nasal_glm   <- plot_importance(glmnet_fit,  "Elastic Net")
nasal_pls   <- plot_importance(pls_fit,     "PLS")

nasal_rf
nasal_knn
nasal_svmr
nasal_svml
nasal_glm
nasal_pls

figure_importance_10 <- cowplot::plot_grid(nasal_rf, nasal_knn,
  nasal_svmr, nasal_svml,
  nasal_glm, nasal_pls,
  ncol = 2, labels = LETTERS[1:6]
)

figure_importance_10

nasal_knn   <- plot_importance(knn_fit,     "kNN and SVM")
figure_importance_10 <- cowplot::plot_grid(
  nasal_rf, nasal_knn,
  nasal_glm, nasal_pls,
  ncol = 2
)

#overlaps
importance_knn <- varImp(knn_fit,)
importance_svmr <- varImp(svmr_fit)
importance_svml <- varImp(svm_fit)
importance_pls <- varImp(pls_fit)
importance_glm <- varImp(glmnet_fit)
importance_rf <- varImp(r_fit)


importance_rf_df <- as.data.frame(importance_rf$importance)
sorted_importance_rf_df <- importance_rf_df[order(importance_rf_df$Overall,
    decreasing = TRUE), , drop = FALSE]


importance_knn_df <- as.data.frame(importance_knn$importance)
sorted_importance_knn_df <- importance_knn_df[order(importance_knn_df$control,
    decreasing = TRUE), , drop = FALSE]


importance_pls_df <- as.data.frame(importance_pls$importance)
sorted_importance_pls_df <- importance_pls_df[order(importance_pls_df$Overall,
    decreasing = TRUE), , drop = FALSE]


importance_glm_df <- as.data.frame(importance_glm$importance)
sorted_importance_glm_df <- importance_glm_df[order(importance_glm_df$Overall,
    decreasing = TRUE), , drop = FALSE]


# which overlap between machine learning methods ####
# top 10 important
top_predictive_genes_knn <- rownames(sorted_importance_knn_df)[1:10]
top_predictive_genes_rf <- rownames(sorted_importance_rf_df)[1:10]
top_predictive_genes_pls <- rownames(sorted_importance_pls_df)[1:10]
top_predictive_genes_glm <- rownames(sorted_importance_glm_df)[1:10]

overlapping_genes_ml <- Reduce(intersect, list(
  top_predictive_genes_knn,
  top_predictive_genes_rf,
  top_predictive_genes_pls,
  top_predictive_genes_glm
))
overlapping_genes_ml


#### Refit nasal models based on common-consensus variable importance genes (NASAL_4) ####

genes <- c("SPIB", "SHISA2", "TESPA1", "CD1B")

dim(vst_mat)
head(vst_mat[, 1:5])
rownames(vst_mat)
colnames(vst_mat)
all(colnames(vst_mat) == rownames(coldata))
sig_genes <- genes

sig_genes_present <- intersect(sig_genes, rownames(vst_mat))
if(length(sig_genes_present) < length(sig_genes)) {
  warning("Some sig_genes not found in VST = matrix")
}

feats <- t(vst_mat[sig_genes_present, , drop = FALSE]) %>% as.data.frame()
head(feats)
feats$status <- factor(coldata[rownames(coldata), "status"])
colnames(feats)
head(feats)
stopifnot(!anyNA(feats$status))
levels(feats$status)
table(feats$status)
genes



#### RANGER ####
set.seed(42)
fit_control <- trainControl(method="repeatedcv",
                            number=10,
                            repeats=10, 
                            summaryFunction=twoClassSummary,
                            savePredictions = TRUE,
                            classProbs = TRUE,
                            verboseIter = TRUE,
                            search = 'random')

set.seed(42)
r_fit <- caret::train(status~.,
                      data = feats,
                      method = 'ranger',
                      metric = 'ROC',
                      tuneLength = 15,
                      trControl = fit_control,
                      importance = 'permutation')

fm_model_r <- evalm(r_fit, gnames='random forest', plots="r", fsize=11)
ggr <- fm_model_r$roc + theme_SL2() + theme(legend.position = "bottom")
dev.off()
fm_model_r


#### ELASTIC NET ####

set.seed(42)
glmnet_fit <- train(
  status ~ .,
  data = feats,
  method = "glmnet",
  metric = "ROC",
  trControl = fit_control,
  preProcess = c("center","scale")
)

fm_model_e <- evalm(glmnet_fit, gnames='elastic net', plots="r", fsize=11)
gge <- fm_model_e$roc + theme_SL2() + theme(legend.position = "bottom")
dev.off()
fm_model_e

#### SVM radial ####

set.seed(42)
svmr_fit <- train(
  status ~ .,
  data = feats,
  method = "svmRadial",
  metric = "ROC",
  trControl = fit_control,
  tuneLength = 15,
  preProcess = c("center","scale")
)

fm_model_svmr <- evalm(svmr_fit, gnames='SVM (radial)', plots="r", fsize=11)
ggsvmr <- fm_model_svmr$roc + theme_SL2() + theme(legend.position = "bottom")
dev.off()
fm_model_svmr

#### SVM linear ####
set.seed(42)
svm_fit <- train(
  status ~ ., data = feats,
  method = "svmLinear",
  metric = "ROC",
  trControl = fit_control,
  tuneLength = 15,
  preProcess = c("center","scale")
)



fm_model_svm <- evalm(svm_fit, gnames='SVM (linear)', plots="r", fsize=11)
ggsvm <- fm_model_svm$roc + theme_SL2() + theme(legend.position = "bottom")
dev.off()
fm_model_svm

#### KNN ####

set.seed(42)
knn_fit <- train(
  status ~ ., data = feats,
  method = "knn",
  metric = "ROC",
  trControl = fit_control,
  tuneLength = 15,
  preProcess = c("center","scale"))


fm_model_knn <- evalm(knn_fit, gnames='kNN', plots="r", fsize=11)
ggk <- fm_model_knn$roc + theme_SL2() + theme(legend.position = "bottom")
dev.off()
fm_model_knn

#### PLS####
set.seed(42)
pls_fit <- train(
  status ~ ., data = feats,
  method   = "pls",
  metric   = "ROC",
  trControl = fit_control,
  preProcess = c("center","scale"))


fm_model_pls <- evalm(pls_fit, gnames='PLS', plots="r", fsize=11)
ggp <- fm_model_pls$roc + theme_SL2() + theme(legend.position = "bottom")
dev.off()
fm_model_pls

#all models
fm_model_r$roc
fm_model_e$roc
fm_model_svm$roc
fm_model_svmr$roc
fm_model_knn$roc
fm_model_pls$roc


#### Figures  ####

all_plots <- list(ggr, gge, ggk, ggsvm, ggsvmr, ggp)

all_plots <- lapply(all_plots, function(p) {
  p + theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6)
  )
})

all_plots <- wrap_plots(all_plots, ncol = 3, nrow = 2)


#### BLOOD ####

#### Modelling dataframe ####
counts_data_blood <- read.csv("data/bloodcounts.csv", header=TRUE, row.names = 1)
coldata_blood <- read.csv("data/bloodcoldata.csv", header=TRUE, row.names = 1)
coldata_blood <- coldata_blood %>%
  dplyr::mutate(
    status = factor(status, levels = c("control", "case")),
    sex = factor(sex, levels = c("male", "female"))
  )
str(coldata_blood)

counts_data_blood <- as.matrix(counts_data_blood)
all(rownames(coldata_blood) %in% colnames(counts_data_blood))
all(rownames(coldata_blood) == colnames(counts_data_blood))
all(colnames(counts_data_blood) %in% rownames(coldata_blood))
all(colnames(counts_data_blood) == rownames(coldata_blood))


batch <- factor(coldata_blood$sex)
sample_group <- factor(coldata_blood$status)
counts_corrected_blood <- ComBat_seq(as.matrix(counts_data_blood),
    batch=batch,
    group=sample_group)
dds <- DESeqDataSetFromMatrix(countData = counts_corrected_blood,
                              colData = coldata_blood,
                              design = ~ status)
levels(dds$status)
dds
mcols(dds)
vsd <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd)
dim(vst_mat)
head(vst_mat[, 1:5])
rownames(vst_mat)
colnames(vst_mat)
all(colnames(vst_mat) == rownames(coldata_blood))

sig_genes_blood <- readRDS("data/sig_genes_blood.rds")
sig_genes_blood_present <- intersect(sig_genes_blood, rownames(vst_mat))
if(length(sig_genes_blood_present) < length(sig_genes_blood)) {
  warning("Some sig_genes not found in VST = matrix")
}

feats <- t(vst_mat[sig_genes_blood_present, , drop = FALSE]) %>% as.data.frame()
head(feats)
feats$status <- factor(coldata_blood[rownames(coldata_blood), "status"])
colnames(feats)
head(feats)
stopifnot(!anyNA(feats$status))
levels(feats$status)
table(feats$status)
feats_blood <- feats

#### Run models ALL DEG blood ####

#### RANGER ####
set.seed(42)
fit_control <- trainControl(method="repeatedcv",
                            number=10,
                            repeats=10, 
                            summaryFunction=twoClassSummary,
                            savePredictions = TRUE,
                            classProbs = TRUE,
                            verboseIter = TRUE,
                            search = 'random')

set.seed(42)
r_fit <- caret::train(status~., data=feats,
                      method='ranger',
                      metric='ROC',
                      tuneLength=15,
                      trControl= fit_control,
                      importance = 'permutation')

fm_model_r <- evalm(r_fit, gnames='random forest', plots="r", fsize=11)
ggr = fm_model_r$roc + theme_SL2() + theme(legend.position = "bottom")
dev.off()
fm_model_r

r_fit_blood238 <- r_fit

#### ELASTIC NET ####

set.seed(42)
glmnet_fit <- train(
  status ~ ., data = feats,
  method = "glmnet",
  metric = "ROC",
  trControl = fit_control,
  preProcess = c("center","scale")
)


fm_model_e <- evalm(glmnet_fit, gnames='elastic net', plots="r", fsize=11)
gge = fm_model_e$roc + theme_SL2() + theme(legend.position = "bottom")
dev.off()
fm_model_e

#### SVM radial ####

set.seed(42)
svmr_fit <- train(
  status ~ ., data = feats,
  method = "svmRadial",
  metric = "ROC",
  trControl = fit_control,
  tuneLength=15,
  preProcess = c("center","scale")
)


fm_model_svmr <- evalm(svmr_fit, gnames='SVM (radial)', plots="r", fsize=11)
ggsvmr <- fm_model_svmr$roc + theme_SL2() + theme(legend.position = "bottom")
dev.off()
fm_model_svmr

#### SVM linear ####
set.seed(42)
svm_fit <- train(
  status ~ ., data = feats,
  method = "svmLinear",
  metric = "ROC",
  trControl = fit_control,
  tuneLength=15,
  preProcess = c("center","scale")
)

fm_model_svm <- evalm(svm_fit, gnames='SVM (linear)', plots="r", fsize=11)
ggsvm <- fm_model_svm$roc + theme_SL2() + theme(legend.position = "bottom")
dev.off()
fm_model_svm

#### KNN ####

set.seed(42)
knn_fit <- train(
  status ~ ., data = feats,
  method = "knn",
  metric = "ROC",
  trControl = fit_control,
  tuneLength=15,
  preProcess = c("center","scale"))


fm_model_knn <- evalm(knn_fit, gnames='kNN', plots="r", fsize=11)
ggk <- fm_model_knn$roc + theme_SL2() + theme(legend.position = "bottom")
dev.off()
fm_model_knn

#### PLS####

set.seed(42)
pls_fit <- train(
  status ~ ., data = feats,
  method   = "pls",                 # PLS-DA via caret
  metric   = "ROC",
  trControl = fit_control,
  preProcess = c("center","scale"),
  tuneGrid = data.frame(ncomp = 1:10)
)


fm_model_pls <- evalm(pls_fit, gnames='PLS', plots="r", fsize=11)
ggp <- fm_model_pls$roc + theme_SL2() + theme(legend.position = "bottom")

dev.off()
fm_model_pls

#all models
fm_model_r$roc
fm_model_e$roc
fm_model_svm$roc
fm_model_svmr$roc
fm_model_knn$roc
fm_model_pls$roc


#### Figures  ####
all_plots <- list(ggr, gge, ggk, ggsvm, ggsvmr, ggp)
all_plots <- lapply(all_plots, function(p) {
  p + theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6)
  )
})

all_plots <- wrap_plots(all_plots, ncol = 3, nrow = 2)



#### added value of combined signatures ####
# nasal_44
r_fit_nasal44
fm_model_r <- evalm(r_fit_nasal44, gnames='random forest', plots="r", fsize=11)
fm_model_r
mle_df <- fm_model_r$probs[[1]]
youden_row <- mle_df %>%
  slice_max(Informedness, n = 1, with_ties = FALSE)
thr_mleval <- youden_row$case
thr_mleval

youden_row %>% select(case, SENS, SPEC, Informedness)
fm_model_r$optres
thr <- thr_mleval
prob_col <- "case"

# extract caret CV predictions with best hyperparameters
preds <- r_fit_nasal44$pred
best <- r_fit_nasal44$bestTune
keep_idx <- Reduce(`&`, lapply(names(best),
    function(col) preds[[col]] == best[[col]]))
preds <- preds[keep_idx, ]
sample_names <- rownames(feats_nasal)
preds$sample_id <- sample_names[preds$rowIndex]

# CV probability per sample and get cm
agg <- preds %>%
  group_by(sample_id) %>%
  summarise(
    true_status = first(obs),
    mean_prob = mean(.data[[prob_col]], na.rm = TRUE)
  )

agg$pred_cv_youden <- ifelse(agg$mean_prob >= thr, "case", "control")
table(agg$pred_cv_youden, agg$true_status)

#samples case/control
samples_pred_case_CV <- agg$sample_id[agg$pred_cv_youden == "case"]
samples_pred_control_CV <- agg$sample_id[agg$pred_cv_youden == "control"]
samples_pred_case_CV
samples_pred_control_CV

case_rows <- coldata[
  rownames(coldata) %in% samples_pred_case_CV,
  c("sample_id", "status", "sample_status")
]

control_rows <- coldata[
  rownames(coldata) %in% samples_pred_control_CV,
  c("sample_id", "status", "sample_status")
]

case_rows_nasal44 <- case_rows
control_rows_nasal44 <- control_rows


#blood 238
r_fit_blood238
fm_model_r <- evalm(r_fit_blood238, gnames='random forest', plots="r", fsize=11)
fm_model_r
mle_df <- fm_model_r$probs[[1]]
youden_row <- mle_df %>%
  slice_max(Informedness, n = 1, with_ties = FALSE)
thr_mleval <- youden_row$case
thr_mleval


youden_row %>% select(case, SENS, SPEC, Informedness)
fm_model_r$optres
thr <- thr_mleval
prob_col <- "case"

# extract caret CV predictions with best hyperparameters
preds <- r_fit_blood238$pred
best <- r_fit_blood238$bestTune
keep_idx <- Reduce(`&`, lapply(names(best),
    function(col) preds[[col]] == best[[col]]))
preds <- preds[keep_idx, ]
sample_names <- rownames(feats_blood)
preds$sample_id <- sample_names[preds$rowIndex]

# CV probability per sample and get cm
agg <- preds %>%
  group_by(sample_id) %>%
  summarise(
    true_status = first(obs),
    mean_prob = mean(.data[[prob_col]], na.rm = TRUE)
  )

agg$pred_cv_youden <- ifelse(agg$mean_prob >= thr, "case", "control")
table(agg$pred_cv_youden, agg$true_status)

# which case/control
samples_pred_case_CV <- agg$sample_id[agg$pred_cv_youden == "case"]
samples_pred_control_CV <- agg$sample_id[agg$pred_cv_youden == "control"]
samples_pred_case_CV
samples_pred_control_CV

case_rows <- coldata_blood[
  rownames(coldata_blood) %in% samples_pred_case_CV,
  c("sample_id", "status", "sample_status")
]

control_rows <- coldata_blood[
  rownames(coldata_blood) %in% samples_pred_control_CV,
  c("sample_id", "status", "sample_status")
]

case_rows_blood238 <- case_rows
control_rows_blood238 <- control_rows

# compare nose and blood missed cases
case_rows_nasal44
control_rows_nasal44
case_rows_blood238
control_rows_blood238

