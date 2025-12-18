# IPA - canonical pathways and upstream regulators 

# This R script imports IPA results and creates barplots for the top 
# canonical pathways and upstream regulators for the nasal and blood
# datasets. 
# Also creates csv results files of IPA canonical pathways and upstream 
# regulators.


#### Load Libraries #### 

library(jamba)
library(readtext)
library(multienrichjam)
library(colorjam)
library(tidyverse)

####  Load Data ####

# path to find IPA outputs
# based on IPA analysis 
path.dir <- file.path("ipa_outs/ipa_export_all")
path <- path.dir
IPA <- list.files(path, pattern = ".txt")
# list of paths to each file
paths <- lapply(IPA, function(x){file.path(path, x)})

# import names of files and get names of comparisons as "tests"
tests <- gsub(".txt", "", IPA)


####   Import IPA results ####  

outs <- lapply(paths, importIPAenrichment,  headerGrep =
                 "(^|\t)((expr.|-log.|)p-value|Pvalue|Score\t|Symbol\t|Ratio\t|Consistency.Score|Master.Regulator\t)",
               ipaNameGrep = c("Pathway", "Regulator$", "Regulators", "Regulator", "Disease",
                               "Toxicity", "Category", "Categories", "Function", "Symbol$", "^ID$",
                               "My.(Lists|Pathways)"),
               geneGrep = c("Molecules in Network", "Molecules"),
               geneCurateFrom = c("[ ]*[(](complex|includes others)[)][ ]*", "^[, ]+|[, ]+$"),
               geneCurateTo = c("", ""),
               method = 1,
               sheet = 1,
               sep = "\t",
               xlsxMultiSheet = TRUE,
               useXlsxSheetNames = FALSE,
               remove_blank_colnames = F,
               verbose = FALSE)

names(outs) <- tests
outs$nasal <- NULL
tests <- c("blood","nasal_filtered")

####   Canonical Pathways ####  

## write entire upstream regulator csv
outs_path <- file.path("outs/ipa/results/core_analysis_outs")
dir.create(outs_path, recursive = T)
dir.create(file.path(outs_path, "canonical_pathways"))
dir.create(file.path(outs_path, "canonical_pathways", "csv_files"))

# extract canonical pathways from IPA output
cp <-  map(outs, ~.[["Canonical Pathways"]] )
names(cp) <- names(outs)

# file path for each file
csv_outs <- file.path(outs_path, "canonical_pathways", "csv_files")        
csvs <- lapply(tests, paste0, "_canonical_pathways.csv") %>% lapply(function(x) file.path(csv_outs, x))

# save files
mapply( write.csv, cp, csvs)

## subset by p-value and z-score
# p value and z-score to dplyr::filter
pval <- 0.05
z_score <- 1.5

# dplyr::filter #
cp_sub <- lapply(cp, dplyr::filter, abs(zScore) >= z_score & `P-value` <= pval) 
cp_sub <- lapply(cp_sub, function(x){arrange(x,  desc(`-log(p-value)`))})
cp_sub <- lapply(cp_sub, function(x) x %>% dplyr::rename(`-log10(pval)`=`-log(p-value)`))

## plot canonical pathways
# only include top 20 at the most
X <- lapply(cp_sub, top_n, n=20, wt=`-log10(pval)`)

g <- lapply(X, function(X){ggplot(data=X,  mapping = aes(x=reorder(`Ingenuity Canonical Pathways`, `-log10(pval)`, decreasing=F), y=`-log10(pval)`,fill=zScore)) +
    geom_bar(stat="identity", width= 0.5) +
    geom_col(position = position_dodge2(width = 0.5, preserve = "single")) +
    scale_fill_gradient2(low="blue", mid="white", high = "red") +
    geom_hline(yintercept=-log10(0.05), linetype=2)+
    coord_flip()+
    xlab("Canonical Pathways")})
lapply(g, print)


####  Upstream Regulators ####  

# extract canonical pathways from IPA output
ur <-  map(outs, ~.[["Upstream Regulators"]] )

names(ur) <- gsub("_all.txt", "", names(ur))
ur$sig_genes_Th1 <- NULL
ur <- lapply(ur, dplyr::mutate, `-log10(pval)`= -log10(`P-value`))

# file path for each file
path <- file.path(outs_path)
dir.create(file.path(path, "upstream_regulators"))
dir.create(file.path(path, "upstream_regulators", "csvs"))

paths_ur <- lapply(tests, paste0, "_upstream_regulators.csv")
paths_ur <- lapply(paths_ur, function(x){file.path(path, "upstream_regulators", "csvs", x)})

# save files
mapply( write.csv, ur, paths_ur)


## dplyr::filtering Upstream Regulators
# p value and z-score to dplyr::filter
pval <- 0.05
z_score <- 1.5

# dplyr::filter #
#ur_sub <- lapply(ur, dplyr::filter, "P-value" <= pval) 
ur <- map(ur, function(df) {
  if ("Activation z-score" %in% names(df)) {
    df <- dplyr::rename(df, zScore = `Activation z-score`)
  }
  df
})

ur <-  map(ur, function(df) {
  if ("neg_log_pval" %in% names(df)) {
    df <- dplyr::rename(df, `-log10(pval)`=neg_log_pval)
  }
  df
})


## plot filtered
z_score <- 1.5

# check if value present first
x <- lapply(ur, function(x){x %>% 
    mutate(zScore = coalesce(zScore, 0)) %>% 
    dplyr::filter(`-log10(pval)` >= -log10(0.05)) %>% 
    dplyr::filter(abs(zScore) >= z_score) %>% 
    arrange( desc(`-log10(pval)`)) %>%                             
    top_n(n=20, wt=`-log10(pval)`)})

g <- lapply(x, function(X){ggplot(data = X, mapping = aes(x=reorder(`Upstream Regulator`, (`-log10(pval)`)), y=`-log10(pval)`,fill=zScore)) +
    geom_bar(stat="identity") +
    geom_hline(yintercept=-log10(0.05), linetype=2)+
    scale_fill_gradient2(low="blue", mid="darkgrey", high = "red") +
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) +
    coord_flip()+
    xlab("Upstream Regulators")})

lapply(g, print)




