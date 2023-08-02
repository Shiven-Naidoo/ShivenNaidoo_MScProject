## DATA RETRIEVAL AND PRE-PROCESSING

## script to normalize, annotate and analyse microarray data sets ()

## check wd

getwd()

## set wd ("D:/Project data/MSc data/!!!example!!!")

## make sure to use /// slashes

#setwd("D:/Project data/MSc data/Microarray_GSE6257")

## download and call libraries needed for data downloading and normalisation

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.17")

BiocManager::install("affy")
BiocManager::install("GEOquery")
BiocManager::install("tidyverse", force = TRUE)

library(affy)
library(GEOquery)
library(tidyverse)
library(dplyr)


## load data
## First step:
## store sample_ID

sample_ID = "GSE2638"

## Now to download the data from GEOquery in the .cel file format and store the data in an object
## There are 3 ways to do this
## 1. download tar file from GEO and untar .cel files in a data/ folder in your wd
## 2. download cel files directly from GEO and store them in your data/ folder in wd
## 3. !!!!INCOMPLETE pipeline in progress :: store .cel file data directly from GEO into an object 

## 1. following method 1...
## download raw data from supp files

getGEOSuppFiles(sample_ID)

## files are compressed, decompress using untar()

untar("GSE2638/GSE2638_RAW.tar", exdir = "data/")

## Following Methods 1. and 2.
## .cel files should be stored in your data/ folder

## ensure that .cel files belong to the same experiment series
## usually only needed if experiment dataset consists of multiple series
## if multiple series of interest are present, apply pipeline to one series at a time

## store data from cel files into an object, make sure to store as a list

raw_data = ReadAffy(celfile.path = "data/")









## (if you are using methods 1 or 2, use the following to get sample info)
## identify useful column names using pData

series_data = getGEO(sample_ID)
series_data = series_data[[1]]

sample_info = pData(series_data)

## after storing pData in sample_info object, trim it down to keep only the important labels

## keep useful info in sample info by using select(sample_info, heading_of_interest_1, heading_of_interest_2)
## these columns amy contain useful information regarding cell type and preparation

sample_info = select(sample_info, title)

sample_info = rename(sample_info, Experiment = title)

## AUTOMATE column name editting









##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

## DATA NORMALIZATION

## normalise data: multiple methods to normalise microarray data, each with its own advantages/disadvantages
## this study will use RMA and the limma package to normalise and process data

## check datasets first to see if normalising methods have already be applied
## data processing should be specified in raw_data$phenoData

exprs(raw_data)

## first, normalise using rma

rma_normalised_data = rma(raw_data)

## rma: fetch normalised expression data and store as a dataframe

rma_normalised_expr = as.data.frame(exprs(rma_normalised_data))

## Rename column names to represent just the samples without the .cel extensions

names(rma_normalised_expr) <- sub(".CEL.gz", "", names(rma_normalised_expr))

summary(rma_normalised_expr)








##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

## DATA ANNOTATION

## attach annotation data to dataframe so that expression data can be linked to their associated gene
## right now dataframe contains probe IDs. map gene symbols to probe IDs (!!!!will be array and study specific)

## store gene symbols in dataframe/matrix to attach to data

GEO_series_data = getGEO(sample_ID, GSEMatrix = TRUE) # gets file for experiment series
GEO_anno_data = GEO_series_data$GSE2638_series_matrix.txt.gz@featureData@data# isolates data file which includes gene symbols in column x

GEO_gene_symbols = GEO_anno_data[, c(1, 11)]  # subsets columns with probe ID and corresponding gene symbols

## merge rma_normalised_expr df with GEO_gene_symbols df

rma_normalised_annotated_expr = rma_normalised_expr %>%
  rownames_to_column(var = "ID") %>%
  inner_join(., GEO_gene_symbols, by = "ID")










## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

## DATA QUALITY CONTROL
## load libraries

library(pheatmap)
library(ggplot2)
library(tidyr)
library(ggrepel)




## BOXPLOT
## visualize data after normalization (check if data is normally distributed)

summary(rma_normalised_expr)

boxplot(rma_normalised_expr, outline = F)





## HEATMAP

## condense normalised expression data into a corMatrix

corMatrix = cor(rma_normalised_expr)

## check heatmap (raw, no labels)

pheatmap(corMatrix)

## rename names in corMatrix to match sample info
## edit columns to have consistent, useful names and add columns for treated vs untreated


sample_info[5, 1] = "HMEC_TNF-treated 1"
sample_info[6, 1] = "HMEC_TNF-treated 2"
sample_info[7, 1] = "HMEC_TNF-treated 3"
sample_info[8, 1] = "HMEC_TNF-treated 4"

sample_info$Treatment = c("Untreated - Control", "Untreated - Control", "Untreated - Control", "Treated - TNF Alpha", "Treated - TNF Alpha", "Treated - TNF Alpha")

## rename column headers if needed using dplyr

# sample_info = rename(sample_info, Experiment = experiment, Treatment = treatment)

## check that rownames in corMatrix match rownames for annotations

rownames(sample_info)
colnames(corMatrix)

## adjust col_names and rownames for cormatrix if needed

colnames(corMatrix) = rownames(sample_info)

rownames(corMatrix) = rownames(sample_info)

## draw heat map

pheatmap(corMatrix, annotation_col = sample_info)






## PCA PLOT

## now to perform a principal component analysis PCA
## first, transpose the expression data using t() and use the prcomp function to store the principal components in an object

pca <- prcomp(t(rma_normalised_expr))

## Join the PCs to the sample information and visualize

cbind(sample_info, pca$x) %>% 
  ggplot(aes(x = PC1, y=PC2, col=Experiment,label=paste("", Experiment))) + geom_point() + geom_text_repel()

## visualise the proportions of each principal component
## first, get summary of PCA and look for columns/rows of interest

pca_sum = summary(pca)

## extract proportion of variance data

pca_proportions = pca_sum$importance["Proportion of Variance", ]

## convert proportions to a percentage out of 100

pca_proportions_percent = pca_proportions * 100

barplot(pca_proportions_percent, main = "Proportion of Variance explained by each Principal Component", xlab = "Principal Component", ylab = "Percentage of Variance")









### Handling outliers

## Option 1: remove outliers identified by PCA
## remove from raw data and redo all steps up to this point

#outliers <- c(outlier_name_1, outlier_name_2)

#gse <- gse[,-outlier_samples]

## Option 2: use the arrayWeights function to give outliers a weaker contribution during DGE
## Make sure to include aw as a parameter when calling the lmFit function

## ONLY USE arrayWeights ON PREPROCESSED DATA (eg., after RMA)

aw = arrayWeights(exprs(rma_normalised_data),design)












## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

## STATISTICAL SIGNIFICANCE AND FOLD CHANGE

## Determine the fold change and significance in expression between treatment and control groups

## for directly comparing expression between samples, isolate PXDN and PXDNL expr data from dataframe
## here GOI stands for Genes of Interest

GOI_expr_data <- rma_normalised_annotated_expr[rma_normalised_annotated_expr$`Gene Symbol` %in% c("PXDN", "PXDNL"), ]

write.table(PXDN_expr_data, file = "GOI_expr_GSE2639", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)


## use limma to visualize DE
## first create a design matrix based on whether a sample is treated or untreated
## design matrix assigns a 0 and 1 values depending on whether the label in sample info is treated/untreated

library(limma)
design <- model.matrix(~0+sample_info$Treatment)
design

## name the columns for the classes you are going to compare. Should match with your design matrix.

colnames(design) = c("TNF_Treated", "Untreated_Control")

## prep data

exprs_summary_table = summary(rma_normalised_expr)

## calculate median expression level, use as cutoff value for now
## The cutoff value chosen will be standardised as the median expression value for now
## This is done to reduce the noise in the data as well as reduce the ""multiple testing burden"
## If PXDN or PXDNL is lowly expressed in any of these samples, it will be filtered out

median_values = exprs_summary_table[3, ]

## PXDN had a relatively high expression as seen earlier in this sample, 
## Since it will not exclude PXDN therefore the median will be used as the cutoff
## PXDNL had a low expression and will be cut out of the data at this point

## use for loop to calculate median from summary(exprs)

# set counts for the number of medians and the average median
average_median = 0
sample_count = 0
for (char_string in median_values){
  # use gsub to remove spaces and characters from median info
  char_string = gsub("[a-zA-Z: ]", "", char_string)
  # only the median value should remain, use as.numeric to convert from string to float
  med_float = as.numeric(char_string)
  average_median = average_median + med_float
  sample_count = sample_count + 1}

cutoff = average_median/sample_count


## store genes considered to be expressed in object


is_expressed = rma_normalised_expr > cutoff

## Identify genes expressed in more than 2 samples

keep <- rowSums(is_expressed) > 2

## check how many genes are removed / retained.
table(keep)

## subset to just those expressed genes
rma_normalised_filtered_data <- rma_normalised_data[keep,]

## Fit data to design matrix

fit <- lmFit(exprs(rma_normalised_filtered_data), design, weights = aw)


head(fit$coefficients)


contrasts <- makeContrasts(TNF_Treated  - Untreated_Control, levels=design)

## can define multiple contrasts
## e.g. makeContrasts(Group1 - Group2, Group2 - Group3,....levels=design)

fit2 <- contrasts.fit(fit, contrasts)

fit2 <- eBayes(fit2)

## store results in an object

results = topTable(fit2, number = Inf)

## Annotation of Differential Gene Expression Analysis


## merge rma_normalised_expr df with GEO_gene_symbols df

annotated_results = results %>%
  rownames_to_column(var = "ID") %>%
  inner_join(., GEO_gene_symbols, by = "ID")

## Volcano plot

## set cutoff values for significant gene expression

p_cutoff = 0.05
fc_cutoff = 0.5

## OPTIONAL set number of differentially expressed genes which you would like to show on the plot

topDEGs = 20

## Use annotated results to create volcano plot

annotated_results %>% 
  mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
  mutate(Rank = 1:n(), Label = ifelse(Rank < topDEGs, `Gene Symbol`, "")) %>% 
  ggplot(aes(x = logFC, y = B, col=Significant,label=Label)) + geom_point() + geom_text_repel(col="black")











## filter for genes of interest

GOI_list = c("PXDN", "PXDNL", "RELA", "CTCF", "TCF21", "RAD21", "EP300", "ERG", "STAG2", "RBPG", "STAG1", "TAl1", "SMC1A", "JUN", "FOS")

filtered_annotated_results <- annotated_results[annotated_results$`Gene Symbol` %in% GOI_list, ]
