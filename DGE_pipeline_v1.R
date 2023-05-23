## data pre-processing

## script to normalize, annotate and analyse microarray data sets ()
## set wd ("D:\Project data\MSc data\!!!example!!!")
## check wd

getwd()

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

## download raw data from supp files

getGEOSuppFiles("GSE2638")

## files are compressed, decompress using untar()

untar("GSE2638/GSE2638_RAW.tar", exdir = "data/")

## store data from cel files into an object, make sure to store as a list

raw_data = ReadAffy(celfile.path = "data/")

## identify useful column names using pData

series_data = getGEO("GSE2638")
series_data = series_data[[1]]

sample_info = pData(series_data)

sample_info = select(sample_info, title)

sample_info = rename(sample_info, experiment = title)






##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

## normalise data: multiple methods to normalise microarray data, each with its own advantages/disadvantages
## this study will use RMA and the limma package to normalise and process data

## first, normalise using rma

rma_normalised_data = rma(raw_data)

## rma: fetch normalised expression data and store as a dataframe

rma_normalised_expr = as.data.frame(exprs(rma_normalised_data))

## rma: attach annotation data to dataframe so that expression data can be linked to their associated gene
## right now dataframe contains probe IDs. map gene symbols to probe IDs (!!!!will be array and study specific)

## store gene symbols in dataframe/matrix to attach to data

GEO_anno_data = getGEO("GSE2638", GSEMatrix = TRUE) # gets file for experiment series
GEO_gene_symbols = GEO_anno_data$GSE2638_series_matrix.txt.gz@featureData@data # isolates data file which includes gene symbols in column x

GEO_gene_symbols = GEO_gene_symbols[, c(1, 11)]  # subsets columns with probe ID and corresponding gene symbols

## merge rma_normalised_expr df with GEO_gene_symbols df

rma_normalised_annotated_expr = rma_normalised_expr %>%
  rownames_to_column(var = "ID") %>%
  inner_join(., GEO_gene_symbols, by = "ID")

## now that the data has been pre-processed, write the final file into a dataframe for use later

write.table(rma_normalised_annotated_expr, file = "rma_normalised_annotated_expr", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)






## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

## data visualization:
## load libraries

library(pheatmap)
library(ggplot2)
library(tidyr)
library(ggrepel)

## visualize data after normalization (check if data is normally distributed)

summary(exprs(rma_normalised_data))

boxplot(exprs(rma_normalised_data), outline = F)

corMatrix = cor(exprs(rma_normalised_data))

pheatmap(corMatrix)

## rename names in corMatrix to match sample info
## edit columns to have consistent, useful names and add columns for treated vs untreated

sample_info[1, 1] = "HMEC_control 1"
sample_info[2, 1] = "HMEC_control 2"
sample_info[3, 1] = "HMEC_control 3"

sample_info[4, 1] = "HMEC_TNF-treated 1"
sample_info[5, 1] = "HMEC_TNF-treated 2"
sample_info[6, 1] = "HMEC_TNF-treated 3"

sample_info$Treatment = c("Untreated - Control", "Untreated - Control", "Untreated - Control", "Treated - TNF Alpha", "Treated - TNF Alpha", "Treated - TNF Alpha")

## rename column headers if needed using dplyr

sample_info = rename(sample_info, Experiment = experiment, Treatment = treatment)

## adjust col_names for cormatrix if needed

colnames(corMatrix) = rownames(sample_info)

## check that rownames in corMatrix match rownames for annotations

rownames(sample_info)
colnames(corMatrix)

rownames(corMatrix) = rownames(sample_info)

## redraw heat map

pheatmap(corMatrix, annotation_col = sample_info)

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






## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

## Determine the fold change and significance in expression between treatment and control groups

## for directly comparing expression between samples, isolate PXDN expr data from dataframe

PXDN_expr_data <- rma_normalised_annotated_expr[rma_normalised_annotated_expr$`Gene Symbol` == "PXDN", ]

## use limma to visualize DE
## first create a design matrix based on whether a sample is treated or untreated
## design matrix assigns a 0 and 1 values depending on whether the label in sample info is treated/untreated

library(limma)
design <- model.matrix(~0+sample_info$Treatment)
design

colnames(design) = c("TNF_Treated", "Untreated_Control")

## prep data

summary(rma_normalised_expr)

## calculate median expression level, use as cutoff value for now
## The cutoff value chosen will be standardised as the median expression value for now
## This is done to reduce the noise in the data as well as reduce the ""multiple testing burden"
## If PXDN or PXDNL is lowly expressed in any of these samples, it will be filtered out

## PXDN had a relatively high expression as seen earlier in this sample, 
## Since it will not exclude PXDN therefore the median will be used as the cutoff
## for this data set the median is 7.5


## AUTOMATE


cutoff <- (7.603 + 7.565 + 7.530 + 7.568 + 7.510 + 7.527)/6

## TRUE or FALSE for whether each gene is "expressed" in each sample

is_expressed <- rma_normalised_expr > cutoff

## Identify genes expressed in more than 2 samples

keep <- rowSums(is_expressed) > 2

## check how many genes are removed / retained.
table(keep)

## subset to just those expressed genes
rma_normalised_filtered_data <- rma_normalised_data[keep,]

## Fit data to design matrix

fit <- lmFit(exprs(rma_normalised_filtered_data), design)


head(fit$coefficients)


contrasts <- makeContrasts(TNF_Treated  - Untreated_Control, levels=design)

## can define multiple contrasts
## e.g. makeContrasts(Group1 - Group2, Group2 - Group3,....levels=design)

fit2 <- contrasts.fit(fit, contrasts)

fit2 <- eBayes(fit2)

## store results in an object

results = topTable(fit2, number = Inf)


## AUTOMATE


## make objects containing identifiers for PXDN based on array probe ID

PXDN_identifier_1 = "212012_at"
PXDN_identifier_2 = "212013_at"

## use identifiers to filter expression data

PXDN_results_1 = results[rownames(results) == PXDN_identifier_1, ]
PXDN_results_2 = results[rownames(results) == PXDN_identifier_2, ]

## call for log fold change and significance values

logFC_1 = PXDN_results_1$logFC
p_val_1 = PXDN_results_1$P.Value

logFC_2 = PXDN_results_2$logFC
p_val_2 = PXDN_results_2$P.Value