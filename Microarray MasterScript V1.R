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

BiocManager::install("org.Hs.eg.db")
BiocManager::install("GEOquery")
BiocManager::install("tidyverse", force = TRUE)

library(affy)
library(GEOquery)
library(tidyverse)
#library(dplyr)
library(limma)

## load data
## First step:
## store sample_ID

sample_ID = "GSE3586"

## Now to download the data from GEOquery in the .cel file format and store the data in an object
## There are 3 ways to do this
## 1. download tar file from GEO and untar .cel files in a data/ folder in your wd
## 2. download cel files directly from GEO and store them in your data/ folder in wd
## 3. !!!!INCOMPLETE pipeline in progress :: store .cel file data directly from GEO into an object 

## 1. following method 1...
## download raw data from supp files

getGEOSuppFiles(sample_ID)

## files are compressed, decompress using untar()

## first create a string which should match file name

sample_string = sample_ID %>%              # sample_ID is input for pipe
  paste(., "_RAW.tar", sep = "") %>%       # join '_RAW.tar' string, standard for tar files downloaded using getGEO
  paste(sample_ID, ., sep = "/")           # join sample using / to mimic directory

## then untar..

untar(sample_string, exdir = "data/")

###Checkpoint!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## Depending on the microarray, you may either have .cel or .gpr files

## .cel files are produced by single-color microarrays analysed with Affymetrix scanners
## .gpr files are produced by double-color microarrays analysed with GenePix scanners

## depending on which data type you are using, you will need to use different tools to load the data

## Following Methods 1. and 2.
## .cel/gpr files should be stored in your data/ folder




### CHOOSE ONLY 1 PATH!!!!!!!!!!!!!!!!!!




###Checkpoint!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## Loading data from .cel files

## ensure that .cel files belong to the same experiment series
## usually only needed if experiment dataset consists of multiple series/platforms
## if multiple series of interest are present, apply pipeline to one series at a time

## store data from cel files into an object, make sure to store as a list

## create a list for all .cel files in "data"

cel_file_names = list.files(path = "data/", pattern = "\\.CEL\\.gz$", full.names = TRUE)

raw_data = read.celfiles(cel_file_names)





###Checkpoint!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## Loading data from .gpr files

## make a list of the gpr file names

gpr_file_names = list.files("data", pattern = "\\.gpr$", full.names = TRUE) #ensure that files are saved with utf-8 encoding

## store data from .gpr files into an object

raw_data = read.maimages(gpr_file_names, source="genepix")













###Checkpoint!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## (if you are using methods 1 or 2, use the following to get sample info)
## identify useful column names using pData

series_data = getGEO(sample_ID)
series_data = series_data[[1]]

sample_info = pData(series_data)

## after storing pData in sample_info object, trim it down to keep only the important labels

## keep useful info in sample info by using select(sample_info, heading_of_interest_1, heading_of_interest_2)
## these columns amy contain useful information regarding cell type and preparation

sample_info = sample_info %>% select(characteristics_ch1, geo_accession)

sample_info = rename(sample_info, Treatment = characteristics_ch1)

## Name according to string automatically

sample_info <- sample_info %>%
  group_by(Treatment) %>%
  mutate(Label = paste(Treatment, row_number(), sep="_"))

sample_info = rename(sample_info, Treatment = Condition)



## ONLY USE BELOW CODE IF DF IS MISSING GEO ACCESSION NUMBERS IN ROWNAMES

sample_info <- sample_info %>% 
  select(Treatment, geo_accession) %>%
  tibble::rownames_to_column() %>%
  column_to_rownames(var = "geo_accession")




##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

## DATA NORMALIZATION FOR .cel FILES

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

plotting_col_names = names(rma_normalised_expr)

## alternative to name plotting columns using sampleinfo

#plotting_col_names = rownames(sample_info)



## DATA NORMALIZATION FOR .gpr FILES

## first perform background correction

data_bc <- backgroundCorrect(raw_data, method = "normexp", offset = 50)

## LOWESS normalization (within-array normalization)
## commonly used in 2 dye arrays. check if this is relevant for your study

norm_data <- normalizeWithinArrays(data_bc, method = "loess")




##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

## DATA ANNOTATION

## attach annotation data to dataframe so that expression data can be linked to their associated gene
## right now dataframe contains probe IDs. map gene symbols to probe IDs (!!!!will be array and study specific)

## store gene symbols in dataframe/matrix to attach to data

GEO_series_data = getGEO(sample_ID, GSEMatrix = TRUE) # gets file for experiment series
GEO_anno_data = GEO_series_data$GSE42955_series_matrix.txt.gz@featureData@data# isolates data file which includes gene symbols in column x

GEO_gene_symbols = GEO_anno_data[, c(1, 10)]  # subsets columns with probe ID and corresponding gene symbols


## WARNING: if using newer microarrays or platforms which include multiple gene names
## use the following code to split and isolate the gene codes you are using
## The following code isolates genesymbols from one such column
## split using '//', then use function to get the second element

gene_symbols = sapply(strsplit(GEO_gene_symbols$gene_assignment, split = "//"), function(x) x[2])

## assign back to df in new column

GEO_gene_symbols$GeneSymbol = gene_symbols





## merge rma_normalised_expr df with GEO_gene_symbols df

### FOR .cel FILES

rma_normalised_annotated_expr = rma_normalised_expr %>%
  rownames_to_column(var = "ID") %>%
  inner_join(., GEO_gene_symbols, by = "ID")

## sometimes the ID value types will not match, in which case use the following code to convert 

GEO_gene_symbols$ID = as.character(GEO_gene_symbols$ID)


### FOR .gpr FILES
norm_data$genes = GEO_gene_symbols$GENE_SYMBOL




## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

## DATA QUALITY CONTROL
## load libraries

library(pheatmap)
library(ggplot2)
library(tidyr)
library(ggrepel)
#library(affy)




## BOXPLOT

## view distribution before rma

boxplot(exprs(raw_data), outline = F, las = 2, names = plotting_col_names)


## visualize data after normalization (check if data is normally distributed)

summary(rma_normalised_expr)

boxplot(rma_normalised_expr, outline = F, las = 2, names = plotting_col_names)




## HISTOGRAMS FOR SAMPLES

hist(raw_data)




## MVA PLOT

MAplot(raw_data,pairs=TRUE,plot.method="smoothScatter", names = plotting_col_names)




## View microarray images to check for damage/artefacts 

image(raw_data)




## RNA degradation plots

deg = AffyRNAdeg(raw_data)

plotAffyRNAdeg(deg)




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

sample_info$Experiment = c("Control_1", "Control_2", "Control_3", "Treated_1", "Treated_2", "Treated_3")

## rename column headers if needed using dplyr

# sample_info = rename(sample_info, Experiment = experiment, Treatment = treatment)

## check that rownames in corMatrix match rownames for annotations

rownames(sample_info)

colnames(corMatrix)
rownames(corMatrix)

## adjust col_names and rownames for cormatrix if needed

colnames(corMatrix) = rownames(sample_info)

rownames(sample_info) = colnames(corMatrix)

## draw heat map

pheatmap(corMatrix, annotation_col = sample_info)






## PCA PLOT

## now to perform a principal component analysis PCA
## first, transpose the expression data using t() and use the prcomp function to store the principal components in an object

pca <- prcomp(t(rma_normalised_expr))

## Join the PCs to the sample information and visualize

cbind(sample_info, pca$x) %>% 
  ggplot(aes(x = PC1, y=PC2, col=Treatment,label=paste("", Experiment))) + geom_point() + geom_text_repel() +
  labs(title = "Example Quality Control PCA") +
  theme(plot.title = element_text(hjust = 0.5))

## visualise the proportions of each principal component
## first, get summary of PCA and look for columns/rows of interest

pca_sum = summary(pca)

## extract proportion of variance data

pca_proportions = pca_sum$importance["Proportion of Variance", ]

## convert proportions to a percentage out of 100

pca_proportions_percent = pca_proportions * 100

barplot(pca_proportions_percent, main = "Proportion of Variance explained by each Principal Component", xlab = "Principal Component", ylab = "Percentage of Variance")




## QC FOR .gpr FILES


## boxplots for raw data (green)
boxplot(raw_data$G, outline = F, las = 2)

## boxplots for raw data (red)
boxplot(raw_data$R, outline = F, las = 2)


## post normalisation boxplot
## boxplot after normalization DGE (M is the log fold change when comparing G/R)
boxplot(norm_data$M, outline = F, las = 2)

## boxplot after normalization QC
boxplot(norm_data$A, outline = F, las = 2)



## print MA plots in 3 by 2 to your wd
plotMA3by2(raw_data)



## print out histogram of average values for normalised data
## used to investigate distribution of intensity, should follow a normal dist
hist(norm_data$A, main = "Histogram of A Values", xlab = "Average Log Intensity", breaks = 50)



## print boxplot of A values
## investigate array quality
boxplot(norm_data$A, main = "Boxplot of A Values", xlab = "Arrays", ylab = "Average Log Intensity")


## Background intensity boxplots
## use to identify arrays with abnormally high background fluorescence for each colour of the array

boxplot(data.frame(log2(raw_data$Gb)),main="Green background")

boxplot(data.frame(log2(raw_data$Rb)),main="Red background")



# define  number of images/samples
num_images = 56


# loop through each image for the red dye
for (i in 1:num_images) {
  # Define the file name
  file_name <- sprintf("appendix_GSE3586_image%d.png", i)
  
  # Start the PNG device
  png(filename = file_name, width = 848, height = 512)
  
  # Create the image plot
  imageplot(log2(raw_data$Rb[, i]), raw_data$printer)
  
  # Close the PNG device
  dev.off()
}

# repeat for the green dye
for (i in 1:num_images) {
  # Define the file name
  file_name <- sprintf("appendix_GSE3586_alt_image%d.png", i)
  
  # Start the PNG device
  png(filename = file_name, width = 848, height = 512)
  
  # Create the image plot
  imageplot(log2(raw_data$Gb[, i]), raw_data$printer)
  
  # Close the PNG device
  dev.off()
}


## Cor Heatmap for .gpr

## calculate correlation matrix
cor_matrix = cor(t(norm_data$M))

## plot graph

pheatmap(cor_matrix, clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         annotation_col = sample_anno$Treatment)




## PCA for .gpr

## first create an annotation df

sample_anno = data.frame(Treatment = rep(sample_info$Treatment, each = 2))

## shoorten for better readability

sample_anno_short = sample_anno %>%
  mutate(Treatment_abbrev = case_when(
    Treatment == "Non failing"  ~ "NF",
    Treatment == "Dilated cardiomyopathy" ~ "DC",
    TRUE              ~ Treatment
  ))

## typically done after normalisation of $M

pca_data = norm_data$M

## use base R for PCA

pca_result = prcomp(t(pca_data), scale. = TRUE)

## plot PCA with % variance of each PC

cbind(sample_anno, pca_result$x) %>% 
  ggplot(aes(x = PC1, y=PC2, col=Treatment,label=paste("", Treatment))) + geom_point() + geom_text_repel() +
  labs(title = "GSE3586 PCA") +
  theme(plot.title = element_text(hjust = 0.5))

plot(pca_result$x[,1], pca_result$x[,2], xlab = "PC1", ylab = "PC2", main = "PCA of GSE3586")
text(pca_result$x[,1], pca_result$x[,2], labels = sample_anno_short$Treatment_abbrev, pos = 4)

## visualise the proportions of each principal component
## first, get summary of PCA and look for columns/rows of interest

pca_sum = summary(pca_result)

## extract proportion of variance data

pca_proportions = pca_sum$importance["Proportion of Variance", ]

## convert proportions to a percentage out of 100

pca_proportions_percent = pca_proportions * 100

barplot(pca_proportions_percent, main = "Proportion of Variance explained by each Principal Component", xlab = "Principal Component", ylab = "Percentage of Variance")










## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

## STATISTICAL SIGNIFICANCE AND FOLD CHANGE

## Determine the fold change and significance in expression between treatment and control groups

## for directly comparing expression between samples, isolate PXDN and PXDNL expr data from dataframe
## here GOI stands for Genes of Interest

GOI_expr_data <- rma_normalised_annotated_expr[rma_normalised_annotated_expr$GeneSymbol %in% c("PXDN", "PXDNL"), ]

write.table(PXDN_expr_data, file = "GOI_expr_GSE2639", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)


## use limma to visualize DE and handle outliers
## first create a design matrix based on whether a sample is treated or untreated
## design matrix assigns a 0 and 1 values depending on whether the label in sample info is treated/untreated

## load limma for arrayweights

library(limma)

## Handling outliers
## Option 1: remove outliers identified by PCA
## remove from raw data and redo all steps up to this point

#outliers <- c(outlier_name_1, outlier_name_2)

#gse <- gse[,-outlier_samples]

## Option 2: use the arrayWeights function to give outliers a weaker contribution during DGE
## Make sure to include aw as a parameter when calling the lmFit function

## ONLY USE arrayWeights ON PREPROCESSED DATA (eg., after RMA)


design <- model.matrix(~0+sample_info$Treatment)
design

aw = arrayWeights(exprs(rma_normalised_data),design)

## name the columns for the classes you are going to compare. Should match with your design matrix.

colnames(design) = c("ICM", "DCM", "Control")

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

#fit <- lmFit(exprs(rma_normalised_filtered_data), design, weights = aw)




## run on unfiltered data (rma_normalised_data) instead of filtered (rma_normalised_filtered_data) to retain all genes

fit <- lmFit(exprs(rma_normalised_data), design, weights = aw)

## fit for .gpr

#fit <- lmFit(norm_data$M, design)


head(fit$coefficients)


contrasts <- makeContrasts(ICM  - DCM, levels=design)

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





## FITTING FOR .gpr

fit <- lmFit(norm_data$M, design)

head(fit$coefficients)

contrasts <- makeContrasts(DCM  - Non_failing, levels=design)

## can define multiple contrasts
## e.g. makeContrasts(Group1 - Group2, Group2 - Group3,....levels=design)

fit2 <- contrasts.fit(fit, contrasts)

fit2 <- eBayes(fit2)

## store results in an object

results = topTable(fit2, number = Inf)

## Annotation of Differential Gene Expression Analysis

results$ID <- raw_data$genes$Name

annotated_results = results %>%
  inner_join(., GEO_gene_symbols, by = "ID")





## Volcano plot

## set cutoff values for significant gene expression

p_cutoff = 0.05
fc_cutoff = 0.5

## OPTIONAL create a list of genes_of_interest to display on volcano plot

## Note, insert your control genes here!!! Use the first list template for older U133 Affy arrays

GOI_list = c("PXDN", "PXDNL", "RELA", "CTCF", "TCF21", "RAD21", "EP300", "ERG", "STAG2", "RBPG", "STAG1", "TAL1", "SMC1A", "JUN", "FOS", 'GAPDH', 'NPPA')

## Use this template for any array you needed to split the gene symbol column for

GOI_list = c(' PXDN ', ' PXDNL ', ' RELA ', ' CTCF ', ' TCF21 ', ' RAD21 ', ' EP300 ', ' ERG ', ' STAG2 ', ' RBPG ', ' STAG1 ', ' TAL1 ', ' SMC1A ', ' JUN ', ' FOS ', ' GAPDH ', ' NPPA ')

annotated_results_final <- annotated_results %>%
  mutate(is_interest = ifelse(annotated_results$GeneSymbol %in% GOI_list, "Interest", "Not interest"))

# check

annotated_results_final$Symbol %in% GOI_list




## Use annotated results to create volcano plot


annotated_results_final %>% 
  mutate(Significant = adj.P.Val < p_cutoff & abs(logFC) > fc_cutoff ) %>% 
  mutate(Label = ifelse(annotated_results_final$is_interest == "Interest", `Gene Symbol`, NA)) %>% 
  ggplot(aes(x = logFC, y = B, col=Significant,label=Label)) + geom_point() + geom_text_repel(col="black") + 
  labs(title="", x="Log2 Fold Change in gene expression", y="B") +
  theme(plot.title = element_text(hjust = 0.5))




## filter for genes of interest

GOI_list = c("PXDN", "PXDNL", "RELA", "CTCF", "TCF21", "RAD21", "EP300", "ERG", "STAG2", "RBPG", "STAG1", "TAL1", "SMC1A", "JUN", "FOS", 'GAPDH', 'NPPA')

filtered_annotated_results <- annotated_results[annotated_results$`Gene Symbol` %in% GOI_list, ]







## Barplots to visualize gene expression levels across samples for each gene

# isolate annotated results for gois in new object

genesOfInterest = rma_normalised_annotated_expr %>%
  filter(`Gene Symbol` %in% GOI_list)


# transform df to long form

long_genesOfInterest = genesOfInterest %>%
  select(., -c('ID'))  %>%    #remove unneeded columns
  pivot_longer(
    cols = -`Gene Symbol`,
    names_to = "Sample",
    values_to = "Expression"
  )

# transform sample_info to use for annotating bar graphs

sample_info_new <- sample_info %>%
  rownames_to_column(var = "Sample")

# merge sample info with long df

long_genesOfInterest_annotated = long_genesOfInterest  %>%
  left_join(sample_info_new, by = 'Sample')

# write results into csv for further analysis/visualisation using external tools

write.csv(long_genesOfInterest_annotated, "GSE2638_GOI_exprs.csv", row.names=FALSE)




## misc code

# write.csv(annotated_results, "vplot data.csv", row.names=FALSE)

#annotated_results %>% 
#  mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
#  mutate(Rank = 1:n(), Label = ifelse(Rank < topDEGs, `Gene Symbol`, "")) %>% 
#  ggplot(aes(x = logFC, y = B, col=Significant,label=Label)) + geom_point() + geom_text_repel(col="black") + 
#  labs(title="Example Volcano plot showing statistically significant changes in gene expression", x="Log2 Fold Change in gene expression", y="Level of statistical significance")



# new vplot, highlights GOIs

annotated_results_final %>%
  mutate(Significant = adj.P.Val < p_cutoff & abs(logFC) > fc_cutoff) %>%
  mutate(Label = ifelse(is_interest == "Interest", GeneSymbol, NA),
         Highlight = is_interest == "Interest") %>%
  ggplot(aes(x = logFC, y = B, col = Significant)) +
  geom_point(aes(size = Highlight)) +  # Control size of points based on Highlight
  scale_size_manual(values = c(1, 2), guide = "none") +  # Adjust the size of points, no legend for size
  geom_point(data = subset(annotated_results_final, is_interest == "Interest"), 
             aes(x = logFC, y = B), color = "black", size = 2) +  # Add black points for GOIs
  geom_text_repel(aes(label = Label), col="black") +  # Ensure Label is used here
  labs(title="Statistically significant changes in gene expression for GSE42955\nDilated cardiomyopathy vs Ischemic cardiomyopathy",
       x="Log2 Fold Change in gene expression",
       y="B") +
  theme(plot.title = element_text(hjust = 0.5))
