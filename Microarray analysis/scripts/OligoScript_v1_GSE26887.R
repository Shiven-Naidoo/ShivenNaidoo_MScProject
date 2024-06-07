### oligo pipeline for newer microarrays

## download and call libraries needed for data downloading and normalisation

#BiocManager::install("pd.hugene.1.0.st.v1")

library(GEOquery)
library(tidyverse)
#library(dplyr)
library(pheatmap)
library(ggplot2)
library(tidyr)
library(ggrepel)
library(limma)
library(oligo)
library(pd.hugene.1.0.st.v1)
library(pheatmap)
library(arrayQualityMetrics)

## check functions built into package

ls("package:oligo")


## DATA RETRIEVAL AND PRE-PROCESSING

## store sample_ID

sample_ID = "GSE26887"

## download raw data from supp files

getGEOSuppFiles(sample_ID)

## files are compressed, decompress using untar()

## create string which should match file name

sample_string = sample_ID %>%              # sample_ID is input for pipe
  paste(., "_RAW.tar", sep = "") %>%       # join '_RAW.tar' string, standard for tar files downloaded using getGEO
  paste(sample_ID, ., sep = "/")           # join sample using / to mimic directory

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

gpr_file_names = list.files("data", pattern = "\\.gpr\\.gz", full.names = TRUE)

## store data from .gpr files into an object

raw_data = read.maimages(gpr_file_names, source="genepix")




###Checkpoint!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

## get sample info

series_data = getGEO(sample_ID)

## identify useful column names using pData

# if there are multiple platforms, use the index to select your dataset of interest

series_data = series_data[[1]]     # [[1]] will extract the first platform etc...         

sample_info = pData(series_data)

## after storing pData in sample_info object, trim it down to keep only the important labels

#View(sample_info)

## keep useful info in sample info by using select(sample_info, heading_of_interest_1, heading_of_interest_2)
## these columns amy contain useful information regarding cell type and preparation

## visualize columns of interest

colnames(sample_info)









## use select from dplr/tidyverse to subset columns of interest by name
## eg, if your columns of interest are 'title' and 'disease state:ch1', you would use...

sample_info = sample_info %>% 
  select('title', 'disease state:ch1')

## rename columns

#sample_info = rename(sample_info, Experiment = title)
#sample_info = rename(sample_info, Treatment = 'disease state:ch1')

## alt rename columns

colnames(sample_info) = c("Experiment", "Treatment")

## treatment renaming

## if you want to shorten the names in your treatment column, use the following code...

## for this dataset, samples fall under 3 categories which need to be counted

treatment_tracker = c("CONTROL" = 0, "DIABETIC, HEART FAILURE" = 0, "NON DIABETIC, HEART FAILURE" = 0)

# create an empty vector which labels will be stored in later

new_labels = character(nrow(sample_info))

# loop through the rows in sample_info
for (i in 1:nrow(sample_info)) {
  # store the current treatment label
  current_treatment = sample_info$Treatment[i]
  
  # update the treatment count
  treatment_tracker[current_treatment] = treatment_tracker[current_treatment] + 1
  
  # create new label which will be added to empty vector later
  if (current_treatment == "CONTROL") {
    new_label <- paste0("C", treatment_tracker[current_treatment], sep = "")
  } else if (current_treatment == "DIABETIC, HEART FAILURE") {
    new_label <- paste0("DHF", treatment_tracker[current_treatment], sep = "")
  } else if (current_treatment == "NON DIABETIC, HEART FAILURE") {
    new_label <- paste0("NDHF", treatment_tracker[current_treatment], sep = "")
  } else {
    new_label <- NA  # for unexpected treatment values
  }
  
  # store new label in empty vec..
  new_labels[i] <- new_label
}


# use new_label vec to add new unique name column or replace an existing column
sample_info$Experiment <- new_labels















##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

## DATA NORMALIZATION

## normalise data: multiple methods to normalise microarray data, each with its own advantages/disadvantages
## this study will use RMA and the limma package to normalise and process data

## check datasets first to see if normalising methods have already be applied
## data processing should be specified in raw_data$phenoData

raw_exprs = exprs(raw_data)

## Rename column names to represent just the samples without the .cel extensions

colnames(raw_exprs) <- sub(".CEL.gz", "", colnames(raw_exprs))

## !!!!!!!!!!!!!!!!!!!!!!!!Optional for when plotting graphs
## Further shorten names

#colnames(raw_exprs) <- sub("GSM.......", "", colnames(raw_exprs))




## first, normalise using rma

rma_normalised_data = rma(raw_data)

## rma: fetch normalised expression data and store as a dataframe

rma_normalised_expr = as.data.frame(exprs(rma_normalised_data))

## Rename column names to represent just the samples without the .cel extensions

names(rma_normalised_expr) <- sub(".CEL.gz", "", names(rma_normalised_expr))

## Further shorten names

names(rma_normalised_expr) <- sub("_.*", "", names(rma_normalised_expr))

## check datafram to see if expr values are normalized

summary(rma_normalised_expr)








##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

## DATA ANNOTATION

## attach annotation data to dataframe so that expression data can be linked to their associated gene
## right now dataframe contains probe IDs. map gene symbols to probe IDs (!!!!will be array and study specific)

## store gene symbols in dataframe/matrix to attach to data

GEO_series_data = getGEO(sample_ID, GSEMatrix = TRUE) # gets file for experiment series
GEO_anno_data = GEO_series_data$GSE26887_series_matrix.txt.gz@featureData@data # isolates data file which includes gene symbols in column x

GEO_gene_symbols = GEO_anno_data[, c(1, 10)]  # subsets columns with probe ID and corresponding gene symbols

## first, convert ID column in GEO_gene_symbols to characters 

GEO_gene_symbols$ID = as.character(GEO_gene_symbols$ID)

## extract gene code from gene assignment column

## split using '//', then use function to get the second element

gene_symbols = sapply(strsplit(GEO_gene_symbols$gene_assignment, split = "//"), function(x) x[2])

## assign back to df in new column

GEO_gene_symbols$GeneSymbol = gene_symbols


## merge rma_normalised_expr df with GEO_gene_symbols df

rma_normalised_annotated_expr = rma_normalised_expr %>%
  rownames_to_column(var = "ID") %>%
  inner_join(., GEO_gene_symbols, by = "ID")












## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

## DATA QUALITY CONTROL

## GENERATE ARRAYQUALITYMETRICS REPORT

arrayQualityMetrics(expressionset = raw_data)

## BOXPLOT
## visualize data after normalization (check if data is normally distributed)

summary(rma_normalised_expr)

## view distribution before rma

boxplot(raw_exprs, outline = F, las = 2)

## view distribution after rma

## get names from raw_exprs

plotting_col_names = colnames(raw_exprs)

boxplot(rma_normalised_expr, outline = F, las = 2, names = plotting_col_names)




colnames(raw_exprs)


## UNDER CONSTRUCTION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## RLE and NUSE plots for array quality






## HEATMAP

## condense normalised expression data into a corMatrix

corMatrix = cor(rma_normalised_expr)

## check heatmap (raw, no labels)

pheatmap(corMatrix)

## rename names in corMatrix to match sample info
## edit columns to have consistent, useful names and add columns for treated vs untreated

#sample_info$Experiment = gsub("Human lymphatic endothelial cells ", "", sample_info$Experiment)

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
  ggplot(aes(x = PC1, y=PC2, col=Treatment,label=paste(""))) + geom_point() + geom_text_repel()

## visualise the proportions of each principal component
## first, get summary of PCA and look for columns/rows of interest

pca_sum = summary(pca)

## extract proportion of variance data

pca_proportions = pca_sum$importance["Proportion of Variance", ]

## convert proportions to a percentage out of 100

pca_proportions_percent = pca_proportions * 100

barplot(pca_proportions_percent, main = "Proportion of Variance explained by each Principal Component", xlab = "Principal Component", ylab = "Percentage of Variance")




## RNA DEGREDATION PLOT

## only load affy when drawing this plot, may interfere with functions from oligo

library

deg = AffyRNAdeg(raw_data)

plotAffyRNAdeg(deg)







### Handling outliers

## Option 1: remove outliers identified by PCA
## remove from raw data and redo all steps up to this point

#outliers <- c("GSM143717", "GSM143907")

#raw_data <- raw_data[,-outliers]

## Option 2: use the arrayWeights function to give outliers a weaker contribution during DGE
## Make sure to include aw as a parameter when calling the lmFit function

## ONLY USE arrayWeights ON PREPROCESSED DATA (eg., after RMA)

## first create a design matrix based on whether a sample is treated or untreated
## design matrix assigns a 0 and 1 values depending on whether the label in sample info is treated/untreated



design <- model.matrix(~0+sample_info$Treatment)
design

aw = arrayWeights(exprs(rma_normalised_data),design)












## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

## STATISTICAL SIGNIFICANCE AND FOLD CHANGE


## name the columns for the classes you are going to compare. Should match with your design matrix.

colnames(design) = c("Control", "DiabeticHF", "NonDiabeticHF")

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

rma_normalised_filtered_data = rma_normalised_data[keep,]

## Fit data to design matrix

fit <- lmFit(exprs(rma_normalised_filtered_data), design, weights = aw)


head(fit$coefficients)


contrasts <- makeContrasts(NonDiabeticHF - Control, levels=design)

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
fc_cutoff = 1


## OPTIONAL create a list of genes_of_interest to display on volcano plot

## when making the gene list, make sure the genes match the entry in you df exactly
## This means space for space. You can check how your genes are listed using manual indexing

# annotated_results$GeneSymbol[11844]

# Here the genes are flanked by sapces, account for this in list or use .strip when iterating through df

GOI_list = c(' PXDN ', ' PXDNL ', ' RELA ', ' CTCF ', ' TCF21 ', ' RAD21 ', ' EP300 ', ' ERG ', ' STAG2 ', ' RBPG ', ' STAG1 ', ' TAL1 ', ' SMC1A ', ' JUN ', ' FOS ')

annotated_results_final <- annotated_results %>%
  mutate(is_interest = ifelse(GeneSymbol %in% GOI_list, "Interest", "Not interest"))

# print results to file for further analysis/recall

write.csv(annotated_results, "control_NonDHF.tsv", row.names=FALSE)


## Use annotated results to create volcano plot

write.csv(annotated_results, "control_NonDHF.tsv", row.names=FALSE)

is_interest = ifelse(annotated_results$GeneSymbol %in% GOI_list, "Interest", "Not interest")

annotated_results_final %>% 
  mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
  mutate(Label = ifelse(annotated_results_final$is_interest == "Interest", `GeneSymbol`, NA)) %>% 
  ggplot(aes(x = logFC, y = B, col=Significant,label=Label)) + geom_point() + geom_text_repel(col="black") + 
  labs(title="Statistically significant changes in gene expression for GSE26887\nControl vs Non diabetic heart failure", x="Log2 Fold Change in gene expression", y="Level of statistical significance") +
  theme(plot.title = element_text(hjust = 0.5))






## Barplots to visualize gene expression levels across samples for each gene

# isolate annotated results for gois in new object

genesOfInterest = rma_normalised_annotated_expr %>%
  filter(GeneSymbol %in% GOI_list)


# transform df to long form

long_genesOfInterest = genesOfInterest %>%
  select(., -c('ID', 'gene_assignment'))  %>%    #remove unneeded columns
  pivot_longer(
    cols = -GeneSymbol,
    names_to = "Sample",
    values_to = "Expression"
  )

# transform sample_info to use for annotating bar graphs

sample_info_new <- sample_info %>%
  rownames_to_column(var = "Sample")

# merge sample info with long df

long_genesOfInterest_annotated = long_genesOfInterest  %>%
  left_join(sample_info_new, by = 'Sample')


# for loop to generate plots for each gene

for (gene_of_interest in GOI_list) {
  p <- ggplot(data = filter(long_genesOfInterest_annotated, GeneSymbol == gene_of_interest), 
              aes(x = Sample, y = Expression, fill = Treatment)) + 
    geom_bar(stat = "identity") +
    labs(title = paste("Expression of", gene_of_interest), 
         y = "Expression Level", 
         x = "Samples") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_brewer(palette = "Set1", name = "Treatment")
  
  print(p)
}


# repeat for raw expression values

# run rma without background correction/normalization etc


raw_data_d = rma(raw_data, normalize=FALSE)

# store raw data in df

raw_df = as.data.frame(exprs(raw_data_d))

# remove sample name extensions

names(raw_df) <- sub("_.*", "", names(raw_df))

## subset to just those expressed genes

raw_df_filtered <- raw_df[keep,]

# Annotate

raw_df_filtered_annotated = raw_df_filtered %>%
  rownames_to_column(var = "ID") %>%
  inner_join(., GEO_gene_symbols, by = "ID")




## Barplots to visualize gene expression levels across samples for each gene

# isolate annotated results for gois in new object

raw_genesOfInterest = raw_df_filtered_annotated %>%
  filter(GeneSymbol %in% GOI_list)


# transform df to long form

long_raw_genesOfInterest = raw_genesOfInterest %>%
  select(., -c('ID', 'gene_assignment'))  %>%    #remove unneeded columns
  pivot_longer(
    cols = -GeneSymbol,
    names_to = "Sample",
    values_to = "Expression"
  )

# merge sample info with long df

long_raw_genesOfInterest_annotated = long_raw_genesOfInterest  %>%
  left_join(sample_info_new, by = 'Sample')


# for loop to generate plots for each gene

for (gene_of_interest in GOI_list) {
  p <- ggplot(data = filter(long_raw_genesOfInterest_annotated, GeneSymbol == gene_of_interest), 
              aes(x = Sample, y = Expression, fill = Treatment)) + 
    geom_bar(stat = "identity") +
    labs(title = paste("Expression of", gene_of_interest), 
         y = "Expression Level", 
         x = "Samples") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_brewer(palette = "Set1", name = "Treatment")
  
  print(p)
}








## Outlier removal

