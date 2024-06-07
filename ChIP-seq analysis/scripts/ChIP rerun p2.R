
##### Data Visualisation
## Install Bioconductor core packages and Load CHIPseeker packages

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install(version = "3.13")
BiocManager::install("ChIPpeakAnno")
BiocManager::install("clusterProfiler")
BiocManager::install("ChIPseeker", force = TRUE)
BiocManager::install("ReactomePA", force = TRUE)
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")

## Load libraries

library(ChIPseeker)
library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

## txdb contains database generated genomic annotation data from the UCSC
## split into 2 objects specific to chromosomes 2 and 8

txdb = TxDb.Hsapiens.UCSC.hg38.knownGene

library(clusterProfiler)

## Load in data

setwd("C:/Users/Shiven/Documents/R/Honours Project/Data/BED Files")




# CHECKPOINT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


PXDNTFs = read.table("PXDNFiltered2.bed", header = TRUE, quote = '"', sep = " ", comment.char = "")
PXDNLTFs = read.table("TFsPXDNL2.bed", header = TRUE, quote = '"', sep = " ", comment.char = "")

PXDNTFs
PXDNLTFs

## Isolate columns for chrom, start, end to make a peak file

peakPXDNTFs = PXDNTFs %>% select(1:3)
peakPXDNLTFs = PXDNLTFs %>% select(1:3)

## Write peak into new bed file
## make sure separators are "\t", and file ends with an emptyline

## !!! access file and make sure there is no header line

write.table(peakPXDNTFs, "peakPXDNTab.bed", sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(peakPXDNLTFs, "peakPXDNLTab.bed", sep = "\t", row.names = FALSE, col.names = FALSE)

## Annotate Peak files in 5kb radius of TSS

peakAnnoPXDN = annotatePeak(
  "peakPXDNTab.bed", 
  tssRegion=c(-5000, 5000), 
  TxDb=txdb, 
  level = "gene", 
  verbose = FALSE)

peakAnnoPXDNL = annotatePeak(
  "peakPXDNLTab.bed", 
  tssRegion=c(-5000, 5000), 
  TxDb=txdb, 
  level = "gene", 
  verbose = FALSE)

## Visualise Peak Data as pie

plotAnnoPie(peakAnnoPXDN, col = blues9)
plotAnnoPie(peakAnnoPXDNL, col = blues9)

## Visualise Peak Data as bar

plotAnnoBar(peakAnnoPXDN)

plotAnnoBar(peakAnnoPXDNL)

## Visualise data relative to TSS

plotDistToTSS(peakAnnoPXDN, title="")

plotDistToTSS(peakAnnoPXDNL, title = "")

## Store annotation data

PXDNanno = as.data.frame(peakAnnoPXDN@anno)

PXDNLanno = as.data.frame(peakAnnoPXDNL@anno)

write.table(PXDNanno, file="PXDNanno.txt", sep = "\t", row.names = FALSE)

write.table(PXDNLanno, file="PXDNLanno.txt", sep = "\t", row.names = FALSE)

## Load in annotated data

PXDNanno <- read.table("PXDNanno.txt", sep = "\t")

PXDNLanno <- read.table("PXDNLanno.txt", sep = "\t")

## Attach annotation data to tf data

library(tidyverse)

## Load missing columns

missingDataPXDN = read.table("PXDNFiltered2.bed", header = FALSE, quote = '"', sep = " ", comment.char = "")
missingDataPXDNL = read.table("TFsPXDNL2.bed", header = FALSE, quote = '"', sep = " ", comment.char = "")

missingDataPXDN
missingDataPXDNL

Macs_SRX_TF_CellType = select(missingDataPXDN, -c(V1, V2, V3))
Macs_SRX_TF_CellTypePXDNL = select(missingDataPXDNL, -c(V1, V2, V3))

Macs_SRX_TF_CellType

fullAnnoPXDN = cbind(PXDNanno, Macs_SRX_TF_CellType)
fullAnnoPXDNL = cbind(PXDNLanno, Macs_SRX_TF_CellTypePXDNL)

fullAnnoPXDN

colnames(fullAnnoPXDN) = c("chr", "start", "end", "width", "strand", "annotation", "geneCh", "geneStart", "geneEnd", "geneLen", "geneStr", "geneID", "distanceToTSS", "MACSqValue", "SRX_ID", "TF", "cell_type")
colnames(fullAnnoPXDNL) = c("chr", "start", "end", "width", "strand", "annotation", "geneCh", "geneStart", "geneEnd", "geneLen", "geneStr", "geneID", "distanceToTSS", "MACSqValue", "SRX_ID", "TF", "cell_type")

write.table(fullAnnoPXDN, file="PXDNFullAnno.bed", sep = "\t", row.names = FALSE)
write.table(fullAnnoPXDNL, file="PXDNLFullAnno.bed", sep = "\t", row.names = FALSE)

## Test for .bed files

peak = readPeakFile("Test2.bed")

peak = readPeakFile("peakPXDNTab.bed")

##### Data Filtering

library(ChIPseeker)
library(dplyr)

## Load in annotated files if needed

#FullAnnoPXDN = read.table("PXDNFullAnno.bed", header = TRUE, quote = "", sep = "\t", comment.char = "")
#FullAnnoPXDNL = read.table("PXDNLFullAnno.bed", header = TRUE, quote = "", sep = "\t", comment.char = "")

## Filter annotated data within 5 kb of TSS (focus on promoter region, first intron and distal intergenic)

### Filter using gene IDs for PXDN and PXDNL
## PXDN: "7837"
## PXDNL: "137902"

filtAnnoPXDN_1 = fullAnnoPXDN %>% 
  filter(geneID == c("7837"))

filtAnnoPXDNL_1 = fullAnnoPXDNL %>% 
  filter(geneID == c("137902"))


## Filter for regions of interest

filtAnnoPXDN = filtAnnoPXDN_1 %>% 
  filter(annotation %in% c("Promoter (1-2kb)", "Promoter (<=1kb)"))

filtAnnoPXDNL = filtAnnoPXDNL_1 %>% 
  filter(annotation %in% c("Promoter (4-5kb)", "Promoter (<=1kb)"))

## Filter for intronic regions

filtAnnoPXDN_intron = filtAnnoPXDN_1 %>% 
  filter(annotation %in% c("Intron (ENST00000252804.9/7837, intron 1 of 22)"))

filtAnnoPXDNL_intron = filtAnnoPXDNL_1 %>% 
  filter(annotation %in% c("Intron (ENST00000356297.5/137902, intron 1 of 22)"))







## Filter for high Macsq values

filtMacsAnnoPXDN = filtAnnoPXDN %>% 
  filter(MACSqValue > 200)

filtMacsAnnoPXDNL = filtAnnoPXDNL %>% 
  filter(MACSqValue > 200)

filtMacsAnnoPXDN_intron = filtAnnoPXDN_intron %>% 
  filter(MACSqValue > 200)

filtMacsAnnoPXDNL_intron = filtAnnoPXDNL_intron %>% 
  filter(MACSqValue > 200)

## remove any unwanted data, like epitope tags

filtMacsAnnoPXDN = filtMacsAnnoPXDN %>% 
  filter(!(TF == "Epitope"))

filtMacsAnnoPXDN_intron = filtMacsAnnoPXDN_intron %>% 
  filter(!(TF == "Epitope"))

## Write into files

write.table(filtMacsAnnoPXDN, file="filtMacsAnnoPXDN.txt", sep = "\t", row.names = FALSE)
write.table(filtMacsAnnoPXDNL, file="filtMacsAnnoPXDNL.txt", sep = "\t", row.names = FALSE)



write.table(filtMacsAnnoPXDN_intron, file="filtMacsAnnoPXDN_intron.txt", sep = "\t", row.names = FALSE)
write.table(filtMacsAnnoPXDNL_intron, file="filtMacsAnnoPXDNL_intron.txt", sep = "\t", row.names = FALSE)




















##### Graph Plotting

library(ggplot2)

# TFs binding to PXDNs promoter region

par(mfrow=c(4,1))

filtMacsAnnoPXDN$n = 1
filtMacsAnnoPXDN$colorCategory <- ifelse(filtMacsAnnoPXDN$MACSqValue > 500, "MACSqscore > 500", "MACSqscore < 500")

TFPXDN = filtMacsAnnoPXDN %>%
  arrange(MACSqValue) %>%
  ggplot(aes(x = TF, n, fill = colorCategory)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("MACSqscore > 500" = "blue", "MACSqscore < 500" = "deepskyblue"))

TFPXDN + 
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Number of studies identifying significant transcription factors \nbinding to PXDN in Cardiovascular Cells") +
  labs(y = "Number of studies", x = "Transcription Factors", fill = "Significance Level") + theme(plot.title = element_text(hjust = 0.5))


# TFs binding to PXDNLs promoter region

filtMacsAnnoPXDNL$n = 1
filtMacsAnnoPXDNL$colorCategory <- ifelse(filtMacsAnnoPXDNL$MACSqValue > 500, "MACSqscore > 500", "MACSqscore < 500")

TFPXDNL = filtMacsAnnoPXDNL %>%
  arrange(MACSqValue) %>%
  ggplot(aes(x = TF, n, fill = colorCategory)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("MACSqscore > 500" = "blue", "MACSqscore < 500" = "deepskyblue"))

TFPXDNL + 
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Number of studies identifying significant transcription factors \nbinding to PXDNL in Cardiovascular Cells") +
  labs(y = "Number of studies", x = "Transcription Factors", fill = "Significance Level") + theme(plot.title = element_text(hjust = 0.5))


# TFs binding to PXDN's first intron

filtMacsAnnoPXDN_intron$n = 1
filtMacsAnnoPXDN_intron$colorCategory <- ifelse(filtMacsAnnoPXDN_intron$MACSqValue > 500, "MACSqscore > 500", "MACSqscore < 500")

TFPXDN = filtMacsAnnoPXDN_intron %>%
  arrange(MACSqValue) %>%
  ggplot(aes(x = TF, n, fill = colorCategory)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("MACSqscore > 500" = "blue", "MACSqscore < 500" = "deepskyblue"))

TFPXDN + 
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Number of studies identifying significant transcription factors \nbinding to PXDN's first Intron in Cardiovascular Cells") +
  labs(y = "Number of studies", x = "Transcription Factors", fill = "Significance Level") + theme(plot.title = element_text(hjust = 0.5))

# TFs in PXDNLs 1st intron


filtMacsAnnoPXDNL_intron$n = 1
filtMacsAnnoPXDNL_intron$colorCategory <- ifelse(filtMacsAnnoPXDNL_intron$MACSqValue > 500, "MACSqscore > 500", "MACSqscore < 500")

TFPXDNL = filtMacsAnnoPXDNL_intron %>%
  arrange(MACSqValue) %>%
  ggplot(aes(x = TF, n, fill = colorCategory)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("MACSqscore > 500" = "blue", "MACSqscore < 500" = "deepskyblue"))

TFPXDNL + 
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Number of studies identifying significant transcription factors \nbinding to PXDNL in Cardiovascular Cells") +
  labs(y = "Number of studies", x = "Transcription Factors", fill = "Significance Level") + theme(plot.title = element_text(hjust = 0.5))


colnames(filtMacsAnnoPXDN)



##### Graph Plotting V2: cell types

library(ggplot2)

# TFs binding to PXDNs promoter region

# first group by cell type

TFPXDN <- filtMacsAnnoPXDN %>%
  group_by(TF, cell_type) %>%
  summarise(Count = n(), .groups = 'drop')

# draw plot, set fill to cell_type
ggplot(TFPXDN, aes(x = TF, y = Count, fill = cell_type)) +
  
  #use geom_bar() to create a stacked chart
  geom_bar(stat = "identity", position = "stack", color = "black") +
  
  #minimal theme makes it easier to edit in pp
  theme_minimal() +
  
  #label axes
  labs(x = "Transcription Factor", y = "Number of Studies", fill = "Cell Type") +
  
  # adjust TF names to slant
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


# TFs binding to PXDNLs promoter region

TFPXDNL <- filtMacsAnnoPXDNL %>%
  group_by(TF, cell_type) %>%
  summarise(Count = n(), .groups = 'drop')

ggplot(TFPXDNL, aes(x = TF, y = Count, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  theme_minimal() +
  labs(x = "Transcription Factor", y = "Number of Studies", fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# TFs binding to PXDN's first intron

TFPXDN_I <- filtMacsAnnoPXDN_intron %>%
  group_by(TF, cell_type) %>%
  summarise(Count = n(), .groups = 'drop')

ggplot(TFPXDN_I, aes(x = TF, y = Count, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  theme_minimal() +
  labs(x = "Transcription Factor", y = "Number of Studies", fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# TFs in PXDNLs 1st intron

TFPXDNL_I <- filtMacsAnnoPXDNL_intron %>%
  group_by(TF, cell_type) %>%
  summarise(Count = n(), .groups = 'drop')

ggplot(TFPXDNL_I, aes(x = TF, y = Count, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  theme_minimal() +
  labs(x = "Transcription Factor", y = "Number of Studies", fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))













## SCRAP CODE FOR FUTURE PROJECTS

# checking annotations in txdb object for 5' UTR
library(AnnotationDbi)
library(GenomicFeatures)




# check valid keytypes of current reference genome

keytypes(txdb)

# check columns of txdb

columns(txdb)

# Get 5' UTRs for all transcripts
five_utrs <- fiveUTRsByTranscript(txdb, use.names = TRUE)

# Check for 3' as well

three_utrs <- threeUTRsByTranscript(txdb, use.names = TRUE)

# Check for the Gene ID type used for UTRs

names(five_utrs)

# create gene_identifiers object using Ensembl gene IDs for PXDN/L

gene_identifiers = c("ENST00000252804.9", "ENST00000356297.5")

# Check if 5' UTRs exist for PXDN and PXDNL

filtered_five_utrs <- five_utrs[names(five_utrs) %in% gene_identifiers]

# repeat for three

filtered_three_utrs <- three_utrs[names(three_utrs) %in% gene_identifiers]


lengths(filtered_five_utrs)

lengths(filtered_three_utrs)








# Find TSS sequence for PXDN and PXDNL

# create object to store gene of interest
gene_numbers = c("7837", "137902")   # gene id for pxdn and pxdnl from NCBI


# get identifiers for PXDN/l

PXDN_name <- select(txdb, keys="7837", keytype="GENEID", columns="GENEID")


all_transcripts <- transcriptsBy(txdb, by="gene")

# Find the TSSs PXDN and PXDNL gene_numbers
# Assuming the gene_id for PXDN is in the returned table and it's the first one (it may return multiple IDs)
all_transcripts
all_transcripts[names(all_transcripts) %in% gene_identifiers]

pxdn_transcripts <- all_transcripts["7837"]

# Get the TSS sites
pxdn_transcripts

# To view the range
pxdn_tss





### comparing number of lower significance peaks to high significance peaks

## first, create a filtered object with only PXDN/L's TFs

PXDN_AllTFs = filtAnnoPXDN_1 %>% 
  filter(!(TF == "Epitope"))

PXDNL_AllTFs = filtAnnoPXDNL_1 %>% 
  filter(!(TF == "Epitope"))

## second, filter to include only the intron and promoter

## Filter for regions of interest

PXDN_pro_int = PXDN_AllTFs %>% 
  filter(annotation %in% c("Promoter (1-2kb)", "Promoter (<=1kb)", "Intron (ENST00000252804.9/7837, intron 1 of 22)"))

PXDNL_pro_int = PXDNL_AllTFs %>% 
  filter(annotation %in% c("Promoter (4-5kb)", "Promoter (<=1kb)", "Intron (ENST00000356297.5/137902, intron 1 of 22)"))

## filter for stringent MACS values

PXDN_pro_int_MACS = PXDN_pro_int %>% 
  filter(MACSqValue > 200)

PXDNL_pro_int_MACS = PXDNL_pro_int %>% 
  filter(MACSqValue > 200)

## view TFs for PXDN

print(unique(PXDN_AllTFs$TF)) #40 TFs, 216 SRX

print(unique(PXDN_pro_int$TF)) #38 TFs, 193 SRX

print(unique(PXDN_pro_int_MACS$TF)) #27 TFs, 161 SRX

## view TFs for PXDNL

print(unique(PXDNL_AllTFs$TF)) #23 TFs, 173 SRX

print(unique(PXDNL_pro_int$TF)) #21 TFs, 152 SRX

print(unique(PXDNL_pro_int_MACS$TF)) #18, 119 SRX

### annotate final results with metadata columns

## get unique SRX IDs from df of interest

PXDN_SRXIDs = unique(PXDN_pro_int_MACS$SRX_ID)

PXDNL_SRXIDs = unique(PXDNL_pro_int_MACS$SRX_ID)

## use the Chrom2/8Cardio.bed files to get treatment info as well as other important metadata

PXDN_TF_metadata = Chrom2Cardio.bed %>%
  
  # filter using SRX IDs
  filter(ID %in% PXDN_SRXIDs) %>%
  
  # select columns of interest from essentialCardio.bed
  select(ID, treatment, TFs, CellTypes) %>%
  
  # group by SRX ID to avoid duplicates
  group_by(ID) %>%
  summarise(
    treatment = first(treatment),  
    TFs = list(unique(TFs)),       
    CellTypes = list(unique(CellTypes))  
  )


# repeat for PXDNL

PXDNL_TF_metadata = Chrom8Cardio.bed %>%
  
  # filter using SRX IDs
  filter(ID %in% PXDNL_SRXIDs) %>%
  
  # select columns of interest from essentialCardio.bed
  select(ID, treatment, TFs, CellTypes) %>%
  
  # group by SRX ID to avoid duplicates
  group_by(ID) %>%
  summarise(
    treatment = first(treatment),  
    TFs = list(unique(TFs)),       
    CellTypes = list(unique(CellTypes))  
  )

# flatten lists so files can be written to csv

PXDN_TF_metadata <- PXDN_TF_metadata %>%
  mutate(across(where(is.list), ~sapply(., toString)))

PXDNL_TF_metadata <- PXDNL_TF_metadata %>%
  mutate(across(where(is.list), ~sapply(., toString)))


# write to tables

write.table(PXDN_TF_metadata, "PXDN_SRX_meta.csv", sep = ";", row.names = FALSE, col.names = TRUE)

write.table(PXDNL_TF_metadata, "PXDNL_SRX_meta.csv", sep = ";", row.names = FALSE, col.names = TRUE)






### Venn diagrams for shared TFs

# load library

library(VennDiagram)

# draw blank grid

grid.newpage()

# first, create lists of the unique TFs for comparisons of interest

unique(filtMacsAnnoPXDN$TF)
PXDN_prom_list = c("RELA",    "SMARCA4", "SMC1A",   "ERG",     "KLF4",    "TAL1",    "RBPJ",    "FLI1",    "CTCF",   
                   "STAG2",   "STAG1",   "RAD21",   "EP300",   "HMGB1",   "TCF21")

unique(filtMacsAnnoPXDN_intron$TF)
PXDN_intron_list = c("RELA",    "SMC1A",   "SMARCA4", "TCF21",   "EP300",   "TAL1",    "KLF4",    "JUN",     "ERG",   
                     "PGR",     "RBPJ",    "NOTCH1",  "MAX",     "MYC",     "CTCF",    "FOXO1",   "FOS",     "NFE2L2", 
                     "JUNB",    "CEBPD",   "FOXC2",   "JUND")

unique(filtMacsAnnoPXDNL$TF)
PXDNL_prom_list = c("IRF1",  "RELA",  "CTCF",  "STAG1")

unique(filtMacsAnnoPXDNL_intron$TF)
PXDNL_intron_list = c("SMC1A",   "CTCF",    "STAG1",   "RAD21",   "STAG2",   "SMARCA4", "RELA",    "FOS",     "TAL1",   
                      "EP300",   "JUN",     "JUNB",    "RBPJ",    "NFE2L2",  "ERG",     "JUND",    "LMNA")  

unique(PXDN_pro_int_MACS$TF)
PXDN_full_list = c("RELA",    "SMC1A",   "SMARCA4", "TCF21",   "EP300",   "TAL1",    "KLF4",    "JUN",     "ERG",    
                   "PGR",     "RBPJ",    "NOTCH1",  "MAX",     "MYC",     "CTCF",    "FOXO1",   "FOS",     "NFE2L2", 
                   "JUNB",    "CEBPD",   "FOXC2",   "JUND",    "FLI1",    "STAG2",   "STAG1",   "RAD21",   "HMGB1")

unique(PXDNL_pro_int_MACS$TF)
PXDNL_full_list = c("SMC1A",   "CTCF",    "STAG1",   "RAD21",   "STAG2",   "SMARCA4", "RELA",    "FOS",     "TAL1",   
                    "EP300",   "JUN",     "JUNB",    "RBPJ",    "NFE2L2",  "ERG",     "JUND",    "LMNA",    "IRF1" )

## show unique for PXDN/L

unique_PXDN = setdiff(PXDN_pro_int_MACS$TF, PXDNL_pro_int_MACS$TF)

unique_PXDNL = setdiff(PXDNL_pro_int_MACS$TF, PXDN_pro_int_MACS$TF)


## show shared

shared_PXDN_PXDNL = c(intersect(PXDN_pro_int_MACS$TF, PXDNL_pro_int_MACS$TF))

## show total number of studies for peaks in promoter/intron

all_SRX_IDs = union(PXDN_pro_int_MACS$SRX_ID, PXDNL_pro_int_MACS$SRX_ID)

## compare PXDN 1st intron and promoter

# create a list of lists to be compared

PXDN_IntProm_venn = list(PXDN_prom = PXDN_prom_list, PXDN_int = PXDN_intron_list)

# store plot in object

venn_plot <- venn.diagram(
  x = PXDN_IntProm_venn,
  filename = NULL,
  category.names = c("PXDN_prom", "PXDN_int"),
  fill = "skyblue2"
)

# draw Venn diagram
grid.draw(venn_plot)

# reset page
grid.newpage()

## compare PXDNL 1st intron and promoter

# create a list of lists to be compared

PXDNL_IntProm_venn = list(PXDNL_prom = PXDNL_prom_list, PXDNL_int = PXDNL_intron_list)

# store plot in object

venn_plot_2 <- venn.diagram(
  x = PXDNL_IntProm_venn,
  filename = NULL,
  category.names = c("PXDNL_prom", "PXDNL_int"),
  fill = "skyblue2"
)

# draw Venn diagram
grid.draw(venn_plot_2)

# reset page
grid.newpage()


## compare PXDN and PXDNL promoter and introns

# create a list of lists to be compared

PXDN_PXDNL_venn = list(PXDN = PXDN_full_list, PXDNL = PXDNL_full_list)

# store plot in object

venn_plot_3 <- venn.diagram(
  x = PXDN_PXDNL_venn,
  filename = NULL,
  category.names = c("PXDN", "PXDNL"),
  fill = "skyblue2"
)

# draw Venn diagram
grid.draw(venn_plot_3)

# reset page
grid.newpage()









## Visualize binding distribution of TFs with MACS 200 in PXDN and PXDNL

## filter for MACS values over 200

filt_PXDN_AllTFs = PXDN_AllTFs %>% 
  filter(MACSqValue > 200)

filt_PXDNL_AllTFs = PXDNL_AllTFs %>% 
  filter(MACSqValue > 200)

## count of unique experiments in PXDN and PXDNL

all_SRX_IDs = union(filt_PXDN_AllTFs$SRX_ID, filt_PXDNL_AllTFs$SRX_ID)


## Isolate columns for chrom, start, end to make a peak file

peakPXDNAllTFs = filt_PXDN_AllTFs %>% select(1:3)
peakPXDNLAllTFs = filt_PXDNL_AllTFs %>% select(1:3)

## Write peak into new bed file
## make sure separators are "\t", and file ends with an emptyline

## !!! access file and make sure there is no header line

write.table(peakPXDNAllTFs, "peakPXDNAllTFs.bed", sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(peakPXDNLAllTFs, "peakPXDNLAllTFs.bed", sep = "\t", row.names = FALSE, col.names = FALSE)

## Annotate Peak files in 5kb radius of TSS


peakAnnoPXDNAll = annotatePeak(
  "peakPXDNAllTFs.bed", 
  tssRegion=c(-10000, 10000), 
  TxDb=txdb, 
  level = "gene", 
  verbose = FALSE)

peakAnnoPXDNLAll = annotatePeak(
  "peakPXDNLAllTFs.bed", 
  tssRegion=c(-10000, 10000), 
  TxDb=txdb, 
  level = "gene", 
  verbose = FALSE)


## Visualise Peak Data as bar

test_PXDN_filt = as.data.frame(peakAnnoPXDNAll)

plotAnnoBar(peakAnnoPXDNAll)

plotAnnoBar(peakAnnoPXDNLAll)

