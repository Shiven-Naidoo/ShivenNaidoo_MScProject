#### Data Cleaning
## Load data cleaning packages 

install.packages("tidyverse")
library(dplyr)  # do not run at the same time as genomicfeatures, as the select() function is used in both
library(tidyr)
library(stringr)

## Set working directory for .bed files

getwd()

## Clean downloaded .bed data. sample metadata is compressed into 1 column. Goal = Extract and sort

cardio.bed = read.table("D:/Project data/MSc data/ChipSeq/Oth.CDV.05.AllAg.AllCell.bed", header = TRUE, quote = "", sep = "\t", comment.char = "")

## Look at data to identify column containing sample metadata

str(cardio.bed)

colnames(cardio.bed) <- c("Chromosome", "Start", "End", "SampleMetadata", "-10Log10(MACS2 Q-value)", "Dot", "Column2", "Column3", "ColorCode")

## Isolate metadata

SampleMetadataCol = cardio.bed[, 4]

## Calculate max number of possible columns if split using ";" as delimiter

colMax = max(stringr::str_count(SampleMetadataCol, ";")) + 1

str(colMax)

## Split metadata into columns

splitMetadata = str_split_fixed(SampleMetadataCol, ";", colMax)

## Assign columns with important metadata to R objects

ID = splitMetadata[, 1]
TF_CellType = splitMetadata[, 2]










## Prototype, treatment identifier
## challenge, many samples are missing treatment information and may need to be manually reviewed

## store treatment info in new column
## isolate treatment by using "treatment=" string, and any character other than ";"

cardio.bed$treatment <- str_extract(cardio.bed$SampleMetadata, "treatment=([^;]+)")





## Remove special characters

specialCharacters = c("%", "20", "(", ")", "@")

Cleaned1_TF_CellType = gsub("%", " ", as.character(TF_CellType))
Cleaned2_TF_CellType = gsub("20", "", as.character(Cleaned1_TF_CellType))
Cleaned3_TF_CellType = gsub("@", "", as.character(Cleaned2_TF_CellType))

## Split TF and Cell type into 2 different columns

splitTF_CellType = str_split_fixed(Cleaned3_TF_CellType, " ", 2)

colnames(splitTF_CellType) <- c("TF_Raw", "CellType")

## Check columns

splitTF_CellType[, 1]

## Clean columns of unwanted characters

Clean_ID = gsub("ID=", "", as.character(ID))

Clean_TFs = gsub("Name=", "", as.character(splitTF_CellType[, 1]))

Clean_CellTypes = gsub("\\(|\\)| ", "", as.character(splitTF_CellType[, 2]))

## Check cleaned columns

Clean_ID

Clean_TFs

Clean_CellTypes

## Remove non-essential columns from cardio.bed

essentialCardio.bed = dplyr::select(cardio.bed, -c(SampleMetadata, Dot, Column2, Column3, ColorCode))

## Check data

essentialCardio.bed

## Combine with processed data for IDs, TFs and CellTypes

essentialCardio.bed[c("ID", "TFs", "CellTypes")] = c(Clean_ID, Clean_TFs, Clean_CellTypes)

## Clean treatment column

specialCharacters = c("%20", "treatment=")

for (char in specialCharacters) {
  essentialCardio.bed$treatment = gsub(char, " ", essentialCardio.bed$treatment)
}


## Check cleaned data

essentialCardio.bed

## Identify unique cell types

unique(essentialCardio.bed$CellTypes)

## 34 unique cell types

## Identify unique TFs

unique(essentialCardio.bed$TFs)

## 63 unique TFs excluding epitope tags

## Store data in separate .bed files for PXDN and PXDNL based on chromosome position
## Select regions 50kb upstream and downstream of 5' and 3' UTRs respectively
## PXDN coords chr2:1629887-1746901
## PXDNl coords chr8:51319577-51809445 

## Isolate PXDN using dplyr function
## Isolate chr2, then filter by upper and lower borders for gene

Chrom2Cardio.bed = filter(essentialCardio.bed, Chromosome == "chr2")
PXDNUpperFilter.bed = filter(Chrom2Cardio.bed, Start >= 1579887)
PXDNFiltered.bed = filter(PXDNUpperFilter.bed, End <= 1796901)

## Save as file outside of R

PXDNFiltered.bed

write.table(PXDNFiltered.bed, "PXDNFiltered2.bed", row.names = FALSE, col.names = FALSE)

## Isolate PXDNL using dplyr function
## Isolate chr8, then filter by upper and lower borders for gene

Chrom8Cardio.bed = filter(essentialCardio.bed, Chromosome == "chr8")
PXDNLUpperFilter.bed = filter(Chrom8Cardio.bed, Start >= 51219577)
PXDNLFiltered.bed = filter(PXDNLUpperFilter.bed, End <= 51864445)

## function for filtering

## function should take chromosome, start and end range as inputs

filter_bed = function(bed, chromosome, start, end) {
  
  # Filter for chromosome
  chrom_filtered_data = filter(bed, Chromosome == chromosome)
  
  # Filter resulting df for the start position
  start_filtered_data = filter(chrom_filtered_data, Start >= start)
  
  # Filter for the end position
  end_filtered_data = filter(start_filtered_data, End <= end)
  
  # Return the filtered data
  return(end_filtered_data)
}

# use to create 10kb files for PXDN and PXDNL

filtered_PXDN_10kb = filter_bed(essentialCardio.bed, "chr2", 1736901, 1756901)

filtered_PXDNL_10kb = filter_bed(essentialCardio.bed, "chr8", 51799445, 51819445)

# view unique TFs

unique(filtered_PXDN_10kb$TFs)

unique(filtered_PXDNL_10kb$TFs)


## Save as file outside of R

PXDNLFiltered.bed

write.table(PXDNLFiltered.bed, "TFsPXDNL2.bed", row.names = FALSE, col.names = FALSE)




