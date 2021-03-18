# Combine GTEx v8 metadata into single metadata file with information relevant to DE w/ limma
# Sample ID, tissue type, sex, age, RIN (RNA integrity number), post-mortem interval
# Subset only the brain samples in individuals older than 55

# Input files
COUNTS <- "/scratch/mjpete11/limma_tutorial/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct"
META1 <- "/scratch/mjpete11/limma_tutorial/data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
META2 <- "/scratch/mjpete11/limma_tutorial/data/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"

# Output files
BRAIN_META <- "/scratch/mjpete11/limma_tutorial/data/brain_metadata.csv"
SUMMARY_STATS <- "/scratch/mjpete11/limma_tutorial/data/metadata_summary_stats.csv" 

# Libraries
library(tidyverse)
library(data.table)
library(stringr)

# Read in files
# Contains: sample ID, tissue type, RIN, post-mortem info
meta1 <- read_tsv(META1, col_names=TRUE) 
# Contains sex info
meta2 <- read_tsv(META2, col_names=TRUE) 
# Count data
counts <- data.frame(fread(COUNTS))

#_______________________________________________________________________________
# metadata preprocessing 
#_______________________________________________________________________________
# Which columns contain relevant info in meta1 df?
meta1[1,] # SAMPID = 1, SMRIN = 5, SMTSD (detailed tissue type) = 7, SMSISCH (post-mort) = 9

# Subset sample ID, tissue type, RIN, and ischemic time
meta1 <- meta1[,c(1,5,7,9)]

# Rename columns
colnames(meta1) <- c("Sample_ID", "RIN", "Tissue", "Ischemic_Time")

# Add individual ID col;
# grep pattern won't work for the leukemia cell line samples, but I will drop all the cell lines
meta1[["Individual_ID"]] <- str_extract(meta1$Sample_ID, "GTEX-[0-9A-Z]+")

# Reformat tissue names to be contain only '_' between words 
meta1$Tissue <- str_replace_all(meta1$Tissue, c("-" = "", " " = "_", "__" = "_"))

# Replace . to - in colnames
colnames(counts) <- str_replace_all(colnames(counts), pattern = "\\.","-")

# Get list of female IDs
fems <- meta2$SUBJID[which(meta2$SEX==2)]

# Make list of sex of each individual
sex <- with(meta1['Individual_ID'], ifelse(Individual_ID %in% fems, "Female", "Male"))

# Add column containing sex
meta1 <- cbind(meta1, "Sex"=sex)

# Add column containing age (only decade intervals are publically accessible)
meta1$Age <- meta2$AGE[match(meta1$Individual_ID, meta2$SUBJID)]

# Rearrange column order
meta <- meta1 %>% select(Individual_ID, Sex, Age, Tissue, Sample_ID, Ischemic_Time, RIN)

# Drop samples in metadata that do not have count data
select_samples <- colnames(counts)[colnames(counts) %in% meta$Sample_ID]
meta <- meta[meta$Sample_ID %in% select_samples, ]

# Subset sample replicates
Reps <- meta %>%
		    group_by(Individual_ID, Tissue) %>%
		    filter(n() > 1) %>%
			ungroup()

# Samples minus replicates and non-brain_meta tissue 
brain_meta <- meta[grepl("Brain", meta$Tissue),]

# Remove the tissue replciates (Cerebellum/Cortex are replicates of Cerebellar/Frontal_Cortex)
brain_meta <- brain_meta[!grepl("Cerebellum|Brain_Cortex", brain_meta$Tissue),]

# Subset samples >= 50; re-do once I get the metadata with the exact ages and not just decade intervals
brain_meta <- brain_meta %>% filter(Age >= 50)

# Number of samples in meta
nrow(brain_meta) # 1,831

# Function to summarise percent missing values
Missing <- function(x){
		   x %>%
		 	 select(everything()) %>%
		  	 summarise_all(list(~sum(is.na(.))))/nrow(x) * 100
}

# How many samples have missing values?
Missing(meta)

# Drop rows missing values
meta <- meta %>% drop_na()

# Check that all samples missing values were removed
Missing(meta)

# Quick view
head(brain_meta)
tail(brain_meta)

# Summary stats: how many tissues per sex
Summary_Stat <- function(x){
	x %>% group_by(Tissue, Sex) %>% tally()
}

# Print summary statistics 
print(Summary_Stat(brain_meta), n=22)

#______________________________________________________________________________________________________
# Write to file 
#______________________________________________________________________________________________________
# summary stats
write.csv((as.data.frame(Summary_Stat(brain_meta))), SUMMARY_STATS)

# metadata files
write.csv(brain_meta, BRAIN_META, row.names=FALSE)
