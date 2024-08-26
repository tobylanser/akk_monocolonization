setwd("~/Partners HealthCare Dropbox/Toby Lanser/luke-RNAseq/Monocolonization Study/tobys_analysis/scripts")


list.of.packages <- c('DESeq2','openxlsx','BiocGenerics','tximport','ggplot2','biomaRt','S4Vectors','apeglm','tidyverse')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)

lapply(list.of.packages, require, character.only = TRUE)
library(BiocFileCache)
library(dbplyr)
library(dplyr)

# devtools::install_version("dbplyr", version = "2.3.4")
# BiocManager::install("biomaRt", force = T)

##---------------------------------FILE PREPARATION----------------------------------##

## useMart is a function used to connect to the selected BioMart database and dataset (ListEnsemblArchives() to list host sites)
# hsapiens_gene_ensembl for humans
# mmusculus_gene_ensembl for mice
ensembl_ms_mart <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
## getBM is a function used to retrieve information from the BioMart database
ensembl_df <- getBM(attributes = c("ensembl_transcript_id_version", "ensembl_gene_id_version", "external_gene_name"), mart = ensembl_ms_mart) %>%
  as_tibble() %>%
  mutate(noVersion = sub('\\..*', '', ensembl_transcript_id_version))


# Build tx2gene table from ensembl_df
tx2gene <- ensembl_df %>%
  dplyr::select(noVersion, external_gene_name) %>%
  rename(c('Transcript' = "noVersion" , "Gene" = "external_gene_name" )) %>%
  data.frame()

## Create a metadata map using your metadata file.
metadata <- read.csv("../metadata_MMCfull.csv")
sample_data <- data.frame(position = metadata$Position, celltype = metadata$CellType,
                          testvar = metadata$TestVar, experiment = metadata$Experiment)

# Outlier filtering  and subsetting(if necessary)
# outliers = c('F11','G01','E10','G04','E11')
# monocyte_outs = c('C04','C01','A10','B07','B06','A08')
# sample_data <- sample_data %>%            # Filter sample_data sheet
#   filter(!position %in% c(outliers,monocyte_outs)) %>%
#   filter(celltype != 'blank')


#################### conditions ###################
sample_data <- subset(sample_data, experiment != "MM04")

outliers = c("A02", 'B01','C01','C11','B09','B08','H04',
             "F04", "F11", "F03", "F02")
sample_data <- sample_data %>%            # Filter sample_data sheet
  #filter(experiment != 'Blank') %>%
  filter(!position %in% c(outliers))

#################### Microglia ###################

# micro_VE_vs_AT
micro_VE_vs_AT <- subset(sample_data, celltype == 'Microglia' & testvar %in% c('VE', 'AT'))

# micro_VE_vs_H3
micro_VE_vs_H3 <- subset(sample_data, celltype == 'Microglia' & testvar %in% c('VE', 'H3'))

# micro_H3_vs_BC
micro_H3_vs_BC <- subset(sample_data, celltype == 'Microglia' & testvar %in% c('H3', 'BC'))

#################### Astrocytes ###################

# astro_VE_vs_AT
astro_VE_vs_AT <- subset(sample_data, celltype == 'Astrocyte' & testvar %in% c('VE', 'AT'))

# astro_VE_vs_H3
astro_VE_vs_AT <- subset(sample_data, celltype == 'Astrocyte' & testvar %in% c('VE', 'H3'))

# astro_H3_vs_BC
astro_H3_vs_BC <- subset(sample_data, celltype == 'Astrocyte' & testvar %in% c('H3', 'BC'))


## List all directories containing data
#*** All the gene quantification files are in the folder "quant_files", make sure metafile sample order matches order of quant files.
all_files <- list.files("../data", full.names = T, pattern=NULL, all.files=FALSE)
quant_files <- file.path(all_files, "quant.sf")

prefix <- micro_VE_vs_AT

position_list <- paste0(prefix$position)
sample_files <- grep(position_list, quant_files, value = T)

sample_files <- file.path("../data", position_list, "/quant.sf")
## QC check - Order of metafile samples must match quant file order!  Fix metadata file if not!
all(file.exists(quant_files))

sample_files

# Import Transcript quantification files
txi<-tximport(files=sample_files,type="salmon", tx2gene = tx2gene, ignoreTxVersion = T)

## Look at the counts and set the counts as a matrix file.
txicounts<-as.matrix(txi$counts)

## Write the counts to file, round them up, and convert to data.frame and remove txicounts file.
count_data <- txi$counts %>%
  round() %>%
  data.frame()
rm(txicounts)


#-----------------------------------ANALYSIS BEGINS HERE---------------------------------#
# Set file prefix for inputs and outputs
file_prefix = 'micro_VE_vs_AT'

## Create DESeq2Dataset object with variables of interest.
dds <- DESeqDataSetFromTximport(txi, colData = prefix, design = ~ testvar)

# Run DESeq (Wald)
dds_run <- DESeq(dds, betaPrior = F)

dds_run <- estimateSizeFactors(dds_run)
# View names of estimated effects
resultsNames(dds_run)

rld <- vst(dds_run, blind=TRUE)

# Plot PCA 
plotPCA(rld, intgroup=c('testvar', "experiment"))

# Add nametags
z <- plotPCA(rld, intgroup=c('testvar'))
PCA <- z + geom_label(aes(label = prefix$testvar))
PCA

PCs <- PCA$data

write.csv(PCs, file=paste0("../review_material/", file_prefix, "_pcs_coordinates.csv")) # extract the PCs coordinates

?prcomp





