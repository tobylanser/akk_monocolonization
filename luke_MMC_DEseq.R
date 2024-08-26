# References: https://hbctraining.github.io/DGE_workshop/lessons/08_DGE_LRT.html

## Gene-level differential expression analysis using DESeq2
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

list.of.packages <- c('DESeq2','openxlsx','BiocGenerics','tximport','ggplot2','biomaRt','S4Vectors','apeglm','tidyverse')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)

lapply(list.of.packages, require, character.only = TRUE)


##---------------------------------FILE PREPARATION----------------------------------##
library(BiocFileCache)
library(dbplyr)
## useMart is a function used to connect to the selected BioMart database and dataset (ListEnsemblArchives() to list host sites)
# hsapiens_gene_ensembl for humans
# mmusculus_gene_ensembl for mice
ensembl_ms_mart <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
## getBM is a function used to retrieve information from the BioMart database
ensembl_df <- getBM(attributes = c("ensembl_transcript_id_version", "ensembl_gene_id_version", "external_gene_name"), mart = ensembl_ms_mart) %>%
  as_tibble() %>%
  mutate(noVersion = sub('\\..*', '', ensembl_transcript_id_version))

library(dplyr)
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

# micro_VE_vs_BC
micro_VE_vs_BC <- subset(sample_data, celltype == 'Microglia' & testvar %in% c('VE', 'BC'))

# micro_AT_vs_H3
micro_AT_vs_H3 <- subset(sample_data, celltype == 'Microglia' & testvar %in% c('AT', 'H3'))

# micro_AT_vs_BC
micro_AT_vs_BC <- subset(sample_data, celltype == 'Microglia' & testvar %in% c('AT', 'BC'))

# micro_H3_vs_BC
micro_H3_vs_BC <- subset(sample_data, celltype == 'Microglia' & testvar %in% c('H3', 'BC'))

#################### Astrocytes ###################

# astro_VE_vs_AT
astro_VE_vs_AT <- subset(sample_data, celltype == 'Astrocyte' & testvar %in% c('VE', 'AT'))

# astro_VE_vs_H3
astro_VE_vs_H3 <- subset(sample_data, celltype == 'Astrocyte' & testvar %in% c('VE', 'H3'))

# astro_VE_vs_BC
astro_VE_vs_BC <- subset(sample_data, celltype == 'Astrocyte' & testvar %in% c('VE', 'BC'))

# astro_AT_vs_H3
astro_AT_vs_H3 <- subset(sample_data, celltype == 'Astrocyte' & testvar %in% c('AT', 'H3'))

# astro_AT_vs_BC
astro_AT_vs_BC <- subset(sample_data, celltype == 'Astrocyte' & testvar %in% c('AT', 'BC'))

# astro_H3_vs_BC
astro_H3_vs_BC <- subset(sample_data, celltype == 'Astrocyte' & testvar %in% c('H3', 'BC'))

astro_data <- subset(sample_data, celltype == 'Astrocyte')

micro_data <- subset(sample_data, celltype == 'Microglia')

## List all directories containing data
#*** All the gene quantification files are in the folder "quant_files", make sure metafile sample order matches order of quant files.
all_files <- list.files("../data", full.names = T, pattern=NULL, all.files=FALSE)
quant_files <- file.path(all_files, "quant.sf")

prefix <- astro_H3_vs_BC

position_list <- paste0(prefix$position)
sample_files <- grep(position_list, quant_files, value = T)

sample_files <- file.path("../data", position_list, "/quant.sf")
## QC check - Order of metafile samples must match quant file order!  Fix metadata file if not!
all(file.exists(quant_files))

sample_files
astro_H3_vs_BC

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
file_prefix = 'astro_H3_vs_BC'

## Create DESeq2Dataset object with variables of interest.
dds <- DESeqDataSetFromTximport(txi, colData = prefix, design = ~ testvar)

#prefiltering on minimum of 5 reads (NOT REQUIRED AS DESEQ2 OPTIMIZES AND FILTERS AUTOMATICALLY)
# dds <- estimateSizeFactors(dds)
# idx <- rowMeans(counts(dds, normalized = TRUE)) >= 5
# dds <- dds[idx,]

## Run DESeq analysis to gather differential expression ../results
# Run DESeq (LRT)
#dds_run <- DESeq(dds, test = 'LRT', reduced = ~ 1)
# Run DESeq (Wald)
dds_run <- DESeq(dds, betaPrior = F)

# View names of estimated effects
resultsNames(dds_run)

# Create table of effects showing log2fold, standard error, test stats, and p-vals
# LRT
#dds_result <- ../results(dds_run, pAdjustMethod = 'BH', independentFiltering = T, cooksCutoff = T)
# Wald
astro_H3_vs_BC
dds_result <- results(dds_run, contrast = c('testvar', 'BC', 'H3'), independentFiltering = T, pAdjustMethod = 'BH')
dds_result <- lfcShrink(dds_run, contrast = c('testvar', 'BC', 'H3'), res = dds_result, type = 'normal')

summary(dds_result, alpha = 0.05)

## View Total number of normalized counts per sample
dds_run <- estimateSizeFactors(dds_run)
DS_norm_counts <- counts(dds_run, normalized = TRUE)
#colSums(DS_norm_counts)

## Plot dispersion estimates
plotDispEsts(dds_run)
#plotCounts(dds_run, 'Fus', intgroup = 'condition')

# Transform counts for data visualization
rld <- vst(dds_run, blind=TRUE)

# Plot PCA 
plotPCA(rld, intgroup=c('testvar', "experiment"))

# Add nametags
z <- plotPCA(rld, intgroup=c('testvar'))
PCA <- z + geom_label(aes(label = prefix$testvar))
PCA
ggsave(paste0("../plots/PCA_for_review/", file_prefix, "_PCA.pdf"), device = "pdf", 
       width = 10, height = 10, units = "in", dpi = 600, plot = PCA)
# Build significant gene table and extract list of sorted DE genes.
sig_res <- dds_result %>%
  data.frame() %>%
  rownames_to_column(var='gene') %>%    # makes a 'gene' column using column1
  as_tibble %>%
  filter(padj < 0.10) %>%
  arrange(padj)
sorted_DEGenes <- sig_res$gene

# Export significant gene count data
name_list <- c('gene', (paste0(dds_run$celltype,"_",dds_run$tissue,'_',dds_run$age,'_',dds_run$microbiome,'_',dds_run$position)))

sig_export <- DS_norm_counts %>%
  data.frame() %>%
  rownames_to_column(var = 'gene_name') %>%
  as_tibble() %>%
  `colnames<-`(name_list) %>%
  semi_join(sig_res) %>%
  slice(match(sorted_DEGenes, gene))

# Write DE genes to file
#write.xlsx(sig_export, file = paste0('../results/', file_prefix, 'DEGene_counts.xlsx'), overwrite = T)
write.xlsx(sig_export, file = paste0('../results/', file_prefix, 'DEGene_counts_pval05.xlsx'), overwrite = T)

# Data subsetting and gene filtering (if necessary)
unwanted_genes = paste(c('^Gm', '^mt-', '^Vmn', '^Rpl', '^Rps', '^Olfr'), collapse = '|')

# Column Names
rownames <- sig_export$gene
sig_export[,1] <- NULL
rownames(sig_export) <- rownames
check_cols <- c('gene', colnames(sig_export))

# Normalized counts
check_counts <- DS_norm_counts %>%
  as.data.frame() %>%
  rownames_to_column('gene')%>% 
  as_tibble %>%
  `colnames<-`(check_cols)

# ../results
check_res <- dds_result %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  filter(!str_detect(gene, unwanted_genes))

# TPM File
my_genes <- rownames(txi$abundance)
my_genes_ann <- ensembl_df[match(my_genes,ensembl_df$external_gene_name),]
test_TPM <- cbind(my_genes_ann$external_gene_name,txi$abundance) %>%
  as.data.frame() %>% 
  as_tibble() %>%
  `colnames<-`(check_cols)


data_final <- subset(check_counts, gene %in% sig_res$gene)


# Output files
write.xlsx(check_res, file = paste0('../results/', file_prefix, 'statistics.xlsx'), overwrite = T)
write.xlsx(check_counts, file = paste0('../results/', file_prefix, 'DS_counts.xlsx'), overwrite = T)
write.xlsx(data_final, file = paste0('../results/', file_prefix, 'DEGene_counts.xlsx'), rowNames = T, overwrite = T)
write.xlsx(sig_res, file = paste0('../results/', file_prefix, 'DEGene_statistics.xlsx'), overwrite = T)
write.xlsx(test_TPM,file = paste0('../results/', file_prefix, 'TPM.xlsx'), overwrite = T)



