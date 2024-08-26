
setwd("/Users/ly794/Partners HealthCare Dropbox/Luke Schwerdtfeger/luke-RNAseq/Monocolonization Study/tobys_analysis/review_material/PC_coordinates/")


install.packages("vegan")
library(vegan)

# IMPORT DATA TABLE and DEFINE GROUPS #
data <- read.table("astro_H3_vs_BC_pcs_coordinates-luke.txt", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
groups <- factor(c("1", "2"))

# RUN ANALYSIS TO DETERMINE NORMALCY OF DATA DISTRIBUTION #
# Shapiro-Wilk Test #
shapiro_results <- apply(data, 2, shapiro.test)
shapiro_results

#if a majority of PCs give shapiro test p value of < 0.05, then use method = manhatten #
#if a majority of PCs give shapiro test p value of > 0.05, then use method = euclidean #


# Perform PERMANOVA using adonis2 #
set.seed(123)
result <- adonis2(data ~ groups, 
                  data = data, 
                  permutations = 999, 
                  #method = "euclidean",
                  method = "manhattan",
                  sqrt.dist = FALSE, 
                  add = FALSE, 
                  by = "terms",
                  parallel = getOption("mc.cores"),
                  na.action = na.fail,
                  strata = NULL)
print(result)

