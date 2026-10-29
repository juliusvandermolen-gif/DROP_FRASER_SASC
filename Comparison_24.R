# %%
suppressPackageStartupMessages({
  library(FRASER)
})
library(data.table)
# %%
dir <- "/exports/sascstudent/jmvandermolen1/Developer/drop_demo_fraser2/Output_24_samples/processed_results/aberrant_splicing/results/GRCh38_104/fraser/"
NMD_results            <- fread(paste0(dir, "Lymphocyte_NMD/results.tsv"))
NMD_results_per_junction <- fread(paste0(dir, "Lymphocyte_NMD/results_per_junction.tsv"))
NMD_inhibition_results <- fread(paste0(dir, "Lymphocyte_NMD_inhibit/results.tsv"))
NMD_inhibition_results_per_junction <- fread(paste0(dir, "Lymphocyte_NMD_inhibit/results_per_junction.tsv"))

dataset_dir     <- "/exports/sascstudent/jmvandermolen1/Developer/drop_demo_fraser2/Output_24_samples/processed_data/aberrant_splicing/datasets"
annotation_name_NMD <- "raw-local-Lymphocyte_NMD"
fds_NMD <- loadFraserDataSet(dir = dataset_dir, name = annotation_name_NMD) #Dataset


# %%
kolom_NMD_results <- colnames(NMD_results)
fds

#Count check
# %% 
summary(NMD_results$totalCounts)
summary(NMD_inhibition_results$totalCounts)

#Number of samples per junction
# %%
summary(NMD_results_per_junction$numSamplesPerJunc)
summary(NMD_inhibition_results_per_junction$numSamplesPerJunc)

#Checking how many junctions(in results) are found in each sample
# %%
table(NMD_results_per_junction$numSamplesPerJunc)
table(NMD_inhibition_results_per_junction$numSamplesPerJunc)

NMD_results[, .N, by = sampleID][order(-N)]
NMD_inhibition_results[, .N, by = sampleID][order(-N)]


# %%
s <- "R1.26.00018"
results_R1_26_00018 <- NMD_results[sampleID == s]

summary(results_R1_26_00018$totalCounts)                         # dit sample
summary(NMD_results[sampleID != s]$totalCounts) # de rest
results_R1_26_00018[, .N, by = type]
results_R1_26_00018[, .N, by = annotatedJunction]
summary(results_R1_26_00018$nonsplitProportion)



# %%
#Aberant plot per sample
