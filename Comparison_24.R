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

dataset_dir     <- "/exports/sascstudent/jmvandermolen1/Developer/drop_demo_fraser2/Output_24_samples/processed_results/aberrant_splicing/datasets"
annotation_name_NMD <- "Lymphocyte_NMD--GRCh38_104"
fds_NMD <- loadFraserDataSet(dir = dataset_dir, name = annotation_name_NMD) #Dataset

dir_plots_24 <- "/exports/sascstudent/jmvandermolen1/Developer/drop_demo_fraser2/Output_24_samples/plots_24_run"


# %%
kolom_NMD_results <- colnames(NMD_results)
assayNames(fds_NMD)

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
#Create directory
dir.create(dir_plots_24, showWarnings = FALSE, recursive = TRUE)

# %%
type_metric <- "jaccard"
#Aberant plot per sample
pdf(file.path(dir_plots_24, "05_aberrant_per_sample_jaccard.pdf"), width = 6, height = 5)
tryCatch({
  p <- plotAberrantPerSample(fds_NMD, type = type_metric, padjCutoff = 0.1)
  print(p)
}, finally = dev.off())


# %%
# R log .drop directory
rds_log <- readRDS("/exports/sascstudent/jmvandermolen1/Developer/drop_demo_fraser2/.drop/tmp/AS/Lymphocyte_NMD_inhibit--GRCh38_104/FRASER_summary.Rds")
