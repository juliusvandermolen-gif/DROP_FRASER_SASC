#!/usr/bin/env Rscript

# Import libraries
# %% cel
suppressPackageStartupMessages({
  library(FRASER)
  library(data.table)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene) 
  library(org.Hs.eg.db)
})

# Parameters: Directory's, dataset
# %%
use_test_data   <- FALSE     # TRUE = testdata, FALSE = real data
dataset_dir     <- NA       # only needed for real data (where DROP saved the FraserDataSet)
annotation_name <- "raw-local" # Name of the annotation. Needs sampleID and BAM files
type            <- "jaccard"
q_dim           <- 3
outdir_plots    <- "fraser_plots_out"
dataset_output_dir <- "ds_output" # Directory where the dataset is stored.


# Data. For now the real data and not the test data
# %%
anno_sample <- fread("path/to/annotation_sample.csv") # Annotation sample file. 
fds <- FraserDataSet(workingDir = dataset_output_dir, colData = anno_sample, name = annotation_name)


# Counting reads in split and non-split reads
# %%
message("Step 1: Counting reads in split and non-split reads")
fds <- countRNASeq(fds)

# Calculate PSI/Jaccard metrics
# %%
message("Step 2: Calculate PSI/Jaccard values")
fds <- calculatePSIValues(fds)

# Filtering junctions based on expression and variability
# %%
message("Step 3: Filtering junctions based on expression and variability")
fds <- filterExpressionAndVariability(fds, minDeltaPsi = 0, minExpressionInOneSample = 10, filter = TRUE)
plotFilterExpression(fds)

# Annotation introns with gene symbols. Check genome version that is compatible with BAM files. For example, if BAM files are aligned to hg38, use the corresponding annotation.
# Change annotation via library(TxDb.Hsapiens.UCSC.hg38.knownGene) and library(org.Hs.eg.db) for hg38.
# %%
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
orgdb <- org.Hs.eg.db
message("Step 4: Annotating ranges with TxDb...")
fds <- annotateRangesWithTxDb(fds, txdb=txdb, orgDb=orgDb)

# Sample co-variation before fitting the model. To find how samples correlate before correction to get rid of non-biological variation. 
# This is important to check if there are any batch effects or other confounding factors that might affect the analysis.
# %%
message("Step 5: Sample co-variation before fitting the model")
plotCountCorHeatmap(fds, type = type, normalized = FALSE, logit = TRUE)