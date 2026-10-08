#!/usr/bin/env Rscript

# Import libraries
# %% cel
suppressPackageStartupMessages({
  library(FRASER)
  library(data.table)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene) # Annotation for gene symbols
  library(org.Hs.eg.db)
})

# Parameters: Directory's, dataset
# %%
use_test_data   <- FALSE     # TRUE = testdata, FALSE = real data
dataset_dir     <- NA       # only needed for real data (where DROP saved the FraserDataSet)
annotation_name <- "raw-local" # Name of the annotation. Needs sampleID and BAM files
type            <- c("jaccard", "psi5", "psi3")
q_dim           <- 3
outdir_plots    <- "fraser_plots_out"
dataset_output_dir <- "ds_output" # Directory where the dataset is stored.
plot_version <- 1


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
fds <- calculatePSIValues(fds, type = type)

# Filtering junctions based on expression and variability
# %%
dir.create(outdir_plots, showWarnings = FALSE, recursive = TRUE) #directory for output plots
message("Step 3: Filtering junctions based on expression and variability")
fds <- filterExpressionAndVariability(fds, minDeltaPsi = 0, minExpressionInOneSample = 10, filter = TRUE)
pdf(file.path(outdir_plots, paste0(plot_version, "_filter_expression.pdf")), width = 6, height = 5)
plotFilterExpression(fds)
dev.off()

# Annotation introns with gene symbols. Check genome version that is compatible with BAM files. For example, if BAM files are aligned to hg38, use the corresponding annotation.
# Change annotation via library(TxDb.Hsapiens.UCSC.hg38.knownGene) and library(org.Hs.eg.db) for hg38.
# %%
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
orgdb <- org.Hs.eg.db
message("Step 4: Annotating ranges with TxDb...")
fds <- annotateRangesWithTxDb(fds, txdb=txdb, orgDb=orgdb)

# Sample co-variation before fitting the model. To find how samples correlate before correction to get rid of non-biological variation. 
# This is important to check if there are any batch effects or other confounding factors that might affect the analysis.
# %%
message("Step 5: Sample co-variation before fitting the model visualization")
for (t in type) {
  pdf(file.path(outdir_plots, paste0(plot_version, "_heatmap_before_correction_", t, ".pdf")), width = 6, height = 5)
  print(plotCountCorHeatmap(fds, type = t, normalized = FALSE, logit = TRUE))
  dev.off()
}

# Optimal dimension selection latent space. Important before fitting model.
# %%
message("Step 6: Optimal dimension selection latent space")
set.seed(42)

q_list_dim <- list()
for (t in type){
  use_oht <- (t == "jaccard") # Use OHT for jaccard, not for psi5 and psi3. Returns True or False with checking jaccard
  fds <- estimateBestQ(fds, type = t, useOHT = use_oht)
  q_dim <- getBestQ(fds, type = t)
  q_list_dim[[t]] <- q_dim
  message("Optimal dimension (q) selected with ", t, ": ", q_dim)
  
  #plotting the optimal dimension selection
  pdf(file.path(outdir_plots, paste0(plot_version, "_Optimal_Dimension_Selection_", t, ".pdf")),
      width = 6, height = 5)
  print(plotEncDimSearch(fds, type = t, plotType = if (use_oht) "sv" else "auc"))
  dev.off()

}

pdf(file.path(outdir_plots, paste0(plot_version, "_Optimal_Dimension_Selection.pdf")), width = 6, height = 5)
plotEncDimSearch(fds, type = type, plotType = "auc")
dev.off()


# Model fitting with FRASER(). Most intens task
# %%
message("Step 7: Model fitting with FRASER()")
fds <- FRASER(fds, type = type, q = c(jaccard = q_dim), implementation = "AE", iterations = 15)

pdf(file.path(outdir_plots, paste0(plot_version, "_Heatmap_after_correction.pdf")), width = 6, height = 5)
plotCountCorHeatmap(fds, type = type, logit = TRUE, normalized = TRUE)
dev.off()


# Results. Usefull for specific cutoffs
# %%
message("Step 8: Results: Extracting results aberrant splicing events")
res <- results(fds, padjCutoff=0.05, deltaPsiCutoff=0.1, minCount = 10) # Cutoffs arguments


# Sample plots
message("Step 10: QC and overview plots")
plotAberrantPerSample(fds, padjCutoff = 0.05, deltaPsiCutoff = 0.3)
plotVolcano(fds, sampleID = "sample1", type = type, aggregate = TRUE)
plotQQ(fds, aggregate = TRUE, global = TRUE)