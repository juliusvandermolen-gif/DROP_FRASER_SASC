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
dataset_dir     <- NA       # only needed for real data (where DROP saved the FraserDataSet)
annotation_name <- "raw-local" # Name of the annotation. Needs sampleID and BAM files
type            <- c("jaccard", "psi5", "psi3")
outdir_plots    <- "fraser_plots_out"
dataset_output_dir <- "ds_output" # Directory where the dataset is stored.
plot_version <- 1
min_cut_expression <- 10 # minimum reads in at least one sample. Paper uses 20


# Data. For now the real data and not the test data
# %%
anno_sample <- fread("path/to/annotation_sample.csv") # Annotation sample file. 
fds <- FraserDataSet(workingDir = dataset_output_dir, colData = anno_sample, name = annotation_name)


# Counting reads in split and non-split reads
# %%
message("Step 1: Counting reads in split and non-split reads")
fds <- countRNAData(fds)

# Calculate PSI/Jaccard metrics
# %%
message("Step 2: Calculate PSI/Jaccard values")
fds <- calculatePSIValues(fds, types = type)

# Filtering junctions based on expression and variability
# %%
dir.create(outdir_plots, showWarnings = FALSE, recursive = TRUE) #directory for output plots
message("Step 3: Filtering junctions based on expression and variability")
fds <- filterExpressionAndVariability(fds, minDeltaPsi = 0, minExpressionInOneSample = min_cut_expression, filter = TRUE)
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
  
  #plotting the optimal dimension selection
  pdf(file.path(outdir_plots, paste0(plot_version, "_Optimal_Dimension_Selection_", t, ".pdf")),
      width = 6, height = 5)
  print(plotEncDimSearch(fds, type = t, plotType = if (use_oht) "sv" else "auc"))
  dev.off()

  q_list_dim[[t]] <- bestQ(fds, type = t)
  message("Optimal dimension (q) for ", t, ": ", q_list_dim[[t]])
}

q_vec <- unlist(q_list_dim) #List change to vector

# Model fitting with FRASER(). Most intens task
# %%
message("Step 7: Model fitting with FRASER()")
fds <- FRASER(fds, type = type, q = q_vec, implementation = "AE", iterations = 15) #Model run

for (t in type) {
  pdf(file.path(outdir_plots, paste0(plot_version, "_heatmap_after_correction_", t, ".pdf")), width = 6, height = 5)
  print(plotCountCorHeatmap(fds, type = t, normalized = TRUE, logit = TRUE))
  dev.off()
}


# Results. Usefull for specific cutoffs
# %%
message("Step 8: Results: Extracting results aberrant splicing events")
res <- results(fds, psiType = type, padjCutoff=0.05, deltaPsiCutoff=0.1, minCount = 10) # Cutoffs arguments


# Sample plots
message("Step 10: QC and overview plots")

plotAberrantPerSample(fds, type = type, padjCutoff = 0.05, deltaPsiCutoff = 0.3)

for (t in type){
  pdf(file.path(outdir_plots, paste0(plot_version, "_Volcano_sample1_", t, ".pdf")),
    width = 6, height = 5)
  print(plotVolcano(fds, sampleID = "sample1", type = t, aggregate = TRUE))
  dev.off()

  pdf(file.path(outdir_plots, paste0(plot_version, "_QQplot_", t, ".pdf")),
      width = 6, height = 5)
  print(plotQQ(fds, type = t, aggregate = TRUE, global = TRUE))
  dev.off()
}