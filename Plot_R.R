#!/usr/bin/env Rscript
# ==============================================================================
# Generate plots for FRASER QC and results.
# ==============================================================================

# Library's
suppressPackageStartupMessages({
  library(FRASER)
})

# ------------------------------------------------------------------------------
# 1. Argumenten / configuratie
#    Voor testdata: laat --use_test_data op TRUE staan.
#    Voor echte data: zet --use_test_data FALSE en geef --dataset_dir en
#    --annotation_name mee (dit is waar DROP je fds-object heeft opgeslagen,
#    typisch onder <drop_output>/processed_data/aberrant_splicing/datasets).

# For test data: keep --use_test_data TRUE.
# For real data: set --use_test_data FALSE and provide --dataset_dir and 
#--annotation_name (this is where DROP has saved your fds object, typically under <drop_output>/processed_data/aberrant_splicing/datasets).
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
use_test_data   <- TRUE     # TRUE = testdata, FALSE = real data
dataset_dir     <- NA       # only needed for real data (where DROP saved the FraserDataSet)
annotation_name <- "raw-local"
type            <- "jaccard"
q_dim           <- 3
outdir          <- "fraser_plots_out"



# ------------------------------------------------------------------------------
# 2. FraserDataSet load (test or real data)
# ------------------------------------------------------------------------------
if (isTRUE(use_test_data)) {
  message("Load FRASER-testdata...")
  fds <- createTestFraserDataSet()
} else {
  if (is.na(dataset_dir)) {
    stop("Give --dataset_dir for real data (where DROP saved the FraserDataSet).")
  }
  message("Load saved FraserDataSet from: ", dataset_dir)
  fds <- loadFraserDataSet(dir = dataset_dir, name = annotation_name)
}

# ------------------------------------------------------------------------------
# 3. PSI/Jaccard calculation (if not yet calculated)
# ------------------------------------------------------------------------------
message("Calculate PSI/Jaccard values...")
fds <- calculatePSIValues(fds)

# ------------------------------------------------------------------------------
# 4. plotFilterExpression -- QC vóór filtering
# ------------------------------------------------------------------------------
message("Make plotFilterExpression...")
fds_unfiltered <- filterExpressionAndVariability(fds, minDeltaPsi = 0, filter = FALSE)

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
pdf(file.path(outdir, "01_filter_expression.pdf"), width = 6, height = 5)
plotFilterExpression(fds_unfiltered, bins = 100)
dev.off()

# Daadwerkelijk filteren voor de vervolgstappen
fds_filtered <- fds_unfiltered[mcols(fds_unfiltered, type = "j")[, "passed"]]
message(sprintf("Introns voor filtering: %d, na filtering: %d",
                 nrow(fds_unfiltered), nrow(fds_filtered)))

# ------------------------------------------------------------------------------
# 5. Correlatie-heatmap vóór correctie (ruw)
# ------------------------------------------------------------------------------
message("Make raw correlation heatmap...")
pdf(file.path(outdir, "02_correlation_raw.pdf"), width = 6, height = 6)
plotCountCorHeatmap(fds_filtered, type = type, logit = TRUE, normalized = FALSE)
dev.off()

# ------------------------------------------------------------------------------
# 6. Model fitten (nodig voor p-waarden / plotQQ)
#    Let op: dit is de rekenintensieve stap. Bij een grote, echte dataset
#    hoort dit als losse SLURM-job te draaien i.p.v. interactief.
# ------------------------------------------------------------------------------
message(sprintf("Fit the model (q=%d, type=%s)...", q_dim, type))
qParam <- setNames(q_dim, type)
fds_filtered <- FRASER(fds_filtered, q = qParam)

# ------------------------------------------------------------------------------
# 7. Correlatie-heatmap ná correctie (genormaliseerd) -- vergelijk met stap 5
# ------------------------------------------------------------------------------
message("Make normalized correlation heatmap...")
pdf(file.path(outdir, "03_correlation_normalized.pdf"), width = 6, height = 6)
plotCountCorHeatmap(fds_filtered, type = type, logit = TRUE, normalized = TRUE)
dev.off()

# ------------------------------------------------------------------------------
# 8. plotQQ -- kalibratiecheck (globaal, gene-level)
# ------------------------------------------------------------------------------
message("Make plotQQ...")
pdf(file.path(outdir, "04_qq_global.pdf"), width = 5, height = 5)
plotQQ(fds_filtered, aggregate = TRUE, global = TRUE)
dev.off()

# ------------------------------------------------------------------------------
# 9. plotAberrantPerSample -- overzicht outliers per patiënt
# ------------------------------------------------------------------------------
message("Make plotAberrantPerSample...")
pdf(file.path(outdir, "05_aberrant_per_sample.pdf"), width = 6, height = 5)
plotAberrantPerSample(fds_filtered)
dev.off()

# ------------------------------------------------------------------------------
# 10. Resultaattabel wegschrijven (handig voor de vergelijking met predicties)
# ------------------------------------------------------------------------------
message("Write results table...")
res <- results(fds_filtered, all = TRUE)
write.table(as.data.frame(res),
            file = file.path(outdir, "results_all.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

message("Done, all plots and the results table are in: ", normalizePath(outdir))