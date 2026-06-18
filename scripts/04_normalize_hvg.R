suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

args       <- commandArgs(trailingOnly = TRUE)
rds_file   <- args[1]
n_features <- as.integer(args[2])   # e.g. 2000 — number of highly variable genes to select per sample
split_by   <- args[3]               # e.g. "experiment" — metadata column to split by for per-sample normalisation

sobj <- readRDS(rds_file)

# split by sample before normalising — required for integration
# each sample is normalised independently to avoid batch effects leaking in before correction
obj_list <- SplitObject(sobj, split.by = split_by)

obj_list <- lapply(obj_list, function(x) {
  x <- NormalizeData(x)
  # vst: variance-stabilising transformation — selects genes with high cell-to-cell variability
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = n_features)
  x
})

# plot HVGs for the first sample as a QC check
top_hvg <- head(VariableFeatures(obj_list[[1]]), 20)
p <- VariableFeaturePlot(obj_list[[1]])
p <- LabelPoints(plot = p, points = top_hvg, repel = TRUE) +
  theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave("hvg_plot.pdf", plot = p, width = 8, height = 5)

# save as a list — integration script (05) reads it directly
saveRDS(obj_list, file = "seurat_list_normalized.rds")
