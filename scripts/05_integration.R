suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
})

args        <- commandArgs(trailingOnly = TRUE)
rds_file    <- args[1]   # seurat_list_normalized.rds from script 04
integration_features  <- as.integer(args[2])   # e.g. 2000 — integration features
method      <- args[3]               # e.g. "cca", "rpca", "harmony"

obj_list <- readRDS(rds_file)

if (method == "harmony") {
  # harmony works on the merged object after PCA — does not use anchors
  # group.by.vars: metadata column that defines batches (e.g. "experiment")
  merged <- merge(obj_list[[1]], y = obj_list[-1])
  merged <- NormalizeData(merged)
  merged <- FindVariableFeatures(merged, nfeatures = integration_features)
  merged <- ScaleData(merged, verbose = FALSE)
  merged <- RunPCA(merged, verbose = FALSE)
  integrated <- RunHarmony(merged, group.by.vars = "orig.ident", verbose = FALSE)

} else {
  # CCA and RPCA both use the anchor-based integration framework
  # CCA: finds shared correlation structure — better for datasets with different cell type compositions
  # RPCA: reciprocal PCA — faster, more conservative, better when datasets are similar
  reduction <- ifelse(method == "rpca", "rpca", "cca")

  features <- SelectIntegrationFeatures(object.list = obj_list,
                                        nfeatures = integration_features)

  if (method == "rpca") {
    # RPCA requires PCA to be run on each sample first
    obj_list <- lapply(obj_list, function(x) {
      x <- ScaleData(x, features = features, verbose = FALSE)
      x <- RunPCA(x, features = features, verbose = FALSE)
      x
    })
  }

  anchors <- FindIntegrationAnchors(object.list = obj_list,
                                    anchor.features = features,
                                    reduction = reduction)

  # creates a new "integrated" assay — use for dim reduction and clustering only
  # use RNA assay for DE analysis
  integrated <- IntegrateData(anchorset = anchors)
}

cat("Integration method:", method, "\n")
cat("Samples integrated:", length(unique(integrated$orig.ident)), "\n")
cat("Integrated object dimensions:", dim(integrated), "\n")

saveRDS(integrated, file = "seurat_integrated.rds")
