suppressPackageStartupMessages({
  library(Seurat)
  library(hdf5r)
})

#Script is ran in the terminal: Rscript 01_create_seurat.R myfile.h5 MySample
#commandArgs(trailingOnly = TRUE) captures everything I type after Rscript 01_create_seurat.R as a vector

args        <- commandArgs(trailingOnly = TRUE)
h5_file     <- args[1] #position 1 after Rscript...
sample_name <- args[2]

#count matrix 
expr <- Read10X_h5(h5_file)
if (is.list(expr)) expr <- expr[["Gene Expression"]]

#create Seurat object
sobj <- CreateSeuratObject(counts = expr, project = sample_name,
                           min.cells = 3, min.features = 200)

#calculate mito genes
mt_pattern <- ifelse(grepl("^[A-Z]", rownames(sobj)[1]), "^MT-", "^mt-")
sobj[["percent.mito"]] <- PercentageFeatureSet(sobj, pattern = mt_pattern) / 100

saveRDS(sobj, file = "seurat_raw.rds")
