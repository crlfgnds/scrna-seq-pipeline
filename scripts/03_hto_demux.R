suppressPackageStartupMessages({
  library(Seurat)
  library(hdf5r)
  library(ggplot2)
})

args             <- commandArgs(trailingOnly = TRUE)
rds_file         <- args[1]
hto_h5           <- args[2]   # same filtered_feature_bc_matrix.h5 — Antibody Capture slot
hashtags         <- strsplit(args[3], ",")[[1]]   # comma-separated e.g. "mmHashtag1,mmHashtag2,mmHashtag3,mmHashtag4,mmHashtag5,mmHashtag6"
hashtag_labels   <- strsplit(args[4], ",")[[1]]   # human-readable labels in same order e.g. "Control_1,Control_2,Control_3,aLy6G_1,aLy6G_2,aLy6G_3"
group_pattern    <- args[5]                        # grep pattern to identify ident1 in experiment column e.g. "Ctl"
ident1           <- args[6]                        # group label for matching cells e.g. "Ctl"
ident2           <- args[7]                        # group label for non-matching cells e.g. "KO"
experiment_col   <- args[8]                        # metadata column to grep for group assignment e.g. "experiment"

sobj <- readRDS(rds_file)

hto_raw    <- Read10X_h5(hto_h5)
hto_subset <- hto_raw[["Antibody Capture"]][hashtags, ]

# IMPORTANT: RNA and HTO dims must match — subset to shared barcodes only
joint.bcs  <- intersect(colnames(sobj), colnames(hto_subset))
hto_subset <- hto_subset[, joint.bcs]
sobj       <- subset(sobj, cells = joint.bcs)

sobj[["HTO"]]       <- CreateAssayObject(counts = hto_subset)
sobj@assays$HTO@key <- "hto_"   # very important line — required for HTODemux to find the assay

# CLR normalisation is standard for antibody/hashtag data (not log-normalise)
sobj <- NormalizeData(sobj, assay = "HTO", normalization.method = "CLR")

# no hashtag should sum to zero after normalisation
if (any(rowSums(sobj@assays$HTO) == 0))
  stop("One or more hashtags have zero counts after normalisation — check hashtag names.")

sobj <- HTODemux(sobj, assay = "HTO", positive.quantile = 0.95)

cat("HTO classification summary:\n")
print(table(sobj$HTO_classification.global))
cat("\nSample assignments:\n")
print(table(sobj$hash.ID))

# map hashtags to human-readable sample labels (e.g. "mmHashtag1" → "Control_1")
sobj$SampleType <- factor(sobj$hash.ID,
                           levels = hashtags,
                           labels = hashtag_labels)

# create Group column from experiment metadata column using grep pattern
# e.g. experiment = "CD45_pos_Ctl" → Group = "Ctl"
sobj$Group <- ifelse(grepl(group_pattern, sobj@meta.data[[experiment_col]]), ident1, ident2)

cat("\nSampleType assignments:\n")
print(table(sobj$SampleType))
cat("\nGroup assignments:\n")
print(table(sobj$Group))

p_ridge <- RidgePlot(sobj, assay = "HTO", features = hashtags, ncol = 2) &
  theme_bw(base_size = 12) &
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave("hto_ridge.pdf", plot = p_ridge, width = 10, height = 6)

saveRDS(sobj, file = "seurat_demuxed.rds")
