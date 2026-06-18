suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

args          <- commandArgs(trailingOnly = TRUE)
rds_file      <- args[1]
mito_thresh   <- as.numeric(args[2])   # e.g. 0.05 
min_features  <- as.numeric(args[3])   # e.g. 200 — removes empty droplets
max_features  <- as.numeric(args[4])   # e.g. 6000 — removes likely doublets
min_counts    <- as.numeric(args[5])   # e.g. 500
max_counts    <- as.numeric(args[6])   # e.g. 30000

sobj <- readRDS(rds_file)

# QC plots before filtering: use to decide cutoffs
p_before <- VlnPlot(sobj,
                    features = c("nFeature_RNA", "nCount_RNA", "percent.mito"),
                    ncol = 3, pt.size = 0.1, group.by = "orig.ident") &
  theme_bw(base_size = 12) &
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x      = element_text(angle = 45, hjust = 1),
        legend.position  = "none")

p_scatter <- FeatureScatter(sobj, feature1 = "nCount_RNA", feature2 = "percent.mito") |
             FeatureScatter(sobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

ggsave("qc_before_filtering.pdf",
       plot  = p_before / p_scatter,
       width = 12, height = 8)

# filter dead cells (high mito) and low-quality cells (too few or too many features/counts)
# high nCount + high nFeature = likely doublets; low = likely empty droplets or dead cells
sobj <- subset(sobj,
               subset = percent.mito < mito_thresh  &
                        nFeature_RNA  > min_features  &
                        nFeature_RNA  < max_features  &
                        nCount_RNA    > min_counts    &
                        nCount_RNA    < max_counts)

p_after <- VlnPlot(sobj,
                   features = c("nFeature_RNA", "nCount_RNA", "percent.mito"),
                   ncol = 3, pt.size = 0, group.by = "orig.ident") &
  theme_bw(base_size = 12) &
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x      = element_text(angle = 45, hjust = 1),
        legend.position  = "none")

ggsave("qc_after_filtering.pdf", plot = p_after, width = 10, height = 4)

cat("Cells after filtering:", ncol(sobj), "\n")

saveRDS(sobj, file = "seurat_filtered.rds")
