suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggsci)
  library(patchwork)
})

args           <- commandArgs(trailingOnly = TRUE)
rds_file       <- args[1]
clusters       <- strsplit(args[2], ",")[[1]]   # e.g. "ILC3s" or "0,1,2,3" — clusters to subset
subset_name    <- args[3]                        # e.g. "ILC3s" — used for output file names
experiment_col <- args[4]                        # e.g. "experiment" — metadata column to filter by
experiments    <- strsplit(args[5], ",")[[1]]   # e.g. "CD45_neg_Ctl,CD45_neg_KO,..." — values to keep
subset_n_dims         <- as.integer(args[6])     # e.g. 8 — PCs for re-clustering (separate from whole-object n_dims)
subset_resolution     <- as.numeric(args[7])     # e.g. 0.3 (separate from whole-object resolution)
subset_palette_arg    <- args[8]                 # "npg","jco","lancet","nejm" or comma-separated hex codes

parse_palette <- function(palette_arg, n) {
  if (grepl("^#", palette_arg)) strsplit(palette_arg, ",")[[1]]
  else switch(palette_arg,
    "npg"    = pal_npg("nrc")(n),
    "jco"    = pal_jco()(n),
    "lancet" = pal_lancet()(n),
    "nejm"   = pal_nejm()(n),
    pal_npg("nrc")(n)
  )
}

sobj <- readRDS(rds_file)

# subset by cluster identity
sub <- subset(sobj, idents = clusters)

# IMPORTANT: also subset by experiment to ensure only relevant cells are kept
# this removes any cells from experiments not relevant to this subset analysis
sub <- subset(sub, cells = colnames(sub)[sub@meta.data[[experiment_col]] %in% experiments])

# re-cluster the subset at higher resolution
sub <- RunPCA(sub, verbose = FALSE)

elbow <- ElbowPlot(sub) +
  theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave(paste0(subset_name, "_elbow.pdf"), plot = elbow, width = 6, height = 4)

sub <- RunUMAP(sub, reduction = "pca", dims = 1:subset_n_dims)
sub <- FindNeighbors(sub, dims = 1:subset_n_dims)

# multi-resolution pass to help decide
resolutions <- c(0.2, 0.4, 0.6, 0.8)
pdf(paste0(subset_name, "_umap_resolutions.pdf"), width = 12, height = 10)
for (res in resolutions) {
  sub <- FindClusters(sub, resolution = res)
  pal <- parse_palette(subset_palette_arg, length(unique(Idents(sub))))
  p   <- DimPlot(sub, reduction = "umap", label = TRUE, repel = TRUE,
                 label.box = FALSE, label.size = 4, pt.size = 0.05,
                 cols = pal, raster = FALSE) +
         NoAxes() + ggtitle(paste0("resolution = ", res))
  print(p)
}
dev.off()

# final clustering at chosen resolution
sub <- FindClusters(sub, resolution = subset_resolution)

pal <- pal_npg("nrc")(length(unique(Idents(sub))))
p_final <- DimPlot(sub, reduction = "umap", label = TRUE, repel = TRUE,
                   label.box = FALSE, label.size = 4, pt.size = 0.05,
                   cols = pal, raster = FALSE) + NoAxes()
ggsave(paste0(subset_name, "_umap_final.pdf"), plot = p_final, width = 7, height = 6)

saveRDS(sub, file = paste0(subset_name, "_subsetted.rds"))
