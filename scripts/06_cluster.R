suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(ggsci)
})

args         <- commandArgs(trailingOnly = TRUE)
rds_file     <- args[1]
n_dims       <- as.integer(args[2])   # e.g. 17 — number of PCs, use elbow plot output to decide
resolution   <- as.numeric(args[3])   # e.g. 0.18 — clustering resolution: higher = more clusters
n_neighbors  <- as.integer(args[4])   # e.g. 150 — UMAP n.neighbors: higher = more global structure
min_dist     <- as.numeric(args[5])   # e.g. 0.4 — UMAP min.dist: lower = tighter clusters
group_by_var <- args[6]               # e.g. "experiment" — metadata column to colour second UMAP by
palette_arg  <- args[7]               # "npg","jco","lancet","nejm" or comma-separated hex codes
                                      # e.g. "#3C5488FF,#E64B35FF,#4DBBD599,#00A087FF"

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

# IMPORTANT: use integrated assay for scaling, PCA, and clustering
# use RNA assay for feature plots and DE
DefaultAssay(sobj) <- "integrated"
sobj <- ScaleData(sobj, verbose = FALSE)
sobj <- RunPCA(sobj, verbose = FALSE)

# use elbow plot to decide number of PCs (n_dims) — look for where curve flattens
elbow <- ElbowPlot(sobj) +
  theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave("elbow_plot.pdf", plot = elbow, width = 6, height = 4)

sobj <- RunUMAP(sobj, reduction = "pca", dims = 1:n_dims,
                n.neighbors = n_neighbors, min.dist = min_dist)
sobj <- FindNeighbors(sobj, reduction = "pca", dims = 1:n_dims)

# first pass: cluster at multiple resolutions to decide visually
# look at the PDFs, pick a resolution, then rerun with --resolution set (yaml)
resolutions <- c(0.2, 0.4, 0.6, 0.8)
pdf("umap_resolutions.pdf", width = 12, height = 10)
for (res in resolutions) {
  sobj <- FindClusters(sobj, resolution = res)
  pal  <- parse_palette(palette_arg, length(unique(Idents(sobj))))
  p    <- DimPlot(sobj, reduction = "umap", label = TRUE, repel = TRUE,
                  label.box = FALSE, label.size = 4, pt.size = 0.05,
                  cols = pal, raster = FALSE) +
          NoAxes() +
          ggtitle(paste0("resolution = ", res))
  print(p)
}
dev.off()

# final clustering at chosen resolution
sobj <- FindClusters(sobj, resolution = resolution)

pal <- parse_palette(palette_arg, length(unique(Idents(sobj))))

p1 <- DimPlot(sobj, reduction = "umap", label = TRUE, repel = TRUE,
              label.box = FALSE, label.size = 4, pt.size = 0.05,
              cols = pal, raster = FALSE) + NoAxes()

# second UMAP coloured by sample/condition to check integration quality
p2 <- DimPlot(sobj, reduction = "umap", group.by = group_by_var,
              pt.size = 0.05, raster = FALSE) + NoAxes()

ggsave("umap_clusters.pdf", plot = p1 | p2, width = 12, height = 5)

# standalone full-size UMAP coloured by group/condition
p_group <- DimPlot(sobj, reduction = "umap", group.by = group_by_var,
                   pt.size = 0.05, raster = FALSE) + NoAxes()
ggsave("umap_by_group.pdf", plot = p_group, width = 6, height = 5)

cat("Number of clusters:", length(unique(Idents(sobj))), "\n")

saveRDS(sobj, file = "seurat_clustered.rds")
