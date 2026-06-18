suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(writexl)
  library(ggsci)
  library(patchwork)
})

args              <- commandArgs(trailingOnly = TRUE)
rds_file          <- args[1]
marker_file       <- args[2]   # optional csv with columns: gene, celltype (pass "NA" to skip)
annotation_file   <- args[3]   # optional csv with columns: cluster, cell_type (pass "NA" to skip)
                                # if provided: renames clusters and stores as cell_type in metadata
                                # if not provided: cell_type is set to cluster numbers (annotate later)
min_pct           <- as.numeric(args[4])   # e.g. 0.25
logfc_thresh      <- as.numeric(args[5])   # e.g. 0.25
top_n_genes       <- as.integer(args[6])   # e.g. 20
palette_arg       <- args[7]               # "npg","jco","lancet","nejm" or comma-separated hex codes

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

# record the resolution used — stored in seurat_clusters column, extract from metadata
resolution_used <- grep("RNA_snn_res|integrated_snn_res", colnames(sobj@meta.data), value = TRUE)
if (length(resolution_used) > 0) {
  sobj@meta.data$resolution_used <- sobj@meta.data[[tail(resolution_used, 1)]]
}

# if annotation file provided: rename cluster identities to cell type names
# annotation_file must be a csv with two columns: cluster (e.g. "0","1"...) and cell_type (e.g. "ILC3s")
if (!is.na(annotation_file) && file.exists(annotation_file)) {
  ann_map <- read.csv(annotation_file)
  # build named vector for RenameIdents: c("0" = "ILC3s", "1" = "NK cells", ...)
  rename_vec <- setNames(ann_map$cell_type, as.character(ann_map$cluster))
  sobj <- RenameIdents(sobj, rename_vec)
  cat("Clusters renamed using annotation file\n")
} else {
  cat("No annotation file provided — cell_type will reflect cluster numbers\n")
}

# IMPORTANT: always store active ident as cell_type in metadata
# downstream scripts (CellChat, composition, DE) rely on this column
sobj@meta.data$cell_type <- as.character(Idents(sobj))

# IMPORTANT: use RNA assay for marker finding and DE — not integrated
DefaultAssay(sobj) <- "RNA"

# find markers for every cluster vs all remaining cells, positive markers only
# min.pct: minimum % of cells in a cluster expressing the gene — filters noise
# logfc.threshold: minimum log fold change — filters weak signals
all_markers <- FindAllMarkers(sobj, only.pos = TRUE,
                               min.pct         = min_pct,
                               logfc.threshold = logfc_thresh)

top_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(n = top_n_genes, order_by = avg_log2FC)

# export full marker table
write_xlsx(all_markers, path = "all_markers.xlsx")

# summary table: top 15 genes per cluster as comma-separated string — useful for manual annotation
top_genes_summary <- all_markers %>%
  group_by(cluster) %>%
  slice_max(n = 15, order_by = avg_log2FC) %>%
  summarise(top_genes = paste(gene, collapse = ","))
write_xlsx(top_genes_summary, path = "top_genes_per_cluster.xlsx")

# IMPORTANT: use scaled RNA assay for heatmap
sobj_scaled <- ScaleData(sobj, verbose = FALSE)
pdf("heatmap_top_markers.pdf", width = 14, height = 10)
DoHeatmap(subset(sobj_scaled, downsample = 100),
          features = top_markers$gene,
          size     = 4,
          disp.max = 2,
          disp.min = -2) +
  scale_fill_gradient2(
    low      = rev(c("#d1e5f0", "#67a9cf", "#2166ac")),
    mid      = "white",
    high     = rev(c("#b2182b", "#ef8a62", "#fddbc7")),
    midpoint = 0,
    guide    = "colourbar"
  )
dev.off()

pal <- parse_palette(palette_arg, length(unique(Idents(sobj))))
p_umap <- DimPlot(sobj, reduction = "umap", label = TRUE, repel = TRUE,
                  label.box = FALSE, label.size = 4, pt.size = 0.05,
                  cols = pal, raster = FALSE) + NoAxes()
ggsave("umap_annotated.pdf", plot = p_umap, width = 8, height = 7)

# dot plot with known canonical markers if provided
# shape = 21 with stroke gives the outlined dot style
if (!is.na(marker_file) && file.exists(marker_file)) {
  known <- read.csv(marker_file)
  p_dot <- DotPlot(sobj, features = unique(known$gene)) +
    scale_color_gradient2(low      = "#9ECAE1",
                          mid      = "white",
                          high     = "red3",
                          midpoint = 0,
                          name     = "Mean expression\nin group") +
    geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5) +
    scale_size(range = c(0.5, 6.5), name = "Fraction of cells\nin group (%)") +
    theme_bw(base_size = 12) +
    theme(panel.grid.major  = element_blank(),
          panel.grid.minor  = element_blank(),
          axis.text.x       = element_text(angle = 45, hjust = 1,
                                            size = 12, face = "italic"),
          axis.text.y       = element_text(size = 12),
          axis.title        = element_blank(),
          legend.position   = "right",
          legend.text       = element_text(size = 10),
          legend.title      = element_text(size = 10)) +
    guides(size = guide_legend(
      override.aes = list(shape = 21, colour = "black", fill = "white")))
  ggsave("dotplot_markers.pdf", plot = p_dot, width = 14, height = 5)

  genes      <- unique(known$gene)
  plot_rows  <- ceiling(length(genes) / 4)

  # feature plots of canonical markers on the UMAP — for visual cluster identification
  p_feat <- FeaturePlot(sobj, features = genes, order = TRUE,
                        cols = c("gray93", "red2"), pt.size = 0.2, ncol = 4) & NoAxes()
  ggsave("featureplot_markers.pdf", plot = p_feat, width = 16, height = 4 * plot_rows)

  # violin plots of canonical markers per cluster
  p_vln <- VlnPlot(sobj, features = genes, pt.size = 0, ncol = 4, cols = pal) & NoLegend()
  ggsave("violin_markers.pdf", plot = p_vln, width = 16, height = 4 * plot_rows)
}

saveRDS(sobj, file = "seurat_annotated.rds")
