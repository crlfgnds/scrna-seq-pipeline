suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(SeuratWrappers)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

args         <- commandArgs(trailingOnly = TRUE)
rds_file     <- args[1]
root_cluster <- args[2]              # e.g. "5" — cluster to use as trajectory root
group_col    <- args[3]              # e.g. "Group" — column to split CT vs KO
group_levels <- strsplit(args[4], ",")[[1]]  # e.g. "Ctl,KO"
palette_arg  <- args[5]              # "npg","jco","lancet","nejm" or comma-separated hex codes

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

# IMPORTANT: Monocle 3 requires UMAP — does not accept TSNE
# convert Seurat object to CellDataSet (cds)
cds <- SeuratWrappers::as.cell_data_set(sobj)
fData(cds)$gene_short_name <- rownames(fData(cds))

# assign partitions — treating all cells as one partition
reacreate.partition <- factor(rep(1, length(cds@colData@rownames)),
                               levels = 1)
names(reacreate.partition)        <- cds@colData@rownames
cds@clusters$UMAP$partitions      <- reacreate.partition

# transfer UMAP coordinates from Seurat to Monocle
cds@int_colData@listData$reducedDims$UMAP <- sobj@reductions$umap@cell.embeddings

# learn trajectory graph
cds <- learn_graph(cds, use_partition = FALSE)

# order cells using the specified root cluster
cds <- order_cells(cds, reduction_method = "UMAP",
                   root_cells = colnames(cds[, clusters(cds) == root_cluster]))

# store pseudotime in metadata
cds$monocle3_pseudotime <- pseudotime(cds)

# plot all cells coloured by pseudotime
p_all <- plot_cells(cds,
                    color_cells_by            = "pseudotime",
                    group_cells_by            = "partition",
                    label_groups_by_cluster   = FALSE,
                    label_branch_points       = FALSE,
                    label_roots               = FALSE,
                    label_leaves              = FALSE,
                    show_trajectory_graph     = TRUE,
                    trajectory_graph_color    = "black",
                    trajectory_graph_segment_size = 0.5,
                    cell_size                 = 0.5,
                    alpha                     = 0.3,
                    scale_to_range            = TRUE)
ggsave("pseudotime_all.pdf", plot = p_all, width = 6, height = 5)

# boxplot: pseudotime distribution per cluster
data.pseudo <- as.data.frame(colData(cds))
n_clusters  <- length(unique(data.pseudo$seurat_clusters))
pal         <- parse_palette(palette_arg, n_clusters)

p_box <- ggplot(data.pseudo,
                aes(monocle3_pseudotime,
                    reorder(seurat_clusters, monocle3_pseudotime, median),
                    fill = seurat_clusters)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = pal) +
  theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position  = "none") +
  xlab("Pseudotime") +
  ylab("Cluster")
ggsave("pseudotime_boxplot.pdf", plot = p_box, width = 6, height = 5)

# plot split by group (e.g. Ctl vs KO)
cds_list <- list()
plots    <- list()

for (grp in group_levels) {
  sub <- subset(sobj, cells = colnames(sobj)[sobj@meta.data[[group_col]] == grp])
  cds_sub <- SeuratWrappers::as.cell_data_set(sub)
  fData(cds_sub)$gene_short_name <- rownames(fData(cds_sub))

  part <- factor(rep(1, length(cds_sub@colData@rownames)), levels = 1)
  names(part) <- cds_sub@colData@rownames
  cds_sub@clusters$UMAP$partitions <- part
  cds_sub@int_colData@listData$reducedDims$UMAP <- sub@reductions$umap@cell.embeddings

  cds_sub <- learn_graph(cds_sub, use_partition = FALSE)
  cds_sub <- order_cells(cds_sub, reduction_method = "UMAP",
                         root_cells = colnames(cds_sub[, clusters(cds_sub) == root_cluster]))

  cds_list[[grp]] <- cds_sub
  plots[[grp]] <- plot_cells(cds_sub,
                              color_cells_by              = "pseudotime",
                              group_cells_by              = "partition",
                              label_groups_by_cluster     = FALSE,
                              label_branch_points         = FALSE,
                              label_roots                 = FALSE,
                              label_leaves                = FALSE,
                              show_trajectory_graph       = TRUE,
                              trajectory_graph_color      = "black",
                              trajectory_graph_segment_size = 0.5,
                              cell_size                   = 0.5,
                              alpha                       = 0.3,
                              scale_to_range              = TRUE) +
    labs(title = grp)
}

p_split <- Reduce("|", plots)
ggsave("pseudotime_split_groups.pdf", plot = p_split,
       width = 5 * length(group_levels), height = 5)

saveRDS(cds,      file = "cds_all.rds")
saveRDS(cds_list, file = "cds_list_groups.rds")
