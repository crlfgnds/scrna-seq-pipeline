suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(writexl)
  library(ComplexHeatmap)
})

args          <- commandArgs(trailingOnly = TRUE)
rds_file      <- args[1]
group_col     <- args[2]              # e.g. "Group" — condition column
group_levels  <- strsplit(args[3], ",")[[1]]   # e.g. "Control,TNF_deltaARE" — must be exactly 2
celltype_col  <- args[4]              # e.g. "cell_type" — must exist in metadata
organism      <- args[5]              # "mouse" or "human"
source_cells  <- strsplit(args[6], ",")[[1]]   # e.g. "Neutrophils,ILC3s" — for bubble/rankNet per sender

sobj <- readRDS(rds_file)

# IMPORTANT: use RNA assay for CellChat
DefaultAssay(sobj) <- "RNA"

# select ligand-receptor database based on organism
CellChatDB <- if (organism == "mouse") CellChatDB.mouse else CellChatDB.human
ppi        <- if (organism == "mouse") PPI.mouse        else PPI.human

# helper: create and run CellChat for one condition
run_cellchat <- function(seurat_sub, db, ppi, celltype_col) {
  meta <- seurat_sub@meta.data
  cc   <- createCellChat(object = seurat_sub, meta = meta, group.by = celltype_col)
  cc   <- addMeta(cc, meta = meta)
  cc   <- setIdent(cc, ident.use = celltype_col)
  cc@DB <- db

  cc <- subsetData(cc)
  cc <- updateCellChat(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  # project data onto PPI network to improve sensitivity
  cc <- projectData(cc, ppi)
  # truncatedMean trims extreme outliers — faster and more robust than default
  cc <- computeCommunProb(cc, type = "truncatedMean", trim = 0.05)
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc <- netAnalysis_computeCentrality(cc, slot.name = "netP")
  cc
}

# run CellChat per condition
Idents(sobj) <- group_col
cc_list <- list()
for (grp in group_levels) {
  sub          <- subset(sobj, idents = grp)
  Idents(sub)  <- celltype_col
  cc_list[[grp]] <- run_cellchat(sub, CellChatDB, ppi, celltype_col)
  saveRDS(cc_list[[grp]], file = paste0("cellchat_", grp, ".rds"))
}

# merge both conditions for comparison
cellchat_merged <- mergeCellChat(cc_list, add.names = names(cc_list))
saveRDS(cellchat_merged, file = "cellchat_merged.rds")

object.list <- cc_list
groupSize   <- as.numeric(table(cc_list[[group_levels[1]]]@idents))

# --- comparison plots ---

# total interactions and interaction strength
gg1 <- compareInteractions(cellchat_merged, show.legend = FALSE, group = c(1, 2))
gg2 <- compareInteractions(cellchat_merged, show.legend = FALSE, group = c(1, 2), measure = "weight")
ggsave("cellchat_compare_interactions.pdf", plot = gg1 + gg2, width = 8, height = 4)

# differential interaction circle plots (count and weight)
pdf("cellchat_diffInteraction_circles.pdf", width = 10, height = 5)
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE)
netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE, measure = "weight")
dev.off()

# differential interaction heatmap (count and weight)
pdf("cellchat_heatmap.pdf", width = 10, height = 8)
gg_h1 <- netVisual_heatmap(cellchat_merged)
gg_h2 <- netVisual_heatmap(cellchat_merged, measure = "weight")
print(gg_h1 + gg_h2)
dev.off()

# circle plots — split by condition, normalised to same scale
weight.max <- getMaxWeight(object.list, attribute = c("idents", "count"))
pdf("cellchat_circle_split.pdf", width = 10, height = 5)
par(mfrow = c(1, 2), xpd = TRUE)
for (i in seq_along(object.list)) {
  netVisual_circle(object.list[[i]]@net$count,
                   weight.scale    = TRUE,
                   label.edge      = FALSE,
                   edge.weight.max = weight.max[2],
                   edge.width.max  = 12,
                   title.name      = paste0("Number of interactions — ", names(object.list)[i]))
}
dev.off()

# per-cell-type circle plots for first condition
mat <- object.list[[group_levels[1]]]@net$weight
pdf("cellchat_circle_per_celltype.pdf", width = 12, height = 9)
par(mfrow = c(3, 4), xpd = TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = TRUE,
                   edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

# scatter: major sources and targets in 2D space — one plot per condition
num.link     <- sapply(object.list, function(x) rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count))
weight.MinMax <- c(min(num.link), max(num.link))
gg_scatter <- list()
for (i in seq_along(object.list)) {
  gg_scatter[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]],
                                                        title = names(object.list)[i],
                                                        weight.MinMax = weight.MinMax)
}
ggsave("cellchat_signalingRole_scatter.pdf",
       plot = patchwork::wrap_plots(gg_scatter, ncol = 2),
       width = 11, height = 5)

# signaling role heatmap — incoming pattern, all pathways union, both conditions
# uses ComplexHeatmap — must be drawn to pdf explicitly
i <- 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i + 1]]@netP$pathways)
ht1 <- netAnalysis_signalingRole_heatmap(object.list[[i]],
                                          pattern  = "incoming",
                                          signaling = pathway.union,
                                          title    = names(object.list)[i],
                                          width    = 5, height = 35)
ht2 <- netAnalysis_signalingRole_heatmap(object.list[[i + 1]],
                                          pattern  = "incoming",
                                          signaling = pathway.union,
                                          title    = names(object.list)[i + 1],
                                          width    = 5, height = 35)
pdf("cellchat_signalingRole_heatmap_incoming.pdf", width = 12, height = 35)
ComplexHeatmap::draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

ht1_out <- netAnalysis_signalingRole_heatmap(object.list[[i]],
                                              pattern  = "outgoing",
                                              signaling = pathway.union,
                                              title    = names(object.list)[i],
                                              width    = 5, height = 35)
ht2_out <- netAnalysis_signalingRole_heatmap(object.list[[i + 1]],
                                              pattern  = "outgoing",
                                              signaling = pathway.union,
                                              title    = names(object.list)[i + 1],
                                              width    = 5, height = 35)
pdf("cellchat_signalingRole_heatmap_outgoing.pdf", width = 12, height = 35)
ComplexHeatmap::draw(ht1_out + ht2_out, ht_gap = unit(0.5, "cm"))
dev.off()

# rankNet — overall information flow, stacked and unstacked
gg_rank1 <- rankNet(cellchat_merged, mode = "comparison", measure = "weight",
                    stacked = TRUE,  do.stat = TRUE)
gg_rank2 <- rankNet(cellchat_merged, mode = "comparison", measure = "count",
                    stacked = FALSE, do.stat = TRUE)
ggsave("cellchat_rankNet_overall.pdf", plot = gg_rank1 + gg_rank2, width = 12, height = 8)

# rankNet + bubble per source cell type
for (src in source_cells) {
  src_safe <- gsub(" ", "_", src)

  # rankNet for this sender
  gg_rn <- tryCatch(
    rankNet(cellchat_merged, mode = "comparison", stacked = FALSE,
            do.stat = TRUE, thresh = 0.01, sources.use = src,
            title = paste("Top signaling pathways from", src)),
    error = function(e) NULL
  )
  if (!is.null(gg_rn))
    ggsave(paste0("cellchat_rankNet_", src_safe, ".pdf"), plot = gg_rn, width = 6, height = 8)

  # bubble plot for this sender — all targets
  p_bubble <- tryCatch(
    netVisual_bubble(cellchat_merged,
                     sources.use = src,
                     targets.use = NULL,
                     comparison  = c(1, 2),
                     angle.x     = 45),
    error = function(e) NULL
  )
  if (!is.null(p_bubble))
    ggsave(paste0("cellchat_bubble_", src_safe, ".pdf"), plot = p_bubble, width = 10, height = 12)

  # extract communication dataframe and export to xlsx
  comm_data <- tryCatch(
    netVisual_bubble(cellchat_merged,
                     sources.use  = src,
                     targets.use  = NULL,
                     comparison   = c(1, 2),
                     return.data  = TRUE)$communication,
    error = function(e) NULL
  )
  if (!is.null(comm_data))
    write_xlsx(comm_data, path = paste0("cellchat_communication_", src_safe, ".xlsx"))
}

cat("CellChat analysis complete\n")
cat("Conditions compared:", paste(group_levels, collapse = " vs "), "\n")
