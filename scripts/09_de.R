suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(writexl)
  library(EnhancedVolcano)
  library(ggplot2)
})

args             <- commandArgs(trailingOnly = TRUE)
rds_file         <- args[1]
group_col        <- args[2]              # e.g. "Group" — metadata column defining conditions
ident1           <- args[3]              # e.g. "KO" — condition to test
ident2           <- args[4]              # e.g. "Ctl" — reference condition
min_pct          <- as.numeric(args[5])  # e.g. 0.1
logfc_thresh     <- as.numeric(args[6])  # e.g. 0.25 — used for FindMarkers filtering only
padj_cutoff      <- as.numeric(args[7])  # e.g. 0.01 — used for output table filtering only
# volcano visual params — separate from filtering thresholds above
# pass "auto" (or omit) on first run; after inspecting the plot set exact values in the yaml and rerun
volcano_xlim     <- if (!is.na(args[8])  && args[8]  != "auto") as.numeric(args[8])  else NULL
volcano_fc       <- if (!is.na(args[9])  && args[9]  != "auto") as.numeric(args[9])  else NULL
volcano_p        <- if (!is.na(args[10]) && args[10] != "auto") as.numeric(args[10]) else NULL

sobj <- readRDS(rds_file)
# input can be the full annotated object (seurat_annotated.rds) or a subsetted object (e.g. ILC3s_subsetted.rds)
# in either case the script loops over whatever Idents are present (cell types or sub-clusters)

# IMPORTANT: use RNA assay for DE
DefaultAssay(sobj) <- "RNA"

# helper: build colCustom keyvals vector — red2 up, dodgerblue3 down, gray93 not DE
# cutoffLineType = 'blank' removes dashed lines; coloring communicates thresholds instead
make_keyvals <- function(df, fc_cut, p_cut) {
  kv <- ifelse(df$avg_log2FC >  fc_cut & df$p_val_adj < p_cut, "red2",
        ifelse(df$avg_log2FC < -fc_cut & df$p_val_adj < p_cut, "dodgerblue3",
               "gray93"))
  kv[is.na(kv)] <- "gray93"
  names(kv)[kv == "red2"]        <- "Upregulated"
  names(kv)[kv == "dodgerblue3"] <- "Downregulated"
  names(kv)[kv == "gray93"]      <- "Not DE"
  kv
}

make_volcano <- function(df, title, xlim_val, fc_cut, p_cut, padj_filter) {
  # auto defaults computed from data — adjust in yaml after first inspection
  if (is.null(xlim_val)) xlim_val <- ceiling(max(abs(df$avg_log2FC), na.rm = TRUE))
  if (is.null(fc_cut))   fc_cut   <- logfc_thresh
  if (is.null(p_cut))    p_cut    <- 1e-5

  keyvals  <- make_keyvals(df, fc_cut, p_cut)
  # label only genes that pass the visual thresholds
  sel_labs <- df$Gene[df$avg_log2FC > fc_cut & df$p_val_adj < p_cut |
                      df$avg_log2FC < -fc_cut & df$p_val_adj < p_cut]
  ylim_val <- ceiling(max(-log10(df$p_val_adj[df$p_val_adj > 0]), na.rm = TRUE)) + 5

  EnhancedVolcano(df,
    lab              = df$Gene,
    selectLab        = sel_labs,
    x                = "avg_log2FC",
    y                = "p_val_adj",
    xlim             = c(-xlim_val, xlim_val),
    ylim             = c(0, ylim_val),
    title            = title,
    pCutoff          = p_cut,
    FCcutoff         = fc_cut,
    colCustom        = keyvals,
    colAlpha         = 0.8,
    pointSize        = 2.8,
    labSize          = 3.5,
    labCol           = "black",
    drawConnectors   = TRUE,
    boxedLabels      = TRUE,
    borderWidth      = 0.5,
    gridlines.major  = FALSE,
    gridlines.minor  = FALSE,
    legendPosition   = "right"
  )
}

# --- whole-object DE: all cells, ident1 vs ident2 regardless of cluster ---
Idents(sobj_whole <- sobj) <- sobj@meta.data[[group_col]]
DE_whole <- tryCatch(
  FindMarkers(sobj_whole, ident.1 = ident1, ident.2 = ident2,
              min.pct = min_pct, logfc.threshold = logfc_thresh),
  error = function(e) NULL
)
if (!is.null(DE_whole) && nrow(DE_whole) > 0) {
  DE_whole$Gene <- rownames(DE_whole)
  rownames(DE_whole) <- NULL
  write_xlsx(DE_whole, path = "DE_whole_object.xlsx")

  p_whole <- make_volcano(DE_whole,
    title      = paste0("All cells: ", ident1, " vs ", ident2),
    xlim_val   = volcano_xlim,
    fc_cut     = volcano_fc,
    p_cut      = volcano_p,
    padj_filter = padj_cutoff)
  ggsave("volcano_whole_object.pdf", plot = p_whole, width = 8, height = 7)
}

cluster_names   <- as.character(sort(unique(Idents(sobj))))
DE_results_list <- list()

# loop through clusters and run FindMarkers per cluster
for (cluster in cluster_names) {
  cluster_sub <- subset(sobj, idents = cluster)
  Idents(cluster_sub) <- cluster_sub@meta.data[[group_col]]

  DE <- tryCatch(
    FindMarkers(cluster_sub, ident.1 = ident1, ident.2 = ident2,
                min.pct = min_pct, logfc.threshold = logfc_thresh),
    error = function(e) NULL
  )

  if (!is.null(DE) && nrow(DE) > 0) {
    DE$Gene    <- rownames(DE)
    DE$Cluster <- cluster
    rownames(DE) <- NULL
    DE_results_list[[paste0("DE_", cluster)]] <- DE
  }
}

# combine all clusters
all_DE <- bind_rows(DE_results_list)

# filter significant results
sign_DE <- all_DE %>%
  filter(p_val_adj < padj_cutoff) %>%
  group_by(Cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE)

write_xlsx(all_DE,   path = "DE_all_results.xlsx")
write_xlsx(sign_DE,  path = "DE_significant_results.xlsx")

# volcano plot per cluster
for (cluster in unique(all_DE$Cluster)) {
  df <- all_DE %>% filter(Cluster == cluster)
  p  <- make_volcano(df,
    title       = paste0(cluster, ": ", ident1, " vs ", ident2),
    xlim_val    = volcano_xlim,
    fc_cut      = volcano_fc,
    p_cut       = volcano_p,
    padj_filter = padj_cutoff)
  ggsave(paste0("volcano_cluster", cluster, ".pdf"), plot = p, width = 8, height = 7)
}

cat("Total DE genes across all clusters:", nrow(all_DE), "\n")
cat("Significant DE genes (padj <", padj_cutoff, "):", nrow(sign_DE), "\n")
