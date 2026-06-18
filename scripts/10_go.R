suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(clusterProfiler)
  library(ggplot2)
  library(writexl)
})

args           <- commandArgs(trailingOnly = TRUE)
de_file        <- args[1]              # DE_all_results.xlsx from script 09 (per-cluster DE)
whole_de_file  <- args[2]              # DE_whole_object.xlsx from script 09 (pass "NA" to skip)
organism       <- args[3]              # "mouse" or "human"
ont            <- args[4]              # e.g. "BP" (biological process), "MF", "CC", or "ALL"
n_top_terms    <- as.integer(args[5])  # e.g. 10 — top GO terms to plot per cluster
min_gs_size    <- as.integer(args[6])  # e.g. 10 — minimum gene set size
max_gs_size    <- as.integer(args[7])  # e.g. 100 — maximum gene set size
score_type     <- args[8]              # e.g. "pos" — gseGO score type
pvalue_cutoff  <- as.numeric(args[9])  # e.g. 0.05 — p-value cutoff for gseGO
direction      <- args[10]             # "up", "down", or "both"

# select OrgDb based on organism
if (organism == "mouse") {
  library(org.Mm.eg.db)
  orgdb <- org.Mm.eg.db
} else {
  library(org.Hs.eg.db)
  orgdb <- org.Hs.eg.db
}

# helper: run gseGO on a gene list and save outputs with a given label and direction
run_gsego <- function(gene_list, label, direction, orgdb, ont, n_top_terms,
                      min_gs_size, max_gs_size, score_type, pvalue_cutoff) {
  if (length(gene_list) < 5) {
    cat("Skipping", label, direction, "— fewer than 5 genes\n")
    return(NULL)
  }

  GO_res <- tryCatch(
    gseGO(
      geneList      = gene_list,
      ont           = ont,
      keyType       = "SYMBOL",
      minGSSize     = min_gs_size,
      maxGSSize     = max_gs_size,
      scoreType     = score_type,
      pvalueCutoff  = pvalue_cutoff,
      verbose       = FALSE,
      OrgDb         = orgdb,
      nPermSimple   = 10000,
      pAdjustMethod = "none"
    ),
    error = function(e) NULL
  )

  if (is.null(GO_res) || nrow(as.data.frame(GO_res)) == 0) {
    cat("No GO results for", label, direction, "\n")
    return(NULL)
  }

  GO_df <- as.data.frame(GO_res) %>% arrange(p.adjust)
  GO_df$enrichmentScore <- round(GO_df$enrichmentScore, 2)

  write_xlsx(GO_df, path = paste0("GO_", label, "_", direction, ".xlsx"))

  top_GO <- head(GO_df, n_top_terms)
  # label each term as "Description (GO:ID)" and rank by significance
  top_GO$NegativeLogPValue <- -log10(top_GO$pvalue)
  top_GO$ID_Description     <- paste0(top_GO$Description, " (", top_GO$ID, ")")

  p <- ggplot(top_GO, aes(x = NegativeLogPValue,
                          y = reorder(ID_Description, NegativeLogPValue))) +
    geom_col(fill = "black") +
    theme_classic(base_size = 12) +
    xlab("-log10(p-value)") +
    ylab(" ") +
    ggtitle(paste0(label, " — GO ", ont, " (", direction, ")")) +
    theme(axis.text    = element_text(color = "black", size = 12),
          plot.margin  = margin(1, 1, 1, 2, "cm"))

  ggsave(paste0("GO_", label, "_", direction, "_barplot.pdf"), plot = p, width = 8, height = 5)
  GO_df
}

# helper: prepare gene lists for up and down, then run both passes
run_both_directions <- function(DE_df, label, orgdb, ont, n_top_terms,
                                min_gs_size, max_gs_size, score_type, pvalue_cutoff, direction) {
  # remove ribosomal and mitochondrial genes — these dominate GO results and are uninformative
  # Rp* = ribosomal (mouse), mt-* = mitochondrial (mouse); same logic applies to human (RP*, MT-)
  DE_df <- DE_df %>% filter(!grepl("^Rp[sl]|^mt-|^RP[SL]|^MT-", Gene))

  # upregulated: filter positive FC, sort decreasingly
  up <- DE_df %>% filter(avg_log2FC > 0) %>% arrange(desc(avg_log2FC))
  gene_list_up <- na.omit(sort(setNames(up$avg_log2FC, up$Gene), decreasing = TRUE))

  # downregulated: filter negative FC, take abs() — required by gseGO which needs positive ranked list
  down <- DE_df %>% filter(avg_log2FC < 0)
  gene_list_down <- na.omit(sort(setNames(abs(down$avg_log2FC), down$Gene), decreasing = TRUE))

  if (direction %in% c("up",   "both"))
    run_gsego(gene_list_up,   label, "upregulated",   orgdb, ont, n_top_terms,
              min_gs_size, max_gs_size, score_type, pvalue_cutoff)
  if (direction %in% c("down", "both"))
    run_gsego(gene_list_down, label, "downregulated", orgdb, ont, n_top_terms,
              min_gs_size, max_gs_size, score_type, pvalue_cutoff)
}

# --- whole-object GO ---
if (!is.na(whole_de_file) && whole_de_file != "NA" && file.exists(whole_de_file)) {
  DE_whole <- read_xlsx(whole_de_file)
  if (!"Gene" %in% colnames(DE_whole)) DE_whole$Gene <- rownames(DE_whole)
  run_both_directions(DE_whole, "whole_object", orgdb, ont, n_top_terms,
                      min_gs_size, max_gs_size, score_type, pvalue_cutoff, direction)
}

# --- per-cluster GO ---
all_DE          <- read_xlsx(de_file)
cluster_names   <- unique(all_DE$Cluster)

for (cluster in cluster_names) {
  DE_cluster <- all_DE %>% filter(Cluster == cluster)
  run_both_directions(DE_cluster, paste0("cluster", cluster), orgdb, ont, n_top_terms,
                      min_gs_size, max_gs_size, score_type, pvalue_cutoff, direction)
}

cat("GO analysis complete\n")
