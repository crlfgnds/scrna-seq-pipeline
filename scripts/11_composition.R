suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggalluvial)
  library(ggsci)
  library(writexl)
})

args         <- commandArgs(trailingOnly = TRUE)
rds_file     <- args[1]
celltype_col <- args[2]   # e.g. "cell_type" — cluster/cell type labels
group_col    <- args[3]   # e.g. "Group" — condition column (Control vs KO)
sample_col   <- args[4]   # e.g. "experiment" — biological replicate column (per mouse/sample)
ident1       <- args[5]   # e.g. "Control" — reference group for t-test
ident2       <- args[6]   # e.g. "a-Ly6G" — comparison group for t-test
palette_arg  <- args[7]   # named ggsci palette ("npg","jco","lancet","nejm") OR comma-separated hex codes
                           # e.g. "npg" or "#3C5488FF,#E64B35FF,#4DBBD599,#00A087FF"
time_col     <- args[8]   # e.g. "time" — timepoint column (pass "none" if no time variable)
time_levels  <- strsplit(args[9], ",")[[1]]   # e.g. "steadyState,12h,24h,48h" or "none"

sobj <- readRDS(rds_file)
md   <- sobj@meta.data

# --- 1. per-sample percentages (when there are hashtags) ---
# this is the statistically meaningful unit — not bulk cell counts
sample_pct <- md %>%
  group_by(.data[[sample_col]], .data[[group_col]], .data[[celltype_col]]) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(.data[[sample_col]], .data[[group_col]]) %>%
  mutate(percentage = n / sum(n) * 100) %>%
  ungroup()

# t-test per cluster/cell type comparing the two groups
sample_pct <- sample_pct %>%
  group_by(.data[[celltype_col]]) %>%
  mutate(p_value = tryCatch(
    t.test(percentage[.data[[group_col]] == ident1],
           percentage[.data[[group_col]] == ident2])$p.value,
    error = function(e) NA_real_
  )) %>%
  ungroup()

write_xlsx(sample_pct, path = "composition_per_sample.xlsx")

# mean percentage per group × celltype — for stacked bar plot
mean_pct <- sample_pct %>%
  group_by(.data[[group_col]], .data[[celltype_col]]) %>%
  summarise(mean_percentage = mean(percentage), .groups = "drop")

n_celltypes <- length(unique(md[[celltype_col]]))
pal <- if (grepl("^#", palette_arg)) {
  strsplit(palette_arg, ",")[[1]]          # custom hex codes passed directly
} else {
  switch(palette_arg,
    "npg"    = pal_npg("nrc")(n_celltypes),
    "jco"    = pal_jco()(n_celltypes),
    "lancet" = pal_lancet()(n_celltypes),
    "nejm"   = pal_nejm()(n_celltypes),
    pal_npg("nrc")(n_celltypes)            # default fallback
  )
}

# --- 2. boxplot + dotplot per cluster across replicates ---
# shows per-sample variation 
p_box <- ggplot(sample_pct,
                aes(x     = .data[[celltype_col]],
                    y     = percentage,
                    fill  = .data[[group_col]])) +
  geom_boxplot(position = position_dodge2(preserve = "single", width = 0.5),
               outlier.shape = NA) +
  geom_dotplot(binaxis  = "y",
               stackdir = "center",
               position = position_dodge(0.75),
               aes(fill = .data[[group_col]]),
               dotsize  = 0.7) +
  labs(x = NULL, y = "Frequency (%)") +
  theme_classic(base_size = 12) +
  theme(axis.text.x  = element_text(size = 12, color = "black", angle = 45, hjust = 1),
        axis.text.y  = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 14),
        axis.line    = element_line(linewidth = 0.7))
ggsave("composition_boxplot.pdf", plot = p_box, width = 8, height = 5)

# --- 3. stacked bar of mean percentages ---
p_bar <- ggplot(mean_pct,
                aes(x    = .data[[group_col]],
                    y    = mean_percentage,
                    fill = .data[[celltype_col]])) +
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  scale_fill_manual(values = pal, name = celltype_col) +
  labs(x = NULL, y = "Frequency") +
  theme_classic(base_size = 12) +
  theme(axis.text.y  = element_text(size = 14, color = "black"),
        axis.text.x  = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 15),
        axis.line    = element_line(linewidth = 0.7),
        legend.position = "right")
ggsave("composition_barplot.pdf", plot = p_bar, width = 5, height = 5)

# --- 4. alluvial — only when there are different timepoints ---
if (time_col != "none") {
  # rebuild with time dimension — proportion per group × time × celltype
  md[[time_col]]       <- factor(md[[time_col]], levels = time_levels)
  md[[group_col]]      <- factor(md[[group_col]])
  md[[celltype_col]]   <- factor(md[[celltype_col]])

  comp_time <- md %>%
    count(.data[[group_col]], .data[[time_col]], .data[[celltype_col]], name = "n_cells") %>%
    complete(.data[[group_col]], .data[[time_col]], .data[[celltype_col]],
             fill = list(n_cells = 0)) %>%
    group_by(.data[[group_col]], .data[[time_col]]) %>%
    mutate(pct = n_cells / sum(n_cells)) %>%
    ungroup()

  pal_named        <- pal
  names(pal_named) <- levels(md[[celltype_col]])

  p_alluvial <- ggplot(comp_time,
    aes(x        = .data[[time_col]],
        y        = pct,
        alluvium = .data[[celltype_col]],
        stratum  = .data[[celltype_col]],
        fill     = .data[[celltype_col]])) +
    ggalluvial::geom_alluvium(alpha = 0.8, knot.pos = 0.4, color = NA) +  # ribbons
    ggalluvial::geom_stratum(width = 0.35, color = "grey20") +             # bars
    facet_wrap(as.formula(paste("~", group_col)), nrow = 1) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_fill_manual(values = pal_named, drop = FALSE, name = celltype_col) +
    labs(x = "Timepoint", y = "Composition (%)") +
    theme_classic(base_size = 12) +
    theme(panel.spacing    = unit(12, "pt"),
          legend.position  = "right",
          axis.text.x      = element_text(hjust = 0.5))
  ggsave("composition_alluvial.pdf", plot = p_alluvial,
         width = 5 * length(unique(md[[group_col]])), height = 5)
}
