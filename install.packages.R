 install.packages(c("Seurat", "ggplot2", "dplyr", "tidyr", "patchwork",
    "ggrepel", "ggsci", "writexl", "readxl", "ggalluvial",
    "hdf5r", "scales", "remotes"))
  
  BiocManager::install(c("EnhancedVolcano", "clusterProfiler", "ComplexHeatmap",
    "org.Mm.eg.db", "org.Hs.eg.db"))

  remotes::install_github(c("satijalab/seurat-wrappers",
    "cole-trapnell-lab/monocle3",
    "jinworks/CellChat",
    "immunogenomics/harmony"))