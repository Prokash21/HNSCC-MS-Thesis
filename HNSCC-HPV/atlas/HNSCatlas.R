library(Seurat)
library(ggplot2)

setwd("~/DWCT/Thesis-Cancer/Thesis/HNSCC-HPV/atlas")
atlas <- readRDS("HNSCC_Atlas.rds")

DimPlot(atlas, reduction = "umap", group.by = "Cell_Type", label = TRUE, repel = TRUE) +
  ggtitle("HNSCC: Annotated Cell Types")
