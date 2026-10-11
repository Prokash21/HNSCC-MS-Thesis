library(Seurat)
library(ggplot2)

setwd("~/DWCT/Thesis-Cancer/Thesis/HNSCC-HPV/atlas")
atlas <- readRDS("HNSCC_Atlas.rds")

DimPlot(atlas, reduction = "umap", group.by = "Cell_Type", label = TRUE, repel = TRUE) +
  ggtitle("HNSCC: Annotated Cell Types")

rn <- table(atlas@meta.data$Source)


unique(atlas@meta.data$Original_Cell_Type)
table(atlas@meta.data$Cell_Type)
