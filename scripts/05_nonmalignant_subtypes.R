## Split the atlas' non-malignant "Epithelial" cells into subtypes, so the
## tumour-vs-normal DGE can use clean squamous epithelium as its reference.
## (The first-pass DGE showed the group mixes squamous cells with salivary /
## mucous gland cells and immune-contaminated cells.)
##
## Input:  data/processed/HNSCatlas/atlas_epithelial.rds (scripts/02)
## Output: data/processed/HNSCatlas/nonmalignant_subtypes.rds (cell -> subtype)
##         HNSCC-HPV/results/02_normal_cleaning/tables/nonmalignant_clusters.csv
##         HNSCC-HPV/results/02_normal_cleaning/Fig3_nonmalignant_cleaning.png/.pdf
## Run from the repo root:  Rscript scripts/05_nonmalignant_subtypes.R

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggrastr)
  library(ggrepel)
  library(patchwork)
})

in_dir  <- "data/processed/HNSCatlas"
nm_file <- file.path(in_dir, "nonmalignant_epithelial.rds")

## small object with only the non-malignant cells (RNA assay), cached
if (!file.exists(nm_file)) {
  message("Extracting non-malignant epithelial cells...")
  epi <- readRDS(file.path(in_dir, "atlas_epithelial.rds"))
  nm <- subset(epi, subset = Compartment == "Non-malignant epithelial")
  rm(epi)
  invisible(gc())
  DefaultAssay(nm) <- "RNA"
  nm <- DietSeurat(nm, assays = "RNA", layers = c("counts", "data"))
  saveRDS(nm, nm_file)
}
nm <- readRDS(nm_file)
message(sprintf("%d non-malignant epithelial cells from %d patients",
                ncol(nm), length(unique(nm$Patient))))

## sub-cluster
set.seed(1)
nm <- NormalizeData(nm, verbose = FALSE)
nm <- FindVariableFeatures(nm, nfeatures = 2000, verbose = FALSE)
nm <- ScaleData(nm, verbose = FALSE)
nm <- RunPCA(nm, npcs = 30, verbose = FALSE)
nm <- FindNeighbors(nm, dims = 1:20, verbose = FALSE)
nm <- FindClusters(nm, resolution = 0.5, random.seed = 1, verbose = FALSE)
nm <- RunUMAP(nm, dims = 1:20, seed.use = 42, verbose = FALSE)

## marker-set scores
sigs <- list(
  Squamous  = c("KRT5", "KRT14", "KRT6A", "KRT13", "KRT4", "TP63", "DSG3", "SPRR3", "CRNN", "S100A2"),
  Glandular = c("MUC5B", "BPIFB1", "PIGR", "SCGB3A1", "ZG16B", "STATH", "PIP", "AZGP1", "LTF", "TFF3", "KRT7"),
  Ciliated  = c("FOXJ1", "CAPS", "TPPP3", "PIFO", "RSPH1", "C20orf85"),
  Immune    = c("PTPRC", "CD74", "HLA-DRA", "LYZ", "CD14", "C1QC", "AIF1", "CD3E", "MS4A1"),
  Stromal   = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "VWF", "PECAM1"))
sigs <- lapply(sigs, intersect, rownames(nm))
nm <- AddModuleScore(nm, features = sigs, name = "sig_", seed = 1)
names(nm@meta.data)[match(paste0("sig_", seq_along(sigs)), names(nm@meta.data))] <- names(sigs)

## per-cluster summary: size, patients, scores, % cells expressing key markers
key <- c("KRT5", "KRT13", "TP63", "MUC5B", "BPIFB1", "PIGR", "FOXJ1", "PTPRC", "COL1A2")
expr <- FetchData(nm, vars = key, layer = "data") > 0
cl <- nm$seurat_clusters
summ <- data.frame(
  cluster  = levels(cl),
  n_cells  = as.integer(table(cl)),
  patients = as.integer(tapply(nm$Patient, cl, function(p) length(unique(p)))),
  top_pt   = tapply(nm$Patient, cl, function(p) sprintf("%s %.0f%%", names(which.max(table(p))), 100 * max(table(p)) / length(p))),
  sapply(names(sigs), function(s) round(tapply(nm@meta.data[[s]], cl, mean), 2)),
  round(100 * apply(expr, 2, function(e) tapply(e, cl, mean))))
print(summ, row.names = FALSE)

message("\nTop markers per cluster:")
mk <- FindAllMarkers(nm, only.pos = TRUE, logfc.threshold = 0.5, min.pct = 0.25,
                     max.cells.per.ident = 300, features = VariableFeatures(nm),
                     verbose = FALSE)
summ$top_markers <- vapply(levels(cl), function(k) {
  paste(head(mk$gene[mk$cluster == k & mk$p_val_adj < 0.05], 8), collapse = ", ")
}, character(1))

## ---- label clusters (rules, applied in this order) ------------------------
## 1. glandular:  >= 30% of cells express PIGR, MUC5B or BPIFB1 (salivary /
##                mucous gland, incl. myoepithelial and ciliated cells)
## 2. likely malignant: one patient makes up >= 80% of the cluster (tumour cells
##                cluster by patient, normal cells by cell type), or >= 25% of
##                cells carry HPV16 or cancer-testis antigen reads (background
##                in gland cells is ~10%)
## 3. squamous:   remaining clusters with a positive squamous score (> 0.25)
counts <- LayerData(nm, assay = "RNA", layer = "counts")
pct_any <- function(genes) {
  genes <- intersect(genes, rownames(counts))
  round(100 * tapply(colSums(counts[genes, , drop = FALSE]) > 0, cl, mean))
}
summ$pct_gland   <- pct_any(c("PIGR", "MUC5B", "BPIFB1"))
summ$pct_viral   <- pct_any(grep("^HPV16", rownames(counts), value = TRUE))
summ$pct_ct_ag   <- pct_any(c("MAGEA1", "MAGEA3", "MAGEA4", "MAGEA6", "MAGEC2", "CSAG1"))
summ$top_share   <- as.integer(tapply(nm$Patient, cl, function(p) round(100 * max(table(p)) / length(p))))
summ$subtype <- with(summ, ifelse(pct_gland >= 30, "Glandular",
                       ifelse(top_share >= 80 | pct_viral >= 25 | pct_ct_ag >= 25, "Likely malignant",
                       ifelse(Squamous > 0.25, "Squamous", "Other"))))
summ$subtype[summ$subtype == "Glandular" & summ$Ciliated > 1] <- "Ciliated"
print(summ[, c("cluster", "n_cells", "patients", "top_share", "Squamous", "Glandular",
               "pct_gland", "pct_viral", "pct_ct_ag", "subtype", "top_markers")], row.names = FALSE)

## per cell: subtype of its cluster; squamous cells that look like immune or
## stromal doublets (score > 0.5) are dropped from the clean reference
nm$subtype <- summ$subtype[match(as.character(cl), summ$cluster)]
nm$squamous_clean <- nm$subtype == "Squamous" & nm$Immune <= 0.5 & nm$Stromal <= 0.5
message(sprintf("\nClean squamous reference: %d of %d cells (%d squamous-cluster cells dropped as doublets)",
                sum(nm$squamous_clean), ncol(nm),
                sum(nm$subtype == "Squamous") - sum(nm$squamous_clean)))
print(table(nm$subtype))

write.csv(summ, "HNSCC-HPV/results/02_normal_cleaning/tables/nonmalignant_clusters.csv", row.names = FALSE)
saveRDS(data.frame(cell = colnames(nm), Patient = nm$Patient, HPV = nm$HPV,
                   cluster = as.character(cl), subtype = nm$subtype,
                   squamous_clean = nm$squamous_clean,
                   nm@meta.data[, names(sigs)], Embeddings(nm, "umap")),
        file.path(in_dir, "nonmalignant_subtypes.rds"))
saveRDS(nm, file.path(in_dir, "nonmalignant_clustered.rds"))

## ---- Figure 3: how the normal reference was cleaned -----------------------
base_size <- 7
ink  <- "#0b0b0b"
ink2 <- "#52514e"
fmt  <- function(n) format(n, big.mark = ",", trim = TRUE)
grp_cols <- c("Squamous (kept)"             = "#4a3aa7",
              "Squamous, doublet (dropped)" = "#c3c2b7",
              "Glandular"                   = "#1baf7a",
              "Ciliated"                    = "#eda100",
              "Likely malignant"            = "#e34948")
grp <- ifelse(nm$squamous_clean, "Squamous (kept)",
              ifelse(nm$subtype == "Squamous", "Squamous, doublet (dropped)", nm$subtype))
df <- data.frame(Embeddings(nm, "umap"), cluster = as.character(cl),
                 group = factor(grp, levels = names(grp_cols)))
names(df)[1:2] <- c("UMAP_1", "UMAP_2")
lab <- aggregate(cbind(UMAP_1, UMAP_2) ~ cluster, df, median)
n_grp <- table(df$group)

pa <- ggplot(df, aes(UMAP_1, UMAP_2, colour = group)) +
  geom_point_rast(size = 0.45, stroke = 0, shape = 16, raster.dpi = 600, dev = "ragg") +
  ## force = 0: number sits on the cluster, the halo keeps it readable
  geom_text_repel(data = lab, aes(label = cluster), colour = ink, size = (base_size - 1) / .pt,
                  bg.color = "white", bg.r = 0.12, force = 0, force_pull = 0,
                  min.segment.length = Inf, max.overlaps = Inf, seed = 1) +
  scale_colour_manual(values = grp_cols, breaks = names(grp_cols), name = NULL,
                      labels = sprintf("%s (%s)", names(grp_cols), fmt(as.integer(n_grp[names(grp_cols)])))) +
  guides(colour = guide_legend(override.aes = list(size = 2.2), ncol = 1)) +
  coord_fixed() +
  labs(title = "Atlas 'non-malignant epithelial' cells, sub-clustered",
       subtitle = sprintf("%s cells, %d patients; numbers = clusters", fmt(ncol(nm)),
                          length(unique(nm$Patient)))) +
  theme_void(base_size = base_size) +
  theme(text = element_text(colour = ink), legend.position = "right",
        legend.key.size = unit(3, "mm"), legend.key.spacing.y = unit(1, "mm"),
        plot.title = element_text(size = base_size + 1, face = "bold"),
        plot.subtitle = element_text(size = base_size - 0.5, colour = ink2),
        plot.tag = element_text(size = base_size + 3, face = "bold"))

## marker dot plot, clusters ordered by their label
feats <- list(Squamous  = c("KRT5", "KRT14", "TP63", "KRT13", "SPRR3"),
              Gland     = c("MUC5B", "BPIFB1", "PIGR", "ACTA2"),
              Cilia     = "FOXJ1",
              Malignant = c("HPV16-E5", "CDKN2A", "MAGEA4", "MKI67"),
              Doublet   = c("PTPRC", "COL1A2", "HBB"))
feats <- lapply(feats, intersect, rownames(nm))
ord <- summ[order(match(summ$subtype, c("Squamous", "Glandular", "Ciliated", "Likely malignant")),
                  as.integer(summ$cluster)), ]
nm$cluster_lab <- factor(paste(cl, nm$subtype, sep = "  "),
                         levels = rev(paste(ord$cluster, ord$subtype, sep = "  ")))
pb <- DotPlot(nm, features = feats, group.by = "cluster_lab", cols = c("#e1e0d9", "#184f95"),
              dot.scale = 3) +
  labs(title = "Marker genes per cluster", x = NULL, y = NULL) +
  theme(text = element_text(size = base_size, colour = ink),
        axis.text.x = element_text(size = base_size - 1, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = base_size - 1),
        strip.text = element_text(size = base_size - 1, face = "bold"),
        legend.text = element_text(size = base_size - 1),
        legend.title = element_text(size = base_size - 1),
        plot.title = element_text(size = base_size + 1, face = "bold"),
        plot.tag = element_text(size = base_size + 3, face = "bold"))

fig3 <- pa / pb + plot_layout(heights = c(0.9, 1.1)) + plot_annotation(tag_levels = "a")
for (ext in c("png", "pdf")) {
  ggsave(file.path("HNSCC-HPV/results/02_normal_cleaning", paste0("Fig3_nonmalignant_cleaning.", ext)), fig3,
         width = 180, height = 175, units = "mm", dpi = 600, bg = "white",
         device = if (ext == "png") ragg::agg_png else "pdf")
}
message("Wrote HNSCC-HPV/results/02_normal_cleaning/Fig3_nonmalignant_cleaning.png/.pdf")
