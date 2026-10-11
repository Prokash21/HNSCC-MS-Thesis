## Publication-style figures for the HNSCatlas epithelial compartment.
##   Figure 1: every cell in the atlas, coloured by annotated cell type
##   Figure 2: the comparison set -- tumour cells + non-malignant epithelial
##             cells (re-embedded on their own), HPV status, cells per patient
##
## Inputs come from scripts/02_atlas_extract_epithelial.R.
## Run from the repo root:  Rscript scripts/03_atlas_figures.R

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggrastr)
  library(ggrepel)
  library(patchwork)
})

in_dir  <- "data/processed/HNSCatlas"
fig_dir <- "HNSCC-HPV/results/01_atlas"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

cells <- readRDS(file.path(in_dir, "atlas_cells_meta_umap.rds"))

## ---- shared style -------------------------------------------------------
base_size <- 7                       # pt, journal body size
ink       <- "#0b0b0b"               # primary text
ink2      <- "#52514e"               # secondary text / axis arrows
bg_grey   <- "#dcdbd5"               # unassigned / background cells
col_tumor <- "#e34948"
col_norm  <- "#4a3aa7"
col_hpv   <- c("HPV+" = "#2a78d6", "HPV-" = "#eb6834")

theme_umap <- function() {
  theme_void(base_size = base_size) +
    theme(text = element_text(colour = ink),
          legend.title = element_text(size = base_size, face = "bold"),
          legend.text = element_text(size = base_size),
          legend.key.size = unit(3, "mm"),
          plot.title = element_text(size = base_size + 1, face = "bold"),
          plot.subtitle = element_text(size = base_size, colour = ink2),
          plot.tag = element_text(size = base_size + 3, face = "bold"),
          plot.margin = margin(4, 4, 4, 4))
}

## short UMAP arrows at the bottom-left instead of full axes
umap_arrows <- function(x, y, frac = 0.12) {
  x0 <- min(x) - 0.03 * diff(range(x))
  y0 <- min(y) - 0.03 * diff(range(y))
  dx <- diff(range(x)) * frac
  dy <- diff(range(y)) * frac
  arr <- arrow(length = unit(1.2, "mm"), type = "closed")
  list(
    annotate("segment", x = x0, y = y0, xend = x0 + dx, yend = y0,
             arrow = arr, linewidth = 0.3, colour = ink2),
    annotate("segment", x = x0, y = y0, xend = x0, yend = y0 + dy,
             arrow = arr, linewidth = 0.3, colour = ink2),
    annotate("text", x = x0, y = y0, label = "UMAP 1", hjust = 0, vjust = 1.5,
             size = (base_size - 1) / .pt, colour = ink2),
    annotate("text", x = x0, y = y0, label = "UMAP 2", hjust = 0, vjust = -0.6,
             angle = 90, size = (base_size - 1) / .pt, colour = ink2)
  )
}

## label anchor = a real cell in the densest grid bin of each group, so the
## label sits on the cluster core rather than between scattered cells
label_pos <- function(df, group, x = "UMAP_1", y = "UMAP_2", bins = 20) {
  bx <- cut(df[[x]], bins, labels = FALSE)
  by <- cut(df[[y]], bins, labels = FALSE)
  do.call(rbind, lapply(split(seq_len(nrow(df)), df[[group]], drop = TRUE), function(i) {
    b <- paste(bx[i], by[i])
    j <- i[b == names(which.max(table(b)))]
    k <- j[which.min((df[[x]][j] - median(df[[x]][j]))^2 + (df[[y]][j] - median(df[[y]][j]))^2)]
    df[k, c(group, x, y)]
  }))
}

## "Name (12,345)" legend labels
n_label <- function(lab, n) {
  sprintf("%s (%s)", lab, format(as.integer(n), big.mark = ",", trim = TRUE))
}

save_fig <- function(p, name, w, h) {
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h,
         units = "mm", dpi = 600, device = ragg::agg_png, bg = "white")
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h,
         units = "mm", bg = "white")
  message("Wrote ", file.path(fig_dir, name), ".png/.pdf")
}

## ---- Figure 1: all cells, by annotated cell type -------------------------
## colours grouped by lineage (hue family), every cluster labelled directly
type_cols <- c(
  Tumor = col_tumor, Epithelial = col_norm,                       # epithelial
  CD4T = "#6da7ec", CD8T = "#184f95", NK = "#1baf7a",             # T / NK
  Multi = "#898781",
  Bcell = "#c98500", PlasmaCell = "#eda100",                      # B lineage
  Macrophage = "#eb6834", Monocyte = "#e87ba4",                   # myeloid
  Neutrophils = "#943b16", panDC = "#9085e9", Mast = "#d55181",
  Fibroblast = "#008300", Endothelial = "#63c55d",                # stroma
  Unassigned = bg_grey)
type_names <- c(
  Tumor = "Tumor", Epithelial = "Non-malignant epithelial",
  CD4T = "CD4+ T", CD8T = "CD8+ T", NK = "NK", Multi = "Multi (ambiguous)",
  Bcell = "B cell", PlasmaCell = "Plasma cell",
  Macrophage = "Macrophage", Monocyte = "Monocyte", Neutrophils = "Neutrophil",
  panDC = "Dendritic cell", Mast = "Mast cell",
  Fibroblast = "Fibroblast", Endothelial = "Endothelial",
  Unassigned = "Unassigned")

cells$type <- factor(ifelse(is.na(cells$Cell_Type), "Unassigned", cells$Cell_Type),
                     levels = names(type_cols))
n_type <- table(cells$type)
legend_labels <- n_label(type_names[names(type_cols)], n_type[names(type_cols)])

## unassigned underneath, most types in random order (none hides another), and
## the rare comparison cells (non-malignant epithelial) on top so they stay visible
set.seed(1)
ord <- c(which(cells$type == "Unassigned"),
         sample(which(!cells$type %in% c("Unassigned", "Epithelial"))),
         which(cells$type == "Epithelial"))
## direct labels for real cell types only (not the ambiguous/unassigned bins)
lab1 <- label_pos(cells[!cells$type %in% c("Unassigned", "Multi"), ], "type")
lab1$label <- type_names[as.character(lab1$type)]

fig1 <- ggplot(cells[ord, ], aes(UMAP_1, UMAP_2, colour = type)) +
  geom_point_rast(size = 0.15, stroke = 0, shape = 16, raster.dpi = 600, dev = "ragg") +
  geom_text_repel(data = lab1, aes(label = label), colour = ink,
                  size = base_size / .pt, bg.color = "white", bg.r = 0.12,
                  min.segment.length = Inf, max.overlaps = Inf, seed = 1) +
  umap_arrows(cells$UMAP_1, cells$UMAP_2) +
  scale_colour_manual(values = type_cols, breaks = names(type_cols),
                      labels = legend_labels, name = NULL, drop = FALSE) +
  guides(colour = guide_legend(override.aes = list(size = 2.2))) +
  coord_fixed() +
  labs(title = "HNSCC single-cell atlas",
       subtitle = sprintf("%s cells, %d patients (HNSCatlas)",
                          format(nrow(cells), big.mark = ","),
                          length(unique(cells$Patient)))) +
  theme_umap()
save_fig(fig1, "Fig1_atlas_all_cells", w = 180, h = 125)

## ---- Figure 2: tumour + non-malignant epithelial cells -------------------
## re-embed the subset on the atlas' integrated assay (batch-corrected)
emb_file <- file.path(in_dir, "epithelial_umap.rds")
if (!file.exists(emb_file)) {
  message("Re-embedding the epithelial subset...")
  epi <- readRDS(file.path(in_dir, "atlas_epithelial.rds"))
  DefaultAssay(epi) <- "integrated"
  feats <- rownames(epi[["integrated"]])
  epi <- ScaleData(epi, features = feats, verbose = FALSE)
  epi <- RunPCA(epi, features = feats, npcs = 30, verbose = FALSE)
  epi <- RunUMAP(epi, dims = 1:30, reduction.name = "umap_epi", seed.use = 42,
                 verbose = FALSE)
  e <- Embeddings(epi, "umap_epi")
  saveRDS(data.frame(cell = rownames(e), sUMAP_1 = e[, 1], sUMAP_2 = e[, 2]),
          emb_file)
  rm(epi); invisible(gc())
}
sub <- merge(cells, readRDS(emb_file), by = "cell")
sub$Compartment <- factor(ifelse(sub$Cell_Type == "Tumor", "Tumor", "Non-malignant epithelial"),
                          levels = c("Tumor", "Non-malignant epithelial"))
comp_cols <- c("Tumor" = col_tumor, "Non-malignant epithelial" = col_norm)
n_comp <- table(sub$Compartment)
comp_labels <- n_label(names(comp_cols), n_comp[names(comp_cols)])
sub$HPV <- factor(sub$HPV, levels = names(col_hpv))

## A: where the selected cells sit in the whole atlas
cells$sel <- factor(ifelse(cells$cell %in% sub$cell,
                           ifelse(cells$Cell_Type == "Tumor", "Tumor", "Non-malignant epithelial"),
                           "Other cells"),
                    levels = c("Other cells", "Tumor", "Non-malignant epithelial"))
pA <- ggplot(cells[order(cells$sel), ], aes(UMAP_1, UMAP_2, colour = sel)) +
  geom_point_rast(size = 0.08, stroke = 0, shape = 16, raster.dpi = 600, dev = "ragg") +
  umap_arrows(cells$UMAP_1, cells$UMAP_2) +
  scale_colour_manual(values = c("Other cells" = bg_grey, comp_cols), guide = "none") +
  coord_fixed() +
  labs(title = "Selected cells in the atlas") +
  theme_umap()

## B: subset re-embedded, by compartment (rarer non-malignant cells on top)
pB <- ggplot(sub[order(sub$Compartment), ], aes(sUMAP_1, sUMAP_2, colour = Compartment)) +
  geom_point_rast(size = 0.2, stroke = 0, shape = 16, raster.dpi = 600, dev = "ragg") +
  umap_arrows(sub$sUMAP_1, sub$sUMAP_2) +
  scale_colour_manual(values = comp_cols, breaks = names(comp_cols),
                      labels = comp_labels, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2.2))) +
  coord_fixed() +
  labs(title = "Epithelial compartment") +
  theme_umap() + theme(legend.position = "bottom", legend.direction = "vertical")

## C: same embedding, by HPV status
set.seed(2)
n_hpv <- table(sub$HPV)
pC <- ggplot(sub[sample(nrow(sub)), ], aes(sUMAP_1, sUMAP_2, colour = HPV)) +
  geom_point_rast(size = 0.2, stroke = 0, shape = 16, raster.dpi = 600, dev = "ragg") +
  umap_arrows(sub$sUMAP_1, sub$sUMAP_2) +
  ## breaks pin each label to its own colour (labels alone are matched by position)
  scale_colour_manual(values = col_hpv, breaks = names(col_hpv), name = NULL,
                      labels = n_label(names(col_hpv), n_hpv[names(col_hpv)])) +
  guides(colour = guide_legend(override.aes = list(size = 2.2))) +
  coord_fixed() +
  labs(title = "HPV status") +
  theme_umap() + theme(legend.position = "bottom", legend.direction = "vertical")

## D: cells per patient -- what a per-patient (pseudobulk) comparison can use
pp <- as.data.frame(table(Patient = sub$Patient, Compartment = sub$Compartment))
pp <- merge(pp, unique(sub[, c("Patient", "HPV")]), by = "Patient")
pp <- pp[pp$Freq > 0, ]                       # zeros cannot sit on a log axis
## order by tumour cells, counting every patient (some have none)
tumor_n <- tapply(sub$Compartment == "Tumor", sub$Patient, sum)
pp$Patient <- factor(pp$Patient, levels = names(sort(tumor_n, decreasing = TRUE)))
pD <- ggplot(pp, aes(Patient, Freq, colour = Compartment)) +
  geom_hline(yintercept = 20, linetype = "22", linewidth = 0.3, colour = ink2) +
  geom_point(size = 1.1, stroke = 0, shape = 16) +
  geom_text(data = data.frame(HPV = factor(names(col_hpv)[1], levels = names(col_hpv))),
            aes(x = -Inf, y = 20, label = "20 cells"), inherit.aes = FALSE,
            hjust = -0.1, vjust = -0.5, size = (base_size - 1) / .pt, colour = ink2) +
  facet_grid(. ~ HPV, scales = "free_x", space = "free_x") +
  scale_y_log10(labels = scales::label_comma()) +
  scale_colour_manual(values = comp_cols, breaks = names(comp_cols), name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  labs(title = "Cells per patient", x = "Patient", y = "Cells (log scale)") +
  theme_classic(base_size = base_size) +
  theme(text = element_text(colour = ink),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = base_size - 2),
        axis.line = element_line(linewidth = 0.3, colour = "#c3c2b7"),
        axis.ticks = element_line(linewidth = 0.3, colour = "#c3c2b7"),
        panel.grid.major.y = element_line(linewidth = 0.2, colour = "#e1e0d9"),
        strip.background = element_blank(),
        strip.text = element_text(size = base_size, face = "bold"),
        legend.position = "inside", legend.position.inside = c(1, 1),
        legend.justification = c(1, 1), legend.background = element_blank(),
        legend.text = element_text(size = base_size), legend.key.size = unit(3, "mm"),
        plot.title = element_text(size = base_size + 1, face = "bold"),
        plot.tag = element_text(size = base_size + 3, face = "bold"))

fig2 <- (pA | pB | pC) / pD +
  plot_layout(heights = c(1.15, 0.85)) +
  plot_annotation(tag_levels = "A")
save_fig(fig2, "Fig2_epithelial_compartment", w = 180, h = 150)
