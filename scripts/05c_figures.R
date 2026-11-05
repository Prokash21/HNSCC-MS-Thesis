## Step 05c (part D): figures
##
## Purpose: Fig47 gene atlas heatmap (ComplexHeatmap, five blocks); Fig48 blood
##   background panel (Cillo PBMC per cell type x group, GTEx); Fig49 LFC
##   concordance scatters (C2 and C1 vs TCGA); Fig50 gene-gene co-expression
##   heatmaps (TCGA all, TCGA oral, single-cell malignant pseudobulks).
## Inputs:  results/step05c_gene_atlas/tables/{A_*, B_gene_atlas_matrix, C_*}.csv,
##          results/step04b_DE/tables/{D_C2_HPVneg_allpatients, D_C1_HPVneg_paired, E_TCGA_HPVneg_paired_tumour_vs_normal}.csv
## Outputs: results/step05c_gene_atlas/Fig47-Fig50 (png, pdf)
## Date: 2026-10-28

suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(ggplot2)
  library(patchwork)
})
source("config/paths.R")

out <- paths$step05c
rd <- function(f) fread(file.path(out, "tables", f))
rd4 <- function(f) fread(file.path(paths$step04b, "tables", f))
base <- 7; ink <- "#0b0b0b"; ink2 <- "#52514e"; rule <- "#c3c2b7"
up_col <- "#c0392b"; down_col <- "#2a6f97"
theme_fig <- function() theme_classic(base_size = base) +
  theme(text = element_text(colour = ink), axis.text = element_text(colour = ink2),
        axis.line = element_line(linewidth = 0.3, colour = rule), axis.ticks = element_line(linewidth = 0.3, colour = rule),
        plot.title = element_text(size = base + 1, face = "bold"), plot.subtitle = element_text(size = base - 0.5, colour = ink2),
        strip.background = element_blank(), strip.text = element_text(face = "bold"), legend.key.size = unit(3, "mm"))
save_fig <- function(p, name, w, h) {
  ggsave(file.path(out, paste0(name, ".png")), p, width = w, height = h, units = "mm", dpi = 300, bg = "white", device = ragg::agg_png)
  ggsave(file.path(out, paste0(name, ".pdf")), p, width = w, height = h, units = "mm", bg = "white")
}
save_ch <- function(draw_fun, name, w, h) {
  ragg::agg_png(file.path(out, paste0(name, ".png")), width = w, height = h, units = "mm", res = 300, background = "white")
  draw_fun(); dev.off()
  pdf(file.path(out, paste0(name, ".pdf")), width = w / 25.4, height = h / 25.4)
  draw_fun(); dev.off()
}
at <- rd("B_gene_atlas_matrix.csv")
grp_lv <- c("identity", "network", "outcome", "down marker", "optional", "reference")
at[, group := factor(group, grp_lv)]
genes <- at$gene
fs <- 6

## =========================================================================
## Fig 47: gene atlas heatmap
## =========================================================================
mk <- function(cols, labs) { m <- as.matrix(at[, ..cols]); dimnames(m) <- list(genes, labs); m }
l2 <- function(x) log2(x + 1)
m1 <- mk(c("lfc_C2", "lfc_C1", "lfc_TCGA", "lfc_TCGA_oral", "lfc_TCGA_larynx", "lfc_C2_HPVpos"),
         c("C2 single cell", "C1 paired", "TCGA all-site", "TCGA oral", "TCGA larynx", "HPV+ C2"))
f1 <- mk(c("fdr_C2", "fdr_C1", "fdr_TCGA", "fdr_TCGA_oral", "fdr_TCGA_larynx", "fdr_C2_HPVpos"), colnames(m1))
sc_cols <- c("sc_malignant", "sc_normal", paste0("sc_", c("CD4T", "CD8T", "NK", "Bcell", "PlasmaCell", "Macrophage", "Monocyte",
                                                         "panDC", "Neutrophils", "Mast", "Fibroblast", "Endothelial")))
m2 <- l2(mk(sc_cols, c("malignant", "normal epith.", "CD4 T", "CD8 T", "NK", "B", "plasma", "macrophage", "monocyte", "DC",
                       "neutrophil", "mast", "fibroblast", "endothelial")))
m3 <- l2(mk(c("TCGA_tumour_median_TPM", "TCGA_normal_median_TPM", "GTEx_esophagus_mucosa", "GTEx_minor_salivary_gland"),
            c("TCGA HPV- tumour", "TCGA normal", "GTEx esophagus mucosa", "GTEx minor salivary")))
m4 <- l2(mk(c("PBMC_HPVneg_max", "PBMC_HD_max", "GTEx_whole_blood"), c("PBMC HNSCC HPV- (max)", "PBMC healthy (max)", "GTEx whole blood")))
m5 <- log(mk(c("hr_TCGA_PFI", "hr_TCGA_oral_PFI", "hr_GSE41613_DSS"), c("TCGA PFI all-site", "TCGA oral PFI", "GSE41613 DSS")))
p5 <- mk(c("p_TCGA_PFI", "p_TCGA_oral_PFI", "p_GSE41613_DSS"), colnames(m5))
dot_fun <- function(sig) function(j, i, x, y, w, h, fill) if (isTRUE(sig[i, j])) grid.points(x, y, pch = 16, size = unit(1.1, "mm"))
na_col <- "grey88"
lim1 <- max(abs(m1), na.rm = TRUE)
col1 <- colorRamp2(c(-lim1, -2, 0, 2, lim1), c("#1d4e72", down_col, "white", up_col, "#7d1f16"))
col2 <- colorRamp2(c(0, 2, 5, 9), c("white", "#fde6b0", "#f08c3a", "#8c2d04"))
col3 <- colorRamp2(c(0, 3, 7, 12), c("white", "#d6e9c6", "#6aa84f", "#1e5631"))
col4 <- colorRamp2(c(0, 1, 3, 6), c("white", "#e3d7ee", "#9a7bc0", "#4b2a7b"))
col5 <- colorRamp2(c(-0.7, 0, 0.7), c(down_col, "white", up_col))
spec_col <- c(tumour_specific = "#c0392b", tumour_enriched = "#f0a17a", shared_with_TME = "#bdbdbd")
net_col <- c(`TF (E2F1/MYBL2/FOXM1)` = "#7a4b94", target = "#c9b3e0", none = "white")
bl_col <- c(`blood-low` = "#2e7d32", `not blood-low` = "white")
ra <- rowAnnotation(
  specificity = at$specificity, network = at$network, `blood-low` = fifelse(at$blood_low %in% TRUE, "blood-low", "not blood-low"),
  col = list(specificity = spec_col, network = net_col, `blood-low` = bl_col),
  gp = gpar(col = "grey70", lwd = 0.3), annotation_name_gp = gpar(fontsize = fs), annotation_name_rot = 60,
  annotation_legend_param = list(specificity = list(title = "tumour specificity", labels_gp = gpar(fontsize = fs), title_gp = gpar(fontsize = fs, fontface = "bold")),
                                 network = list(title = "network", labels_gp = gpar(fontsize = fs), title_gp = gpar(fontsize = fs, fontface = "bold")),
                                 `blood-low` = list(title = "blood background", labels_gp = gpar(fontsize = fs), title_gp = gpar(fontsize = fs, fontface = "bold"))))
hp <- function(m, col, title, name, cf = NULL, ...) Heatmap(m, col = col, name = name, na_col = na_col, cluster_rows = FALSE, cluster_columns = FALSE,
  row_split = at$group, row_title_gp = gpar(fontsize = fs, fontface = "bold"), row_title_rot = 0, row_names_gp = gpar(fontsize = fs),
  column_names_gp = gpar(fontsize = fs), column_names_rot = 50, column_title = title, column_title_gp = gpar(fontsize = fs + 0.5, fontface = "bold"),
  rect_gp = gpar(col = "white", lwd = 0.5), cell_fun = cf,
  heatmap_legend_param = list(title_gp = gpar(fontsize = fs, fontface = "bold"), labels_gp = gpar(fontsize = fs), legend_height = unit(18, "mm")), ...)
ht <- hp(m1, col1, "log2FC\ntumour vs normal", "log2FC", dot_fun(f1 < 0.05), width = unit(6 * 4.2, "mm")) +
  hp(m2, col2, "single cell (HPV-)\nlog2(median CPM + 1)", "sc log2CPM", width = unit(14 * 4.2, "mm")) +
  hp(m3, col3, "bulk\ntissue", "bulk log2TPM", width = unit(4 * 4.2, "mm")) +
  hp(m4, col4, "blood\n", "blood log2", width = unit(3 * 4.2, "mm")) +
  hp(m5, col5, "outcome\n", "log HR", dot_fun(p5 < 0.05), width = unit(3 * 4.2, "mm"),
     right_annotation = ra)
save_ch(function() { draw(ht, padding = unit(c(14, 2, 2, 2), "mm"), ht_gap = unit(2.5, "mm"), merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right",
                        column_title = "Gene atlas: 24 genes (HPV- main; reference genes as QC, log2FC near 0)",
                        column_title_gp = gpar(fontsize = fs + 2, fontface = "bold"))
                   grid.text(paste("log2FC: dot = FDR < 0.05; grey = not tested.  Bulk: log2(median TPM + 1).  Blood: log2(CPM or TPM + 1), PBMC = max over cell types.",
                                   "Outcome: log HR per SD (TCGA PFI multivariable, step 07a; TCGA oral PFI, step 07c; GSE41613 DSS, step 07b); dot = p < 0.05; grey = not tested.",
                                   "Specificity: step 04b part F rule. Network: TF or target of E2F1/MYBL2/FOXM1 (pySCENIC or CollecTRI). blood-low: practical pre-specified flag.", sep = "\n"),
                             x = unit(4, "mm"), y = unit(3, "mm"), just = c("left", "bottom"), gp = gpar(fontsize = fs - 0.5, col = ink2))
        },
        "Fig47_gene_atlas_heatmap", 250, 160)

## =========================================================================
## Fig 48: blood background panel
## =========================================================================
pc <- rd("A_cillo_pseudobulk_cpm.csv")
bl <- rd("A_blood_background.csv")
gord <- genes
pc[, `:=`(gene = factor(gene, rev(gord)), group = factor(group, c("HNSCC HPV-", "HNSCC HPV+", "healthy donor")),
          lab = fifelse(used_for_max, "", "x"))]
pc[, xlab := paste0(cell_type, " (", n_cells, ")")]
ct_ord <- c("CD4 T", "CD8 T", "NK", "B", "plasma", "CD14 monocyte", "CD16 monocyte", "cDC", "pDC", "platelet")
pc[, cell_type := factor(cell_type, ct_ord)]
setorder(pc, group, cell_type)
pc[, xlab := factor(paste(group, xlab, sep = "__"), unique(paste(group, xlab, sep = "__")))]
pa <- ggplot(pc, aes(xlab, gene, fill = log2(CPM + 1))) + geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = lab), size = 1.6, colour = ink2) +
  scale_x_discrete(labels = function(x) sub(".*__", "", x)) +
  facet_grid(. ~ group, scales = "free_x", space = "free_x") + scale_fill_gradientn(colours = c("white", "#e3d7ee", "#9a7bc0", "#4b2a7b"), values = c(0, 1 / 8, 3 / 8, 1), limits = c(0, 8),
                                               name = "log2(CPM + 1)") +
  labs(x = NULL, y = NULL, title = "Cillo 2020 (GSE139324) CD45+ PBMC: pseudobulk CPM per cell type",
       subtitle = "cell types labelled with standard markers (no author annotation in GEO); (n) = cells; x = < 50 cells (not used for the flag)") +
  theme_fig() + theme(axis.text.x = element_text(angle = 50, hjust = 1), axis.line = element_blank())
gx <- melt(bl[, .(gene, `GTEx whole blood` = GTEx_whole_blood, `GTEx esophagus mucosa` = GTEx_esophagus_mucosa,
                  `GTEx minor salivary` = GTEx_minor_salivary_gland)], id.vars = "gene", variable.name = "tissue", value.name = "TPM")
gx[, gene := factor(gene, rev(gord))]
pb <- ggplot(gx, aes(tissue, gene, fill = log2(TPM + 1))) + geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = fifelse(TPM < 1, "<1", sprintf("%.0f", TPM))), size = 1.6, colour = ink) +
  scale_fill_gradientn(colours = c("white", "#d6e9c6", "#6aa84f", "#1e5631"), limits = c(0, 16), name = "log2(TPM + 1)") +
  labs(x = NULL, y = NULL, title = "GTEx v8 median TPM") + theme_fig() + theme(axis.text.x = element_text(angle = 50, hjust = 1), axis.line = element_blank(),
                                                                               axis.text.y = element_blank())
fl <- bl[, .(gene = factor(gene, rev(gord)), flag = fifelse(blood_low, "blood-low", "not blood-low"),
             post = fifelse(blood_low_lib1M_posthoc, "yes", "no"), umi = sprintf("%.0f", PBMC_max_UMI))]
pf <- ggplot(fl, aes(x = 1, y = gene)) + geom_tile(aes(fill = flag), colour = "white") + geom_text(aes(x = 2, label = post), size = 1.8) +
  geom_text(aes(x = 3, label = umi), size = 1.8, colour = ink2) +
  scale_fill_manual(values = c(`blood-low` = "#2e7d32", `not blood-low` = "grey90"), name = NULL) +
  scale_x_continuous(breaks = 1:3, labels = c("flag", "post hoc", "UMIs"), limits = c(0.5, 3.5)) +
  labs(x = NULL, y = NULL, title = "Flag", subtitle = "post hoc: >= 1 M-count\npseudobulks only;\nUMIs at PBMC max") + theme_fig() + theme(axis.text.y = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
                                                                legend.position = "bottom", axis.text.x = element_text(angle = 50, hjust = 1))
f48 <- (pa | pb | pf) + plot_layout(widths = c(5, 1.2, 1.1)) +
  plot_annotation(caption = paste("blood-low (practical, pre-specified, not literature-based): max PBMC cell-type CPM < 1 AND GTEx whole blood TPM < 1.",
                                  "CD45+ PBMC show immune-cell background only; they cannot show whether tumour RNA or tumour cells are detectable in blood.",
                                  sep = "\n"), theme = theme(plot.caption = element_text(size = base - 1, colour = ink2, hjust = 0)))
save_fig(f48, "Fig48_blood_background", 240, 130)

## =========================================================================
## Fig 49: LFC concordance scatters
## =========================================================================
tt <- rd4("E_TCGA_HPVneg_paired_tumour_vs_normal.csv")[, .(gene, tcga = logFC)]
cc <- rd("C_lfc_concordance.csv")
sc_plot <- function(f, nm) {
  d <- merge(rd4(f)[, .(gene, sc = logFC)], tt, by = "gene")
  d[, hl := gene %in% genes]
  d[hl == TRUE, grp := at$group[match(gene, at$gene)]]
  s <- cc[single_cell == nm]
  ggplot(d, aes(sc, tcga)) + geom_hline(yintercept = 0, colour = rule, linewidth = 0.3) + geom_vline(xintercept = 0, colour = rule, linewidth = 0.3) +
    geom_point(data = d[hl == FALSE], colour = "grey75", size = 0.3, alpha = 0.5) +
    geom_point(data = d[hl == TRUE], aes(colour = grp), size = 1.6) +
    ggrepel::geom_text_repel(data = d[hl == TRUE], aes(label = gene, colour = grp), size = 2, max.overlaps = 50, segment.size = 0.2,
                             min.segment.length = 0, show.legend = FALSE) +
    scale_colour_manual(values = c(identity = "#c0392b", network = "#7a4b94", outcome = "#e69f00", `down marker` = "#2a6f97",
                                   optional = "#2e7d32", reference = "#555555"), limits = grp_lv, drop = FALSE, name = NULL) +
    labs(x = sprintf("%s log2FC (single cell, HPV-)", nm), y = "TCGA paired log2FC (40 HPV- pairs)",
         title = sprintf("%s vs TCGA", nm),
         subtitle = sprintf("all %d genes tested in both: Pearson r %.2f, Spearman rho %.2f\natlas genes tested in both (n = %d): r %.2f, rho %.2f",
                            s$n, s$pearson, s$spearman, s$n_24, s$pearson_24, s$spearman_24)) +
    theme_fig()
}
f49 <- (sc_plot("D_C2_HPVneg_allpatients.csv", "C2") | (sc_plot("D_C1_HPVneg_paired.csv", "C1") + guides(colour = "none"))) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
save_fig(f49, "Fig49_lfc_concordance", 200, 110)

## =========================================================================
## Fig 50: co-expression heatmaps
## =========================================================================
co <- rd("C_coexpression.csv")
net_block <- c("E2F1", "MYBL2", "FOXM1", "KIF23", "HMGB3", "HMMR")
titles <- c(TCGA_all = "TCGA HPV- tumours", TCGA_oral = "TCGA HPV- oral cavity", sc_malignant = "single-cell malignant pseudobulks")
colr <- colorRamp2(c(-1, 0, 1), c(down_col, "white", up_col))
hts <- lapply(names(titles), function(k) {
  d <- co[dataset == k]
  g <- sort(unique(c(d$gene1, d$gene2)))
  R <- diag(1, length(g)); P <- matrix(0, length(g), length(g)); dimnames(R) <- dimnames(P) <- list(g, g)
  R[cbind(d$gene1, d$gene2)] <- R[cbind(d$gene2, d$gene1)] <- d$rho
  P[cbind(d$gene1, d$gene2)] <- P[cbind(d$gene2, d$gene1)] <- d$p
  hc <- hclust(as.dist(1 - R), method = "average")
  nc <- fifelse(g %in% net_block, "#7a4b94", ink)
  Heatmap(R, col = colr, name = paste0("rho_", k), cluster_rows = hc, cluster_columns = hc, show_heatmap_legend = FALSE,
          heatmap_legend_param = list(title = "Spearman rho", title_gp = gpar(fontsize = fs, fontface = "bold"), labels_gp = gpar(fontsize = fs)),
          row_names_gp = gpar(fontsize = fs, col = nc, fontface = fifelse(g %in% net_block, "bold", "plain")),
          column_names_gp = gpar(fontsize = fs, col = nc, fontface = fifelse(g %in% net_block, "bold", "plain")),
          row_dend_width = unit(6, "mm"), column_dend_height = unit(6, "mm"), rect_gp = gpar(col = "white", lwd = 0.4),
          column_title = sprintf("%s (n = %d)", titles[[k]], d$n[1]), column_title_gp = gpar(fontsize = fs + 1, fontface = "bold"),
          cell_fun = function(j, i, x, y, w, h, fill) {
            if (i != j && P[i, j] < 0.05) grid.points(x, y, pch = 16, size = unit(0.7, "mm"), gp = gpar(col = ink))
            if (rownames(R)[i] %in% net_block && colnames(R)[j] %in% net_block) grid.rect(x, y, w, h, gp = gpar(col = "#7a4b94", fill = NA, lwd = 0.8))
          })
})
draw50 <- function() {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 4, heights = unit(c(14, 1), c("mm", "null")), widths = unit(c(1, 1, 1, 22), c("null", "null", "null", "mm")))))
  grid.text("Co-expression of the 19 target genes (Spearman; dot = nominal p < 0.05; purple boxes = network block E2F1, MYBL2, FOXM1, KIF23, HMGB3, HMMR)",
            y = unit(0.65, "npc"), vp = viewport(layout.pos.row = 1, layout.pos.col = 1:4), gp = gpar(fontsize = fs + 1.5, fontface = "bold"))
  grid.text("single-cell n = 27 pseudobulks is small; clustering on 1 - rho, average linkage, separate per panel",
            y = unit(0.2, "npc"), vp = viewport(layout.pos.row = 1, layout.pos.col = 1:4), gp = gpar(fontsize = fs, col = ink2))
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 4))
  draw(Legend(col_fun = colr, title = "Spearman\nrho", title_gp = gpar(fontsize = fs, fontface = "bold"), labels_gp = gpar(fontsize = fs),
              legend_height = unit(30, "mm")), x = unit(2, "mm"), just = "left")
  upViewport()
  for (i in 1:3) {
    pushViewport(viewport(layout.pos.row = 2, layout.pos.col = i))
    draw(hts[[i]], newpage = FALSE)
    upViewport()
  }
}
save_ch(draw50, "Fig50_coexpression", 300, 115)
