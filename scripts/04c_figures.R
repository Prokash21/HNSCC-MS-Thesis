## Step 04c (part C): figures
##
## Purpose: Fig51 log2FC scatter edgeR vs DESeq2 (C2, C1, TCGA; HPV-); Fig52
##   overlap of significant genes (FDR < 0.05) by method, up and down; Fig53
##   the 24 atlas genes, edgeR vs DESeq2 side by side.
## Inputs:  results/step04c_deseq2/tables/{A_deseq2_vs_edgeR_*, B1_*, B2_*, B4_atlas_genes}.csv
## Outputs: results/step04c_deseq2/Fig51-Fig53 (png, pdf)
## Date: 2026-10-28

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})
source("config/paths.R")

out <- paths$step04c
rd <- function(f) fread(file.path(out, "tables", f))
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
labs_k <- c(C2_HPVneg = "C2 single cell (~ study + state)", C1_HPVneg = "C1 paired (~ patient + state)",
            TCGA_HPVneg = "TCGA paired, 40 pairs (~ patient + state)")
atlas <- unique(rd("B4_atlas_genes.csv")$gene)
b1 <- rd("B1_lfc_agreement.csv")

## Fig 51: log2FC scatter
f51 <- wrap_plots(lapply(names(labs_k), function(k) {
  d <- rd(paste0("A_deseq2_vs_edgeR_", k, ".csv"))
  s <- b1[contrast == k & genes == "all tested"]
  ggplot(d, aes(edgeR_logFC, log2FC)) +
    geom_abline(slope = 1, intercept = 0, colour = rule, linewidth = 0.3) +
    geom_point(colour = "grey70", size = 0.25, alpha = 0.5) +
    geom_point(data = d[gene %in% atlas], colour = up_col, size = 1.1) +
    ggrepel::geom_text_repel(data = d[gene %in% atlas], aes(label = gene), size = 1.6, colour = up_col, max.overlaps = 30,
                             segment.size = 0.15, min.segment.length = 0) +
    labs(x = "edgeR QL log2FC", y = "DESeq2 log2FC (unshrunken)", title = labs_k[[k]],
         subtitle = sprintf("%d genes: Pearson %.3f, Spearman %.3f", s$n, s$pearson, s$spearman)) +
    theme_fig()
}), nrow = 1) + plot_annotation(title = "edgeR QL vs DESeq2 log2FC, HPV- tumour vs normal",
                                subtitle = "same genes, counts and designs as step 04b; red = 24 atlas genes; grey line = identity",
                                theme = theme(plot.title = element_text(size = base + 1, face = "bold"),
                                              plot.subtitle = element_text(size = base - 0.5, colour = ink2)))
save_fig(f51, "Fig51_lfc_edgeR_vs_DESeq2", 240, 90)

## Fig 52: overlap of significant genes (FDR < 0.05), up and down
b2 <- rd("B2_significance_agreement.csv")
b2[, contrast := factor(contrast, c("C2_HPVneg", "C1_HPVneg", "TCGA_HPVneg", "C2_HPVpos", "C1_HPVpos"))]
circ <- function(x0, r = 1, n = 200) data.table(x = x0 + r * cos(seq(0, 2 * pi, length.out = n)), y = r * sin(seq(0, 2 * pi, length.out = n)))
venn_one <- function(k, d) {
  v <- b2[contrast == k & direction == d]; col <- if (d == "up") up_col else down_col
  ggplot() +
    geom_polygon(data = circ(-0.6), aes(x, y), fill = col, alpha = 0.18, colour = col, linewidth = 0.3) +
    geom_polygon(data = circ(0.6), aes(x, y), fill = "grey60", alpha = 0.15, colour = "grey40", linewidth = 0.3) +
    annotate("text", x = c(-1.05, 0, 1.05), y = 0, label = c(v$edgeR_only, v$both, v$DESeq2_only), size = 2.6) +
    annotate("text", x = c(-0.9, 0.9), y = 1.25, label = c(sprintf("edgeR (%d)", v$edgeR_total), sprintf("DESeq2 (%d)", v$DESeq2_total)),
             size = 2, colour = ink2) +
    annotate("text", x = 0, y = -1.3, label = sprintf("Jaccard %.2f", v$jaccard), size = 2, colour = ink2) +
    coord_equal(xlim = c(-2, 2), ylim = c(-1.45, 1.45)) + labs(title = sprintf("%s, %s", sub("_", " ", k), d)) +
    theme_void(base_size = base) + theme(plot.title = element_text(face = "bold", hjust = 0.5, size = base))
}
ks <- levels(b2$contrast)
f52 <- wrap_plots(c(lapply(ks, venn_one, d = "up"), lapply(ks, venn_one, d = "down")), nrow = 2) +
  plot_annotation(title = "Genes with FDR < 0.05: edgeR QL vs DESeq2 Wald",
                  subtitle = "no fold-change threshold; direction by each method's own log2FC; HPV+ is the secondary arm",
                  theme = theme(plot.title = element_text(size = base + 1, face = "bold"),
                                plot.subtitle = element_text(size = base - 0.5, colour = ink2)))
save_fig(f52, "Fig52_significant_overlap", 240, 105)

## Fig 53: the 24 atlas genes side by side
b4 <- rd("B4_atlas_genes.csv")
gl <- fread(file.path(paths$step08, "tables", "A_final_gene_list.csv"))
b4[, `:=`(gene = factor(gene, rev(gl$gene)), contrast = factor(contrast, names(labs_k), c("C2", "C1", "TCGA")),
          group = factor(group, unique(gl$group)))]
lg <- melt(b4, id.vars = c("gene", "group", "contrast", "flag"), measure.vars = list(c("edgeR_logFC", "DESeq2_log2FC"), c("edgeR_FDR", "DESeq2_padj")),
           variable.name = "method", value.name = c("log2FC", "FDR"))
lg[, method := factor(c("edgeR QL", "DESeq2")[method], c("edgeR QL", "DESeq2"))]
lg[, sig := !is.na(FDR) & FDR < 0.05]
f53 <- ggplot(lg, aes(log2FC, gene)) +
  geom_vline(xintercept = 0, colour = rule, linewidth = 0.3) +
  geom_line(aes(group = gene), colour = "grey70", linewidth = 0.3) +
  geom_point(aes(colour = method, shape = sig), size = 1.5, stroke = 0.5) +
  geom_text(data = unique(lg[flag != "", .(gene, group, contrast, flag)]), aes(x = Inf, y = gene, label = "*"), hjust = 1.2, size = 3,
            colour = ink, inherit.aes = FALSE) +
  facet_grid(group ~ contrast, scales = "free", space = "free_y") +
  scale_colour_manual(values = c(`edgeR QL` = "#2a6f97", DESeq2 = "#e69f00"), name = NULL) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "FDR < 0.05", `FALSE` = "FDR >= 0.05"), name = NULL) +
  labs(x = "log2FC tumour vs normal (edgeR QL; DESeq2 unshrunken MLE)", y = NULL, title = "The 24 atlas genes: edgeR QL vs DESeq2",
       subtitle = "* = direction or significance (FDR < 0.05) differs between methods; C1 not tested for CASC9, IGF2BP3, CA9 (filtered in step 04b)") +
  theme_fig() + theme(legend.position = "bottom", strip.text.y = element_text(angle = 0, hjust = 0))
save_fig(f53, "Fig53_atlas_genes_methods", 200, 170)
