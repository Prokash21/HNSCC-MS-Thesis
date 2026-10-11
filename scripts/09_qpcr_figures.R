## Figures for the qRT-PCR panel screen (scripts/08_qpcr_panel.R).
##   Fig7: where each named candidate is expressed across atlas cell types, and
##         the screen of all HPV+ candidates (tumour share of the whole-tissue
##         signal vs expression level)
##   Fig8: what to expect in tissue qPCR (whole-tissue expression by HPV for
##         the panel and the TME markers) and the reference-gene pre-screen
##
## Run from the repo root:  Rscript scripts/09_qpcr_figures.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

res     <- readRDS("data/processed/HNSCatlas/qpcr_screen.rds")
fig_dir <- "HNSCC-HPV/results/04_qPCR_panel"
th      <- res$thresholds

## ---- shared style (as in scripts/07) --------------------------------------
base_size <- 7
ink       <- "#0b0b0b"
ink2      <- "#52514e"
rule      <- "#c3c2b7"
grid_col  <- "#e1e0d9"
hpv_cols  <- c("HPV+" = "#2a78d6", "HPV-" = "#eb6834")
nov_cols  <- c(Established = "#898781", Partial = "#c98500", Open = "#1baf7a",
               "Not assessed" = "#d6d5cf")

theme_fig <- function() {
  theme_classic(base_size = base_size) +
    theme(text = element_text(colour = ink),
          axis.text = element_text(colour = ink2),
          axis.line = element_line(linewidth = 0.3, colour = rule),
          axis.ticks = element_line(linewidth = 0.3, colour = rule),
          plot.title = element_text(size = base_size + 1, face = "bold"),
          plot.subtitle = element_text(size = base_size - 0.5, colour = ink2),
          plot.tag = element_text(size = base_size + 3, face = "bold"),
          strip.background = element_blank(),
          strip.text = element_text(size = base_size - 0.5, face = "bold"),
          legend.text = element_text(size = base_size - 0.5),
          legend.title = element_text(size = base_size - 0.5),
          legend.key.size = unit(3, "mm"),
          plot.margin = margin(4, 6, 4, 4))
}
save_fig <- function(p, name, w, h) {
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h,
         units = "mm", dpi = 600, device = ragg::agg_png, bg = "white")
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h,
         units = "mm", bg = "white")
  message("Wrote ", file.path(fig_dir, name), ".png/.pdf")
}
p_fmt <- function(p) if (p < 0.001) formatC(p, format = "e", digits = 0) else sprintf("%.3f", p)

scr   <- res$scr
panel <- res$panel[!is.na(res$panel$gene), ]

## ---- Figure 7a: expression across atlas cell types ------------------------
groups <- c(Tumour = "Tumour", "Normal squamous" = "Normal squamous",
            "Epithelial, other" = "Other epithelial", CD4T = "CD4+ T", CD8T = "CD8+ T",
            NK = "NK", Bcell = "B cell", PlasmaCell = "Plasma cell",
            Macrophage = "Macrophage", Monocyte = "Monocyte", panDC = "Dendritic",
            Neutrophils = "Neutrophil", Mast = "Mast", Fibroblast = "Fibroblast",
            Endothelial = "Endothelial")
named  <- setdiff(intersect(names(res$novelty), scr$gene), panel$gene)
dot_genes <- c(panel$gene, sort(named), res$markers)
dot_set   <- c(rep("Suggested panel", nrow(panel)), rep("Other named candidates", length(named)),
               rep("TME markers", length(res$markers)))

## median CPM per cell group across patients with enough cells of that group
med_cpm <- function(group) {
  min_n <- if (group == "Normal squamous") th$min_normal else th$min_cells
  s <- res$smp$Group == group & res$smp$n_cells >= min_n
  apply(res$cpm_g[dot_genes, s, drop = FALSE], 1, median)
}
cpm_grp <- sapply(names(groups), med_cpm)
pct_grp <- 100 * sweep(res$detected[dot_genes, names(groups)], 2, res$n_group[names(groups)], "/")
z <- t(scale(t(log2(cpm_grp + 1))))
dd <- data.frame(gene  = rep(dot_genes, times = length(groups)),
                 set   = rep(dot_set, times = length(groups)),
                 group = rep(names(groups), each = length(dot_genes)),
                 z     = as.vector(z), pct = as.vector(pct_grp))
dd$gene  <- factor(dd$gene, levels = rev(dot_genes))
dd$group <- factor(dd$group, levels = names(groups))
dd$set   <- factor(dd$set, levels = c("Suggested panel", "Other named candidates", "TME markers"))

p7a <- ggplot(dd, aes(group, gene)) +
  geom_point(aes(size = pct, colour = z)) +
  facet_grid(set ~ ., scales = "free_y", space = "free_y") +
  scale_x_discrete(labels = groups) +
  scale_colour_gradient(low = "#e1e0d9", high = "#184f95", name = "Relative\nexpression") +
  scale_size_area(max_size = 3.2, name = "% cells\nexpressing") +
  labs(title = "Where the candidates are expressed",
       subtitle = "atlas cell types; colour = expression relative to each gene's own range",
       x = NULL, y = NULL) +
  theme_fig() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(colour = ink),
        panel.grid.major = element_line(colour = grid_col, linewidth = 0.2),
        strip.text.y = element_text(angle = 0, hjust = 0))

## ---- Figure 7b: the screen of all HPV+ candidates ----------------------
scr$share_pct <- 100 * scr$tissue_tumour_share
scr$status <- ifelse(scr$passes_all, "Passes all filters", "Fails at least one")
lab <- scr[scr$gene %in% c(names(res$novelty), panel$gene), ]
lab$face <- ifelse(lab$gene %in% panel$gene, "bold", "plain")
p7b <- ggplot(scr, aes(share_pct, tumour_HPVpos_CPM)) +
  annotate("rect", xmin = 100 * th$min_share, xmax = Inf, ymin = th$min_cpm, ymax = Inf,
           fill = "#f0efec") +
  geom_vline(xintercept = 100 * th$min_share, linetype = "22", linewidth = 0.25, colour = ink2) +
  geom_hline(yintercept = th$min_cpm, linetype = "22", linewidth = 0.25, colour = ink2) +
  geom_point(aes(colour = novelty, shape = status), size = 1.8, stroke = 0.5) +
  geom_text_repel(data = lab, aes(label = gene, fontface = face), size = (base_size - 1.5) / .pt,
                  colour = ink, min.segment.length = 0, segment.size = 0.2,
                  segment.colour = ink2, box.padding = 0.4, force = 3, max.time = 3,
                  max.iter = 1e5, max.overlaps = Inf, bg.color = "white", bg.r = 0.1,
                  seed = 1) +
  scale_y_log10(labels = scales::label_comma()) +
  scale_colour_manual(values = nov_cols, name = "Novelty\n(strategy note)") +
  scale_shape_manual(values = c("Passes all filters" = 16, "Fails at least one" = 1), name = NULL) +
  labs(title = sprintf("qPCR screen of the %d HPV+ candidates", nrow(scr)),
       subtitle = sprintf(paste0("shaded: >= %d CPM and >= %d%% of tissue signal from tumour\n",
                                 "filled: passes all filters; bold: suggested panel"),
                          th$min_cpm, round(100 * th$min_share)),
       x = "Share of whole-tissue signal from tumour cells (%)",
       y = "Expression in HPV+ tumour cells (CPM, log scale)") +
  theme_fig()

fig7 <- (p7a | p7b) + plot_layout(widths = c(0.95, 1.05)) + plot_annotation(tag_levels = "a")
save_fig(fig7, "Fig7_qPCR_candidate_screen", w = 180, h = 140)

## ---- Figure 8a: whole-tissue expression by HPV (what tissue qPCR may show) -
wt_genes <- c(panel$gene, res$markers)
lc <- res$lc_wt[wt_genes, , drop = FALSE]
stat_of <- function(g) {
  r <- res$wt[res$wt$gene == g, ]
  if (nrow(r) == 0) return(g)
  sprintf("%s\nlog2FC %.1f\nFDR %s", g, r$logFC, p_fmt(r$FDR))
}
df <- data.frame(gene = rep(wt_genes, times = ncol(lc)),
                 HPV  = factor(rep(res$wt_meta$HPV, each = length(wt_genes)), levels = names(hpv_cols)),
                 y    = as.vector(lc))
df$facet <- factor(vapply(df$gene, stat_of, character(1)),
                   levels = vapply(wt_genes, stat_of, character(1)))
p8a <- ggplot(df, aes(HPV, y)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, colour = ink2, linewidth = 0.3, fill = NA) +
  geom_point(aes(colour = HPV), position = position_jitter(width = 0.15, seed = 1),
             size = 0.9, stroke = 0, shape = 16) +
  facet_wrap(~ facet, nrow = 2, scales = "free_y") +
  scale_colour_manual(values = hpv_cols, guide = "none") +
  labs(title = "What tissue qPCR may show: HPV+ vs HPV- tumours",
       subtitle = sprintf(paste0("all cells of each sample summed ('virtual bulk'); %d patients from ",
                                 "the unsorted studies; study-adjusted edgeR"), nrow(res$wt_meta)),
       x = NULL, y = "log2 CPM (whole sample)") +
  theme_fig()

## ---- Figure 8b: reference-gene pre-screen ---------------------------------
rs <- res$rs
rs$gene  <- factor(rs$gene, levels = rev(rs$gene))
rs$class <- ifelse(rs$stability_rank <= 3, "Top 3 (test these)",
                   ifelse(rs$gene %in% c("GAPDH", "ACTB"), "Common choice", "Other"))
p8b <- ggplot(rs, aes(sd_wholetissue, gene)) +
  geom_segment(aes(x = 0, xend = sd_wholetissue, yend = gene), colour = grid_col, linewidth = 0.5) +
  geom_point(aes(colour = class), size = 1.8) +
  scale_colour_manual(values = c("Top 3 (test these)" = "#184f95", "Common choice" = "#eb6834",
                                 Other = "#898781"), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(title = "Reference genes",
       subtitle = "ranked by overall stability;\nconfirm with geNorm/NormFinder",
       x = "SD of log2 CPM across samples", y = NULL) +
  theme_fig() +
  theme(axis.text.y = element_text(colour = ink), legend.position = "bottom",
        legend.direction = "vertical")

fig8 <- (p8a | p8b) + plot_layout(widths = c(2.3, 1)) + plot_annotation(tag_levels = "a")
save_fig(fig8, "Fig8_qPCR_panel_expectations", w = 180, h = 110)
