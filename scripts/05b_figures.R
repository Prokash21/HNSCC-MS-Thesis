## Step 05b (part G): figures
##
## Purpose: Hallmark GSEA dot plot (C2, C1, TCGA; HPV+ secondary); Venn of HPV-
##   up and down genes (single-cell C2 candidates vs TCGA paired FDR < 0.05);
##   PROGENy heatmap; network diagram (E2F1/MYBL2/FOXM1/NFIC -> eligible and
##   qPCR genes); story diagram (pathways -> TFs -> panel genes).
##   The GO cluster heatmap (Fig43) is drawn by 05b_ora.R.
## Inputs:  results/step05b_functional/tables/ (B_*, D_*, E_*, F_*), results/step04b_DE/tables/{D_candidates_HPVneg.csv,
##          E_TCGA_HPVneg_paired_tumour_vs_normal.csv}
## Outputs: results/step05b_functional/Fig41, Fig42, Fig44, Fig45, Fig46 (png, pdf)
## Date: 2026-10-28

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})
source("config/paths.R")

out <- paths$step05b
rd <- function(f) fread(file.path(out, "tables", f))
fdr_cut <- 0.05
base <- 7; ink <- "#0b0b0b"; ink2 <- "#52514e"; rule <- "#c3c2b7"
up_col <- "#c0392b"; down_col <- "#2a6f97"
theme_fig <- function() theme_classic(base_size = base) +
  theme(text = element_text(colour = ink), axis.text = element_text(colour = ink2),
        axis.line = element_line(linewidth = 0.3, colour = rule), axis.ticks = element_line(linewidth = 0.3, colour = rule),
        plot.title = element_text(size = base + 1, face = "bold"), plot.subtitle = element_text(size = base - 0.5, colour = ink2),
        strip.background = element_blank(), strip.text = element_text(face = "bold"), legend.key.size = unit(3, "mm"))
save_fig <- function(p, name, w, h) {
  ggsave(file.path(out, paste0(name, ".png")), p, width = w, height = h, units = "mm", dpi = 300, bg = "white",
         device = ragg::agg_png)
  ggsave(file.path(out, paste0(name, ".pdf")), p, width = w, height = h, units = "mm", bg = "white")
}
pretty_set <- function(s) gsub("_", " ", sub("^HALLMARK_", "", s))

## =========================================================================
## Fig 41: Hallmark GSEA dot plot
## =========================================================================
g <- rd("B_gsea_all.csv")[collection == "H"]
w <- rd("B_hpvneg_wide.csv")[collection == "H"]
keep <- w[FDR_C2_HPVneg < fdr_cut | FDR_TCGA_HPVneg < fdr_cut, set]
ord <- w[set %in% keep][order(-replicated, NES_C2_HPVneg), set]
g <- g[set %in% keep]
g[, `:=`(set_lab = factor(paste0(pretty_set(set), ifelse(set %in% w[replicated == TRUE, set], "  *", "")),
                          paste0(pretty_set(ord), ifelse(ord %in% w[replicated == TRUE, set], "  *", ""))),
         contrast = factor(contrast, c("C2_HPVneg", "C1_HPVneg", "TCGA_HPVneg", "C2_HPVpos", "C1_HPVpos"),
                           c("C2\nsingle cell", "C1\npaired", "TCGA\n40 pairs", "C2\nsingle cell", "C1\npaired")),
         arm = factor(fifelse(grepl("HPVneg", contrast), "HPV- (main)", "HPV+ (secondary)"), c("HPV- (main)", "HPV+ (secondary)")),
         sig = FDR < fdr_cut, mlog = pmin(-log10(FDR), 10))]
f41 <- ggplot(g, aes(contrast, set_lab)) +
  geom_point(aes(size = mlog, fill = NES, colour = sig), shape = 21, stroke = 0.4) +
  scale_fill_gradient2(low = down_col, mid = "white", high = up_col, midpoint = 0, name = "NES") +
  scale_colour_manual(values = c(`TRUE` = ink, `FALSE` = rule), labels = c(`TRUE` = "FDR < 0.05", `FALSE` = "FDR >= 0.05"), name = NULL) +
  scale_size_continuous(range = c(0.6, 4), name = "-log10 FDR\n(capped at 10)") +
  facet_grid(. ~ arm, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = NULL, title = "Hallmark GSEA: tumour vs normal",
       subtitle = paste("fgsea on sign(logFC) x sqrt(F); sets with FDR < 0.05 in HPV- C2 or TCGA;",
                        "\n* replicated = same direction, FDR < 0.05 in C2 and TCGA")) +
  theme_fig() + theme(legend.position = "right")
save_fig(f41, "Fig41_hallmark_gsea", 175, 140)

## =========================================================================
## Fig 42: Venn, HPV- up and down: single-cell C2 candidates vs TCGA paired FDR < 0.05
## =========================================================================
cand <- fread(file.path(paths$step04b, "tables", "D_candidates_HPVneg.csv"))[candidate_list == TRUE]
tcga <- fread(file.path(paths$step04b, "tables", "E_TCGA_HPVneg_paired_tumour_vs_normal.csv"))
venn <- rbindlist(lapply(c("up", "down"), function(d) {
  a <- cand[candidate == d, gene]
  b <- tcga[FDR < fdr_cut & sign(logFC) == fifelse(d == "up", 1, -1), gene]
  data.table(direction = d, sc_only = length(setdiff(a, b)), both = length(intersect(a, b)), tcga_only = length(setdiff(b, a)),
             sc_not_in_tcga = sum(!a %in% tcga$gene), sc_total = length(a), tcga_total = length(b))
}))
fwrite(venn, file.path(out, "tables", "G_venn_counts.csv"))
circ <- function(x0, r = 1, n = 200) data.table(x = x0 + r * cos(seq(0, 2 * pi, length.out = n)), y = r * sin(seq(0, 2 * pi, length.out = n)))
venn_plot <- function(d) {
  v <- venn[direction == d]; col <- if (d == "up") up_col else down_col
  ggplot() +
    geom_polygon(data = circ(-0.6), aes(x, y), fill = col, alpha = 0.18, colour = col, linewidth = 0.3) +
    geom_polygon(data = circ(0.6), aes(x, y), fill = "grey60", alpha = 0.15, colour = "grey40", linewidth = 0.3) +
    annotate("text", x = c(-1.05, 0, 1.05), y = 0, label = c(v$sc_only, v$both, v$tcga_only), size = 3, colour = ink) +
    annotate("text", x = c(-1.15, 1.15), y = 1.25, size = 2.1, colour = ink2,
             label = c(sprintf("single-cell C2 candidates\n%s (n = %d)", d, v$sc_total), sprintf("TCGA paired FDR < 0.05\n%s (n = %d)", d, v$tcga_total))) +
    annotate("text", x = 0, y = -1.25, size = 2, colour = ink2,
             label = sprintf("%d of %d candidates not tested in TCGA; %.0f%% of candidates in TCGA set", v$sc_not_in_tcga, v$sc_total,
                             100 * v$both / v$sc_total)) +
    coord_equal(xlim = c(-2, 2), ylim = c(-1.4, 1.55)) + labs(title = paste("HPV-", d, "in tumour")) +
    theme_void(base_size = base) + theme(plot.title = element_text(face = "bold", hjust = 0.5))
}
f42 <- (venn_plot("up") | venn_plot("down")) +
  plot_annotation(title = "HPV- tumour vs normal: single-cell candidates and TCGA paired DE",
                  subtitle = "candidates = step 04b C2 significant with C1 in the same direction; TCGA = 40 HPV- pairs, FDR < 0.05, same direction (no fold-change cut-off)",
                  theme = theme(plot.title = element_text(size = base + 1, face = "bold"),
                                plot.subtitle = element_text(size = base - 1, colour = ink2)))
save_fig(f42, "Fig42_venn_sc_vs_tcga", 175, 80)

## =========================================================================
## Fig 44: PROGENy heatmap (effect = activity difference, malignant/tumour - normal)
## =========================================================================
ps <- rd("D_progeny_single_cell.csv")
tr <- rd("D_progeny_hpvneg_summary.csv")
hm <- rbind(ps[arm == "HPV-", .(pathway, test = "HPV- single cell\nlm ~ study + state", effect = all_effect, FDR = all_FDR)],
            ps[arm == "HPV-", .(pathway, test = "HPV- single cell\npaired (7)", effect = paired_diff, FDR = paired_FDR)],
            tr[, .(pathway, test = "HPV- TCGA\npaired (40)", effect = tcga_paired_diff, FDR = tcga_paired_FDR)],
            ps[arm == "HPV+", .(pathway, test = "HPV+ single cell\nlm ~ study + state", effect = all_effect, FDR = all_FDR)],
            ps[arm == "HPV+", .(pathway, test = "HPV+ single cell\npaired (8)", effect = paired_diff, FDR = paired_FDR)])
hm[, `:=`(test = factor(test, unique(test)), star = fcase(FDR < 0.001, "***", FDR < 0.01, "**", FDR < 0.05, "*", default = ""),
          pathway = factor(pathway, tr[order(tcga_paired_diff), pathway]))]
lim <- max(abs(hm$effect))
f44 <- ggplot(hm, aes(test, pathway, fill = effect)) + geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = star), size = 2.6, colour = ink, vjust = 0.75) +
  scale_fill_gradient2(low = down_col, mid = "white", high = up_col, limits = c(-lim, lim), name = "activity\ndifference") +
  labs(x = NULL, y = NULL, title = "PROGENy pathway activity: tumour vs normal",
       subtitle = "decoupleR multivariate linear model, top 500 genes per pathway; BH FDR per test: * < 0.05, ** < 0.01, *** < 0.001") +
  theme_fig() + theme(axis.line = element_blank(), axis.ticks = element_blank())
save_fig(f44, "Fig44_progeny_heatmap", 160, 95)

## =========================================================================
## Fig 45: network diagram, core TFs -> eligible and qPCR genes
## =========================================================================
core <- c("E2F1", "MYBL2", "FOXM1", "NFIC")
e1 <- rd("E_core_tf_edges.csv")
mem <- rd("E_gene_target_membership.csv")
cols <- grep("_(pySCENIC|CollecTRI)$", names(mem), value = TRUE)
ed <- melt(mem[, c("gene", "group", cols), with = FALSE], id.vars = c("gene", "group"))[value == TRUE]
ed[, `:=`(TF = sub("_.*", "", variable), source = sub(".*_", "", variable))]
ed <- ed[, .(source = paste(sort(unique(source)), collapse = "+")), by = .(TF, gene, group)]
## edges to genes that are themselves core TFs (FOXM1, NFIC) are drawn as TF-TF edges
ed <- ed[gene != TF & !gene %in% core]
nsz <- rd("E_target_set_overlap.csv")[set == "HALLMARK_G2M_CHECKPOINT", .(TF, source, targets)]
tf_lab <- sapply(core, function(t) {
  x <- nsz[TF == t]
  if (!nrow(x)) return(t)
  sprintf("%s\n(%s)", t, paste(sprintf("%s %d", x$source, x$targets), collapse = ", "))
})
tgt <- unique(ed[, .(gene, group)])
tgt[, grp := fifelse(group == "07a eligible", "07a eligible gene", "qPCR panel gene")]
setorder(tgt, grp, gene)
nodes <- rbind(data.table(id = core, x = seq(1, by = 2.2, length.out = 4), y = 2, type = "core TF", label = tf_lab[core]),
               data.table(id = tgt$gene, x = seq(0, 7.6, length.out = nrow(tgt)), y = 0, type = tgt$grp, label = tgt$gene))
ed2 <- merge(ed, nodes[, .(gene = id, xe = x, ye = y)], by = "gene")
ed2 <- merge(ed2, nodes[type == "core TF", .(TF = id, xs = x, ys = y)], by = "TF")
tt <- merge(e1, nodes[, .(TF = id, xs = x, ys = y)], by = "TF")
tt <- merge(tt, nodes[, .(target = id, xe = x, ye = y)], by = "target")
src_col <- c(pySCENIC = "#7a4b94", CollecTRI = "#e69f00", `CollecTRI+pySCENIC` = ink)
f45 <- ggplot() +
  geom_segment(data = ed2, aes(x = xs, y = ys - 0.18, xend = xe, yend = ye + 0.12, colour = source), linewidth = 0.35, alpha = 0.8) +
  geom_curve(data = tt[xe > xs], aes(x = xs, y = ys + 0.25, xend = xe, yend = ye + 0.25, colour = sources, linetype = factor(collectri_sign)),
             curvature = -0.35, linewidth = 0.45, arrow = arrow(length = unit(1.4, "mm"), type = "closed")) +
  geom_curve(data = tt[xe < xs], aes(x = xs, y = ys + 0.25, xend = xe, yend = ye + 0.25, colour = sources, linetype = factor(collectri_sign)),
             curvature = 0.35, linewidth = 0.45, arrow = arrow(length = unit(1.4, "mm"), type = "closed")) +
  geom_label(data = nodes[type == "core TF"], aes(x, y, label = label), size = 2.1, fill = "#f3efe6", label.size = 0.2, lineheight = 0.9) +
  geom_point(data = nodes[type != "core TF"], aes(x, y, fill = type), shape = 21, size = 3, stroke = 0.3) +
  geom_text(data = nodes[type != "core TF"], aes(x, y - 0.22, label = label), angle = 60, hjust = 1, size = 2.1, colour = ink) +
  scale_colour_manual(values = src_col, name = "edge source") +
  scale_linetype_manual(values = c(`1` = "solid", `-1` = "dashed"), labels = c(`1` = "activating", `-1` = "repressing"),
                        name = "TF-TF edge (CollecTRI sign)") +
  scale_fill_manual(values = c(`07a eligible gene` = "#c0392b", `qPCR panel gene` = "#2a6f97"), name = NULL) +
  coord_cartesian(xlim = c(-0.6, 8.2), ylim = c(-1.3, 3.1), clip = "off") +
  labs(title = "E2F1, MYBL2, FOXM1 and NFIC: edges among them and to panel / eligible genes",
       subtitle = paste("Regulons are inferred (pySCENIC: co-expression + motif) or curated (CollecTRI), not experimentally proven here.",
                        "TF labels: total targets per source.", "FOXM1 and NFIC are also panel genes; their incoming edges are the TF-TF arcs. Genes without an edge to these TFs are not shown.", sep = "\n")) +
  theme_void(base_size = base) + theme(plot.title = element_text(face = "bold"), plot.subtitle = element_text(colour = ink2, size = base - 1),
                                       legend.position = "bottom", legend.box = "vertical")
save_fig(f45, "Fig45_network_core_tfs", 180, 125)

## =========================================================================
## Fig 46: story diagram, pathways -> TFs -> panel genes
## =========================================================================
st <- rd("F_story_table.csv")[block != "qPCR reference"]
e2 <- rd("E_target_set_overlap.csv")[FDR < fdr_cut & grepl("^HALLMARK", set)]
rep <- rd("B_replicated_sets.csv")[collection == "H"]
## gene -> replicated Hallmark set (shared leading edge)
pg <- st[replicated_hallmark_LE != "", .(set = trimws(unlist(strsplit(replicated_hallmark_LE, ";")))), by = .(gene, block)]
## TF -> gene: core TFs (either source) and pySCENIC consensus regulators
tg <- rbind(st[core_tf_targeting != "", .(TF = unique(sub("_.*", "", trimws(unlist(strsplit(core_tf_targeting, ";")))))), by = .(gene, block)],
            st[regulators_pySCENIC != "", .(TF = unlist(strsplit(regulators_pySCENIC, ";"))), by = .(gene, block)])
tg <- unique(tg)[TF != gene]
## pathway -> TF: core TF targets enriched in the set (E2, FDR < 0.05, either source)
pt <- unique(e2[, .(set = sub("^HALLMARK_", "", set), TF)])
genes <- unique(c(pg$gene, tg$gene))
lv_p <- unique(pg$set); lv_t <- unique(c(intersect(core, tg$TF), sort(setdiff(tg$TF, core))))
lv_g <- st[gene %in% genes][order(block, gene), gene]
pos <- function(v, x) data.table(id = v, x = x, y = if (length(v) > 1) seq(1, 0, length.out = length(v)) else 0.5)
np <- pos(lv_p, 0); nt <- pos(lv_t, 1); ng <- pos(lv_g, 2)
ng[, block := st$block[match(id, st$gene)]]
seg <- function(d, a, b, na, nb) { d <- merge(d, na[, .(k = id, x0 = x, y0 = y)], by.x = a, by.y = "k")
                                   merge(d, nb[, .(k = id, x1 = x, y1 = y)], by.x = b, by.y = "k") }
s_pt <- seg(pt[set %in% lv_p & TF %in% lv_t], "set", "TF", np, nt)
s_tg <- seg(tg, "TF", "gene", nt, ng)
s_pg <- seg(pg, "set", "gene", np, ng)
outc <- st[gene %in% lv_g, .(id = gene, gse = as.numeric(sub(" .*", "", GSE41613_DSS_HR)), tcga = as.numeric(sub(" .*", "", TCGA_allsite_PFI_HR)))]
ng <- merge(ng, outc, by = "id")
f46 <- ggplot() +
  geom_segment(data = s_pg, aes(x0, y0, xend = x1, yend = y1), colour = "grey80", linewidth = 0.25, linetype = "dotted") +
  geom_segment(data = s_pt, aes(x0, y0, xend = x1, yend = y1), colour = "#7a4b94", linewidth = 0.5, alpha = 0.7) +
  geom_segment(data = s_tg, aes(x0, y0, xend = x1, yend = y1), colour = "#e69f00", linewidth = 0.35, alpha = 0.8) +
  geom_label(data = np, aes(x, y, label = gsub("_", " ", id)), hjust = 1, size = 2, label.size = 0.15, fill = "#fbeaea") +
  geom_label(data = nt, aes(x, y, label = id), size = 2, label.size = 0.15,
             fill = fifelse(nt$id %in% core, "#f3efe6", "white"), fontface = fifelse(nt$id %in% core, "bold", "plain")) +
  geom_point(data = ng, aes(x, y, fill = block), shape = 21, size = 2.6, stroke = 0.3) +
  geom_text(data = ng, aes(x + 0.05, y, label = sprintf("%s   TCGA %.2f | GSE41613 %.2f", id, tcga, gse)), hjust = 0, size = 2, colour = ink) +
  scale_fill_manual(values = c(`07a eligible` = "#c0392b", `qPCR core` = "#2a6f97", `qPCR reserve` = "#8fb8d6",
                               `qPCR optional (Arm B)` = "#e69f00"), name = NULL) +
  annotate("text", x = c(0, 1, 2), y = 1.08, label = c("replicated Hallmark set", "regulating TF", "gene (HR per SD: TCGA PFI | GSE41613 DSS)"),
           fontface = "bold", size = 2.4, hjust = c(1, 0.5, 0)) +
  coord_cartesian(xlim = c(-0.9, 3), ylim = c(-0.05, 1.12), clip = "off") +
  labs(title = "Story: pathways -> TFs -> panel and eligible genes (HPV-)",
       subtitle = paste("dotted: gene in the shared leading edge (C2 and TCGA) of a replicated Hallmark set; purple: core-TF targets enriched in the set",
                        "(hypergeometric FDR < 0.05); orange: TF -> gene edge (core TFs from either source; other TFs = pySCENIC consensus regulons).",
                        "Regulons are inferred, not experimentally proven. Genes without a replicated pathway or regulator are not shown.", sep = "\n")) +
  theme_void(base_size = base) + theme(plot.title = element_text(face = "bold"), plot.subtitle = element_text(colour = ink2, size = base - 1.5),
                                       legend.position = "bottom")
save_fig(f46, "Fig46_story_diagram", 200, 150)
