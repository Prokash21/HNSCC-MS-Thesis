## Figures for the pseudobulk DGE (scripts/06_dge_pseudobulk.R).
##   Fig4: volcano plots (A1-A3, B) + number of DE genes per comparison
##   Fig5: overlap of up-regulated genes between comparisons (UpSet + Venn)
##   Fig6: built-in checks (CDKN2A by HPV; HPV16-E6 and MAGEA4 tumour vs
##         normal squamous, paired by patient)
##
## Run from the repo root:  Rscript scripts/07_dge_figures.R

suppressPackageStartupMessages({
  library(edgeR)
  library(ggplot2)
  library(ggrastr)
  library(ggrepel)
  library(patchwork)
})

res     <- readRDS("data/processed/HNSCatlas/dge_results.rds")
fig_dir <- "HNSCC-HPV/results/03_DGE"

## ---- shared style (as in scripts/03) --------------------------------------
base_size <- 7
ink       <- "#0b0b0b"
ink2      <- "#52514e"
rule      <- "#c3c2b7"
grid_col  <- "#e1e0d9"
ns_col    <- "#d6d5cf"
col_tumor <- "#e34948"
col_norm  <- "#4a3aa7"
col_hpvp  <- "#2a78d6"
col_hpvn  <- "#eb6834"
col_third <- "#1baf7a"

theme_fig <- function() {
  theme_classic(base_size = base_size) +
    theme(text = element_text(colour = ink),
          axis.text = element_text(colour = ink2),
          axis.line = element_line(linewidth = 0.3, colour = rule),
          axis.ticks = element_line(linewidth = 0.3, colour = rule),
          plot.title = element_text(size = base_size + 1, face = "bold"),
          plot.subtitle = element_text(size = base_size - 0.5, colour = ink2),
          plot.tag = element_text(size = base_size + 3, face = "bold"),
          legend.text = element_text(size = base_size - 0.5),
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

fmt   <- function(n) format(n, big.mark = ",", trim = TRUE)
p_fmt <- function(p) if (p < 0.001) formatC(p, format = "e", digits = 1) else sprintf("%.3f", p)
up    <- function(t) t$gene[t$significant & t$logFC > 0]
dn    <- function(t) t$gene[t$significant & t$logFC < 0]

n_pat <- function(s, hpv = NULL) {
  length(unique(s$Patient[is.null(hpv) | s$HPV %in% hpv]))
}

## ---- Figure 4: volcano plots ----------------------------------------------
volcano <- function(tab, title, pos, neg, highlight = character(0), n_lab = 8) {
  tab$y <- -log10(pmax(tab$FDR, 1e-300))
  tab$dir <- ifelse(!tab$significant, "ns", ifelse(tab$logFC > 0, "pos", "neg"))
  tab$score <- abs(tab$logFC) * tab$y
  top <- function(d) head(d[order(-d$score), ], n_lab)
  lab <- rbind(top(tab[tab$dir == "pos", ]), top(tab[tab$dir == "neg", ]),
               tab[tab$gene %in% highlight, ])
  lab <- lab[!duplicated(lab$gene), ]
  ggplot(tab, aes(logFC, y)) +
    geom_point_rast(data = tab[tab$dir == "ns", ], colour = ns_col, size = 0.35,
                    shape = 16, raster.dpi = 600, dev = "ragg") +
    geom_hline(yintercept = -log10(0.05), linetype = "22", linewidth = 0.25, colour = ink2) +
    geom_vline(xintercept = c(-1, 1), linetype = "22", linewidth = 0.25, colour = ink2) +
    geom_point(data = tab[tab$dir != "ns", ], aes(colour = dir), size = 0.9,
               shape = 16, stroke = 0) +
    geom_text_repel(data = lab, aes(label = gene), size = (base_size - 1.5) / .pt,
                    colour = ink, max.overlaps = Inf, min.segment.length = 0,
                    segment.size = 0.2, segment.colour = ink2, box.padding = 0.3,
                    force = 2, max.time = 3, max.iter = 1e5,
                    bg.color = "white", bg.r = 0.1, seed = 1) +
    scale_colour_manual(values = c(pos = pos$col, neg = neg$col), guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
    labs(title = title,
         subtitle = sprintf("%s higher in %s  |  %s higher in %s",
                            fmt(sum(tab$dir == "pos")), pos$lab,
                            fmt(sum(tab$dir == "neg")), neg$lab),
         x = "log2 fold change", y = "-log10 FDR") +
    theme_fig()
}

tum  <- list(col = col_tumor, lab = "tumour")
norm <- list(col = col_norm,  lab = "normal squamous")
pA1 <- volcano(res$A1, sprintf("A1  All paired patients (n = %d)", n_pat(res$sA)),
               tum, norm, c("MAGEA4", "HPV16-E6"))
pA2 <- volcano(res$A2, sprintf("A2  HPV+ patients (n = %d)", n_pat(res$sA, "HPV+")),
               tum, norm, "HPV16-E6")
pA3 <- volcano(res$A3, sprintf("A3  HPV- patients (n = %d)", n_pat(res$sA, "HPV-")),
               tum, norm, "MAGEA4")
pB  <- volcano(res$B, sprintf("B  HPV+ vs HPV- tumour cells (%d vs %d patients, study-adjusted)",
                              n_pat(res$sB, "HPV+"), n_pat(res$sB, "HPV-")),
               list(col = col_hpvp, lab = "HPV+"), list(col = col_hpvn, lab = "HPV-"),
               "CDKN2A", n_lab = 10)

## number of DE genes per comparison (right = up, left = down)
tests <- c(A1 = "A1  all", A2 = "A2  HPV+", A3 = "A3  HPV-",
           A4 = "A4  HPV+ vs HPV- change", B = "B  HPV+ vs HPV-")
cnt <- do.call(rbind, lapply(names(tests), function(nm) {
  t <- res[[nm]]
  data.frame(test = tests[[nm]], dir = c("pos", "neg"),
             n = c(sum(t$significant & t$logFC > 0), -sum(t$significant & t$logFC < 0)))
}))
cnt$test <- factor(cnt$test, levels = rev(tests))
cnt$fill <- ifelse(grepl("^A[123]", cnt$test),
                   ifelse(cnt$dir == "pos", col_tumor, col_norm),
                   ifelse(cnt$dir == "pos", col_hpvp, col_hpvn))
pCnt <- ggplot(cnt, aes(n, test)) +
  geom_col(aes(fill = fill), width = 0.65) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = rule) +
  geom_text(aes(label = fmt(abs(n)), hjust = ifelse(dir == "pos", -0.15, 1.15)),
            size = (base_size - 1) / .pt, colour = ink2) +
  scale_fill_identity() +
  scale_x_continuous(labels = function(x) fmt(abs(x)), expand = expansion(mult = 0.25)) +
  labs(title = "DE genes per comparison",
       subtitle = "right = higher in tumour / HPV+\nleft = higher in normal squamous / HPV-",
       x = "Genes (FDR < 0.05, |log2FC| >= 1)", y = NULL) +
  theme_fig()

fig4 <- (pA1 | pA2 | pA3) / ((pB | pCnt) + plot_layout(widths = c(2, 1))) +
  plot_annotation(tag_levels = "a")
save_fig(fig4, "Fig4_DGE_volcano", w = 180, h = 150)

## ---- Figure 5: overlaps (UpSet + Venn) ------------------------------------
sets <- list("A1 tumour up (all)"      = up(res$A1),
             "A2 tumour up (HPV+)"     = up(res$A2),
             "A3 tumour up (HPV-)"     = up(res$A3),
             "B higher in HPV+ tumour" = up(res$B),
             "B higher in HPV- tumour" = dn(res$B))
genes <- unique(unlist(sets))
code  <- apply(sapply(sets, function(s) genes %in% s), 1, paste, collapse = "")
code  <- gsub("TRUE", "1", gsub("FALSE", "0", code))
inter <- head(sort(table(code), decreasing = TRUE), 14)
ids   <- names(inter)
has   <- function(id, k) substr(id, k, k) == "1"

bars <- data.frame(id = factor(ids, levels = ids), n = as.integer(inter),
                   cls = ifelse(has(ids, 2) & has(ids, 4), "HPV+ candidates",
                                ifelse(has(ids, 3) & has(ids, 5), "HPV- candidates", "other")))
dots <- expand.grid(id = ids, k = seq_along(sets), stringsAsFactors = FALSE)
dots$on  <- has(dots$id, dots$k)
dots$set <- factor(names(sets)[dots$k], levels = rev(names(sets)))
dots$id  <- factor(dots$id, levels = ids)
sizes <- data.frame(set = factor(names(sets), levels = rev(names(sets))), n = lengths(sets))

p_bar <- ggplot(bars, aes(id, n, fill = cls)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = fmt(n)), vjust = -0.4, size = (base_size - 1) / .pt, colour = ink2) +
  scale_fill_manual(values = c("HPV+ candidates" = col_hpvp, "HPV- candidates" = col_hpvn,
                               other = "#898781"),
                    breaks = c("HPV+ candidates", "HPV- candidates"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Overlap of up-regulated genes between comparisons",
       x = NULL, y = "Genes in intersection") +
  theme_fig() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(), legend.position = "inside",
        legend.position.inside = c(1, 1), legend.justification = c(1, 1))
p_mat <- ggplot(dots, aes(id, set)) +
  geom_point(data = dots[!dots$on, ], colour = grid_col, size = 2.2) +
  geom_line(data = dots[dots$on, ], aes(group = id), colour = ink, linewidth = 0.4) +
  geom_point(data = dots[dots$on, ], colour = ink, size = 2.2) +
  labs(x = NULL, y = NULL) +
  theme_fig() +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(),
        axis.line = element_blank(), axis.text.y = element_text(colour = ink))
p_set <- ggplot(sizes, aes(n, set)) +
  geom_col(fill = "#898781", width = 0.6) +
  geom_text(aes(label = fmt(n)), hjust = 1.15, size = (base_size - 1) / .pt, colour = ink2) +
  scale_x_reverse(expand = expansion(mult = c(0.45, 0))) +   # room for labels on the left
  labs(x = "Set size", y = NULL) +
  theme_fig() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank())
upset <- plot_spacer() + p_bar + p_set + p_mat +
  plot_layout(ncol = 2, widths = c(1, 3.2), heights = c(2.2, 1.3))

## 3-set Venn drawn directly (equal circles; numbers = genes per region)
venn3 <- function(sets, cols, title, subtitle, bold_codes) {
  g <- unique(unlist(sets))
  cd <- apply(sapply(sets, function(s) g %in% s), 1, paste, collapse = "")
  cd <- gsub("TRUE", "1", gsub("FALSE", "0", cd))
  reg <- data.frame(code = c("100", "010", "001", "110", "101", "011", "111"),
                    x = c(0, -1.1, 1.1, -0.65, 0.65, 0, 0),
                    y = c(1.15, -0.75, -0.75, 0.3, 0.3, -0.95, -0.1))
  reg$n <- vapply(reg$code, function(k) sum(cd == k), integer(1))
  reg$face <- ifelse(reg$code %in% bold_codes, "bold", "plain")
  ctr <- data.frame(cx = c(0, -0.6, 0.6), cy = c(0.55, -0.45, -0.45))
  th <- seq(0, 2 * pi, length.out = 240)
  circ <- do.call(rbind, lapply(1:3, function(i) {
    data.frame(set = names(sets)[i], x = ctr$cx[i] + cos(th), y = ctr$cy[i] + sin(th))
  }))
  circ$set <- factor(circ$set, levels = names(sets))
  lab <- data.frame(x = c(0, -1.3, 1.3), y = c(1.8, -1.75, -1.75),
                    label = sprintf("%s\n(%s)", names(sets), fmt(lengths(sets))))
  ggplot() +
    geom_polygon(data = circ, aes(x, y, group = set, fill = set, colour = set),
                 alpha = 0.12, linewidth = 0.5) +
    geom_text(data = reg, aes(x, y, label = fmt(n), fontface = face),
              size = base_size / .pt, colour = ink) +
    geom_text(data = lab, aes(x, y, label = label), size = (base_size - 0.5) / .pt,
              lineheight = 0.9, colour = ink) +
    scale_fill_manual(values = setNames(cols, names(sets)), guide = "none") +
    scale_colour_manual(values = setNames(cols, names(sets)), guide = "none") +
    coord_fixed(xlim = c(-2.3, 2.3), ylim = c(-2.15, 2.15), clip = "off") +
    labs(title = title, subtitle = subtitle) +
    theme_void(base_size = base_size) +
    theme(plot.title = element_text(size = base_size + 1, face = "bold", colour = ink),
          plot.subtitle = element_text(size = base_size - 0.5, colour = ink2),
          plot.tag = element_text(size = base_size + 3, face = "bold"))
}
venn_cols <- c(col_hpvp, col_hpvn, col_third)
n_genes <- function(n) sprintf("%d %s", n, if (n == 1) "gene" else "genes")
v_pos <- venn3(sets[c(2, 3, 4)], venn_cols, "HPV+ candidates",
               sprintf("bold = tumour-up in HPV+ AND higher in HPV+ tumours (%s)",
                       n_genes(nrow(res$cand_pos))), c("101", "111"))
v_neg <- venn3(sets[c(2, 3, 5)], venn_cols, "HPV- candidates",
               sprintf("bold = tumour-up in HPV- AND higher in HPV- tumours (%s)",
                       n_genes(nrow(res$cand_neg))), c("011", "111"))

fig5 <- wrap_elements(full = upset) / (v_pos | v_neg) +
  plot_layout(heights = c(1.1, 1)) +
  plot_annotation(tag_levels = "a")
save_fig(fig5, "Fig5_DGE_overlaps", w = 180, h = 170)

## ---- Figure 6: built-in checks --------------------------------------------
## expression from the unfiltered pseudobulk: filterByExpr drops genes such
## as the viral transcripts that only a few samples express
pb_counts <- readRDS("data/processed/HNSCatlas/pseudobulk.rds")$counts
lc_all <- cpm(DGEList(pb_counts), log = TRUE, prior.count = 2)
stat <- function(tab, g) {
  r <- tab[tab$gene == g, ]
  if (nrow(r) == 0) return("not tested (too few samples express it)")
  sprintf("log2FC %.2f, FDR %s", r$logFC, p_fmt(r$FDR))
}

hpv_cols <- c("HPV+" = col_hpvp, "HPV-" = col_hpvn)
## dot plot by HPV label with manual jitter, so a label can sit on its point
by_hpv <- function(df, title, subtitle, ylab, label_ids = character(0)) {
  set.seed(1)
  df$x <- as.numeric(df$HPV) + runif(nrow(df), -0.15, 0.15)
  ggplot(df, aes(x, y)) +
    geom_boxplot(aes(x = as.numeric(HPV), group = HPV), outlier.shape = NA, width = 0.5,
                 colour = ink2, linewidth = 0.3, fill = NA) +
    geom_point(aes(colour = HPV), size = 1.2, stroke = 0, shape = 16) +
    geom_text_repel(data = df[df$Patient %in% label_ids, ], aes(label = Patient),
                    size = (base_size - 1.5) / .pt, colour = ink, min.segment.length = 0,
                    segment.size = 0.2, nudge_x = 0.3, seed = 1) +
    scale_x_continuous(breaks = 1:2, labels = levels(df$HPV), expand = expansion(add = 0.5)) +
    scale_colour_manual(values = hpv_cols, guide = "none") +
    labs(title = title, subtitle = subtitle, x = NULL, y = ylab) +
    theme_fig()
}

## a: p16 in tumour cells; the low HPV+ outlier is labelled
dfa <- data.frame(Patient = res$sB$Patient, y = lc_all["CDKN2A", res$sB$sample],
                  HPV = factor(res$sB$HPV, levels = names(hpv_cols)))
low <- dfa$Patient[dfa$HPV == "HPV+"][which.min(dfa$y[dfa$HPV == "HPV+"])]
pa <- by_hpv(dfa, "CDKN2A (p16), tumour cells", paste("B:", stat(res$B, "CDKN2A")),
             "log2 CPM (pseudobulk per patient)", low)

## b: HPV16 reads in tumour cells vs clinical label (only Puram 2023 mapped viral
## genes)
vir <- grep("^HPV", rownames(pb_counts), value = TRUE)
puram <- res$sB$Study == "Puram 2023 (GSE182227)"
cntB <- pb_counts[, res$sB$sample[puram]]
dfb <- data.frame(Patient = res$sB$Patient[puram],
                  y = colSums(cntB[vir, , drop = FALSE]) / colSums(cntB) * 1e4,
                  HPV = factor(res$sB$HPV[puram], levels = names(hpv_cols)))
pb <- by_hpv(dfb, "HPV16 reads, tumour cells",
             "Puram 2023 only (the other studies did not map viral reads)",
             "Viral reads per 10,000 reads")

paired_plot <- function(gene, hpv, tab, test, studies = unique(res$sA$Study)) {
  s <- res$sA[res$sA$HPV %in% hpv & res$sA$Study %in% studies, ]
  df <- data.frame(Patient = s$Patient, y = lc_all[gene, s$sample],
                   Compartment = factor(ifelse(s$Compartment == "Tumor", "Tumour", "Normal squamous"),
                                        levels = c("Normal squamous", "Tumour")))
  ggplot(df, aes(Compartment, y)) +
    geom_line(aes(group = Patient), colour = ns_col, linewidth = 0.3) +
    geom_point(aes(colour = Compartment), size = 1.2, stroke = 0, shape = 16) +
    scale_colour_manual(values = c("Normal squamous" = col_norm, "Tumour" = col_tumor),
                        guide = "none") +
    labs(title = sprintf("%s, tumour vs normal squamous", gene),
         subtitle = paste0(test, ": ", stat(tab, gene)),
         x = NULL, y = "log2 CPM (pseudobulk per patient)") +
    theme_fig()
}
pc <- paired_plot("HPV16-E6", "HPV+", res$A2, "A2 (HPV+, Puram 2023 patients shown)",
                  "Puram 2023 (GSE182227)")
pd <- paired_plot("MAGEA4", c("HPV+", "HPV-"), res$A1, "A1 (all paired patients)")

fig6 <- (pa | pb) / (pc | pd) + plot_annotation(tag_levels = "a")
save_fig(fig6, "Fig6_DGE_checks", w = 180, h = 140)
