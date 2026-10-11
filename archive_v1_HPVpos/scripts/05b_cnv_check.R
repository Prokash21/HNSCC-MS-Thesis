## Independent check of the normal-reference cleaning (scripts/05) with copy-
## number inference, so the malignant/normal call does not rest on the marker
## genes used for cleaning (viral reads, MAGE antigens).
##
## Method: the inferCNV approach of Tirosh et al. 2016 / Puram et al. 2017,
## implemented directly (the infercnv package needs the JAGS system library).
## Per patient: log2(CPM/10 + 1), relative to that patient's immune/stromal
## reference cells, clipped to +-3, 101-gene moving average along each
## chromosome, centred per cell and on the reference. Per cell:
##   CNV signal      = mean square of the smoothed profile
##   CNV correlation = correlation with the patient's average tumour profile
## A cell has a tumour-like CNV profile when its signal is above the 99th
## percentile of the patient's reference cells and its correlation > 0.4.
## Only patients whose own tumour cells mostly pass are informative (HPV+
## tumours often have few CNVs).
##
## Inputs:  HNSCC-HPV/atlas/HNSCC_Atlas.rds (only to build the cache)
##          data/processed/HNSCatlas/nonmalignant_subtypes.rds   (scripts/05)
##          HNSCC-HPV/results/01_atlas/tables/patient_metadata.csv (scripts/04)
##          data/reference/ensembl/*.cdna.all.fa.gz (gene positions; scripts/10)
## Outputs: HNSCC-HPV/results/02_normal_cleaning/ (Fig3b + tables/)
##          data/processed/HNSCatlas/normal_reference_cells.rds (final normal
##          reference: clean squamous cells minus CNV-positive ones)
##          data/processed/HNSCatlas/cnv_input.rds (cache)
## Run from the repo root:  Rscript scripts/05b_cnv_check.R

suppressPackageStartupMessages({
  library(Matrix)
  library(ggplot2)
  library(ggrastr)
  library(patchwork)
})

in_dir  <- "data/processed/HNSCatlas"
out_dir <- "HNSCC-HPV/results/02_normal_cleaning"
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
set.seed(1)
max_tumour <- 300     # tumour cells sampled per patient (enough for its profile)
max_ref    <- 300     # immune/stromal reference cells sampled per patient
min_ref    <- 30
n_genes    <- 7000    # most highly expressed genes with a known position
window     <- 101
ref_types  <- c("CD4T", "CD8T", "NK", "Bcell", "PlasmaCell", "Macrophage", "Monocyte",
                "panDC", "Mast", "Fibroblast", "Endothelial")

## ---- cache: counts for the cells we need ----------------------------------
cache_file <- file.path(in_dir, "cnv_input.rds")
if (!file.exists(cache_file)) {
  message("Selecting cells from the full atlas (slow)...")
  suppressPackageStartupMessages(library(Seurat))
  atlas  <- readRDS("HNSCC-HPV/atlas/HNSCC_Atlas.rds")
  counts <- LayerData(atlas, assay = "RNA", layer = "counts")
  md     <- atlas@meta.data[, c("Patient", "Source", "Cell_Type")]
  rm(atlas)
  invisible(gc())
  sq <- readRDS(file.path(in_dir, "nonmalignant_subtypes.rds"))
  group <- rep(NA_character_, nrow(md))
  group[md$Cell_Type %in% ref_types] <- "Reference"
  group[md$Cell_Type %in% "Tumor" & md$Source != "CD45p"] <- "Tumour"
  nm_group <- ifelse(sq$squamous_clean, "Normal squamous",
                     ifelse(sq$subtype == "Squamous", "Squamous doublet", sq$subtype))
  group[match(sq$cell, rownames(md))] <- nm_group
  pick <- function(g, n) {
    idx <- which(group == g)
    unlist(lapply(split(idx, md$Patient[idx]), function(i) {
      if (length(i) > n) sample(i, n) else i
    }))
  }
  keep <- sort(c(pick("Tumour", max_tumour), pick("Reference", max_ref),
                 which(group %in% c("Normal squamous", "Likely malignant", "Glandular",
                                    "Ciliated", "Squamous doublet"))))
  saveRDS(list(counts = counts[, keep],
               meta = data.frame(cell = rownames(md)[keep], Patient = md$Patient[keep],
                                 Cell_Type = md$Cell_Type[keep], group = group[keep])),
          cache_file)
  rm(counts)
  invisible(gc())
}
cnv_in <- readRDS(cache_file)
meta <- cnv_in$meta
counts <- cnv_in$counts
pm <- read.csv("HNSCC-HPV/results/01_atlas/tables/patient_metadata.csv")
meta$HPV <- pm$HPV[match(meta$Patient, pm$Patient)]

## ---- gene positions from the Ensembl cDNA FASTA headers -------------------
fa <- list.files("data/reference/ensembl", pattern = "cdna.all.fa.gz$", full.names = TRUE)[1]
hdr <- system(sprintf("gzcat '%s' | grep '^>'", fa), intern = TRUE)
hdr <- hdr[grepl("chromosome:GRCh38:", hdr) & grepl("gene_symbol:", hdr)]
loc <- regmatches(hdr, regexpr("chromosome:GRCh38:[^ ]+", hdr))
sym <- sub("^gene_symbol:", "", regmatches(hdr, regexpr("gene_symbol:[^ ]+", hdr)))
parts <- do.call(rbind, strsplit(sub("chromosome:GRCh38:", "", loc), ":"))
pos <- data.frame(gene = sym, chr = parts[, 1], start = as.numeric(parts[, 2]))
pos <- pos[pos$chr %in% as.character(1:22), ]
pos <- aggregate(start ~ gene + chr, pos, min)
pos <- pos[!duplicated(pos$gene), ]

## genes: positioned autosomal genes, the n_genes most highly expressed
lib <- colSums(counts)
genes <- intersect(rownames(counts), pos$gene)
## log2(CPM/10 + 1) kept sparse (zeros stay zero), mean per gene
norm <- counts[genes, ] %*% Diagonal(x = 1e5 / lib)
norm@x <- log2(norm@x + 1)
mean_expr <- setNames(Matrix::rowMeans(norm), genes)
rm(norm)
genes <- names(sort(mean_expr, decreasing = TRUE))[seq_len(min(n_genes, length(genes)))]
pos <- pos[match(genes, pos$gene), ]
ord <- order(as.integer(pos$chr), pos$start)
genes <- genes[ord]
pos <- pos[ord, ]
chr_runs <- split(seq_along(genes), factor(pos$chr, levels = as.character(1:22)))
message(sprintf("%d genes on chromosomes 1-22; %d cells", length(genes), ncol(counts)))

## moving average within each chromosome (window truncated at the ends)
smooth <- function(m) {
  half <- (window - 1) %/% 2
  out <- m
  for (idx in chr_runs) {
    if (length(idx) < 2) next
    cs <- rbind(0, apply(m[idx, , drop = FALSE], 2, cumsum))
    n <- length(idx)
    lo <- pmax(1, seq_len(n) - half)
    hi <- pmin(n, seq_len(n) + half)
    out[idx, ] <- (cs[hi + 1, , drop = FALSE] - cs[lo, , drop = FALSE]) / (hi - lo + 1)
  }
  out
}

## ---- per patient ----------------------------------------------------------
patients <- names(which(table(meta$Patient[meta$group == "Reference"]) >= min_ref))
patients <- intersect(patients, meta$Patient[meta$group == "Tumour"])
example_profiles <- list()
scores <- do.call(rbind, lapply(patients, function(p) {
  cells <- which(meta$Patient == p)
  e <- log2(t(t(as.matrix(counts[genes, cells])) / lib[cells]) * 1e5 + 1)   # log2(CPM/10 + 1)
  is_ref <- meta$group[cells] == "Reference"
  r <- pmin(pmax(e - rowMeans(e[, is_ref, drop = FALSE]), -3), 3)
  s <- smooth(r)
  s <- sweep(s, 2, apply(s, 2, median))
  s <- s - rowMeans(s[, is_ref, drop = FALSE])
  tumour_profile <- rowMeans(s[, meta$group[cells] == "Tumour", drop = FALSE])
  signal <- colMeans(s^2)
  corr <- as.vector(cor(s, tumour_profile))
  thr <- quantile(signal[is_ref], 0.99)
  if (p %in% c("BHN63", "P8")) example_profiles[[p]] <<- list(s = s, group = meta$group[cells])
  data.frame(cell = meta$cell[cells], Patient = p, HPV = meta$HPV[cells][1],
             group = meta$group[cells], signal = signal, signal_rel = signal / thr,
             correlation = corr, tumour_like = signal > thr & corr > 0.4)
}))

## informative patients: at least half of their own tumour cells are tumour-like
inf <- tapply(scores$tumour_like[scores$group == "Tumour"], scores$Patient[scores$group == "Tumour"], mean)
informative <- names(inf)[inf >= 0.5]
scores$informative <- scores$Patient %in% informative
write.csv(scores, file.path(out_dir, "tables", "cnv_scores_per_cell.csv"), row.names = FALSE)

## final normal reference for the DGE (scripts/06) and qPCR screen (scripts/08):
## the clean squamous cells, minus those with a tumour-like CNV profile in
## patients where the test works (elsewhere it cannot be checked)
sq <- readRDS(file.path(in_dir, "nonmalignant_subtypes.rds"))
cnv_pos <- scores$cell[scores$group == "Normal squamous" & scores$informative & scores$tumour_like]
normal_ref <- setdiff(sq$cell[sq$squamous_clean], cnv_pos)
saveRDS(normal_ref, file.path(in_dir, "normal_reference_cells.rds"))
message(sprintf("Normal reference: %d cells kept, %d removed as CNV-positive",
                length(normal_ref), length(cnv_pos)))

grp_levels <- c("Reference", "Normal squamous", "Glandular", "Likely malignant", "Tumour")
summ <- do.call(rbind, lapply(grp_levels, function(g) {
  d <- scores[scores$group == g & scores$informative, ]
  data.frame(group = g, cells = nrow(d), patients = length(unique(d$Patient)),
             pct_tumour_like = round(100 * mean(d$tumour_like), 1))
}))
write.csv(summ, file.path(out_dir, "tables", "cnv_summary.csv"), row.names = FALSE)
message(sprintf("Informative patients: %d of %d (HPV+ %d, HPV- %d)", length(informative), length(patients),
                sum(pm$HPV[match(informative, pm$Patient)] == "HPV+"),
                sum(pm$HPV[match(informative, pm$Patient)] == "HPV-")))
print(summ, row.names = FALSE)

## ---- Figure 3b ------------------------------------------------------------
base_size <- 7
ink <- "#0b0b0b"; ink2 <- "#52514e"; rule <- "#c3c2b7"; grid_col <- "#e1e0d9"
grp_cols <- c(Reference = "#898781", "Normal squamous" = "#4a3aa7", Glandular = "#1baf7a",
              "Likely malignant" = "#c98500", Tumour = "#e34948")
theme_fig <- function() {
  theme_classic(base_size = base_size) +
    theme(text = element_text(colour = ink), axis.text = element_text(colour = ink2),
          axis.line = element_line(linewidth = 0.3, colour = rule),
          axis.ticks = element_line(linewidth = 0.3, colour = rule),
          plot.title = element_text(size = base_size + 1, face = "bold"),
          plot.subtitle = element_text(size = base_size - 0.5, colour = ink2),
          plot.tag = element_text(size = base_size + 3, face = "bold"),
          strip.background = element_blank(),
          strip.text = element_text(size = base_size - 0.5, face = "bold"),
          legend.text = element_text(size = base_size - 0.5), legend.key.size = unit(3, "mm"),
          plot.margin = margin(4, 6, 4, 4))
}

## a: CNV heatmaps for two example patients
heat <- function(p) {
  ex <- example_profiles[[p]]
  g <- factor(ex$group, levels = grp_levels)
  keep <- unlist(lapply(split(seq_along(g), g), function(i) if (length(i) > 80) sample(i, 80) else i))
  keep <- keep[order(g[keep])]
  m <- ex$s[, keep]
  df <- data.frame(gene = rep(seq_len(nrow(m)), times = ncol(m)),
                   cell = rep(seq_len(ncol(m)), each = nrow(m)),
                   value = pmax(pmin(as.vector(m), 0.4), -0.4))
  chr_mid <- vapply(chr_runs, function(i) mean(i), numeric(1))
  chr_end <- vapply(chr_runs, max, numeric(1))
  grp_n <- c(table(droplevels(g[keep])))     # plain named vector, not a table
  grp_end <- cumsum(grp_n)
  ggplot(df, aes(gene, cell, fill = value)) +
    geom_raster() +
    geom_vline(xintercept = chr_end[-length(chr_end)] + 0.5, linewidth = 0.15, colour = ink2) +
    geom_hline(yintercept = grp_end[-length(grp_end)] + 0.5, linewidth = 0.3, colour = ink) +
    scale_fill_gradient2(low = "#2a78d6", mid = "#f0efec", high = "#e34948", midpoint = 0,
                         limits = c(-0.4, 0.4), name = "Inferred\nCNV") +
    scale_x_continuous(breaks = chr_mid[c(TRUE, FALSE)], labels = names(chr_runs)[c(TRUE, FALSE)],
                       expand = c(0, 0)) +
    scale_y_continuous(breaks = grp_end - grp_n / 2, labels = names(grp_n), expand = c(0, 0)) +
    labs(title = sprintf("Patient %s (%s)", p, pm$HPV[pm$Patient == p]),
         x = "Chromosome", y = NULL) +
    theme_fig() + theme(axis.line = element_blank(), axis.ticks.y = element_blank(),
                        axis.text.y = element_text(colour = ink))
}
pa <- wrap_plots(lapply(intersect(c("BHN63", "P8"), names(example_profiles)), heat), ncol = 1)

## b: signal vs correlation, informative patients, one panel per group
d <- scores[scores$informative, ]
d$group <- factor(d$group, levels = grp_levels)
d <- d[!is.na(d$group), ]
pb <- ggplot(d, aes(pmin(signal_rel, 6), correlation)) +
  geom_point_rast(aes(colour = group), size = 0.35, stroke = 0, shape = 16, raster.dpi = 600,
                  dev = "ragg") +
  annotate("rect", xmin = 1, xmax = Inf, ymin = 0.4, ymax = Inf, fill = NA, colour = ink2,
           linewidth = 0.25, linetype = "22") +
  facet_wrap(~ group, nrow = 1) +
  scale_colour_manual(values = grp_cols, guide = "none") +
  labs(title = "CNV signal vs similarity to the patient's tumour",
       subtitle = sprintf("%d informative patients; dashed box = tumour-like CNV profile", length(informative)),
       x = "CNV signal (1 = 99th percentile of reference cells)", y = "Correlation with tumour CNV") +
  theme_fig()

## c: % tumour-like per group
summ$group <- factor(summ$group, levels = grp_levels)
pc <- ggplot(summ, aes(group, pct_tumour_like, fill = group)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = sprintf("%.0f%%\n(%s cells)", pct_tumour_like, format(cells, big.mark = ",", trim = TRUE))),
            vjust = -0.3, size = (base_size - 1.5) / .pt, colour = ink2, lineheight = 0.9) +
  scale_fill_manual(values = grp_cols, guide = "none") +
  scale_y_continuous(limits = c(0, 115), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  labs(title = "Cells with a tumour-like CNV profile", x = NULL, y = "% of cells") +
  theme_fig() + theme(axis.text.x = element_text(angle = 30, hjust = 1, colour = ink))

fig3b <- (pa | pc) / pb + plot_layout(heights = c(1.3, 1)) + plot_annotation(tag_levels = "a")
for (ext in c("png", "pdf")) {
  ggsave(file.path(out_dir, paste0("Fig3b_CNV_check.", ext)), fig3b, width = 180, height = 170,
         units = "mm", dpi = 600, bg = "white", device = if (ext == "png") ragg::agg_png else "pdf")
}
message("Wrote ", file.path(out_dir, "Fig3b_CNV_check.png/.pdf"))
