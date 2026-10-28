## Step 05b (part A-B): ranking statistics and GSEA
##
## Purpose: read MSigDB v2024.1 Hallmark, GO BP and Reactome sets (HGNC-mapped),
##   rank genes by sign(logFC) * sqrt(F) for the step 04b contrasts (HPV- C2,
##   C1, TCGA paired; HPV+ C2, C1), run fgsea per contrast and collection (plus
##   the step 07a programme sets), call replicated sets (C2 and TCGA HPV-) and
##   mark qPCR and eligible genes in the leading edges.
## Inputs:  data/reference/msigdb/{h.all,c5.go.bp,c2.cp.reactome}.v2024.1.Hs.symbols.gmt, HGNC complete set,
##          results/step04b_DE/tables/{D_C2_*, D_C1_*, E_TCGA_HPVneg_paired_tumour_vs_normal}.csv,
##          results/step07a_recurrence/tables/{A_gene_sets.csv, D_recurrence_ranked.csv},
##          results/step06_qPCR/tables/D_panel_final.csv
## Outputs: results/step05b_functional/tables/B_*.csv, data/processed/v2/step05b_gene_sets.rds
## Date: 2026-10-28

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
})
source("config/paths.R")
source("scripts/lib_de.R")

tdir <- file.path(paths$step05b, "tables")
dir.create(tdir, recursive = TRUE, showWarnings = FALSE)
t04 <- file.path(paths$step04b, "tables")
fdr_cut <- 0.05
min_size <- 15
max_size <- 500
hg <- read_hgnc(paths$hgnc)

## =========================================================================
## A1. Gene sets (MSigDB v2024.1), symbols mapped to HGNC as the DE tables
## =========================================================================
read_gmt <- function(f) {
  stopifnot(file.exists(f))
  l <- strsplit(readLines(f), "\t")
  setNames(lapply(l, function(x) x[-(1:2)]), sapply(l, `[`, 1))
}
raw <- list(H = read_gmt(paths$msigdb_h), GOBP = read_gmt(paths$msigdb_gobp), REACTOME = read_gmt(paths$msigdb_reactome))
syms <- unique(unlist(raw))
mp <- hgnc_map(syms, hg)
sym_stats <- mp[, .N, by = how][order(-N)]
fwrite(sym_stats, file.path(tdir, "A_geneset_symbol_mapping.csv"))
map_set <- function(s) unique(na.omit(mp$mapped[match(s, mp$original)]))
sets <- lapply(raw, function(col) lapply(col, map_set))
gs07 <- fread(file.path(paths$step07a, "tables", "A_gene_sets.csv"))
sets$RECUR <- split(gs07$gene, gs07$set)
saveRDS(sets, paths$step05b_sets)
print(sym_stats)

## gene groups for marking
panel <- fread(file.path(paths$step06, "tables", "D_panel_final.csv"))
qpcr <- panel[set %in% c("core", "reserve", "optional (Arm B)"), gene]
eligible <- fread(file.path(paths$step07a, "tables", "D_recurrence_ranked.csv"))$gene
stopifnot(length(qpcr) == 20, length(eligible) == 7)

## =========================================================================
## A2. Ranking statistics: sign(logFC) * sqrt(F)
## =========================================================================
contrasts <- c(C2_HPVneg = "D_C2_HPVneg_allpatients.csv", C1_HPVneg = "D_C1_HPVneg_paired.csv",
               TCGA_HPVneg = "E_TCGA_HPVneg_paired_tumour_vs_normal.csv",
               C2_HPVpos = "D_C2_HPVpos_allpatients.csv", C1_HPVpos = "D_C1_HPVpos_paired.csv")
stats <- lapply(contrasts, function(f) {
  d <- fread(file.path(t04, f))
  stopifnot(!anyDuplicated(d$gene))
  s <- setNames(sign(d$logFC) * sqrt(pmax(d$F, 0)), d$gene)
  sort(s[is.finite(s)], decreasing = TRUE)
})
fwrite(data.table(contrast = names(stats), genes = lengths(stats),
                  ties = sapply(stats, function(s) sum(duplicated(s)))), file.path(tdir, "A_ranking_stats.csv"))

## =========================================================================
## B1. fgsea per contrast and collection
## =========================================================================
set.seed(1)
res <- rbindlist(lapply(names(stats), function(k) rbindlist(lapply(names(sets), function(col) {
  set.seed(1)
  r <- fgseaMultilevel(sets[[col]], stats[[k]], minSize = min_size, maxSize = max_size, eps = 0, nproc = 1)
  r <- as.data.table(r)
  if (!nrow(r)) return(NULL)
  r[, `:=`(contrast = k, collection = col,
           leadingEdge_genes = sapply(leadingEdge, paste, collapse = ";"),
           qpcr_in_LE = sapply(leadingEdge, function(x) paste(intersect(x, qpcr), collapse = ";")),
           eligible_in_LE = sapply(leadingEdge, function(x) paste(intersect(x, eligible), collapse = ";")))]
  r[, leadingEdge := NULL]
  message(sprintf("B1. %s %s: %d sets tested, %d FDR < %.2f", k, col, nrow(r), sum(r$padj < fdr_cut, na.rm = TRUE), fdr_cut))
  r
}))))
setnames(res, c("pathway", "padj"), c("set", "FDR"))
setcolorder(res, c("contrast", "collection", "set", "size", "NES", "ES", "pval", "FDR", "log2err"))
fwrite(res, file.path(tdir, "B_gsea_all.csv"))
fwrite(res[, .(tested = .N, sig = sum(FDR < fdr_cut, na.rm = TRUE), up = sum(FDR < fdr_cut & NES > 0, na.rm = TRUE),
               down = sum(FDR < fdr_cut & NES < 0, na.rm = TRUE)), by = .(contrast, collection)],
       file.path(tdir, "B_gsea_counts.csv"))

## =========================================================================
## B2. Replicated sets (HPV-): same direction and FDR < 0.05 in C2 and TCGA
## =========================================================================
w <- dcast(res[contrast %in% c("C2_HPVneg", "C1_HPVneg", "TCGA_HPVneg")], collection + set ~ contrast,
           value.var = c("NES", "FDR", "size"))
w[, replicated := (FDR_C2_HPVneg < fdr_cut & FDR_TCGA_HPVneg < fdr_cut & sign(NES_C2_HPVneg) == sign(NES_TCGA_HPVneg)) %in% TRUE]
w[, direction := fifelse(replicated, fifelse(NES_C2_HPVneg > 0, "up in tumour", "down in tumour"), "")]
w[, C1_same_direction := sign(NES_C1_HPVneg) == sign(NES_C2_HPVneg)]
le <- res[contrast %in% c("C2_HPVneg", "TCGA_HPVneg"), .(collection, set, contrast, leadingEdge_genes, qpcr_in_LE, eligible_in_LE)]
le_w <- dcast(le, collection + set ~ contrast, value.var = c("leadingEdge_genes", "qpcr_in_LE", "eligible_in_LE"))
w <- merge(w, le_w, by = c("collection", "set"), all.x = TRUE)
## leading-edge genes shared by C2 and TCGA
w[, `:=`(LE_shared = mapply(function(a, b) paste(intersect(strsplit(a, ";")[[1]], strsplit(b, ";")[[1]]), collapse = ";"),
                           fifelse(is.na(leadingEdge_genes_C2_HPVneg), "", leadingEdge_genes_C2_HPVneg),
                           fifelse(is.na(leadingEdge_genes_TCGA_HPVneg), "", leadingEdge_genes_TCGA_HPVneg)))]
w[, `:=`(qpcr_in_LE_shared = sapply(strsplit(LE_shared, ";"), function(x) paste(intersect(x, qpcr), collapse = ";")),
         eligible_in_LE_shared = sapply(strsplit(LE_shared, ";"), function(x) paste(intersect(x, eligible), collapse = ";")))]
setorder(w, collection, -replicated, FDR_C2_HPVneg)
fwrite(w, file.path(tdir, "B_hpvneg_wide.csv"))
fwrite(w[replicated == TRUE], file.path(tdir, "B_replicated_sets.csv"))
print(w[replicated == TRUE, .N, by = .(collection, direction)])
print(w[collection == "H" & replicated == TRUE, .(set, NES_C2 = round(NES_C2_HPVneg, 2), FDR_C2 = signif(FDR_C2_HPVneg, 2),
                                                  NES_TCGA = round(NES_TCGA_HPVneg, 2), FDR_TCGA = signif(FDR_TCGA_HPVneg, 2),
                                                  NES_C1 = round(NES_C1_HPVneg, 2), qpcr_in_LE_shared)])
print(w[collection == "RECUR", .(set, NES_C2 = round(NES_C2_HPVneg, 2), FDR_C2 = signif(FDR_C2_HPVneg, 2),
                                 NES_TCGA = round(NES_TCGA_HPVneg, 2), FDR_TCGA = signif(FDR_TCGA_HPVneg, 2),
                                 NES_C1 = round(NES_C1_HPVneg, 2), FDR_C1 = signif(FDR_C1_HPVneg, 2))])
