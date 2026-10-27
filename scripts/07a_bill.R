## Step 07a (part C): Bill 2023 HPV- relapse vs primary (malignant cells)
##
## Purpose: score the programme gene sets in Bill 2023 HPV- malignant cells
##   (Tirosh 2016 single-cell score, patient mean; pseudobulk mean z-score),
##   test relapse vs primary per set (Wilcoxon), and run malignant-cell
##   pseudobulk DE relapse vs primary (edgeR QL, ~ site + group).
## Inputs:  data/processed/v2/step05_export/ (HPV- primary malignant cells, raw counts, HGNC symbols)
##          data/processed/v2/{bill_recurrence_epithelial_labelled.rds, cell_labels_v2.csv.gz}, HGNC complete set
##          results/step07a_recurrence/tables/A_gene_sets.csv
## Outputs: results/step07a_recurrence/tables/C_*.csv
##          data/processed/v2/step07a_bill_cells.rds (cache: Bill malignant cell counts and cell table)
## Date: 2026-10-28

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(edgeR)
})
source("config/paths.R")
source("scripts/lib_de.R")

tdir <- file.path(paths$step07a, "tables")
min_cells <- 20
n_bins <- 30
n_ctrl <- 100
set.seed(1)
hg <- read_hgnc(paths$hgnc)
gs <- fread(file.path(tdir, "A_gene_sets.csv"))
free <- function() invisible(gc(full = TRUE))

## =========================================================================
## C1. Bill malignant cells: primaries (step 05 export) and relapses
## =========================================================================
if (!file.exists(paths$step07a_bill)) {
  ex <- paths$step05_export
  g <- fread(file.path(ex, "genes.csv"))$gene
  ec <- fread(file.path(ex, "cells.csv"))
  x <- readMM(file.path(ex, "matrix.mtx"))
  dimnames(x) <- list(g, ec$cell)
  prim <- ec[state == "malignant" & study == "Bill2023"]
  xp <- as(x[, prim$cell], "CsparseMatrix")
  rm(x); free()
  lab <- fread(paths$cell_labels_v2, select = c("cell", "patient_id_unique", "role", "hpv", "site", "sample_type",
                                                "cell_state_4b", "primary_4b", "unit"))
  rel <- lab[role == "recurrence" & hpv == "HPV-" & primary_4b == TRUE & cell_state_4b == "malignant"]
  b <- readRDS(paths$bill_labelled)
  xb <- map_rows(b$counts[, rel$cell], hgnc_map(rownames(b$counts), hg))
  rm(b); free()
  gi <- intersect(g, rownames(xb))
  xb2 <- Matrix(0, length(g), ncol(xb), sparse = TRUE, dimnames = list(g, colnames(xb)))
  xb2[gi, ] <- xb[gi, ]
  cells <- rbind(prim[, .(cell, patient, group = "primary", sample_type = "primary", sample = patient)],
                 rel[, .(cell, patient = patient_id_unique, group = "relapse", sample_type, sample = unit)])
  cells[, site := lab$site[match(cell, lab$cell)]]
  counts <- cbind(xp, xb2)
  stopifnot(identical(colnames(counts), cells$cell))
  saveRDS(list(counts = counts, cells = cells, genes_missing_in_bill_relapse = setdiff(g, rownames(xb))), paths$step07a_bill)
  rm(xp, xb, xb2, counts, lab); free()
}
bc <- readRDS(paths$step07a_bill)
cells <- bc$cells
## patient-level table; HN23 and HN45 (two recurrences of one patient) are pooled as BHN23
pt <- cells[, .(n_cells = .N, samples = paste(sort(unique(sample)), collapse = "+"),
                sample_type = paste(sort(unique(sample_type)), collapse = "+"),
                site = paste(sort(unique(na.omit(site[site != ""]))), collapse = "+")), by = .(patient, group)]
pt[, included := n_cells >= min_cells]
pt[, site2 := factor(fifelse(site == "Oral cavity", "oral", "other_or_unknown"), c("oral", "other_or_unknown"))]
pt[, locoregional := group == "primary" | sample_type == "recurrence"]
fwrite(pt[order(group, -n_cells)], file.path(tdir, "C_bill_patients.csv"))
pt <- pt[included == TRUE]
pt[, group := factor(group, c("primary", "relapse"))]
message(sprintf("C1. Bill patients with >= %d malignant cells: %d primary, %d relapse (%d locoregional recurrence, %d metastasis)",
                min_cells, sum(pt$group == "primary"), sum(pt$group == "relapse"),
                sum(pt$group == "relapse" & pt$locoregional), sum(pt$group == "relapse" & !pt$locoregional)))
keep_cells <- cells$patient %in% pt$patient
counts <- bc$counts[, keep_cells]
cells <- cells[keep_cells]

## =========================================================================
## C2. Programme scores
## =========================================================================
## single cells (Tirosh 2016): relative expression = log2(CPM/10 + 1) minus the
## gene mean over scored cells; score = set mean minus control-set mean
## (100 controls per set gene from the same of 30 expression bins); computed
## on the sparse matrix: mean(x - mu) over genes = mean(x) - mean(mu)
ln <- counts %*% Diagonal(x = 1e6 / colSums(counts))
ln@x <- log2(ln@x / 10 + 1)
dimnames(ln) <- dimnames(counts)
mu <- rowMeans(ln)
expd <- names(mu)[mu > 0]
bins <- setNames(cut(rank(mu[expd], ties.method = "first"), n_bins, labels = FALSE), expd)
set_score <- function(s) colMeans(ln[s, , drop = FALSE]) - mean(mu[s])
sets <- split(gs$gene, gs$set)
sc_cell <- sapply(names(sets), function(k) {
  s <- intersect(sets[[k]], expd)
  ctrl <- unique(unlist(lapply(s, function(g) sample(names(bins)[bins == bins[g]], n_ctrl, replace = TRUE))))
  set_score(s) - set_score(ctrl)
})
cov_sc <- data.table(set = names(sets), genes = lengths(sets), scored_single_cell = sapply(sets, function(s) length(intersect(s, expd))))
sc_pt <- as.data.table(sc_cell)[, patient := cells$patient][, lapply(.SD, mean), by = patient]

## pseudobulk: mean z-score of log CPM (TMM, prior count 2) across included patients
pbc <- pseudobulk(counts, cells[, .(cell, key = patient)])
y <- DGEList(as.matrix(pbc$counts))
y <- normLibSizes(y[filterByExpr(y, group = pt$group[match(colnames(y), pt$patient)]), , keep.lib.sizes = FALSE])
lc <- edgeR::cpm(y, log = TRUE, prior.count = 2)
zz <- (lc - rowMeans(lc)) / matrixStats::rowSds(lc)
sc_pb <- sapply(sets, function(s) colMeans(zz[intersect(s, rownames(zz)), , drop = FALSE]))
cov_sc[, scored_pseudobulk := sapply(sets, function(s) length(intersect(s, rownames(zz))))]

scores <- merge(melt(sc_pt, id.vars = "patient", variable.name = "set", value.name = "score_single_cell"),
                melt(data.table(patient = rownames(sc_pb), sc_pb), id.vars = "patient", variable.name = "set",
                     value.name = "score_pseudobulk"), by = c("patient", "set"))
scores <- merge(scores, pt[, .(patient, group, sample_type, site2, locoregional, n_cells)], by = "patient")
scores <- merge(scores, unique(gs[, .(set, programme, role)]), by = "set")
fwrite(scores, file.path(tdir, "C_bill_programme_scores.csv"))
fwrite(cov_sc, file.path(tdir, "C_bill_set_coverage.csv"))

wtest <- function(d, label) d[, {
  r <- score_single_cell[group == "relapse"]; p <- score_single_cell[group == "primary"]
  rb <- score_pseudobulk[group == "relapse"]; pb_ <- score_pseudobulk[group == "primary"]
  .(analysis = label, n_primary = length(p), n_relapse = length(r),
    median_primary = median(p), median_relapse = median(r), diff_median = median(r) - median(p),
    wilcox_p = wilcox.test(r, p, exact = FALSE)$p.value,
    pb_diff_median = median(rb) - median(pb_), pb_wilcox_p = wilcox.test(rb, pb_, exact = FALSE)$p.value)
}, by = .(set, programme, role)][, FDR := p.adjust(wilcox_p, "BH"), by = role][]
tests <- rbind(wtest(scores, "all relapse"), wtest(scores[locoregional == TRUE], "locoregional only"))
fwrite(tests, file.path(tdir, "C_bill_programme_tests.csv"))
print(tests[, .(analysis, set, role, n_primary, n_relapse, diff = round(diff_median, 3), p = signif(wilcox_p, 2),
                FDR = signif(FDR, 2), pb_diff = round(pb_diff_median, 3), pb_p = signif(pb_wilcox_p, 2))])

## =========================================================================
## C3. Pseudobulk DE relapse vs primary: edgeR QL, ~ site + group (raw counts)
## =========================================================================
de_run <- function(d) {
  smp <- d[match(colnames(pbc$counts)[colnames(pbc$counts) %in% d$patient], patient)]
  smp[, state := group]
  de_qlf(pbc$counts[, smp$patient], smp, ~ site2 + group, coef = "grouprelapse", fdr = 0.05, lfc = 1)
}
de_all <- de_run(pt)
de_loc <- de_run(pt[locoregional == TRUE])
de <- merge(de_all$table[, .(gene, logFC, logCPM, PValue, FDR, direction)],
            de_loc$table[, .(gene, logFC_locoregional = logFC, FDR_locoregional = FDR)], by = "gene", all.x = TRUE)
de[, missing_in_bill_relapse_matrix := gene %in% bc$genes_missing_in_bill_relapse]
setorder(de, PValue)
fwrite(de, file.path(tdir, "C_bill_DE.csv"))
de_sum <- data.table(analysis = c("all relapse", "locoregional only"), samples = c(de_all$n_samples, de_loc$n_samples),
                     genes_tested = c(de_all$n_genes_tested, de_loc$n_genes_tested),
                     up = c(sum(de_all$table$direction == "up"), sum(de_loc$table$direction == "up")),
                     down = c(sum(de_all$table$direction == "down"), sum(de_loc$table$direction == "down")),
                     min_FDR = c(min(de_all$table$FDR), min(de_loc$table$FDR)))
fwrite(de_sum, file.path(tdir, "C_bill_DE_counts.csv"))
print(de_sum)
