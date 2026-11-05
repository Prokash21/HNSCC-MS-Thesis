## Step 07a (part B): TCGA HPV- expression, programme scores, gene annotation
##
## Purpose: build the TCGA HPV- primary tumour expression matrix (log CPM, TMM)
##   and TPM; score every programme gene set (mean z-score); annotate every set
##   and panel gene with malignant pseudobulk CPM and tumour-cell share
##   (step 04b part F method). No outcome data are used here.
## Inputs:  data/raw/step07_bulk/TCGA-HNSC.star_counts.tsv.gz (GDC STAR counts, UCSC Xena GDC hub)
##          data/raw/tcga_hnsc_star/ (GDC STAR files: gene lengths for TPM), cBioPortal patient table
##          data/processed/v2/{pseudobulk_step04b.rds, celltype_pseudobulk_v2.rds, cell_labels_v2.csv.gz}
##          results/step04b_DE/tables/{D_pseudobulk_samples.csv, E_tcga_manifest.csv, F_specificity_annotation.csv}
##          results/step07a_recurrence/tables/A_gene_sets.csv, results/step06_qPCR/tables/D_panel_final.csv
## Outputs: results/step07a_recurrence/tables/B_*.csv
##          data/processed/v2/step07a_tcga_hpvneg.rds (cache: log CPM, TPM, samples)
## Date: 2026-10-28

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(edgeR)
})
source("config/paths.R")
source("scripts/lib_de.R")

tdir <- file.path(paths$step07a, "tables")
min_type <- 20
min_mal  <- 20
tme_types <- c("CD4T", "CD8T", "NK", "Bcell", "PlasmaCell", "Macrophage", "Monocyte",
               "panDC", "Neutrophils", "Mast", "Fibroblast", "Endothelial")
hg <- read_hgnc(paths$hgnc)
gs <- fread(file.path(tdir, "A_gene_sets.csv"))
panel <- fread(file.path(paths$step06, "tables", "D_panel_final.csv"))[set %in% c("core", "reserve", "optional (Arm B)")]

## =========================================================================
## B1. TCGA HPV- primary tumours: counts, TPM, log CPM
## =========================================================================
if (!file.exists(paths$step07a_tcga)) {
  x <- fread(paths$tcga_all_counts)
  ids <- x[[1]]; x[, 1 := NULL]
  smp <- data.table(sample = names(x), case = substr(names(x), 1, 12), type = substr(names(x), 14, 15))
  flow <- data.table(step = "Xena GDC TCGA-HNSC samples", n = nrow(smp))
  smp <- smp[type == "01"]
  flow <- rbind(flow, data.table(step = "primary tumour (-01) samples", n = nrow(smp)))
  smp <- smp[order(case, sample)][, .SD[1], by = case]
  flow <- rbind(flow, data.table(step = "patients (one primary tumour sample each)", n = nrow(smp)))
  cb <- fread(paths$cbio_pancan, select = c("patientId", "SUBTYPE"))
  smp[, hpv := cb$SUBTYPE[match(case, cb$patientId)]]
  smp <- smp[hpv %in% "HNSC_HPV-"]
  flow <- rbind(flow, data.table(step = "HPV- (cBioPortal SUBTYPE)", n = nrow(smp)))
  cnt <- round(2^as.matrix(x[, smp$sample, with = FALSE]) - 1)
  rownames(cnt) <- ids
  rm(x); invisible(gc())

  ## gene lengths from the GDC STAR files: count / tpm is proportional to length
  ## within a sample; per-sample scale removed by the median log ratio
  man <- fread(file.path(paths$step04b, "tables", "E_tcga_manifest.csv"))
  rl <- lapply(man$file_name, function(fn) {
    r <- fread(file.path(paths$tcga_dir, fn), skip = 1)[grepl("^ENSG", gene_id), .(gene_id, unstranded, tpm_unstranded)]
    v <- setNames(ifelse(r$unstranded > 0 & r$tpm_unstranded > 0, log(r$unstranded / r$tpm_unstranded), NA), r$gene_id)
    v - median(v, na.rm = TRUE)
  })
  lm_ <- do.call(cbind, lapply(rl, function(v) v[ids]))
  len <- exp(apply(lm_, 1, median, na.rm = TRUE))
  ok <- is.finite(len)
  tpm_of <- function(m) { r <- m[ok, , drop = FALSE] / len[ok]; t(t(r) / colSums(r)) * 1e6 }
  ## check against the GDC TPM for the tumour samples present in both
  chk <- man[sample %in% smp$sample & sample_type == "Primary Tumor"]
  tpm_chk <- rbindlist(lapply(seq_len(nrow(chk)), function(i) {
    r <- fread(file.path(paths$tcga_dir, chk$file_name[i]), skip = 1)[grepl("^ENSG", gene_id)]
    mine <- tpm_of(cnt[, chk$sample[i], drop = FALSE])[, 1]
    gdc <- r$tpm_unstranded[match(names(mine), r$gene_id)]
    e <- gdc > 1
    data.table(sample = chk$sample[i], genes_compared = sum(e), pearson_log = cor(log(mine[e]), log(gdc[e])),
               median_rel_diff = median(abs(mine[e] - gdc[e]) / gdc[e]),
               frac_within_5pct = mean(abs(mine[e] - gdc[e]) / gdc[e] < 0.05),
               max_rel_diff = max(abs(mine[e] - gdc[e]) / gdc[e]))
  }))
  tpm <- tpm_of(cnt)
  ens <- sub("\\..*$", "", ids)
  sym <- hg$symbol[match(ens, hg$ensembl_gene_id)]
  keep <- !is.na(sym) & !grepl("_PAR_Y$", ids)
  cnt_s <- rowsum(cnt[keep, ], sym[keep])
  keep_t <- !is.na(sym[ok]) & !grepl("_PAR_Y$", ids[ok])
  tpm_s <- rowsum(tpm[keep_t, ], sym[ok][keep_t])
  y <- DGEList(cnt_s)
  fk <- filterByExpr(y)
  y <- normLibSizes(y[fk, , keep.lib.sizes = FALSE])
  lcpm <- edgeR::cpm(y, log = TRUE, prior.count = 1)
  saveRDS(list(lcpm = lcpm, tpm_median = apply(tpm_s, 1, median), samples = smp, flow = flow, tpm_check = tpm_chk,
               genes_before_filter = nrow(cnt_s)), paths$step07a_tcga)
  rm(cnt, cnt_s, tpm, tpm_s, y); invisible(gc())
}
tc <- readRDS(paths$step07a_tcga)
fwrite(tc$flow, file.path(tdir, "B_tcga_expression_flow.csv"))
fwrite(tc$tpm_check, file.path(tdir, "B_tcga_tpm_check.csv"))
message(sprintf("B1. TCGA HPV-: %d tumours, %d genes after filterByExpr (of %d)", ncol(tc$lcpm), nrow(tc$lcpm),
                tc$genes_before_filter))
print(tc$tpm_check[, .(samples = .N, min_pearson_log = min(pearson_log), max_median_rel_diff = max(median_rel_diff),
                       min_frac_within_5pct = min(frac_within_5pct), max_rel_diff = max(max_rel_diff))])

## =========================================================================
## B2. Programme scores in TCGA: mean z-score of log CPM
## =========================================================================
z <- tc$lcpm
z <- (z - rowMeans(z)) / matrixStats::rowSds(z)
sc <- sapply(split(gs$gene, gs$set), function(s) colMeans(z[intersect(s, rownames(z)), , drop = FALSE]))
cov <- gs[, .(genes = .N, in_tcga = sum(gene %in% rownames(z))), by = .(set, programme, role)]
cov[, coverage := round(in_tcga / genes, 3)]
fwrite(data.table(sample = rownames(sc), case = substr(rownames(sc), 1, 12), sc), file.path(tdir, "B_tcga_programme_scores.csv"))
fwrite(cov, file.path(tdir, "B_tcga_set_coverage.csv"))
print(cov)

## =========================================================================
## B3. Gene annotation: malignant pseudobulk CPM, tumour-cell share, TCGA TPM
## =========================================================================
genes <- unique(c(gs$gene, panel$gene))
## malignant pseudobulk CPM: median over step 04b HPV- malignant pseudobulks (>= 20 cells), as in step 06
pb <- readRDS(paths$pseudobulk_v2b)
smp4 <- fread(file.path(paths$step04b, "tables", "D_pseudobulk_samples.csv"))
mal_k <- smp4[group == "malignant" & hpv == "HPV-" & primary_4b == TRUE & n_cells >= min_mal]
m <- pb$full$counts
g_in <- intersect(genes, rownames(m))
pb_cpm <- setNames(rep(NA_real_, length(genes)), genes)
pb_cpm[g_in] <- apply(as.matrix(m[g_in, mal_k$key, drop = FALSE]) %*% diag(1e6 / colSums(m[, mal_k$key, drop = FALSE]),
                                                                       nrow = nrow(mal_k)), 1, median)
rm(pb, m); invisible(gc())
## tumour-cell share: step 04b part F (median over HPV- patients of malignant / all cell-type counts)
ct <- readRDS(paths$celltype_pb_v2)
ctc <- map_rows(ct$counts, hgnc_map(rownames(ct$counts), hg))
lab_pt <- fread(paths$cell_labels_v2, select = c("patient_id_unique", "hpv", "role", "flag_excluded", "primary_4b"))
neg_pts <- unique(lab_pt[role == "discovery" & !flag_excluded & hpv == "HPV-" & primary_4b == TRUE, patient_id_unique])
cs <- ct$samples[patient %in% neg_pts]
g_ct <- intersect(genes, rownames(ctc))
lib_ct <- colSums(ctc)
cpm_ct <- as.matrix(ctc[g_ct, cs$key, drop = FALSE]) %*% diag(1e6 / lib_ct[cs$key], nrow = nrow(cs))
colnames(cpm_ct) <- cs$key
med_type <- function(g) { j <- cs[group == g & n_cells >= min_type, key]; apply(cpm_ct[, j, drop = FALSE], 1, median) }
tum <- med_type("malignant")
tme <- sapply(tme_types, med_type)
share_pts <- cs[group == "malignant" & n_cells >= min_type, patient]
share <- sapply(share_pts, function(p) {
  tot <- Matrix::rowSums(ctc[g_ct, cs[patient == p, key], drop = FALSE])
  ifelse(tot > 0, ctc[g_ct, cs[patient == p & group == "malignant", key]] / tot, NA_real_)
})
ann <- data.table(gene = genes, malignant_pb_CPM = pb_cpm[genes])
ann[, `:=`(malignant_median_CPM = tum[gene], max_TME_CPM = apply(tme, 1, max)[gene],
           max_TME_type = tme_types[apply(tme, 1, which.max)][match(gene, g_ct)],
           tumour_share = apply(share, 1, median, na.rm = TRUE)[gene],
           tcga_median_TPM = tc$tpm_median[gene], in_tcga_lcpm = gene %in% rownames(tc$lcpm),
           panel_set = panel$set[match(gene, panel$gene)])]
ann <- merge(ann, gs[, .(sets = paste(set, collapse = ";"), programmes = paste(unique(programme), collapse = ";")), by = gene],
             by = "gene", all.x = TRUE)
## check against step 04b part F for shared genes
f4 <- fread(file.path(paths$step04b, "tables", "F_specificity_annotation.csv"))
chk <- merge(ann[, .(gene, tumour_share)], f4[, .(gene, f_share = tumour_share)], by = "gene")
message(sprintf("B3. tumour share vs step 04b part F: %d shared genes, max abs difference %.2g", nrow(chk),
                max(abs(chk$tumour_share - chk$f_share), na.rm = TRUE)))
fwrite(ann, file.path(tdir, "B_gene_annotation.csv"))
fwrite(data.table(malignant_pseudobulks = nrow(mal_k), share_patients = length(share_pts), celltype_patients = uniqueN(cs$patient)),
       file.path(tdir, "B_annotation_patients.csv"))
