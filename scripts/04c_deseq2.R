## Step 04c (part A-B): DESeq2 sensitivity analysis of the step 04b pseudobulk DE
##
## Purpose: rebuild the step 04b C1, C2 (HPV- and HPV+) and TCGA HPV- paired
##   count sets exactly, confirm edgeR QL reproduces the saved step 04b
##   tables, fit DESeq2 (Wald, BH; unshrunken and apeglm-shrunken log2FC) on
##   the same genes and designs, and compare with edgeR: log2FC agreement,
##   significance overlap, the candidate rule re-applied, the 24 atlas genes.
## Inputs:  data/processed/v2/{pseudobulk_step04b, tcga_hnsc_paired_counts}.rds,
##          results/step04b_DE/tables/{D_pseudobulk_samples, D_C1_*, D_C2_*, E_TCGA_HPVneg_paired_tumour_vs_normal,
##          D_candidates_HPVneg}.csv, data/raw/step01_metadata/cbioportal_hnsc_tcga_pan_can_atlas_2018_patient.csv,
##          results/step08_primers/tables/A_final_gene_list.csv
## Outputs: results/step04c_deseq2/tables/{A_*, B_*}.csv
## Date: 2026-10-28

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(edgeR)
  library(DESeq2)
})
source("config/paths.R")
source("scripts/lib_de.R")

tdir <- file.path(paths$step04c, "tables")
dir.create(tdir, recursive = TRUE, showWarnings = FALSE)
rd4 <- function(f) fread(file.path(paths$step04b, "tables", f))
min_mal <- 20; min_norm <- 10; fdr_cut <- 0.05; lfc_cut <- 1; excl_normal <- "P9"
free <- function() invisible(gc(full = TRUE))

## =========================================================================
## A1. Rebuild the step 04b sets (make_set logic) and TCGA pairs
## =========================================================================
pb <- readRDS(paths$pseudobulk_v2b)
sm_sh <- rd4("D_pseudobulk_samples.csv")
stopifnot(setequal(sm_sh$key, colnames(pb$shared$counts)))
make_set <- function(hpv_arm, paired) {
  s <- sm_sh[hpv == hpv_arm & primary_4b == TRUE]
  s <- s[!(group == "normal" & patient == excl_normal)]
  agg <- function(g, state, mn) s[group == g, .(keys = list(key), n_cells = sum(n_cells)), by = .(patient, study)][n_cells >= mn][, state := state]
  ss <- rbind(agg("malignant", "malignant", min_mal), agg("normal", "normal", min_norm))
  if (paired) ss <- ss[patient %in% intersect(ss[state == "malignant", patient], ss[state == "normal", patient])]
  ss[, sample := paste(patient, state, sep = "|")]
  cnt <- do.call(cbind, lapply(ss$keys, function(k) Matrix::rowSums(pb$shared$counts[, k, drop = FALSE])))
  colnames(cnt) <- ss$sample
  ss[, `:=`(state = factor(state, levels = c("normal", "malignant")), patient = factor(patient), study = factor(study))]
  list(smp = ss[, !"keys"], counts = cnt)
}
## study_zero: zero counts in every malignant pseudobulk (>= 20 cells) of a study (full matrix)
fk <- data.table(pbkey = colnames(pb$full$counts))
fk[, `:=`(patient = sub("\\|.*$", "", pbkey), group = sub("^.*\\|", "", pbkey), n_cells = as.integer(pb$full$n_cells[pbkey]))]
fk[, study := sm_sh$study[match(patient, sm_sh$patient)]]
mz <- fk[group == "malignant" & n_cells >= min_mal & !is.na(study)]
st_sum <- sapply(split(mz$pbkey, mz$study), function(k) Matrix::rowSums(pb$full$counts[, k, drop = FALSE]))
study_zero <- apply(st_sum == 0, 1, function(z) paste(colnames(st_sum)[z], collapse = " "))

sets <- list(C1_HPVneg = make_set("HPV-", TRUE), C2_HPVneg = make_set("HPV-", FALSE),
             C1_HPVpos = make_set("HPV+", TRUE), C2_HPVpos = make_set("HPV+", FALSE))
rm(pb); free()
tc <- readRDS(paths$tcga_counts)
cb <- fread(paths$cbio_pancan, select = c("patientId", "SUBTYPE"))
ts <- copy(tc$samples)[, `:=`(hpv = cb$SUBTYPE[match(case, cb$patientId)],
                              state = factor(fifelse(sample_type == "Primary Tumor", "tumour", "normal"), levels = c("normal", "tumour")))]
ts <- ts[hpv == "HNSC_HPV-"][, patient := factor(case)]
sets$TCGA_HPVneg <- list(smp = ts, counts = tc$counts[, ts$sample])
rm(tc); free()
design <- list(C1_HPVneg = ~ patient + state, C2_HPVneg = ~ study + state, C1_HPVpos = ~ patient + state,
               C2_HPVpos = ~ study + state, TCGA_HPVneg = ~ patient + state)
coef_e <- c(C1_HPVneg = "statemalignant", C2_HPVneg = "statemalignant", C1_HPVpos = "statemalignant",
            C2_HPVpos = "statemalignant", TCGA_HPVneg = "statetumour")
saved <- c(C1_HPVneg = "D_C1_HPVneg_paired.csv", C2_HPVneg = "D_C2_HPVneg_allpatients.csv", C1_HPVpos = "D_C1_HPVpos_paired.csv",
           C2_HPVpos = "D_C2_HPVpos_allpatients.csv", TCGA_HPVneg = "E_TCGA_HPVneg_paired_tumour_vs_normal.csv")

## =========================================================================
## A2. Reconstruction check: edgeR QL on the rebuilt sets == saved step 04b tables
## =========================================================================
edg <- list(); chk <- list()
for (k in names(sets)) {
  e <- de_qlf(sets[[k]]$counts, sets[[k]]$smp, design[[k]], coef_e[[k]], fdr_cut, lfc_cut)$table
  s <- rd4(saved[[k]])
  m <- merge(s[, .(gene, logFC, FDR)], e[, .(gene, logFC_re = logFC, FDR_re = FDR)], by = "gene")
  chk[[k]] <- data.table(contrast = k, samples = ncol(sets[[k]]$counts), genes_saved = nrow(s), genes_rebuilt = nrow(e),
                         max_abs_diff_logFC = max(abs(m$logFC - m$logFC_re)), max_abs_diff_FDR = max(abs(m$FDR - m$FDR_re)))
  edg[[k]] <- s
}
chk <- rbindlist(chk)
fwrite(chk, file.path(tdir, "A_reconstruction_check.csv"))
print(chk)
stopifnot(all(chk$genes_saved == chk$genes_rebuilt), all(chk$max_abs_diff_logFC < 1e-6), all(chk$max_abs_diff_FDR < 1e-6))

## =========================================================================
## A3. DESeq2 on the same genes and designs
## =========================================================================
dres <- list(); filt <- list()
for (k in names(sets)) {
  g <- edg[[k]]$gene
  cnt <- round(as.matrix(sets[[k]]$counts[g, , drop = FALSE]))
  storage.mode(cnt) <- "integer"
  cd <- as.data.frame(sets[[k]]$smp); rownames(cd) <- colnames(cnt)
  dds <- DESeqDataSetFromMatrix(cnt, cd, design[[k]])
  dds <- DESeq(dds, quiet = TRUE)
  cn <- grep("^state_", resultsNames(dds), value = TRUE)
  r <- results(dds, name = cn, alpha = fdr_cut)
  sh <- lfcShrink(dds, coef = cn, type = "apeglm", quiet = TRUE)
  rep_genes <- if ("replace" %in% names(mcols(dds))) sum(mcols(dds)$replace %in% TRUE) else 0L
  d <- data.table(gene = rownames(r), baseMean = r$baseMean, log2FC = r$log2FoldChange, lfcSE = r$lfcSE, stat = r$stat,
                  pvalue = r$pvalue, padj = r$padj, log2FC_apeglm = sh$log2FoldChange,
                  maxCooks = if (!is.null(assays(dds)[["cooks"]])) apply(assays(dds)[["cooks"]], 1, max) else NA_real_)
  d[, `:=`(cooks_outlier = is.na(pvalue) & baseMean > 0, indep_filtered = !is.na(pvalue) & is.na(padj))]
  dres[[k]] <- d
  filt[[k]] <- data.table(contrast = k, genes = nrow(d), cooks_outlier_pvalue_NA = sum(d$cooks_outlier),
                          independent_filtering_padj_NA = sum(d$indep_filtered), cooks_replaced_genes = rep_genes,
                          coefficient = cn)
  message(sprintf("A3. %s: %d genes, Cook's NA %d, filtered %d, replaced %d", k, nrow(d), sum(d$cooks_outlier), sum(d$indep_filtered), rep_genes))
}
fwrite(rbindlist(filt), file.path(tdir, "A_deseq2_filtering_cooks.csv"))
both <- lapply(names(sets), function(k) {
  m <- merge(edg[[k]][, .(gene, edgeR_logFC = logFC, edgeR_logCPM = logCPM, edgeR_FDR = FDR)], dres[[k]], by = "gene")
  m[, contrast := k]
  fwrite(m, file.path(tdir, paste0("A_deseq2_vs_edgeR_", k, ".csv")))
  m
})
names(both) <- names(sets)

## =========================================================================
## B1. log2FC agreement; B2. significance agreement
## =========================================================================
cand <- rd4("D_candidates_HPVneg.csv")[candidate_list == TRUE]
cand_genes <- cand$gene
b1 <- rbindlist(lapply(names(both), function(k) {
  m <- both[[k]]
  f <- function(x) data.table(n = nrow(x), pearson = cor(x$edgeR_logFC, x$log2FC), spearman = cor(x$edgeR_logFC, x$log2FC, method = "spearman"))
  rbind(cbind(contrast = k, genes = "all tested", f(m)),
        if (grepl("HPVneg", k)) cbind(contrast = k, genes = "step 04b candidates", f(m[gene %in% cand_genes])))
}))
fwrite(b1, file.path(tdir, "B1_lfc_agreement.csv"))
print(b1)
b2 <- rbindlist(lapply(names(both), function(k) {
  m <- both[[k]]
  rbindlist(lapply(c("up", "down"), function(d) {
    s <- if (d == "up") 1 else -1
    e <- m[edgeR_FDR < fdr_cut & sign(edgeR_logFC) == s, gene]
    q <- m[padj < fdr_cut & !is.na(padj) & sign(log2FC) == s, gene]
    data.table(contrast = k, direction = d, edgeR_only = length(setdiff(e, q)), both = length(intersect(e, q)),
               DESeq2_only = length(setdiff(q, e)), edgeR_total = length(e), DESeq2_total = length(q),
               jaccard = if (length(union(e, q))) length(intersect(e, q)) / length(union(e, q)) else NA_real_)
  }))
}))
fwrite(b2, file.path(tdir, "B2_significance_agreement.csv"))
print(b2)

## =========================================================================
## B3. Candidate rule re-applied with DESeq2 (HPV-)
## =========================================================================
c2 <- both$C2_HPVneg; c1 <- both$C1_HPVneg
dc <- merge(c2[, .(gene, C2_log2FC = log2FC, C2_padj = padj, C2_baseMean = baseMean, C2_cooks = cooks_outlier, C2_filtered = indep_filtered,
                   C2_edgeR_logFC = edgeR_logFC, C2_edgeR_FDR = edgeR_FDR)],
            c1[, .(gene, C1_log2FC = log2FC, C1_edgeR_logFC = edgeR_logFC)], by = "gene", all.x = TRUE)
dc[, C2_dir := fifelse(!is.na(C2_padj) & C2_padj < fdr_cut & abs(C2_log2FC) >= lfc_cut, fifelse(C2_log2FC > 0, "up", "down"), "ns")]
dc[, candidate := fcase(C2_dir == "up" & C1_log2FC > 0, "up", C2_dir == "down" & C1_log2FC < 0, "down", default = "")]
dc[is.na(candidate), candidate := ""]
dc[, study_zero := study_zero[gene]][is.na(study_zero), study_zero := ""]
dc[, candidate_list := candidate != "" & study_zero == ""]
cmp <- merge(cand[, .(gene, edgeR_candidate = candidate)], dc[candidate_list == TRUE, .(gene, DESeq2_candidate = candidate)], by = "gene", all = TRUE)
cmp[, status := fcase(is.na(DESeq2_candidate), "lost in DESeq2", is.na(edgeR_candidate), "gained in DESeq2",
                      edgeR_candidate == DESeq2_candidate, "recovered", default = "direction differs")]
cmp <- merge(cmp, dc[, .(gene, C2_log2FC, C2_padj, C2_baseMean, C2_cooks, C2_filtered, C1_log2FC, C2_edgeR_logFC, C2_edgeR_FDR, C1_edgeR_logFC)],
             by = "gene", all.x = TRUE)
cmp[status == "lost in DESeq2", reason := fcase(C2_filtered %in% TRUE, "C2 removed by DESeq2 independent filtering (padj NA)",
                                                C2_cooks %in% TRUE, "C2 Cook's outlier (pvalue NA)",
                                                !(C2_padj < fdr_cut) %in% TRUE & abs(C2_log2FC) >= lfc_cut, "C2 padj >= 0.05",
                                                C2_padj < fdr_cut & abs(C2_log2FC) < lfc_cut, "C2 |log2FC| < 1 in DESeq2",
                                                !(C2_padj < fdr_cut) %in% TRUE & abs(C2_log2FC) < lfc_cut, "C2 padj >= 0.05 and |log2FC| < 1",
                                                default = "C1 sign differs")]
cmp[status == "gained in DESeq2", reason := fcase(!(C2_edgeR_FDR < fdr_cut), "edgeR C2 FDR >= 0.05",
                                                  abs(C2_edgeR_logFC) < lfc_cut, "edgeR C2 |logFC| < 1", default = "other")]
fwrite(cmp[order(status, gene)], file.path(tdir, "B3_candidates_compared.csv"))
fwrite(cmp[status != "recovered", .N, by = .(status, direction = fifelse(is.na(edgeR_candidate), DESeq2_candidate, edgeR_candidate), reason)][order(status, direction)],
       file.path(tdir, "B3_reasons.csv"))
b3 <- rbind(cmp[!is.na(edgeR_candidate), .(edgeR = .N, recovered = sum(status == "recovered"), lost = sum(status == "lost in DESeq2")), by = .(direction = edgeR_candidate)],
            fill = TRUE)
b3 <- merge(b3, cmp[status == "gained in DESeq2", .(gained = .N), by = .(direction = DESeq2_candidate)], by = "direction", all = TRUE)
fwrite(b3, file.path(tdir, "B3_candidate_counts.csv"))
print(b3)
print(fread(file.path(tdir, "B3_reasons.csv")))
print(cmp[status == "lost in DESeq2", .(gene, status, reason, C2_edgeR_logFC = round(C2_edgeR_logFC, 2), C2_edgeR_FDR = signif(C2_edgeR_FDR, 2),
                                    C2_log2FC = round(C2_log2FC, 2), C2_padj = signif(C2_padj, 2), baseMean = round(C2_baseMean, 1),
                                    cooks = C2_cooks, filtered = C2_filtered)])

## =========================================================================
## B4-B5. The 24 atlas genes (step 08 list) side by side
## =========================================================================
gl <- fread(file.path(paths$step08, "tables", "A_final_gene_list.csv"))[, .(gene, group)]
b4 <- rbindlist(lapply(c("C2_HPVneg", "C1_HPVneg", "TCGA_HPVneg"), function(k) {
  m <- both[[k]][gene %in% gl$gene]
  m[, .(contrast = k, gene, edgeR_logFC, edgeR_FDR, DESeq2_log2FC = log2FC, DESeq2_log2FC_apeglm = log2FC_apeglm, DESeq2_padj = padj,
        baseMean, cooks_outlier, indep_filtered)]
}))
b4 <- merge(gl, b4, by = "gene", all.x = TRUE, sort = FALSE)
b4[, `:=`(edgeR_sig = edgeR_FDR < fdr_cut, DESeq2_sig = !is.na(DESeq2_padj) & DESeq2_padj < fdr_cut)]
b4[, `:=`(direction_changes = sign(edgeR_logFC) != sign(DESeq2_log2FC), significance_changes = edgeR_sig != DESeq2_sig)]
b4[, flag := fcase(direction_changes & significance_changes, "direction and significance change",
                   direction_changes, "direction changes", significance_changes, "significance changes", default = "")]
fwrite(b4, file.path(tdir, "B4_atlas_genes.csv"))
print(b4[, .(gene, contrast, edgeR = round(edgeR_logFC, 2), eFDR = signif(edgeR_FDR, 2), DESeq2 = round(DESeq2_log2FC, 2),
             dFDR = signif(DESeq2_padj, 2), flag, group)][flag != "" | group == "pilot reference"], nrows = 100)
