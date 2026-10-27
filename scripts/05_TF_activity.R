## Step 05 (R parts): final candidate list, export for pySCENIC, TF activity
## (decoupleR + CollecTRI), convergence of evidence per TF, figures
##
## Purpose: finalise the step 04b candidate tiers (adds tier 1b) and flag TFs;
##   export HPV- malignant and squamous normal cells (raw counts, HGNC
##   symbols) for pySCENIC; infer TF activity with ULM at contrast and
##   patient level; combine DE, TCGA, TF activity, regulon and specificity
##   evidence per TF.
## Inputs:  results/step04b_DE/tables/ (DE, TCGA, candidate tables)
##          data/processed/v2/{cell_labels_v2.csv.gz, choi_nl_lp_reannotated.csv.gz,
##          atlas_epithelial_labelled.rds, choi_nl_lp_allcells.rds,
##          pseudobulk_step04b.rds, celltype_pseudobulk_v2.rds}
##          data/reference/collectri/collectri_omnipath.tsv, HGNC complete set,
##          Lambert et al. 2018 human TF list (downloaded)
##          data/processed/v2/step05_scenic/ (outputs of scripts/05_pyscenic.py)
## Outputs: results/step05_TF/ (tables, figures)
##          data/processed/v2/step05_export/ (matrix.mtx, genes.csv, cells.csv)
## Date: 2026-10-27
## Run order (from the repo root):
##   1. Rscript scripts/05_TF_activity.R      (parts A-C and the export)
##   2. envs/py310/bin/python scripts/05_pyscenic.py
##   3. Rscript scripts/05_TF_activity.R      (again: part E and figures)

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(edgeR)
  library(decoupleR)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})
source("config/paths.R")
source("scripts/lib_de.R")
out <- paths$step05
tdir <- file.path(out, "tables")
dir.create(tdir, recursive = TRUE, showWarnings = FALSE)
free <- function() invisible(gc(full = TRUE))
set.seed(1)
t04 <- file.path(paths$step04b, "tables")
rd <- function(f) fread(file.path(t04, f))

fdr_cut  <- 0.05
lfc_cut  <- 1
min_cpm  <- 10
min_type <- 20
min_mal  <- 20
min_norm <- 10
excl_normal <- "P9"
not_ranked  <- c("TRIB3", "MTERF3")   # candidate genes, but not DNA-binding TFs: left out of the TF ranking
tme_types <- c("CD4T", "CD8T", "NK", "Bcell", "PlasmaCell", "Macrophage", "Monocyte",
               "panDC", "Neutrophils", "Mast", "Fibroblast", "Endothelial")
pos_up   <- c("MAGEA4", "PRAME", "MMP1", "LAMC2", "PDPN")
pos_down <- c("KRT4", "KRT13", "CRNN")
check_tfs <- c("HOXD11", "HOXB7", "HOXA10", "HOXC13", "HOXC9", "TP63", "SOX2", "KLF5", "TFAP2A",
               "FOSL1", "HIF1A", "EPAS1", "MYC", "E2F1", "FOXM1")

hg <- read_hgnc(paths$hgnc)

## =========================================================================
## A. Final candidate list (tiers 1, 1b, 2) and TF flags
## =========================================================================
## TF lists: CollecTRI sources (complex members split) + Lambert et al. 2018
if (!file.exists(paths$lambert_tfs)) {
  dir.create(dirname(paths$lambert_tfs), recursive = TRUE, showWarnings = FALSE)
  download.file("http://humantfs.ccbr.utoronto.ca/download/v_1.01/TF_names_v_1.01.txt", paths$lambert_tfs, quiet = TRUE)
}
lambert <- unique(na.omit(hgnc_map(readLines(paths$lambert_tfs), hg)$mapped))
raw_net <- fread(paths$collectri)
collectri_tfs <- unique(na.omit(hgnc_map(unique(unlist(strsplit(raw_net$source_genesymbol, "_"))), hg)$mapped))
tf_all <- union(lambert, collectri_tfs)

## specificity and detectability for any gene (method of step 04b part F)
pb <- readRDS(paths$pseudobulk_v2b)
smp <- rd("D_pseudobulk_samples.csv")
map_atlas <- hgnc_map(rownames(readRDS(paths$atlas_labelled)$counts), hg)
free()
ct <- readRDS(paths$celltype_pb_v2)
ctc <- map_rows(ct$counts, map_atlas)
## HPV- primary patients as in step 04b (all discovery patients in the labels)
lab_pt <- fread(paths$cell_labels_v2, select = c("patient_id_unique", "hpv", "role", "flag_excluded", "primary_4b"))
neg_pts <- unique(lab_pt[role == "discovery" & !flag_excluded & hpv == "HPV-" & primary_4b == TRUE, patient_id_unique])
rm(lab_pt)
cs <- ct$samples[patient %in% neg_pts]
lib_ct <- colSums(ctc)
specificity_of <- function(genes) {
  genes <- intersect(genes, rownames(ctc))
  cp <- as.matrix(ctc[genes, cs$key, drop = FALSE]) %*% diag(1e6 / lib_ct[cs$key], nrow = nrow(cs))
  colnames(cp) <- cs$key
  med <- function(g) { j <- cs[group == g & n_cells >= min_type, key]; apply(cp[, j, drop = FALSE], 1, median) }
  tum <- med("malignant"); tme <- sapply(tme_types, med)
  mx <- apply(tme, 1, max); mxt <- tme_types[apply(tme, 1, which.max)]
  sp <- sapply(cs[group == "malignant" & n_cells >= min_type, patient], function(p) {
    tot <- Matrix::rowSums(ctc[genes, cs[patient == p, key], drop = FALSE])
    ifelse(tot > 0, ctc[genes, cs[patient == p & group == "malignant", key]] / tot, NA_real_)
  })
  mal_k <- smp[group == "malignant" & hpv == "HPV-" & primary_4b == TRUE & n_cells >= min_mal]
  nor_k <- smp[group == "normal" & hpv == "HPV-" & primary_4b == TRUE & n_cells >= min_norm & patient != excl_normal]
  med_pb <- function(k, m) apply(as.matrix(m[genes, k$key, drop = FALSE]) %*% diag(1e6 / colSums(m[, k$key, drop = FALSE]), nrow = nrow(k)), 1, median)
  g_sh <- intersect(genes, rownames(pb$shared$counts))
  nor_med <- setNames(rep(NA_real_, length(genes)), genes)
  nor_med[g_sh] <- apply(as.matrix(pb$shared$counts[g_sh, nor_k$key, drop = FALSE]) %*%
                           diag(1e6 / colSums(pb$shared$counts[, nor_k$key, drop = FALSE]), nrow = nrow(nor_k)), 1, median)
  d <- data.table(gene = genes, malignant_median_CPM = tum, max_TME_CPM = mx, max_TME_type = mxt,
                  fold_vs_max_TME = tum / mx, tumour_share = apply(sp, 1, median, na.rm = TRUE),
                  malignant_pb_median_CPM = med_pb(mal_k, pb$full$counts), normal_pb_median_CPM = nor_med[genes])
  d[, specificity := fcase(malignant_median_CPM >= 2 * max_TME_CPM & tumour_share >= 0.5, "tumour_specific",
                           malignant_median_CPM > max_TME_CPM, "tumour_enriched", default = "shared_with_TME")]
  d[]
}

cand <- rd("D_candidates_HPVneg.csv")
tcga <- rd("E_TCGA_HPVneg_paired_tumour_vs_normal.csv")[, .(gene, TCGA_logFC = logFC, TCGA_FDR = FDR)]
f04 <- rd("F_final_candidates_HPVneg.csv")
cf <- merge(cand, tcga, by = "gene", all.x = TRUE)
cf <- cf[candidate != "" | C2_sig_not_tested_C1 == TRUE]
cf[, direction := fifelse(C2_logFC > 0, "up", "down")]
cf[, TCGA_replicated := !is.na(TCGA_FDR) & TCGA_FDR < fdr_cut & sign(TCGA_logFC) == sign(C2_logFC)]
cf <- merge(cf, specificity_of(cf$gene), by = "gene", all.x = TRUE)
## detectable: up genes in malignant pseudobulks, down genes in normal pseudobulks
cf[, detectable := fifelse(direction == "up", malignant_pb_median_CPM >= min_cpm, normal_pb_median_CPM >= min_cpm)]
cf[is.na(detectable), detectable := FALSE]
cf[, tier := fcase(
  candidate == "up" & candidate_list & TCGA_replicated & specificity == "tumour_specific" & detectable, "1",
  C2_sig_not_tested_C1 & TCGA_replicated & detectable & study_zero == "", "1b",
  candidate == "up" & candidate_list & TCGA_replicated & specificity == "tumour_enriched" & detectable, "2",
  candidate == "down" & candidate_list & TCGA_replicated, "down_replicated",
  default = "")]
## the tier 1/2 genes must equal step 04b
stopifnot(setequal(cf[tier == "1", gene], f04[tier == "1_specific", gene]),
          setequal(cf[tier == "2", gene], f04[tier == "2_enriched", gene]))
cf[, `:=`(is_TF = gene %in% tf_all, TF_list = fcase(gene %in% lambert & gene %in% collectri_tfs, "Lambert+CollecTRI",
                                                    gene %in% lambert, "Lambert", gene %in% collectri_tfs, "CollecTRI",
                                                    default = ""))]
## positive controls: rule unchanged; record their C2 direction and TCGA status
pcs <- merge(data.table(gene = c(pos_up, pos_down), expected = rep(c("up", "down"), c(5, 3))),
             merge(rd("D_C2_HPVneg_allpatients.csv")[, .(gene, C2_logFC = logFC, C2_FDR = FDR)], tcga, by = "gene", all = TRUE),
             by = "gene", all.x = TRUE)
pcs[, note := sprintf("positive control: C2 %s (log2FC %.2f, FDR %.2g); TCGA %s; not a candidate under the unchanged rule (7 normals)",
                      fifelse(sign(C2_logFC) == fifelse(expected == "up", 1, -1), "right direction", "WRONG direction"),
                      C2_logFC, C2_FDR,
                      fifelse(!is.na(TCGA_FDR) & TCGA_FDR < fdr_cut & sign(TCGA_logFC) == fifelse(expected == "up", 1, -1),
                              sprintf("replicated (log2FC %.2f, FDR %.2g)", TCGA_logFC, TCGA_FDR), "not tested / not replicated"))]
cf[, positive_control_note := pcs$note[match(gene, pcs$gene)]]
pc_extra <- pcs[!gene %in% cf$gene, .(gene, C2_logFC, C2_FDR, TCGA_logFC, TCGA_FDR, positive_control_note = note)]
cf <- rbind(cf, pc_extra, fill = TRUE)
cf[is.na(tier), tier := ""]
setorder(cf, tier, C2_FDR)
cf <- rbind(cf[tier %in% c("1", "1b", "2")][order(factor(tier, levels = c("1", "1b", "2")), C2_FDR)], cf[!tier %in% c("1", "1b", "2")])
fwrite(cf, file.path(tdir, "A_candidates_final.csv"))
tier_sum <- cf[tier != "", .(genes = .N, TFs = sum(is_TF %in% TRUE)), by = tier]
fwrite(tier_sum, file.path(tdir, "A_tier_summary.csv"))
cand_genes <- cf[tier %in% c("1", "1b", "2"), gene]

## =========================================================================
## B2. Export for pySCENIC (raw counts, HGNC symbols, shared atlas/Choi genes)
## =========================================================================
ex_dir <- paths$step05_export
if (!file.exists(file.path(ex_dir, "matrix.mtx"))) {
  message("B. Export for pySCENIC")
  dir.create(ex_dir, recursive = TRUE, showWarnings = FALSE)
  lab <- fread(paths$cell_labels_v2, select = c("cell", "dataset", "role", "hpv", "patient_id_unique", "study",
                                                "site", "cell_state_4b", "primary_4b"))
  cm <- fread(paths$choi_reannotated, select = c("cell", "patient", "hpv", "tissue", "cell_state_4b"))
  base <- lab[role == "discovery" & hpv == "HPV-" & primary_4b == TRUE]
  mal <- base[dataset == "HNSCatlas" & cell_state_4b == "malignant"]
  nor_a <- base[dataset == "HNSCatlas" & cell_state_4b == "normal_epithelial"]
  nor_c <- cm[tissue == "NL" & hpv == "HPV-" & cell_state_4b == "normal_epithelial" & patient != excl_normal]
  atl <- readRDS(paths$atlas_labelled)
  xa <- map_rows(atl$counts[, c(mal$cell, nor_a$cell)], map_atlas)
  rm(atl); free()
  ch <- readRDS(paths$choi_nl_lp_all)
  xc <- map_rows(ch$counts[, nor_c$cell, drop = FALSE], hgnc_map(rownames(ch$counts), hg))
  rm(ch); free()
  g <- intersect(rownames(xa), rownames(xc))
  x <- cbind(xa[g, ], xc[g, ])
  cells <- rbind(mal[, .(cell, patient = patient_id_unique, state = "malignant", study, source = "atlas")],
                 nor_a[, .(cell, patient = patient_id_unique, state = "normal", study, source = "atlas")],
                 nor_c[, .(cell, patient, state = "normal", study = "Choi2023", source = "Choi NL")])
  stopifnot(identical(colnames(x), cells$cell))
  cells[, n_counts := colSums(x)]
  writeMM(as(x, "CsparseMatrix"), file.path(ex_dir, "matrix.mtx"))
  fwrite(data.table(gene = g), file.path(ex_dir, "genes.csv"))
  fwrite(cells, file.path(ex_dir, "cells.csv"))
  fwrite(data.table(gene = cand_genes), file.path(ex_dir, "candidate_genes.csv"))
  rm(xa, xc, x); free()
}
ex_cells <- fread(file.path(ex_dir, "cells.csv"))
ex_sum <- ex_cells[, .(cells = .N, patients = uniqueN(patient)), by = state]
fwrite(ex_cells[, .N, by = .(state, patient, study)][order(state, -N)], file.path(tdir, "B_export_cells_per_patient.csv"))

## =========================================================================
## C. TF activity: decoupleR ULM with CollecTRI
## =========================================================================
## same processing as decoupleR::get_collectri(split_complexes = FALSE) (v1)
src <- raw_net$source_genesymbol
cplx <- grepl("_", src)
src[cplx & grepl("NFKB|REL", src)] <- "NFKB"
src[cplx & grepl("JUN|FOS|ATF", src)] <- "AP1"
inhibits <- as.logical(raw_net$is_inhibition) & !as.logical(raw_net$is_stimulation)
net <- data.table(source = src, target = raw_net$target_genesymbol, mor = fifelse(inhibits, -1, 1))
net <- net[!grepl("_", source)]
## HGNC symbols for TFs and targets (AP1/NFKB are not HGNC symbols: kept)
nm <- hgnc_map(unique(c(net$source, net$target)), hg)
net[, `:=`(source = nm$mapped[match(source, nm$original)], target = nm$mapped[match(target, nm$original)])]
net <- unique(net[!is.na(source) & !is.na(target)], by = c("source", "target"))
fwrite(data.table(TFs = uniqueN(net$source), interactions = nrow(net)), file.path(tdir, "C_collectri_size.csv"))
ulm <- function(mat) as.data.table(run_ulm(mat = mat, net = net, .source = "source", .target = "target",
                                           .mor = "mor", minsize = 5))
n_targets <- net[, .N, by = source]

## C1. contrast level: signed QL statistic sign(log2FC) * sqrt(F)
stat_of <- function(tab) setNames(sign(tab$logFC) * sqrt(pmax(tab$F, 0)), tab$gene)
contrast_ulm <- function(files) rbindlist(lapply(names(files), function(nm) {
  s <- stat_of(rd(files[[nm]]))
  r <- ulm(matrix(s, ncol = 1, dimnames = list(names(s), nm)))
  r[, FDR := p.adjust(p_value, "BH")]
  r[, .(contrast = condition, TF = source, score, p_value, FDR)]
}))
ca <- contrast_ulm(c(C2_HPVneg = "D_C2_HPVneg_allpatients.csv", C1_HPVneg = "D_C1_HPVneg_paired.csv",
                     C2_HPVpos = "D_C2_HPVpos_allpatients.csv", C1_HPVpos = "D_C1_HPVpos_paired.csv"))
fwrite(ca, file.path(tdir, "C_contrast_activity_all.csv"))
wide_c <- function(arm) {
  w <- dcast(ca[grepl(arm, contrast)], TF ~ contrast, value.var = c("score", "FDR"))
  setnames(w, sub(paste0("_", arm), "", names(w)))
  w[, `:=`(C2_sig = FDR_C2 < fdr_cut, C1_sig = FDR_C1 < fdr_cut,
           agree = sign(score_C2) == sign(score_C1))]
  w[, n_targets := n_targets$N[match(TF, n_targets$source)]]
  w[order(FDR_C2)]
}
cneg <- wide_c("HPVneg"); cpos <- wide_c("HPVpos")
fwrite(cneg, file.path(tdir, "C_contrast_activity_HPVneg.csv"))
fwrite(cpos, file.path(tdir, "C_contrast_activity_HPVpos.csv"))

## C1b. bulk support: ULM on the paired TCGA HPV- tumour vs normal statistic (40 pairs, step 04b);
## bulk tissue includes the microenvironment, so this is support, not a tumour-cell measure
s_tcga <- stat_of(rd("E_TCGA_HPVneg_paired_tumour_vs_normal.csv"))
tca <- ulm(matrix(s_tcga, ncol = 1, dimnames = list(names(s_tcga), "TCGA_HPVneg")))
tca[, FDR := p.adjust(p_value, "BH")]
tca <- tca[, .(TF = source, act_TCGA = score, act_TCGA_p = p_value, act_TCGA_FDR = FDR)]
fwrite(tca[order(act_TCGA_p)], file.path(tdir, "C_TCGA_activity_HPVneg.csv"))

## C2. patient level: ULM on gene z-scores of log CPM per pseudobulk
ps <- smp[hpv == "HPV-" & primary_4b == TRUE &
            ((group == "malignant" & n_cells >= min_mal) | (group == "normal" & n_cells >= min_norm & patient != excl_normal))]
ps[, state := factor(fifelse(group == "malignant", "malignant", "normal"), levels = c("normal", "malignant"))]
y <- DGEList(as.matrix(pb$shared$counts[, ps$key]))
y <- normLibSizes(y[filterByExpr(y, group = ps$state), , keep.lib.sizes = FALSE])
z <- t(scale(t(edgeR::cpm(y, log = TRUE, prior.count = 2))))
act <- ulm(z)
am <- dcast(act, source ~ condition, value.var = "score")
amat <- as.matrix(am[, -1]); rownames(amat) <- am$source
amat <- amat[, ps$key]
pairs <- intersect(ps[state == "malignant", patient], ps[state == "normal", patient])
pk_m <- ps[state == "malignant"][match(pairs, patient), key]
pk_n <- ps[state == "normal"][match(pairs, patient), key]
study <- factor(ps$study)
pl <- rbindlist(lapply(rownames(amat), function(tf) {
  a <- amat[tf, ]
  d <- a[pk_m] - a[pk_n]
  tt <- t.test(d)
  co <- summary(lm(a ~ study + ps$state))$coefficients
  data.table(TF = tf, paired_diff = mean(d), paired_p = tt$p.value,
             all_effect = co["ps$statemalignant", "Estimate"], all_p = co["ps$statemalignant", "Pr(>|t|)"])
}))
pl[, `:=`(paired_FDR = p.adjust(paired_p, "BH"), all_FDR = p.adjust(all_p, "BH"))]
fwrite(pl[order(all_p)], file.path(tdir, "C_patient_activity_HPVneg.csv"))
fwrite(data.table(pseudobulks_malignant = sum(ps$state == "malignant"), pseudobulks_normal = sum(ps$state == "normal"),
                  pairs = length(pairs), pair_patients = paste(pairs, collapse = " "), genes = nrow(y)),
       file.path(tdir, "C_patient_level_design.csv"))

## C3. robust TFs: C2 significant and same direction at patient level (FDR < 0.1)
tfc <- merge(cneg, pl, by = "TF", all = TRUE)
tfc[, robust := C2_sig %in% TRUE & all_FDR < 0.1 & sign(all_effect) == sign(score_C2)]
tfc[is.na(robust), robust := FALSE]
fwrite(tfc[order(-robust, FDR_C2)], file.path(tdir, "C_TF_activity_HPVneg_combined.csv"))

## C4. specific checks: activity and expression
de2 <- rd("D_C2_HPVneg_allpatients.csv"); de1 <- rd("D_C1_HPVneg_paired.csv")
chk <- data.table(TF = unique(c(check_tfs, grep("^HOX[ABCD][0-9]+$", tf_all, value = TRUE))))
chk <- merge(chk, tfc[, .(TF, n_targets, score_C2, FDR_C2, score_C1, FDR_C1, all_effect, all_FDR, paired_FDR, robust)],
             by = "TF", all.x = TRUE)
chk <- merge(chk, de2[, .(TF = gene, expr_C2_logFC = logFC, expr_C2_FDR = FDR)], by = "TF", all.x = TRUE)
chk <- merge(chk, de1[, .(TF = gene, expr_C1_logFC = logFC, expr_C1_FDR = FDR)], by = "TF", all.x = TRUE)
chk <- merge(chk, tcga[, .(TF = gene, TCGA_logFC, TCGA_FDR)], by = "TF", all.x = TRUE)
chk <- merge(chk, tca[, .(TF, act_TCGA, act_TCGA_FDR)], by = "TF", all.x = TRUE)
chk[, in_check_list := TF %in% check_tfs]
chk <- chk[in_check_list | !is.na(expr_C2_logFC) | !is.na(score_C2)]
setorder(chk, -in_check_list, TF)
fwrite(chk, file.path(tdir, "C_specific_checks.csv"))

## =========================================================================
## E. Convergence of evidence per TF (needs scripts/05_pyscenic.py outputs)
## =========================================================================
sc_dir <- paths$step05_scenic
sc_ok <- file.exists(file.path(sc_dir, "regulons_consensus.csv")) && file.exists(file.path(tdir, "D_regulon_activity.csv"))
if (!sc_ok) {
  message("pySCENIC outputs not found: run envs/py310/bin/python scripts/05_pyscenic.py, then this script again")
  quit(save = "no")
}
reg <- fread(file.path(sc_dir, "regulons_consensus.csv"))     # regulon, TF, target, weight (mean), target_runs
ra <- fread(file.path(tdir, "D_regulon_activity.csv"))         # per regulon statistics
stab <- fread(file.path(tdir, "D_regulon_stability.csv"))
pa <- fread(file.path(tdir, "D_regulon_auc_per_patient.csv"))   # patient, state, regulon, mean_auc
## regulon AUC per patient ~ study + state, the same model as the patient-level activity test (C2)
pa[, `:=`(study = factor(study_of(patient)), state = factor(state, levels = c("normal", "malignant")))]
rs <- pa[, {
  co <- summary(lm(mean_auc ~ study + state))$coefficients
  .(study_diff = co["statemalignant", "Estimate"], study_p = co["statemalignant", "Pr(>|t|)"])
}, by = regulon]
rs[, study_FDR := p.adjust(study_p, "BH")]
ra <- merge(ra, rs, by = "regulon", all.x = TRUE)
fwrite(ra[order(study_p), .(regulon, TF, n_targets, runs_found, study_diff, study_p, study_FDR, all_diff, all_FDR,
                            paired_diff, paired_FDR, candidate_targets)],
       file.path(tdir, "E_regulon_study_adjusted.csv"))
tf_univ <- unique(c(tf_all, ra$TF))
ev <- data.table(TF = tf_univ)
ev <- merge(ev, de2[, .(TF = gene, C2_logFC = logFC, C2_FDR = FDR)], by = "TF", all.x = TRUE)
ev <- merge(ev, de1[, .(TF = gene, C1_logFC = logFC, C1_FDR = FDR)], by = "TF", all.x = TRUE)
ev <- merge(ev, tcga[, .(TF = gene, TCGA_logFC, TCGA_FDR)], by = "TF", all.x = TRUE)
ev <- merge(ev, tfc[, .(TF, act_C2 = score_C2, act_C2_FDR = FDR_C2, act_C1 = score_C1, act_patient = all_effect,
                        act_patient_FDR = all_FDR)], by = "TF", all.x = TRUE)
ev <- merge(ev, tca[, .(TF, act_TCGA, act_TCGA_FDR)], by = "TF", all.x = TRUE)
ra_best <- ra[order(TF, study_p)][, .SD[1], by = TF]
ev <- merge(ev, ra_best[, .(TF, regulon, regulon_targets = n_targets, regulon_diff = study_diff, regulon_FDR = study_FDR,
                            regulon_MW_FDR = all_FDR, regulon_paired_FDR = paired_FDR, regulon_runs = runs_found)],
            by = "TF", all.x = TRUE)
sp_tf <- specificity_of(ev$TF)
ev <- merge(ev, sp_tf, by.x = "TF", by.y = "gene", all.x = TRUE)
ev[, detectable := malignant_pb_median_CPM >= min_cpm]
## lines of evidence for a TF activated in malignant cells
ev[, `:=`(
  e_DE = C2_FDR < fdr_cut & C2_logFC >= lfc_cut & C1_logFC > 0,
  e_TCGA = TCGA_FDR < fdr_cut & TCGA_logFC > 0,
  e_act_contrast = act_C2_FDR < fdr_cut & act_C2 > 0,
  e_act_patient = act_patient_FDR < 0.1 & act_patient > 0,
  e_act_TCGA = act_TCGA_FDR < fdr_cut & act_TCGA > 0,
  e_regulon = regulon_FDR < fdr_cut & regulon_diff > 0 & regulon_runs >= 2,
  e_expression = specificity %in% c("tumour_specific", "tumour_enriched") & detectable %in% TRUE)]
ecols <- c("e_DE", "e_TCGA", "e_act_contrast", "e_act_patient", "e_act_TCGA", "e_regulon", "e_expression")
for (e in ecols) set(ev, which(is.na(ev[[e]])), e, FALSE)
ev[, score := rowSums(.SD), .SDcols = ecols]
## the mirror score for TFs lower / less active in malignant cells
ev[, score_down := (C2_FDR < fdr_cut & C2_logFC <= -lfc_cut & C1_logFC < 0) %in% TRUE + (TCGA_FDR < fdr_cut & TCGA_logFC < 0) %in% TRUE +
                   (act_C2_FDR < fdr_cut & act_C2 < 0) %in% TRUE + (act_patient_FDR < 0.1 & act_patient < 0) %in% TRUE +
                   (act_TCGA_FDR < fdr_cut & act_TCGA < 0) %in% TRUE +
                   (regulon_FDR < fdr_cut & regulon_diff < 0 & regulon_runs >= 2) %in% TRUE]
## targets that are tier 1/1b genes (CollecTRI and pySCENIC regulon)
t1 <- cf[tier %in% c("1", "1b"), gene]
ct_t <- net[target %in% t1, .(collectri_tier1_targets = paste(sort(unique(target)), collapse = " ")), by = source]
rg_t <- reg[target %in% t1, .(regulon_tier1_targets = paste(sort(unique(target)), collapse = " ")), by = TF]
ev <- merge(ev, ct_t, by.x = "TF", by.y = "source", all.x = TRUE)
ev <- merge(ev, rg_t, by = "TF", all.x = TRUE)
ev[, `:=`(is_TF_lambert = TF %in% lambert, tier = cf$tier[match(TF, cf$gene)])]
## not DNA-binding TFs: kept as candidate genes (part A), left out of the TF ranking
ev <- ev[!TF %in% not_ranked]
setorder(ev, -score, C2_FDR, act_C2_FDR, na.last = TRUE)
fwrite(ev, file.path(tdir, "E_convergence_TF.csv"))
shortlist <- ev[score >= 4]
fwrite(shortlist, file.path(tdir, "E_master_regulator_shortlist.csv"))
fwrite(ev[order(-score_down, act_C2_FDR)][1:20], file.path(tdir, "E_convergence_TF_down_top20.csv"))

## =========================================================================
## F. Figures
## =========================================================================
base <- 7; ink <- "#0b0b0b"; ink2 <- "#52514e"; rule <- "#c3c2b7"
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

## Fig 22: TF activity, C2 (top 15 up and 15 down by FDR)
fb <- cneg[FDR_C2 < fdr_cut]
fb <- rbind(fb[score_C2 > 0][order(FDR_C2)][1:min(.N, 15)], fb[score_C2 < 0][order(FDR_C2)][1:min(.N, 15)])
fb[, `:=`(TF = factor(TF, levels = TF[order(score_C2)]),
          status = fifelse(TF %in% tfc[robust == TRUE, TF], "robust (patient level agrees)", "contrast level only"))]
f22 <- ggplot(fb, aes(score_C2, TF, fill = score_C2 > 0, alpha = status)) + geom_col(width = 0.7) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = ink2) +
  scale_fill_manual(values = c(`TRUE` = "#e34948", `FALSE` = "#2a78d6"), guide = "none") +
  scale_alpha_manual(values = c(`robust (patient level agrees)` = 1, `contrast level only` = 0.4), name = NULL) +
  labs(title = "TF activity, HPV- malignant vs normal epithelium (C2)",
       subtitle = sprintf("decoupleR ULM on signed QL statistics, CollecTRI; top TFs with FDR < %.2f; faded = not robust at patient level", fdr_cut),
       x = "ULM activity score (t)", y = NULL) + theme_fig()
save_fig(f22, "Fig22_TF_activity_C2", 130, 120)

## Fig 23: regulon AUC per patient (top regulons, malignant vs normal)
top_reg <- c(ra[study_FDR < fdr_cut & study_diff > 0 & runs_found >= 2][order(study_p)][1:min(.N, 20), regulon],
             ra[study_FDR < fdr_cut & study_diff < 0 & runs_found >= 2][order(study_p)][1:min(.N, 10), regulon])
if (length(top_reg)) {
  hm <- pa[regulon %in% top_reg]
  hm[, z := as.vector(scale(mean_auc)), by = regulon]
  ordc <- unique(hm[order(state, patient), paste(state, patient, sep = "|")])
  hm[, `:=`(col = factor(paste(state, patient, sep = "|"), levels = ordc), regulon = factor(regulon, levels = rev(top_reg)))]
  f23 <- ggplot(hm, aes(col, regulon, fill = pmax(pmin(z, 2.5), -2.5))) + geom_tile() +
    scale_fill_gradient2(low = "#2a78d6", mid = "#f0efec", high = "#e34948", name = "z (mean AUC)") +
    scale_x_discrete(labels = function(k) sub("^.*[|]", "", k)) +
    facet_grid(~ state, scales = "free_x", space = "free_x") +
    labs(title = "pySCENIC regulon activity per patient (HPV-)",
         subtitle = "mean AUCell per patient; regulons significant malignant vs normal (~ study + state, FDR < 0.05), found in >= 2 of 3 runs",
         x = NULL, y = NULL) +
    theme_fig() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = base - 2),
                        axis.line = element_blank(), axis.ticks = element_blank())
  save_fig(f23, "Fig23_regulon_auc_heatmap", 180, 30 + 3.2 * length(top_reg))
}

## Fig 24: convergence dot plot (top 25 TFs)
labs_e <- c(e_DE = "DE up (C2, C1)", e_TCGA = "TCGA up", e_act_contrast = "Activity (contrast)",
            e_act_patient = "Activity (patient)", e_act_TCGA = "Activity (TCGA bulk)", e_regulon = "Regulon up (stable)",
            e_expression = "Tumour-enriched, detectable")
top_tf <- ev[score > 0][1:min(.N, 25), TF]
dp <- melt(ev[TF %in% top_tf, c("TF", "score", ecols), with = FALSE], id.vars = c("TF", "score"), variable.name = "evidence", value.name = "met")
dp[, `:=`(TF = factor(TF, levels = rev(top_tf)), evidence = factor(labs_e[as.character(evidence)], levels = labs_e))]
f24 <- ggplot(dp, aes(evidence, TF)) + geom_point(aes(colour = met, size = met)) +
  geom_text(data = unique(dp[, .(TF, score)]), aes(x = length(labs_e) + 0.8, y = TF, label = score), size = 2.2, colour = ink, inherit.aes = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "#e34948", `FALSE` = "#e6e5df"), guide = "none") +
  scale_size_manual(values = c(`TRUE` = 2.6, `FALSE` = 1.2), guide = "none") +
  coord_cartesian(clip = "off") +
  labs(title = "Convergence of evidence: TFs activated in HPV- malignant cells",
       subtitle = "top 25 TFs by number of lines of evidence (number at right)", x = NULL, y = NULL) +
  theme_fig() + theme(axis.text.x = element_text(angle = 35, hjust = 1), plot.margin = margin(5, 20, 5, 5))
save_fig(f24, "Fig24_convergence_dotplot", 130, 130)

## Fig 25: regulon network of the top TF with a stable regulon
top_net <- ev[score >= 1 & !is.na(regulon) & regulon_runs >= 2 & regulon_diff > 0][1, TF]
if (length(top_net) && !is.na(top_net)) {
  rr <- reg[TF == top_net][order(-weight)]
  rr <- rr[1:min(.N, 40)]
  n <- nrow(rr); ang <- seq(0, 2 * pi, length.out = n + 1)[-1]
  nd <- data.table(target = rr$target, x = cos(ang), y = sin(ang), weight = rr$weight,
                   cls = fcase(rr$target %in% cf[tier %in% c("1", "1b"), gene], "tier 1 / 1b",
                               rr$target %in% cf[tier == "2", gene], "tier 2",
                               rr$target %in% tf_all, "TF", default = "other target"))
  f25 <- ggplot() + geom_segment(data = nd, aes(x = 0, y = 0, xend = x, yend = y, linewidth = weight), colour = rule) +
    geom_point(data = nd, aes(x, y, colour = cls), size = 2.2) +
    geom_text_repel(data = nd, aes(x, y, label = target), size = 2, colour = ink, segment.size = 0.2, max.overlaps = 60) +
    annotate("point", x = 0, y = 0, size = 12, colour = "#e34948") +
    annotate("text", x = 0, y = 0, label = top_net, size = 2.4, fontface = "bold", colour = "white") +
    scale_colour_manual(values = c(`tier 1 / 1b` = "#e34948", `tier 2` = "#eda100", TF = "#2a78d6", `other target` = "#898781"), name = NULL) +
    scale_linewidth(range = c(0.2, 1.2), guide = "none") + coord_equal() +
    labs(title = sprintf("%s regulon (pySCENIC consensus of 3 runs)", top_net),
         subtitle = sprintf("%d targets found in >= 2 of 3 runs; top %d by mean GRNBoost2 importance shown", reg[TF == top_net, uniqueN(target)], n)) +
    theme_void(base_size = base) + theme(plot.title = element_text(face = "bold"), plot.subtitle = element_text(colour = ink2))
  save_fig(f25, "Fig25_regulon_network_top_TF", 140, 130)
}
print(tier_sum); print(ex_sum); print(cneg[C2_sig == TRUE][1:20]); print(tfc[robust == TRUE]); print(chk[in_check_list == TRUE])
print(ev[1:25, c("TF", "score", ecols), with = FALSE])
message("Done.")
