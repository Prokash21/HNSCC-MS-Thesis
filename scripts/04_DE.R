## Step 04: tumour-cell-specific genes (pseudobulk differential expression)
##
## Purpose: apply the step 03 decisions to the cell labels (CNV call at the
##   95th percentile of reference cells, strict 99th-percentile call kept,
##   label_only cells and prior head-and-neck primaries moved to sensitivity,
##   unique patient IDs, leukoplakia split by CNV); build patient-level
##   pseudobulks of malignant and normal epithelial cells; run edgeR
##   quasi-likelihood DE of malignant vs normal epithelium (paired and
##   all-patient designs, HPV-negative primary, HPV-positive secondary); test
##   tumour-cell specificity against the microenvironment and qRT-PCR
##   detectability; sensitivity re-runs.
## Inputs:  data/processed/v2/cell_labels_v2.csv.gz (step 03)
##          data/processed/v2/cnv/{cnv_units.rds, cnv_genes.csv, doublets.csv.gz}
##          data/processed/v2/atlas_epithelial_labelled.rds,
##          data/processed/v2/choi_nl_lp_epithelial_labelled.rds
##          data/processed/HNSCatlas/celltype_pseudobulk.rds (v1 cell-type pseudobulk)
##          reference-cell inputs of step 03 (only to rebuild the reference quantiles)
## Outputs: results/step04_DE/ (tables, figures)
##          data/processed/v2/cell_labels_v2.csv.gz (step 04 columns added)
##          data/processed/v2/cnv/ref_signal_quantiles.csv
##          data/processed/v2/pseudobulk_step04.rds
## Date: 2026-10-24
## Run from the repo root:  Rscript scripts/04_DE.R

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(edgeR)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})
source("config/paths.R")
source("scripts/lib_cnv.R")
source("scripts/lib_de.R")
out <- paths$step04
tdir <- file.path(out, "tables")
dir.create(tdir, recursive = TRUE, showWarnings = FALSE)
free <- function() invisible(gc(full = TRUE))
set.seed(1)

min_cells     <- 20
fdr_cut       <- 0.05
lfc_cut       <- 1
min_cpm       <- 10      # qRT-PCR detectability: median CPM in malignant pseudobulk
min_share     <- 0.5     # tumour share of the whole-tissue signal
min_fold      <- 2       # malignant CPM >= 2x every microenvironment type
prior_primary <- c("BHN39", "BHN46", "BHN49")
tumour_types  <- c("primary_tumour", "recurrence", "metastasis", "lymph_node")
normal_types  <- c("normal_mucosa", "adjacent_normal")
pos_up        <- c("MAGEA4", "PRAME", "MMP1", "LAMC2", "PDPN")
pos_down      <- c("KRT4", "KRT13", "CRNN")
squamous_krt  <- c("KRT5", "KRT14", "KRT13", "KRT4", "KRT6A", "KRT15", "KRT17", "SFN")
min_squamous  <- 1000    # squamous keratin CPM a normal pseudobulk must reach (check)

## =========================================================================
## A. Relabelling (cached CNV scores; reference quantiles rebuilt, verified)
## =========================================================================
if (!file.exists(paths$ref_quantiles)) {
  message("A. Rebuilding reference-cell signal quantiles")
  units <- readRDS(file.path(paths$cnv_cache, "cnv_units.rds"))$units
  fwrite(ref_signal_quantiles(paths, units), paths$ref_quantiles)
  rm(units); free()
}
rq <- fread(paths$ref_quantiles)
if (!all(rq$match_99)) stop("rebuilt 99th percentile differs from the step 03 threshold in: ",
                            paste(rq[match_99 == FALSE, unit], collapse = ", "))

lab <- fread(paths$cell_labels_v2)
s04_cols <- c("thr_p95", "thr_p99", "cnv_call_p95", "informative_p95", "label_only_p95",
              "cell_state_primary", "cell_state_strict", "patient_id_unique", "study",
              "primary_step04", "sensitivity_reason")
if (length(intersect(s04_cols, names(lab)))) lab[, intersect(s04_cols, names(lab)) := NULL]
for (v in c("unit", "cnv_call", "dbl_class", "bill_patient_id")) set(lab, which(lab[[v]] == ""), v, NA)
lab[, `:=`(thr_p95 = rq$q95[match(unit, rq$unit)], thr_p99 = rq$q99[match(unit, rq$unit)])]
cnv_call_at <- function(signal, corr, thr) fcase(
  is.na(signal), NA_character_,
  is.na(corr), fifelse(signal > thr, "uncertain", "non_malignant"),
  signal > thr & corr > cnv_corr_cut, "malignant",
  signal <= thr & corr <= cnv_corr_cut, "non_malignant",
  default = "uncertain")
## the rebuilt 99th percentile must reproduce every step 03 call
stopifnot(identical(cnv_call_at(lab$signal, lab$corr, lab$thr_p99), fifelse(is.na(lab$signal), NA_character_, lab$cnv_call)))
lab[, cnv_call_p95 := cnv_call_at(signal, corr, thr_p95)]

## informative sample, step 03 rule
informative_rule <- function(call, is_lab) {
  n_mal <- sum(call == "malignant")
  pct_lab <- if (sum(is_lab) >= 20) 100 * mean(call[is_lab] == "malignant") else NA_real_
  n_mal >= 20 & (if (is.na(pct_lab)) 100 * n_mal / length(call) >= 20 else pct_lab >= 50)
}
usum <- lab[!is.na(unit), {
  is_lab <- compartment == "tumour_label" & !suspect_label & sample_type %in% tumour_types
  .(dataset = dataset[1], hpv = hpv[1], scored = .N, thr_p99 = thr_p99[1], thr_p95 = thr_p95[1],
    malignant_p99 = sum(cnv_call == "malignant"), malignant_p95 = sum(cnv_call_p95 == "malignant"),
    uncertain_p99 = sum(cnv_call == "uncertain"), uncertain_p95 = sum(cnv_call_p95 == "uncertain"),
    informative_p99 = informative_rule(cnv_call, is_lab), informative_p95 = informative_rule(cnv_call_p95, is_lab),
    informative_step03 = informative[1])
}, by = unit]
stopifnot(identical(usum$informative_p99, usum$informative_step03))
usum[, informative_step03 := NULL]
usum[, `:=`(pct_malignant_p99 = round(100 * malignant_p99 / scored, 1),
            pct_malignant_p95 = round(100 * malignant_p95 / scored, 1))]
setorder(usum, dataset, hpv, -scored)
fwrite(usum, file.path(tdir, "A_relabel_per_unit.csv"))
lab[, informative_p95 := usum$informative_p95[match(unit, usum$unit)]]

## cell state: step 03 rule with label_only cells and leukoplakia split out
assign_state <- function(call, informative) with(lab, fcase(
  flag_excluded, "excluded",
  dbl_class %in% "doublet", "doublet",
  sample_type == "leukoplakia" & call %in% "malignant", "LP_CNVpos",
  sample_type == "leukoplakia" & call %in% "non_malignant", "LP_CNVneg",
  sample_type == "leukoplakia", "LP_uncertain",
  call %in% "malignant" & sample_type %in% normal_types, "uncertain",
  call %in% "malignant", "malignant",
  compartment == "tumour_label" & !suspect_label & informative %in% FALSE & sample_type %in% tumour_types,
    "malignant_label_only",
  call %in% "non_malignant" & !suspect_label & (sample_type %in% normal_types | squamous_clean), "normal_epithelial",
  default = "uncertain"))
lab[, cell_state_primary := assign_state(cnv_call_p95, informative_p95)]
lab[, cell_state_strict := assign_state(cnv_call, informative)]
lab[, label_only_p95 := cell_state_primary == "malignant_label_only"]
## strict labels must equal step 03 cell_state_v2 apart from the two splits
chk <- lab[, .N, by = .(cell_state_v2, cell_state_strict)]
stopifnot(all(chk[cell_state_v2 != cell_state_strict,
                  (cell_state_v2 == "malignant" & cell_state_strict == "malignant_label_only") |
                  (cell_state_v2 == "leukoplakia_epithelial" & grepl("^LP_", cell_state_strict))]))

## unique patient: Bill Table S1 patient IDs (HN45 = second recurrence of HN23)
lab[, patient_id_unique := fifelse(!is.na(bill_patient_id) & bill_patient_id != "",
                                   paste0("B", bill_patient_id), patient)]
merged_ids <- unique(lab[patient_id_unique != patient, .(patient, patient_id_unique)])
lab[, study := study_of(patient_id_unique)]
lab[dataset == "GSE172577_Peng", study := "Peng"]
lab[, sensitivity_reason := fcase(flag_excluded, "excluded: nasopharynx",
                                  flag_treated, "previously treated (OP8, OP13)",
                                  patient_id_unique %in% prior_primary, "prior head-and-neck primary",
                                  default = "")]
lab[, primary_step04 := sensitivity_reason == ""]
fwrite(lab, paths$cell_labels_v2)

## A2. safety check: CNV-malignant cells in normal samples at the new threshold
safety <- lab[sample_type %in% normal_types & !is.na(cnv_call_p95),
              .(epithelial = .N, malignant_p99 = sum(cnv_call == "malignant"),
                malignant_p95 = sum(cnv_call_p95 == "malignant"),
                pct_malignant_p95 = round(100 * mean(cnv_call_p95 == "malignant"), 1),
                pct_uncertain_p95 = round(100 * mean(cnv_call_p95 == "uncertain"), 1)),
              by = .(dataset, sample, patient, sample_type, hpv)][order(dataset, sample)]
fwrite(safety, file.path(tdir, "A_safety_normal_samples.csv"))
if (any(safety[sample_type == "normal_mucosa", pct_malignant_p95] > 5))
  stop("Choi normal mucosa sample(s) above 5 % CNV-malignant at the 95th percentile: ",
       paste(safety[sample_type == "normal_mucosa" & pct_malignant_p95 > 5, sample], collapse = ", "))

state_tab <- merge(lab[, .(cells = .N), by = .(dataset, hpv, state = cell_state_primary)],
                   lab[, .(cells_strict = .N), by = .(dataset, hpv, state = cell_state_strict)],
                   by = c("dataset", "hpv", "state"), all = TRUE)
state_tab[is.na(state_tab)] <- 0L
setnames(state_tab, "cells", "cells_primary_p95")
fwrite(state_tab[order(dataset, hpv, state)], file.path(tdir, "A_cell_states_primary_vs_strict.csv"))
lp <- lab[sample_type == "leukoplakia", .N, by = .(sample, cell_state_primary, cell_state_strict)][order(sample)]
fwrite(lp, file.path(tdir, "A_leukoplakia_split.csv"))

## A7. discovery patient table (cells per patient and group)
disc <- lab[role == "discovery" & !flag_excluded]
pt <- disc[, .(study = study[1], site = site[site != ""][1],
               malignant = sum(cell_state_primary == "malignant"),
               malignant_label_only = sum(cell_state_primary == "malignant_label_only"),
               malignant_strict = sum(cell_state_strict == "malignant"),
               normal = sum(cell_state_primary == "normal_epithelial"),
               normal_strict = sum(cell_state_strict == "normal_epithelial"),
               primary = primary_step04[1], sensitivity_reason = sensitivity_reason[1]),
           by = .(hpv, patient_id_unique)]
pt[, `:=`(tumour_pb = malignant >= min_cells, normal_pb = normal >= min_cells)]
pt[, paired := tumour_pb & normal_pb]
setorder(pt, hpv, -primary, patient_id_unique)
fwrite(pt, file.path(tdir, "A_discovery_patients.csv"))

## =========================================================================
## B. Pseudobulk (raw counts; atlas 45,075 genes, Choi 20,000 genes)
## =========================================================================
groups <- list(malignant = quote(cell_state_primary == "malignant"),
               malignant_label_only = quote(cell_state_primary == "malignant_label_only"),
               malignant_strict = quote(cell_state_strict == "malignant"),
               normal = quote(cell_state_primary == "normal_epithelial"),
               normal_strict = quote(cell_state_strict == "normal_epithelial"),
               LP_CNVpos = quote(cell_state_primary == "LP_CNVpos"),
               LP_CNVneg = quote(cell_state_primary == "LP_CNVneg"))
key_cells <- function(d) rbindlist(lapply(names(groups), function(g)
  d[eval(groups[[g]]), .(cell, key = paste(patient_id_unique, g, sep = "|"))]))

if (!file.exists(paths$pseudobulk_v2)) {
  message("B. Pseudobulk: atlas")
  obj <- readRDS(paths$atlas_labelled)
  pb_atl <- pseudobulk(obj$counts, key_cells(disc[dataset == "HNSCatlas"]))
  rm(obj); free()
  message("B. Pseudobulk: Choi normal mucosa / leukoplakia")
  obj <- readRDS(paths$choi_labelled)
  pb_choi <- pseudobulk(obj$counts, key_cells(disc[dataset == "GSE181919_Choi"]))
  rm(obj); free()
  ## shared genes: atlas columns restricted to Choi's genes, summed by key
  g_sh <- intersect(rownames(pb_choi$counts), rownames(pb_atl$counts))
  keys <- union(colnames(pb_atl$counts), colnames(pb_choi$counts))
  sh <- Matrix(0, length(g_sh), length(keys), dimnames = list(g_sh, keys), sparse = TRUE)
  sh[, colnames(pb_atl$counts)] <- pb_atl$counts[g_sh, ]
  sh[, colnames(pb_choi$counts)] <- sh[, colnames(pb_choi$counts)] + pb_choi$counts[g_sh, ]
  n_sh <- setNames(numeric(length(keys)), keys)
  n_sh[names(pb_atl$n_cells)] <- pb_atl$n_cells
  n_sh[names(pb_choi$n_cells)] <- n_sh[names(pb_choi$n_cells)] + pb_choi$n_cells
  saveRDS(list(full = pb_atl, shared = list(counts = sh, n_cells = n_sh),
               source = list(atlas = colnames(pb_atl$counts), choi = colnames(pb_choi$counts))),
          paths$pseudobulk_v2)
  rm(pb_atl, pb_choi, sh); free()
}
pb <- readRDS(paths$pseudobulk_v2)
ptinfo <- disc[, .(hpv = hpv[1], study = study[1], site = site[site != ""][1], primary_step04 = primary_step04[1]),
               by = patient_id_unique]
pb_meta <- function(m, n_cells) {
  s <- data.table(id = colnames(m), patient = sub("\\|.*$", "", colnames(m)), group = sub("^.*\\|", "", colnames(m)),
                  n_cells = as.integer(n_cells[colnames(m)]), lib = colSums(m))
  setnames(s, "id", "key")
  merge(s, ptinfo, by.x = "patient", by.y = "patient_id_unique", all.x = TRUE, sort = FALSE)
}
sm_full <- pb_meta(pb$full$counts, pb$full$n_cells)
sm_sh <- pb_meta(pb$shared$counts, pb$shared$n_cells)
sm_sh[, source := fifelse(key %in% pb$source$choi & key %in% pb$source$atlas, "atlas+Choi",
                          fifelse(key %in% pb$source$choi, "Choi", "atlas"))]
fwrite(sm_sh[order(hpv, group, patient)], file.path(tdir, "B_pseudobulk_samples.csv"))
dropped <- sm_sh[group %in% c("malignant", "normal") & n_cells < min_cells & primary_step04,
                 .(hpv, group, patient, n_cells)][order(hpv, group, -n_cells)]
fwrite(dropped, file.path(tdir, "B_pseudobulks_dropped_lt20.csv"))
## check: do normal pseudobulks look like squamous epithelium? (summed CPM of
## squamous keratins; Choi's "Epithelial.cells" label in some normal mucosa
## samples covers cells without keratin expression)
sq_cpm <- function(keys) colSums(as.matrix(pb$shared$counts[squamous_krt, keys, drop = FALSE])) /
  colSums(pb$shared$counts[, keys, drop = FALSE]) * 1e6
sq <- sm_sh[group %in% c("normal", "malignant") & n_cells >= min_cells & primary_step04,
            .(key, patient, group, hpv, source, n_cells)]
sq[, squamous_keratin_CPM := round(sq_cpm(key))]
sq[, squamous_ok := squamous_keratin_CPM >= min_squamous]
setorder(sq, group, squamous_keratin_CPM)
fwrite(sq, file.path(tdir, "B_squamous_keratin_check.csv"))
nonsq_normal <- sq[group == "normal" & !squamous_ok, patient]

## =========================================================================
## C. Malignant vs normal epithelium (discovery only)
## =========================================================================
## build a contrast from pseudobulk groups; `mal` / `nor` are group names,
## a vector sums groups per patient (e.g. malignant + label_only)
make_set <- function(hpv_arm, mal = "malignant", nor = "normal", pts_ok = NULL, paired = TRUE, nor_drop = NULL) {
  s <- sm_sh[hpv == hpv_arm & (primary_step04 | patient %in% pts_ok)]
  s <- s[!(group %in% nor & patient %in% nor_drop)]
  agg <- function(gs, state) {
    x <- s[group %in% gs, .(keys = list(key), n_cells = sum(n_cells)), by = .(patient, study, site)]
    x[n_cells >= min_cells][, state := state]
  }
  ss <- rbind(agg(mal, "malignant"), agg(nor, "normal"))
  if (paired) ss <- ss[patient %in% intersect(ss[state == "malignant", patient], ss[state == "normal", patient])]
  ss[, sample := paste(patient, state, sep = "|")]
  cnt <- do.call(cbind, lapply(ss$keys, function(k) Matrix::rowSums(pb$shared$counts[, k, drop = FALSE])))
  colnames(cnt) <- ss$sample
  ss[, `:=`(state = factor(state, levels = c("normal", "malignant")), patient = factor(patient), study = factor(study))]
  list(smp = ss[, !"keys"], counts = cnt)
}
run_c1 <- function(set) de_qlf(set$counts, set$smp, ~ patient + state, "statemalignant", fdr_cut, lfc_cut)
save_de <- function(res, name) {
  fwrite(res$table, file.path(tdir, paste0(name, ".csv")))
  res$table
}
counts_line <- function(name, res, set) data.table(
  contrast = name, patients = uniqueN(set$smp$patient), malignant_pb = sum(set$smp$state == "malignant"),
  normal_pb = sum(set$smp$state == "normal"), genes_tested = res$n_genes_tested,
  up = sum(res$table$direction == "up"), down = sum(res$table$direction == "down"),
  gene_set = sprintf("shared atlas/Choi (%d genes)", nrow(pb$shared$counts)))

message("C1. HPV- paired")
set_c1 <- make_set("HPV-")
c1 <- run_c1(set_c1); C1 <- save_de(c1, "C1_HPVneg_paired_malignant_vs_normal")
message("C2. HPV- all patients")
set_c2 <- make_set("HPV-", paired = FALSE)
study_state <- dcast(set_c2$smp, study ~ state, value.var = "sample", fun.aggregate = length)
fwrite(study_state, file.path(tdir, "C2_study_by_state.csv"))
d2 <- model.matrix(~ study + state, set_c2$smp)
c2_estimable <- qr(d2)$rank == ncol(d2) && any(study_state$normal > 0 & study_state$malignant > 0)
c2 <- if (c2_estimable) de_qlf(set_c2$counts, set_c2$smp, ~ study + state, "statemalignant", fdr_cut, lfc_cut) else
  de_qlf(set_c2$counts, set_c2$smp, ~ state, "statemalignant", fdr_cut, lfc_cut)
C2 <- save_de(c2, "C2_HPVneg_allpatients_malignant_vs_normal")
## C2 rows from patients that contribute both states are not independent;
## report the within-patient correlation the paired data imply
message("C1b. HPV+ paired")
set_c1p <- make_set("HPV+")
c1p <- run_c1(set_c1p); C1p <- save_de(c1p, "C1_HPVpos_paired_malignant_vs_normal")

de_counts <- rbind(counts_line("C1 HPV- paired (~ patient + state)", c1, set_c1),
                   counts_line(sprintf("C2 HPV- all (~ %s)", if (c2_estimable) "study + state" else "state"), c2, set_c2),
                   counts_line("C1 HPV+ paired (~ patient + state)", c1p, set_c1p))
## replication of C1 hits in C2
rep_tab <- merge(C1[sig == TRUE, .(gene, direction, logFC_C1 = logFC)],
                 C2[, .(gene, logFC_C2 = logFC, FDR_C2 = FDR)], by = "gene", all.x = TRUE)
rep_tab[, replicated := !is.na(FDR_C2) & FDR_C2 < fdr_cut & sign(logFC_C2) == sign(logFC_C1)]
rep_sum <- rep_tab[, .(C1_sig = .N, tested_in_C2 = sum(!is.na(FDR_C2)), replicated_C2 = sum(replicated),
                       pct = round(100 * mean(replicated), 1)), by = direction]
fwrite(rep_sum, file.path(tdir, "C_C1_replication_in_C2.csv"))
fwrite(de_counts, file.path(tdir, "C_de_counts.csv"))

## positive controls
pc <- rbindlist(lapply(c(pos_up, pos_down), function(g) {
  r1 <- C1[gene == g]; r2 <- C2[gene == g]
  data.table(gene = g, expected = if (g %in% pos_up) "up" else "down",
             C1_logFC = if (nrow(r1)) round(r1$logFC, 2) else NA, C1_FDR = if (nrow(r1)) signif(r1$FDR, 3) else NA,
             C2_logFC = if (nrow(r2)) round(r2$logFC, 2) else NA, C2_FDR = if (nrow(r2)) signif(r2$FDR, 3) else NA)
}))
pc[, C1_pass := !is.na(C1_FDR) & C1_FDR < fdr_cut & sign(C1_logFC) == fifelse(expected == "up", 1, -1) & abs(C1_logFC) >= lfc_cut]
pc[, C2_pass := !is.na(C2_FDR) & C2_FDR < fdr_cut & sign(C2_logFC) == fifelse(expected == "up", 1, -1) & abs(C2_logFC) >= lfc_cut]
fwrite(pc, file.path(tdir, "C_positive_controls.csv"))
pc_cpm <- cpm(DGEList(set_c1$counts))[c(pos_up, pos_down), , drop = FALSE]
pc_cpm <- data.table(sample = colnames(pc_cpm), round(t(pc_cpm), 1))
pc_cpm[, normal_source := fifelse(grepl("\\|normal$", sample),
                                  sm_sh[match(sub("\\|normal$", "|normal", sample), key), source], "")]
fwrite(pc_cpm[order(grepl("normal", sample), sample)], file.path(tdir, "C_positive_controls_CPM_C1_samples.csv"))

## =========================================================================
## D. Tumour-cell specificity (v1 cell-type pseudobulk, unsorted studies)
## =========================================================================
ctp <- readRDS(paths$celltype_pseudobulk)
cs <- as.data.table(ctp$samples)
cs[, study := study_of(Patient)]
d_pts <- pt[hpv == "HPV-" & primary, patient_id_unique]
cs <- cs[study != "Kurten2021" & Patient %in% d_pts]   # Kurten 2021: CD45-sorted
tme_types <- c("CD4T", "CD8T", "NK", "Bcell", "PlasmaCell", "Macrophage", "Monocyte",
               "panDC", "Neutrophils", "Mast", "Fibroblast", "Endothelial")
up_genes <- intersect(union(C1[direction == "up", gene], C2[direction == "up", gene]), rownames(ctp$counts))
lib_ct <- colSums(ctp$counts)
cpm_ct <- as.matrix(ctp$counts[up_genes, cs$sample, drop = FALSE]) %*% Diagonal(x = 1e6 / lib_ct[cs$sample])
cpm_ct <- as.matrix(cpm_ct); colnames(cpm_ct) <- cs$sample
mean_cpm <- function(g) {
  j <- cs[Group == g & n_cells >= min_cells, sample]
  if (!length(j)) return(rep(NA_real_, length(up_genes)))
  rowMeans(cpm_ct[, j, drop = FALSE])
}
tum_cpm <- mean_cpm("Tumour")
tme <- sapply(tme_types, mean_cpm)
max_tme <- apply(tme, 1, max, na.rm = TRUE)
max_tme_type <- tme_types[apply(tme, 1, which.max)]
## share of whole-tissue counts from tumour cells per patient, median (v1 scripts/08)
share_pts <- cs[Group == "Tumour" & n_cells >= min_cells, Patient]
share <- sapply(share_pts, function(p) {
  j <- cs[Patient == p, sample]
  tot <- Matrix::rowSums(ctp$counts[up_genes, j, drop = FALSE])
  tum <- ctp$counts[up_genes, cs[Patient == p & Group == "Tumour", sample]]
  ifelse(tot > 0, tum / tot, NA_real_)
})
tissue_share <- apply(share, 1, median, na.rm = TRUE)
## detectability: median CPM over HPV- primary malignant pseudobulks (full genes)
mal_full <- sm_full[group == "malignant" & hpv == "HPV-" & primary_step04 & n_cells >= min_cells]
cpm_mal <- as.matrix(pb$full$counts[up_genes, mal_full$key, drop = FALSE]) %*% Diagonal(x = 1e6 / mal_full$lib)
med_mal <- apply(as.matrix(cpm_mal), 1, median)
nor_full <- sm_full[group == "normal" & hpv == "HPV-" & primary_step04 & n_cells >= min_cells]

ann <- data.table(gene = up_genes, C1_logFC = C1[match(up_genes, gene), logFC], C1_FDR = C1[match(up_genes, gene), FDR],
                  C2_logFC = C2[match(up_genes, gene), logFC], C2_FDR = C2[match(up_genes, gene), FDR],
                  tumour_mean_CPM = tum_cpm, max_TME_CPM = max_tme, max_TME_type = max_tme_type,
                  fold_vs_max_TME = tum_cpm / max_tme, tumour_share = tissue_share,
                  malignant_pb_median_CPM = med_mal)
ann <- cbind(ann, setNames(as.data.table(round(tme, 2)), paste0("CPM_", tme_types)))
ann[, `:=`(in_C1_up = gene %in% C1[direction == "up", gene], in_C2_up = gene %in% C2[direction == "up", gene])]
ann[, specificity := fcase(tumour_mean_CPM >= min_fold * max_TME_CPM & tumour_share >= min_share, "tumour_specific",
                           tumour_mean_CPM > max_TME_CPM, "tumour_enriched",
                           default = "shared_with_TME")]
ann[, detectable := malignant_pb_median_CPM >= min_cpm]
## genes with zero counts in every malignant pseudobulk of a study (gene
## absent from that study's annotation, not low expression)
mal_any <- sm_full[group == "malignant" & n_cells >= min_cells]
st_sum <- sapply(split(mal_any$key, mal_any$study), function(k) Matrix::rowSums(pb$full$counts[ann$gene, k, drop = FALSE]))
ann[, zero_in_studies := apply(st_sum == 0, 1, function(z) paste(colnames(st_sum)[z], collapse = " "))]
ann[, C2_replicated := !is.na(C2_FDR) & C2_FDR < fdr_cut & C2_logFC > 0]
## C1 (primary) genes first by C1 FDR, then C2-only genes by C2 FDR
ann[, rank_key := fifelse(in_C1_up, C1_FDR, 1 + C2_FDR)]
setorder(ann, rank_key)
ann[, rank_key := NULL]
fwrite(ann, file.path(tdir, "D_gene_annotation_up_genes.csv"))
spec_sum <- rbind(ann[in_C1_up == TRUE, .(set = "C1 up (primary)", .N), by = .(specificity, detectable)],
                  ann[in_C2_up == TRUE, .(set = "C2 up", .N), by = .(specificity, detectable)],
                  ann[in_C2_up == TRUE & in_C1_up == FALSE, .(set = "C2 up only (extended)", .N), by = .(specificity, detectable)])
setorder(spec_sum, set, specificity, -detectable)
fwrite(spec_sum, file.path(tdir, "D_specificity_summary.csv"))
d_n <- c(patients = length(unique(cs$Patient)), tumour_share_patients = length(share_pts),
         malignant_pb = nrow(mal_full))
fwrite(as.data.table(as.list(d_n)), file.path(tdir, "D_patients_used.csv"))

## =========================================================================
## E. Sensitivity re-runs of C1 (HPV-)
## =========================================================================
sig_up <- function(t) t[direction == "up", gene]
overlap <- function(name, set, res, ref = sig_up(C1), ref_set = set_c1) {
  a <- sig_up(res$table)
  data.table(analysis = name, patients = uniqueN(set$smp$patient),
             patients_added = paste(setdiff(levels(set$smp$patient), levels(ref_set$smp$patient)), collapse = " "),
             patients_lost = paste(setdiff(levels(ref_set$smp$patient), levels(set$smp$patient)), collapse = " "),
             normal_pb_lost = paste(setdiff(ref_set$smp[state == "normal", as.character(patient)],
                                            set$smp[state == "normal", as.character(patient)]), collapse = " "),
             sig_up = length(a), shared_with_primary = length(intersect(a, ref)),
             jaccard = round(length(intersect(a, ref)) / length(union(a, ref)), 3),
             pct_primary_retained = round(100 * length(intersect(a, ref)) / length(ref), 1))
}
message("E. Sensitivity")
add_back <- c(prior_primary, "OP8")
set_a <- make_set("HPV-", mal = "malignant_strict", nor = "normal_strict")
set_b <- make_set("HPV-", mal = c("malignant", "malignant_label_only"))
set_c <- make_set("HPV-", pts_ok = add_back)
e_a <- run_c1(set_a); save_de(e_a, "E_a_C1_strict_labels")
e_b <- run_c1(set_b); save_de(e_b, "E_b_C1_with_label_only")
e_c <- run_c1(set_c); save_de(e_c, "E_c_C1_with_prior_primaries_OP8")
## the added patients enter C1 only if they have a normal pseudobulk; the
## same additions are also tested in the all-patient design
set_c2c <- make_set("HPV-", pts_ok = add_back, paired = FALSE)
e_c2 <- de_qlf(set_c2c$counts, set_c2c$smp, if (c2_estimable) ~ study + state else ~ state, "statemalignant", fdr_cut, lfc_cut)
save_de(e_c2, "E_c_C2_with_prior_primaries_OP8")
## (d) additional: pairs whose normal pseudobulk is Choi normal mucosa only
nl_pts <- sm_sh[group == "normal" & source == "Choi", patient]
set_d <- list(smp = droplevels(set_c1$smp[patient %in% nl_pts]))
set_d$counts <- set_c1$counts[, set_d$smp$sample]
e_d <- run_c1(set_d); save_de(e_d, "E_d_C1_Choi_normal_mucosa_pairs")
## (e) additional: normal pseudobulks without squamous keratin expression removed
set_e <- make_set("HPV-", nor_drop = nonsq_normal)
e_e <- run_c1(set_e); save_de(e_e, "E_e_C1_squamous_normals_only")
set_e2 <- make_set("HPV-", paired = FALSE, nor_drop = nonsq_normal)
e_e2 <- de_qlf(set_e2$counts, set_e2$smp, if (c2_estimable) ~ study + state else ~ state, "statemalignant", fdr_cut, lfc_cut)
save_de(e_e2, "E_e_C2_squamous_normals_only")
set_ep <- make_set("HPV+", nor_drop = nonsq_normal)
e_ep <- run_c1(set_ep); save_de(e_ep, "E_e_C1_HPVpos_squamous_normals_only")
sens <- rbind(overlap("primary C1", set_c1, c1),
              overlap("(a) strict CNV labels (99th pct)", set_a, e_a),
              overlap("(b) label_only cells added", set_b, e_b),
              overlap("(c) BHN39/46/49 + OP8 added (C1)", set_c, e_c),
              overlap("(d) additional: Choi normal-mucosa pairs only", set_d, e_d),
              overlap("(e) additional: non-squamous normal pseudobulks removed (C1)", set_e, e_e),
              overlap("(e') additional: non-squamous normals removed (C2 design; vs primary C2)", set_e2, e_e2, ref = sig_up(C2), ref_set = set_c2),
              overlap("(e'') additional: HPV+ C1, non-squamous normals removed (vs HPV+ C1)", set_ep, e_ep, ref = sig_up(C1p), ref_set = set_c1p),
              overlap("(c') BHN39/46/49 + OP8 added (C2 design; vs primary C2)", set_c2c, e_c2, ref = sig_up(C2), ref_set = set_c2))
fwrite(sens, file.path(tdir, "E_sensitivity_overlap.csv"))
## flag genes that stay up when non-squamous normal pseudobulks are removed
ann[, `:=`(up_C1_squamous_normals = gene %in% sig_up(e_e$table), up_C2_squamous_normals = gene %in% sig_up(e_e2$table))]
fwrite(ann, file.path(tdir, "D_gene_annotation_up_genes.csv"))

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
dir_cols <- c(up = "#e34948", down = "#2a78d6", ns = "#c3c2b7")
spec_cols <- c(tumour_specific = "#e34948", tumour_enriched = "#eda100", shared_with_TME = "#898781")

## Fig 14: volcano C1
v <- copy(C1)[, `:=`(nlp = -log10(PValue), face = fifelse(gene %in% c(pos_up, pos_down), "bold", "plain"))]
lab_g <- unique(c(v[direction == "up"][order(PValue)][1:15, gene], v[direction == "down"][order(PValue)][1:10, gene],
                  intersect(c(pos_up, pos_down), v$gene)))
f14 <- ggplot(v, aes(logFC, nlp, colour = direction)) +
  ggrastr::geom_point_rast(size = 0.4, stroke = 0, shape = 16, raster.dpi = 300, dev = "ragg") +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "22", linewidth = 0.3, colour = ink2) +
  geom_text_repel(data = v[gene %in% lab_g], aes(label = gene, fontface = face), size = 2, colour = ink,
                  max.overlaps = 40, segment.size = 0.2, min.segment.length = 0) +
  scale_colour_manual(values = dir_cols, name = NULL,
                      labels = c(up = sprintf("up (%d)", sum(v$direction == "up")),
                                 down = sprintf("down (%d)", sum(v$direction == "down")), ns = "n.s.")) +
  labs(title = sprintf("C1: HPV- malignant vs normal epithelium, paired (%d patients)", uniqueN(set_c1$smp$patient)),
       subtitle = sprintf("edgeR QL, ~ patient + state; FDR < %.2f and |log2FC| >= %g; bold = positive controls", fdr_cut, lfc_cut),
       x = "log2 fold change (malignant / normal)", y = "-log10 P") + theme_fig()
save_fig(f14, "Fig14_volcano_C1", 150, 110)

## Fig 15: specificity scatter
sp <- copy(ann)[, `:=`(x = log10(max_TME_CPM + 1), y = log10(tumour_mean_CPM + 1),
                       set = factor(fifelse(in_C1_up, sprintf("C1 up, primary (%d)", sum(in_C1_up)),
                                            sprintf("C2 up only, extended (%d)", sum(!in_C1_up)))))]
top_spec <- ann[specificity == "tumour_specific" & detectable == TRUE][1:min(.N, 30), gene]
f15 <- ggplot(sp, aes(x, y, colour = specificity)) +
  geom_abline(intercept = 0, slope = 1, linewidth = 0.3, colour = ink2) +
  geom_point(aes(shape = detectable), size = 0.9, alpha = 0.8) +
  geom_text_repel(data = sp[gene %in% top_spec[1:min(20, length(top_spec))]], aes(label = gene), size = 2, colour = ink,
                  segment.size = 0.2, max.overlaps = 30) +
  facet_wrap(~ set) +
  scale_colour_manual(values = spec_cols, name = NULL) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), name = sprintf("median CPM >= %d", min_cpm)) +
  labs(title = "Tumour-cell specificity of up-regulated genes (C1 primary; C2-only genes as extension)",
       subtitle = sprintf("mean CPM over HPV- patients (unsorted studies, >= %d cells per type); line: equal expression", min_cells),
       x = "log10(CPM + 1), highest microenvironment cell type", y = "log10(CPM + 1), tumour cells") + theme_fig()
save_fig(f15, "Fig15_specificity_scatter", 180, 95)

## Fig 16: heatmap of the top 30 tumour-specific, detectable genes
if (length(top_spec)) {
  hs <- sm_sh[group %in% c("malignant", "normal") & hpv == "HPV-" & primary_step04 & n_cells >= min_cells]
  lc <- cpm(DGEList(as.matrix(pb$shared$counts[, hs$key])), log = TRUE, prior.count = 2)[top_spec, , drop = FALSE]
  z <- t(scale(t(lc)))
  ord <- hs[, .(key, group, study, source, patient)][order(group, source, study, patient)]
  hd <- data.table(gene = rep(rownames(z), ncol(z)), id = rep(colnames(z), each = nrow(z)), z = pmax(pmin(as.vector(z), 2.5), -2.5))
  setnames(hd, "id", "key")
  hd[, panel := factor(fifelse(grepl("\\|malignant$", key), "malignant", "normal epithelium"),
                       levels = c("malignant", "normal epithelium"))]
  hd[, `:=`(gene = factor(gene, levels = rev(top_spec)), key = factor(key, levels = ord$key))]
  f16 <- ggplot(hd, aes(key, gene, fill = z)) + geom_tile() +
    scale_fill_gradient2(low = "#2a78d6", mid = "#f0efec", high = "#e34948", name = "z (log CPM)") +
    scale_x_discrete(labels = function(k) paste0(sub("\\|.*$", "", k),
                                                 fifelse(sub("\\|.*$", "", k) %in% nonsq_normal & grepl("normal$", k), "*", ""))) +
    facet_grid(~ panel, scales = "free_x", space = "free_x") +
    labs(title = sprintf("Top %d tumour-specific, detectable genes (C1 genes first by C1 FDR, then C2-only by C2 FDR)", length(top_spec)),
         subtitle = sprintf("HPV- primary pseudobulks with >= %d cells, %d genes shared with Choi; * = normal pseudobulk without squamous keratins",
                            min_cells, nrow(pb$shared$counts)), x = NULL, y = NULL) +
    theme_fig() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = base - 2),
                        axis.line = element_blank(), axis.ticks = element_blank())
  save_fig(f16, "Fig16_heatmap_top30_tumour_specific", 180, 120)
}

## console summary
print(usum[, .(units = .N, malignant_p99 = sum(malignant_p99), malignant_p95 = sum(malignant_p95),
               informative_p99 = sum(informative_p99), informative_p95 = sum(informative_p95)), by = dataset])
print(safety); print(merged_ids); print(de_counts); print(rep_sum); print(pc); print(spec_sum); print(sens)
message("Done.")
