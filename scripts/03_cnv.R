## Step 03: doublets, CNV-based malignant cell calling, final normal reference
##
## Purpose: flag doublets (scDblFinder) in Peng and Bill; infer CNV per sample
##   for every epithelial cell (atlas, Choi NL/LP, Bill recurrences/metastases,
##   Peng) against same-sample immune/stromal reference cells; check calls
##   with copyKAT on a subset; add Bill Table S1 information; assign final
##   cell states and the normal reference.
## Inputs:  data/processed/v2/{atlas_epithelial_counts, choi_nl_lp_allcells,
##          bill_recurrence_epithelial, peng_allcells_qc, peng_epithelial}.rds
##          data/processed/v2/epithelial_v2_meta.csv.gz (step 02)
##          data/processed/HNSCatlas/cnv_input.rds (atlas reference cells)
##          data/processed/HNSCatlas/nonmalignant_subtypes.rds (squamous_clean)
##          data/raw/GSE234933/counts/*.rds, GSE234933 cell/sample annotation
##          data/reference/ensembl (gene positions)
##          Bill et al. 2023 Table S1 (docs/references/), transcribed below
## Outputs: results/step03_cnv/ (tables, figures)
##          data/processed/v2/cnv/ (caches), data/processed/v2/*_labelled.rds,
##          data/processed/v2/cell_labels_v2.csv.gz
## Date: 2026-10-12
## Run from the repo root:  Rscript scripts/03_cnv.R

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
  library(ggrastr)
  library(patchwork)
})
source("config/paths.R")
source("scripts/lib_cnv.R")
out <- paths$step03
cdir <- paths$cnv_cache
dir.create(file.path(out, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(cdir, recursive = TRUE, showWarnings = FALSE)
free <- function() invisible(gc(full = TRUE))
set.seed(1)
n_genes_cnv <- 7000
max_ref     <- 300
nasopharynx <- c("BHN77", "BHN34", "BHN61", "BHN73")
heat_units  <- c("BHN70", "Bill_HN23", "P15")

meta2 <- fread(paths$epithelial_v2_meta)
setkey(meta2, cell)

## =========================================================================
## A. Doublets (scDblFinder, per sample; flag only)
## =========================================================================
dbl_file <- file.path(cdir, "doublets.csv.gz")
if (!file.exists(dbl_file)) {
  suppressPackageStartupMessages({ library(scDblFinder); library(SingleCellExperiment) })
  run_dbl <- function(m) {
    set.seed(1)
    sce <- scDblFinder(SingleCellExperiment(list(counts = m)), verbose = FALSE)
    data.table(cell = colnames(m), dbl_class = sce$scDblFinder.class, dbl_score = sce$scDblFinder.score)
  }
  message("A. Doublets: Peng")
  pa <- readRDS(paths$peng_all_v2)
  smp <- sub("^Peng_([^_]+)_.*$", "\\1", colnames(pa$counts))
  d_peng <- rbindlist(lapply(split(seq_along(smp), smp), function(i) run_dbl(pa$counts[, i])))
  d_peng[, dataset := "GSE172577_Peng"]
  rm(pa); free()
  message("A. Doublets: Bill (annotated cells; doublet removal by the authors is UNCLEAR)")
  bc <- fread(paths$bill_cells)
  d_bill <- rbindlist(lapply(list.files(file.path(paths$bill_dir, "counts"), full.names = TRUE), function(f) {
    m <- readRDS(f)
    s <- sub("_gex_raw_counts.rds$", "", basename(f))
    m <- m[, intersect(bc[sample == s, sample_barcode], colnames(m))]
    r <- run_dbl(m); r[, cell := paste0("Bill_", cell)]; r
  }))
  d_bill[, dataset := "GSE234933_Bill"]
  fwrite(rbind(d_peng, d_bill), dbl_file)
  free()
}
dbl <- fread(dbl_file)
dbl[, sample := fifelse(dataset == "GSE172577_Peng", sub("^Peng_([^_]+)_.*$", "\\1", cell),
                        sub("^Bill_([^_]+)_.*$", "\\1", cell))]
dbl_tab <- dbl[, .(cells = .N, doublets = sum(dbl_class == "doublet"),
                   pct = round(100 * mean(dbl_class == "doublet"), 1)), by = .(dataset, sample)]
fwrite(dbl_tab, file.path(out, "tables", "doublets_per_sample.csv"))

## =========================================================================
## B. CNV gene list: positioned autosomal genes present in all datasets,
##    top 7,000 by mean expression in atlas epithelial + reference cells
## =========================================================================
pos_file <- file.path(cdir, "cnv_genes.csv")
if (!file.exists(pos_file)) {
  pos <- gene_positions(paths$ensembl_dir)
  g_all <- Reduce(intersect, list(rownames(readRDS(paths$choi_epi_v2)$counts),
                                  rownames(readRDS(paths$bill_epi_v2)$counts),
                                  rownames(readRDS(paths$peng_epi_v2)$counts)))
  atl <- readRDS(paths$atlas_epi_v2)$counts
  g <- intersect(intersect(g_all, rownames(atl)), pos$gene)
  idx <- sample(ncol(atl), min(20000, ncol(atl)))
  x <- atl[g, idx] %*% Diagonal(x = 1e5 / colSums(atl[, idx]))
  x@x <- log2(x@x + 1)
  top <- names(sort(Matrix::rowMeans(x), decreasing = TRUE))[seq_len(n_genes_cnv)]
  pos <- pos[gene %in% top]
  fwrite(pos, pos_file)
  rm(atl, x); free()
}
pos <- fread(pos_file)
setorder(pos, chr_n, start)

## =========================================================================
## C. CNV per sample
## =========================================================================
keep_heat <- function(res, unit, groups) {
  ## subsample the smoothed matrix for heatmap examples (<= 80 cells per group)
  if (!unit %in% heat_units || is.null(res$s)) return(NULL)
  keep <- unlist(lapply(split(seq_along(groups), groups), function(i) if (length(i) > 80) sample(i, 80) else i))
  list(s = res$s[, keep], group = groups[keep], chr = res$chr)
}
units <- list(); heat <- list()
summarise_unit <- function(unit, dataset, res, ep_meta, prof_cells) {
  r <- res$cells[is_epi == TRUE]
  r[, `:=`(unit = unit, dataset = dataset)]
  units[[unit]] <<- list(cells = r, profile = res$profile, mal_profile = res$mal_profile,
                         genes = res$genes, thr = res$thr, n_ref = res$n_ref, n_profile = res$n_profile)
}

cnv_file <- file.path(cdir, "cnv_units.rds")
if (!file.exists(cnv_file)) {
  ## ---- C1. atlas patients (Choi patients in C2) ----------------------------
  message("C1. CNV: atlas patients")
  atl <- readRDS(paths$atlas_epi_v2)
  ref <- readRDS(paths$cnv_input)
  ref_meta <- as.data.table(ref$meta)[group == "Reference"]
  am <- meta2[dataset == "HNSCatlas"]
  for (p in setdiff(unique(am$patient), nasopharynx)) {
    if (grepl("^P[0-9]", p)) next
    ep <- am[patient == p]
    rc <- ref_meta[Patient == p, cell]
    if (length(rc) < 30 || nrow(ep) < 1) next
    cnt <- cbind(atl$counts[, ep$cell, drop = FALSE], ref$counts[, rc, drop = FALSE])
    is_ref <- colnames(cnt) %in% rc
    in_prof <- colnames(cnt) %in% ep[sample_type != "adjacent_normal" & !suspect_label, cell]
    res <- run_cnv(cnt, colSums(cnt), is_ref, !is_ref, in_prof, pos, keep_s = p %in% heat_units)
    grp <- fifelse(is_ref, "Reference", paste0("Epithelial: ", res$cells$cnv_call))
    heat[[p]] <- keep_heat(res, p, grp)
    summarise_unit(p, "HNSCatlas", res)
    rm(cnt, res); free()
  }
  ## ---- C2. Choi patients: CA (atlas) + NL/LP (GEO), Choi's 20k genes -------
  message("C2. CNV: Choi patients (CA + NL + LP)")
  ch <- readRDS(paths$choi_nl_lp_all)
  chm <- ch$meta
  g_choi <- intersect(rownames(ch$counts), rownames(atl$counts))
  choi_pts <- union(unique(am[grepl("^P[0-9]", patient), patient]), unique(chm$patient))
  for (p in choi_pts) {
    ca <- am[patient == p, cell]
    ca_ref <- ref_meta[Patient == p, cell]
    nl_epi <- chm[patient == p & celltype == "Epithelial.cells", cell]
    nl_ref <- chm[patient == p & celltype != "Epithelial.cells",
                  .(cell = if (.N > max_ref) sample(cell, max_ref) else cell), by = sample]$cell
    cnt <- cbind(atl$counts[g_choi, ca, drop = FALSE], ref$counts[g_choi, ca_ref, drop = FALSE],
                 ch$counts[g_choi, c(nl_epi, nl_ref), drop = FALSE])
    is_ref <- colnames(cnt) %in% c(ca_ref, nl_ref)
    is_epi <- !is_ref
    ## profile from the tumour (CA) sample; LP-only patients use their LP cells
    in_prof <- if (length(ca)) colnames(cnt) %in% am[patient == p & !suspect_label, cell] else is_epi
    if (sum(is_ref) < 30) next
    res <- run_cnv(cnt, colSums(cnt), is_ref, is_epi, in_prof, pos, keep_s = p %in% heat_units)
    src <- fifelse(colnames(cnt) %in% ca, "CA", sub("^Choi_([A-Z]+)[0-9]+_.*$", "\\1", colnames(cnt)))
    grp <- fifelse(is_ref, "Reference", paste0(src, ": ", res$cells$cnv_call))
    heat[[p]] <- keep_heat(res, p, grp)
    summarise_unit(p, "Choi (CA atlas + NL/LP)", res)
    rm(cnt, res); free()
  }
  rm(atl, ref, ch); free()

  ## ---- C3. Bill recurrences / metastases ------------------------------------
  message("C3. CNV: Bill recurrences/metastases")
  bill <- readRDS(paths$bill_epi_v2)
  bc <- fread(paths$bill_cells)
  bm <- meta2[dataset == "GSE234933_Bill"]
  for (s in unique(sub("^Bill_([^_]+)_.*$", "\\1", bm$cell))) {
    u <- paste0("Bill_", s)
    if (paste0("B", s) %in% nasopharynx) next
    raw <- readRDS(file.path(paths$bill_dir, "counts", paste0(s, "_gex_raw_counts.rds")))
    nonepi <- bc[sample == s & subclustering != "Tumor cells", sample_barcode]
    nonepi <- setdiff(nonepi, sub("^Bill_", "", dbl[dbl_class == "doublet", cell]))
    if (length(nonepi) > max_ref) nonepi <- sample(nonepi, max_ref)
    ep <- bm[startsWith(cell, paste0(u, "_"))]
    cnt <- cbind(bill$counts[, ep$cell, drop = FALSE], raw[rownames(bill$counts), nonepi, drop = FALSE])
    lib <- c(colSums(raw[, sub("^Bill_", "", ep$cell)]), colSums(raw[, nonepi, drop = FALSE]))
    colnames(cnt)[-(seq_len(nrow(ep)))] <- paste0("Bill_", nonepi)
    is_ref <- !colnames(cnt) %in% ep$cell
    in_prof <- colnames(cnt) %in% ep[compartment != "nonmalignant_label", cell]
    if (sum(is_ref) < 30) next
    res <- run_cnv(cnt, lib, is_ref, !is_ref, in_prof, pos, keep_s = u %in% heat_units)
    grp <- fifelse(is_ref, "Reference", paste0("Epithelial: ", res$cells$cnv_call))
    heat[[u]] <- keep_heat(res, u, grp)
    summarise_unit(u, "GSE234933_Bill", res)
    rm(raw, cnt, res); free()
  }
  rm(bill); free()

  ## ---- C4. Peng ---------------------------------------------------------------
  message("C4. CNV: Peng")
  pa <- readRDS(paths$peng_all_v2)
  cl <- readRDS(paths$peng_epi_v2)$clusters
  pa$meta[, lineage := cl$lineage[match(cluster, cl$cluster)]]
  pm <- meta2[dataset == "GSE172577_Peng"]
  for (s in unique(pa$meta$sample)) {
    u <- paste0("Peng_", s)
    ep <- pm[patient == u, cell]
    rc <- pa$meta[sample == s & lineage %in% c("immune", "fibroblast", "endothelial") &
                  !cell %in% dbl[dbl_class == "doublet", cell], cell]
    if (length(rc) > max_ref) rc <- sample(rc, max_ref)
    cnt <- pa$counts[, c(ep, rc)]
    is_ref <- colnames(cnt) %in% rc
    in_prof <- !is_ref & !colnames(cnt) %in% dbl[dbl_class == "doublet", cell]
    res <- run_cnv(cnt, colSums(cnt), is_ref, !is_ref, in_prof, pos)
    summarise_unit(u, "GSE172577_Peng", res)
    rm(cnt, res); free()
  }
  rm(pa); free()
  saveRDS(list(units = units, heat = heat), cnv_file)
}
cnv <- readRDS(cnv_file)
units <- cnv$units
cells <- rbindlist(lapply(units, `[[`, "cells"))
cells <- merge(cells, meta2[, .(cell, patient, sample, sample_type, compartment, hpv, site, suspect_label,
                                flag_excluded, flag_treated)], by = "cell", all.x = TRUE)

## informative: >= 20 malignant cells and >= 50% of non-suspect tumour-labelled
## cells malignant; samples without tumour labels (Peng, LP-only): >= 20% of
## epithelial cells malignant
usum <- cells[, {
  lab <- compartment == "tumour_label" & !suspect_label & sample_type %in% c("primary_tumour", "recurrence", "metastasis", "lymph_node")
  n_mal <- sum(cnv_call == "malignant")
  pct_lab <- if (sum(lab) >= 20) 100 * mean(cnv_call[lab] == "malignant") else NA_real_
  .(dataset = dataset[1], hpv = hpv[1], site = site[1], epithelial_scored = .N, malignant = n_mal,
    pct_malignant = round(100 * n_mal / .N, 1), pct_uncertain = round(100 * mean(cnv_call == "uncertain"), 1),
    tumour_labelled = sum(lab), pct_tumour_labelled_malignant = round(pct_lab, 1))
}, by = unit]
usum[, n_ref := vapply(unit, function(u) units[[u]]$n_ref, 0L)]
usum[, informative := malignant >= 20 & fifelse(is.na(pct_tumour_labelled_malignant),
                                                pct_malignant >= 20, pct_tumour_labelled_malignant >= 50)]
setorder(usum, dataset, hpv, -epithelial_scored)
fwrite(usum, file.path(out, "tables", "cnv_per_sample.csv"))
cells[, informative := usum$informative[match(unit, usum$unit)]]

## =========================================================================
## D. Bill Table S1 (Bill et al. 2023 Science, Table S1, PMC10755760)
## =========================================================================
## Table S1 columns: patient ID, sex, age, HPV, smoking, original anatomic site,
## anatomic location of scRNA-seq specimen. It has NO treatment or timing
## columns. Entries that differ from or add to the GEO sample annotation:
tbl_s1 <- fread(text = "
sample,specimen_table_s1,ebv_note
HN38,Locoregional recurrence (persistent),
HN39,Primary (second),
HN45,Locoregional recurrence (second after HN23),
HN46,Primary (third),
HN49,Primary (second),
HN50,Primary (second),
HN34,Distant metastasis (Lung),HPV-positive EBV-negative
HN61,Distant metastasis (Liver),HPV-negative EBV-positive
HN73,Locoregional recurrence,HPV-negative EBV-negative
HN77,Primary,HPV-positive EBV-negative")
bs <- fread(paths$bill_samples)
setnames(bs, c("sample", "sex", "age", "hpv", "smoking", "site", "specimen_geo"))
bs <- merge(bs, tbl_s1, by = "sample", all.x = TRUE)
bs[is.na(specimen_table_s1), specimen_table_s1 := specimen_geo]
bs[, patient_id := fifelse(sample == "HN45", "HN23", sample)]
bs[, atlas_id := paste0("B", sample)]
bs[, in_atlas := atlas_id %in% meta2[dataset == "HNSCatlas", patient]]
bs[, prior_hn_primary := grepl("second|third", specimen_table_s1) & grepl("^Primary", specimen_table_s1)]
bs[, `:=`(prior_surgery = "UNCLEAR", prior_radiation = "UNCLEAR", prior_chemotherapy = "UNCLEAR",
          prior_immunotherapy = "UNCLEAR",
          treatment_note = fcase(grepl("recurrence|metastasis", specimen_table_s1, ignore.case = TRUE),
                                 "recurrent/metastatic: previously treated by implication; type not reported in Table S1",
                                 prior_hn_primary, "second/third primary: earlier HNC, treatment not reported",
                                 default = "not reported in Table S1"),
          time_to_recurrence = "not reported in Table S1",
          excluded = atlas_id %in% nasopharynx)]
setorder(bs, -in_atlas, sample)
fwrite(bs[, .(sample, atlas_id, patient_id, in_atlas, sex, age, hpv, site, specimen_geo, specimen_table_s1,
              prior_hn_primary, prior_surgery, prior_radiation, prior_chemotherapy, prior_immunotherapy,
              treatment_note, time_to_recurrence, ebv_note, excluded)],
       file.path(out, "bill_treatment_pairing.csv"))

## =========================================================================
## E. Specific checks
## =========================================================================
prof <- function(u) units[[u]]$mal_profile
cor_prof <- function(a, b) {
  pa <- prof(a); pb <- prof(b)
  if (is.null(pa) || is.null(pb)) return(NA_real_)
  g <- intersect(names(pa), names(pb))
  if (!length(g)) g <- seq_along(pa)
  round(cor(pa[g], pb[g]), 3)
}
## names on profiles: gene order of each unit
for (u in names(units)) {
  if (!is.null(units[[u]]$mal_profile)) names(units[[u]]$mal_profile) <- units[[u]]$genes
}
bg <- c("BHN70", "BHN64", "OP12", "HN01", "BHN60")
pairs <- rbind(
  data.table(a = "Bill_HN13", b = "BHN46", type = "test: HN13 vs BHN46"),
  data.table(a = "Bill_HN45", b = "Bill_HN23", type = "positive control: HN23 vs HN45 (same patient, Table S1)"),
  data.table(a = "Bill_HN13", b = bg, type = "background: HN13 vs unrelated HPV- primary"),
  data.table(a = "BHN46", b = bg, type = "background: BHN46 vs unrelated HPV- primary"))
pairs[, r := mapply(cor_prof, a, b)]
fwrite(pairs, file.path(out, "tables", "check_hn13_bhn46_profiles.csv"))

nl_lp <- cells[sample_type %in% c("normal_mucosa", "leukoplakia"),
               .(epithelial = .N, cnv_malignant = sum(cnv_call == "malignant"),
                 pct_malignant = round(100 * mean(cnv_call == "malignant"), 1),
                 pct_uncertain = round(100 * mean(cnv_call == "uncertain"), 1),
                 unit_informative = informative[1]),
               by = .(sample, patient, sample_type, hpv)]
fwrite(nl_lp, file.path(out, "tables", "check_choi_normal_leukoplakia.csv"))

kurten <- cells[grepl("^HN[0-9]", patient) & compartment == "tumour_label",
                .(tumour_labelled = .N, suspect = sum(suspect_label),
                  suspect_malignant = sum(suspect_label & cnv_call == "malignant"),
                  suspect_uncertain = sum(suspect_label & cnv_call == "uncertain"),
                  nonsuspect_malignant = sum(!suspect_label & cnv_call == "malignant"),
                  informative = informative[1]), by = .(patient, hpv)]
setorder(kurten, -suspect)
fwrite(kurten, file.path(out, "tables", "check_kurten_suspect.csv"))

op34 <- cells[sample == "Puram_OP34norm", .N, by = .(compartment, cnv_call)]
fwrite(op34, file.path(out, "tables", "check_op34norm.csv"))

## =========================================================================
## F. Final labels
## =========================================================================
lab <- copy(meta2)
lab <- merge(lab, cells[, .(cell, unit, signal, signal_rel, corr, cnv_call, informative)], by = "cell", all.x = TRUE)
lab <- merge(lab, dbl[, .(cell, dbl_class)], by = "cell", all.x = TRUE)
sq <- as.data.table(readRDS(paths$nonmalignant_subtypes))
lab[, squamous_clean := cell %in% sq[squamous_clean == TRUE, cell]]
lab[, label_only := FALSE]
lab[, cell_state_v2 := fcase(
  flag_excluded, "excluded",
  dbl_class %in% "doublet", "doublet",
  sample_type == "leukoplakia", "leukoplakia_epithelial",
  cnv_call %in% "malignant" & sample_type %in% c("normal_mucosa", "adjacent_normal"), "uncertain",
  cnv_call %in% "malignant", "malignant",
  compartment == "tumour_label" & !suspect_label & informative %in% FALSE &
    sample_type %in% c("primary_tumour", "recurrence", "metastasis", "lymph_node"), "malignant",
  cnv_call %in% "non_malignant" & !suspect_label &
    (sample_type %in% c("normal_mucosa", "adjacent_normal") | squamous_clean), "normal_epithelial",
  default = "uncertain")]
lab[cell_state_v2 == "malignant" & !cnv_call %in% "malignant", label_only := TRUE]
lab[, state_reason := fcase(
  cell_state_v2 != "uncertain", "",
  is.na(cnv_call), "not scored (sample without enough reference cells)",
  cnv_call == "malignant", "CNV-malignant cell in a normal sample",
  cnv_call == "uncertain", "CNV signal and correlation disagree",
  suspect_label, "suspect label, not CNV-malignant",
  compartment == "tumour_label", "tumour label, CNV non-malignant in an informative sample",
  compartment == "nonmalignant_label", "CNV non-malignant, not clean squamous (gland/ciliated/other or non-HN epithelium)",
  compartment == "unannotated", "no lineage label, CNV non-malignant",
  default = "other")]
## Bill Table S1 columns for Bill patients (atlas BHN* and extra samples)
bcols <- bs[, .(patient = atlas_id, bill_patient_id = patient_id, bill_specimen_table_s1 = specimen_table_s1,
                bill_prior_hn_primary = prior_hn_primary, bill_treatment_note = treatment_note)]
lab <- merge(lab, bcols, by = "patient", all.x = TRUE)
lab[, primary_analysis := !flag_excluded & !flag_treated]
fwrite(lab, paths$cell_labels_v2)

state_tab <- dcast(lab, dataset + hpv ~ cell_state_v2, value.var = "cell", fun.aggregate = length)
fwrite(state_tab, file.path(out, "tables", "cell_state_by_dataset.csv"))
reason_tab <- lab[cell_state_v2 == "uncertain", .N, by = .(dataset, state_reason)][order(dataset, -N)]
fwrite(reason_tab, file.path(out, "tables", "uncertain_reasons.csv"))

## patients with >= 50 malignant cells (primary analysis set)
mal_pt <- lab[primary_analysis & cell_state_v2 == "malignant",
              .(malignant = .N, malignant_cnv_only = sum(!label_only)), by = .(role, dataset, hpv, patient, sample_type)]
fwrite(mal_pt, file.path(out, "tables", "malignant_cells_per_patient.csv"))
## one row per patient (OP35 has tumour + lymph-node cells)
mal_pp <- mal_pt[, .(malignant = sum(malignant), malignant_cnv_only = sum(malignant_cnv_only)), by = .(role, hpv, patient)]
mal_cnt <- mal_pp[, .(with_label_only = sum(malignant >= 50), cnv_only = sum(malignant_cnv_only >= 50)),
                  by = .(role, hpv)][order(role, hpv)]
fwrite(mal_cnt, file.path(out, "tables", "patients_50plus_malignant.csv"))

## final normal reference
normref <- lab[cell_state_v2 == "normal_epithelial",
               .(cells = .N), by = .(hpv, source = fcase(sample_type == "normal_mucosa", "Choi normal mucosa",
                                                          sample_type == "adjacent_normal", "Puram adjacent normal",
                                                          default = "atlas clean squamous"),
                                     patient, informative_unit = informative, primary_analysis)]
setorder(normref, hpv, source, -cells)
fwrite(normref, file.path(out, "tables", "normal_reference_cells.csv"))
fwrite(normref[, .(patients = uniqueN(patient), cells = sum(cells)), by = .(hpv, source)],
       file.path(out, "tables", "normal_reference_summary.csv"))

## per-dataset objects with labels
lab_cols <- c("cell", "unit", "signal", "signal_rel", "corr", "cnv_call", "informative", "dbl_class",
              "cell_state_v2", "label_only", "state_reason")
for (d in list(c("HNSCatlas", "atlas_epi_v2", "atlas_labelled"), c("GSE181919_Choi", "choi_epi_v2", "choi_labelled"),
               c("GSE234933_Bill", "bill_epi_v2", "bill_labelled"), c("GSE172577_Peng", "peng_epi_v2", "peng_labelled"))) {
  f <- paths[[d[3]]]
  if (file.exists(f)) next
  obj <- readRDS(paths[[d[2]]])
  obj$meta <- merge(as.data.table(obj$meta), lab[dataset == d[1], ..lab_cols], by = "cell", all.x = TRUE)[match(colnames(obj$counts), cell)]
  saveRDS(obj, f)
  rm(obj); free()
}

## =========================================================================
## G. Figures
## =========================================================================
base <- 7; ink <- "#0b0b0b"; ink2 <- "#52514e"; rule <- "#c3c2b7"
call_cols <- c(malignant = "#e34948", uncertain = "#eda100", non_malignant = "#4a3aa7")
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
## Fig 1: CNV heatmaps
heat_plot <- function(h, title) {
  g <- factor(h$group, levels = sort(unique(h$group)))
  o <- order(g)
  m <- h$s[, o]
  df <- data.table(gene = rep(seq_len(nrow(m)), ncol(m)), cell = rep(seq_len(ncol(m)), each = nrow(m)),
                   v = pmax(pmin(as.vector(m), 0.4), -0.4))
  runs <- split(seq_along(h$chr), factor(h$chr, levels = as.character(1:22)))
  ends <- vapply(runs, function(i) if (length(i)) max(i) else NA_real_, 0)
  mids <- vapply(runs, function(i) if (length(i)) mean(i) else NA_real_, 0)
  gn <- c(table(droplevels(g[o])))
  ggplot(df, aes(gene, cell, fill = v)) + geom_raster() +
    geom_vline(xintercept = na.omit(ends) + 0.5, linewidth = 0.1, colour = ink2) +
    geom_hline(yintercept = cumsum(gn)[-length(gn)] + 0.5, linewidth = 0.3) +
    scale_fill_gradient2(low = "#2a78d6", mid = "#f0efec", high = "#e34948", limits = c(-0.4, 0.4), name = "CNV") +
    scale_x_continuous(breaks = mids[c(TRUE, FALSE)], labels = names(runs)[c(TRUE, FALSE)], expand = c(0, 0)) +
    scale_y_continuous(breaks = cumsum(gn) - gn / 2, labels = names(gn), expand = c(0, 0)) +
    labs(title = title, x = "Chromosome", y = NULL) + theme_fig() +
    theme(axis.line = element_blank(), axis.ticks.y = element_blank())
}
hp <- Filter(Negate(is.null), cnv$heat)
titles <- c(BHN70 = "Atlas BHN70 (HPV-, larynx/hypopharynx)", Bill_HN23 = "Bill HN23 (HPV-, recurrence)",
            P15 = "Choi P15 (HPV-, oral cavity): tumour, normal mucosa, leukoplakia")
f1 <- wrap_plots(lapply(names(hp), function(u) heat_plot(hp[[u]], titles[[u]])), ncol = 1)
save_fig(f1, "Fig11_cnv_heatmaps", 180, 200)
## Fig 2: signal vs correlation
pd <- cells[!is.na(corr)][, .SD[sample(.N, min(.N, 15000))], by = dataset]
f2 <- ggplot(pd, aes(pmin(signal_rel, 6), corr, colour = cnv_call)) +
  geom_point_rast(size = 0.25, stroke = 0, shape = 16, raster.dpi = 300, dev = "ragg") +
  geom_vline(xintercept = 1, linetype = "22", linewidth = 0.3) + geom_hline(yintercept = cnv_corr_cut, linetype = "22", linewidth = 0.3) +
  facet_wrap(~ dataset, nrow = 1) + scale_colour_manual(values = call_cols, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  labs(title = "CNV signal vs correlation with the sample's tumour profile",
       subtitle = "dashed: 99th percentile of reference cells (x = 1) and correlation 0.4; up to 15,000 cells per dataset",
       x = "CNV signal (relative to reference 99th percentile, capped at 6)", y = "Correlation") + theme_fig()
save_fig(f2, "Fig12_cnv_signal_correlation", 180, 75)
## Fig 3: % malignant per sample
pu <- copy(usum)[, unit := factor(unit, levels = unit[order(dataset, -pct_malignant)])]
pu[, informative := factor(as.character(informative), levels = c("TRUE", "FALSE"))]
f3 <- ggplot(pu, aes(unit, pct_malignant, fill = hpv, alpha = informative)) + geom_col(width = 0.75) +
  facet_grid(~ dataset, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c("HPV-" = "#eb6834", "HPV+" = "#2a78d6", UNCLEAR = "#898781"), name = "HPV") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35), name = "Informative") +
  labs(title = "CNV-malignant epithelial cells per sample", x = NULL, y = "% of scored epithelial cells") +
  theme_fig() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = base - 2.5))
save_fig(f3, "Fig13_pct_malignant_per_sample", 180, 90)

## =========================================================================
## H. copyKAT on a subset (Gao et al. 2021) and concordance
## =========================================================================
ck_units <- c("BHN70", "BHN64", "BHN60", "OP12", "OP10", "HN01", "Bill_HN23", "Bill_HN43", "Bill_HN35",
              "Peng_SYSMH4", "Peng_SYSMH5")
ck_dir <- file.path(cdir, "copykat")
dir.create(ck_dir, showWarnings = FALSE)
ck_res <- list()
for (u in ck_units) {
  f <- file.path(ck_dir, paste0(u, "_pred.csv"))
  if (!file.exists(f)) {
    message("H. copyKAT: ", u)
    suppressPackageStartupMessages(library(copykat))
    ds <- units[[u]]$cells$dataset[1]
    ep <- units[[u]]$cells[!is.na(corr), cell]
    if (length(ep) > 600) ep <- sample(ep, 600)
    if (ds == "HNSCatlas") {
      a <- readRDS(paths$atlas_epi_v2)$counts; r <- readRDS(paths$cnv_input)
      rc <- as.data.table(r$meta)[group == "Reference" & Patient == u, cell]
      cnt <- cbind(a[, ep], r$counts[, rc]); rm(a, r)
    } else if (ds == "GSE234933_Bill") {
      s <- sub("^Bill_", "", u)
      raw <- readRDS(file.path(paths$bill_dir, "counts", paste0(s, "_gex_raw_counts.rds")))
      bc <- fread(paths$bill_cells)
      rc <- setdiff(bc[sample == s & subclustering != "Tumor cells", sample_barcode],
                    sub("^Bill_", "", dbl[dbl_class == "doublet", cell]))
      cnt <- raw[, c(sub("^Bill_", "", ep), rc)]
      colnames(cnt) <- paste0("Bill_", colnames(cnt)); rc <- paste0("Bill_", rc); rm(raw)
    } else {
      pa <- readRDS(paths$peng_all_v2); cl <- readRDS(paths$peng_epi_v2)$clusters
      pa$meta[, lineage := cl$lineage[match(cluster, cl$cluster)]]
      rc <- pa$meta[sample == sub("^Peng_", "", u) & lineage %in% c("immune", "fibroblast", "endothelial") &
                      !cell %in% dbl[dbl_class == "doublet", cell], cell]
      cnt <- pa$counts[, c(ep, rc)]; rm(pa)
    }
    if (length(rc) > 200) rc <- sample(rc, 200)
    cnt <- cnt[, c(ep, rc)]
    free()
    wd <- getwd(); dir.create(file.path(ck_dir, u), showWarnings = FALSE); setwd(file.path(ck_dir, u))
    ck <- tryCatch(copykat(rawmat = as.matrix(cnt), id.type = "S", ngene.chr = 5, win.size = 25, KS.cut = 0.1,
                           sam.name = u, distance = "euclidean", norm.cell.names = rc, plot.genes = "FALSE",
                           output.seg = "FALSE", genome = "hg20", n.cores = 4),
                   error = function(e) { message("copyKAT failed for ", u, ": ", conditionMessage(e)); NULL })
    setwd(wd)
    ## keep the prediction and heatmap only; the per-cell result tables are ~0.4 GB
    unlink(list.files(file.path(ck_dir, u), "_results_.*\\.txt$", full.names = TRUE))
    pred <- if (is.null(ck)) data.table(cell = character(), copykat_pred = character()) else
      data.table(cell = ck$prediction$cell.names, copykat_pred = ck$prediction$copykat.pred)
    fwrite(pred, f)
    rm(cnt, ck); free()
  }
  ck_res[[u]] <- fread(f)[, unit := u]
}
ck <- rbindlist(ck_res, fill = TRUE)
ck <- merge(ck, cells[, .(cell, cnv_call)], by = "cell")
conc <- ck[, .(epithelial_compared = .N,
               ours_malignant = sum(cnv_call == "malignant"), copykat_aneuploid = sum(copykat_pred == "aneuploid"),
               both_malignant = sum(cnv_call == "malignant" & copykat_pred == "aneuploid"),
               both_non = sum(cnv_call == "non_malignant" & copykat_pred == "diploid"),
               ours_mal_ck_diploid = sum(cnv_call == "malignant" & copykat_pred == "diploid"),
               ours_non_ck_aneuploid = sum(cnv_call == "non_malignant" & copykat_pred == "aneuploid"),
               ours_uncertain = sum(cnv_call == "uncertain"),
               ours_uncertain_ck_aneuploid = sum(cnv_call == "uncertain" & copykat_pred == "aneuploid"),
               ck_not_defined = sum(!copykat_pred %in% c("aneuploid", "diploid"))),
           by = unit]
conc[, agreement_pct := round(100 * (both_malignant + both_non) /
                                (both_malignant + both_non + ours_mal_ck_diploid + ours_non_ck_aneuploid), 1)]
fwrite(conc, file.path(out, "tables", "copykat_concordance.csv"))
print(conc)
message("Done.")
