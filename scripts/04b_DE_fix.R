## Step 04b: normal reference fix, gene symbol harmonisation, DE re-run with
## independent bulk replication
##
## Purpose: re-annotate all Choi normal mucosa / leukoplakia cells and gate
##   squamous epithelium per cell (also atlas normal cells); score new
##   epithelial cells for CNV against the step 03 references; map gene symbols
##   to current HGNC symbols; rebuild the per-patient cell-type pseudobulk with
##   v2 CNV labels; re-run pseudobulk DE (paired and all-patient designs);
##   replicate candidates in paired TCGA-HNSC tumour/normal RNA-seq; annotate
##   tumour specificity, detectability and proliferation; sensitivity re-runs.
## Inputs:  data/atlas/HNSCC_Atlas.rds (Part C only)
##          data/processed/HNSCatlas/atlas_cells_meta_umap.rds
##          data/processed/v2/{cell_labels_v2.csv.gz, atlas_epithelial_labelled.rds,
##          choi_nl_lp_allcells.rds}, data/processed/v2/cnv/ (step 03/04 caches)
##          step 03 reference-cell inputs (replayed for the Choi units)
##          HGNC complete set (downloaded), GDC TCGA-HNSC STAR counts (downloaded)
##          data/raw/step01_metadata/cbioportal_hnsc_tcga_pan_can_atlas_2018_patient.csv
## Outputs: results/step04b_DE/ (tables, figures)
##          data/processed/v2/{celltype_pseudobulk_v2.rds, choi_nl_lp_reannotated.csv.gz,
##          pseudobulk_step04b.rds, tcga_hnsc_paired_counts.rds}
##          data/processed/v2/cell_labels_v2.csv.gz (step 04b columns added)
## Date: 2026-10-25
## Run from the repo root:  Rscript scripts/04b_DE_fix.R
##   Rscript scripts/04b_DE_fix.R celltype   builds only the Part C cache
##   (loads the full atlas; run it alone first)

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
out <- paths$step04b
tdir <- file.path(out, "tables")
dir.create(tdir, recursive = TRUE, showWarnings = FALSE)
free <- function() invisible(gc(full = TRUE))
set.seed(1)

min_mal       <- 20      # cells per malignant pseudobulk
min_norm      <- 10      # cells per normal pseudobulk
min_type      <- 20      # cells per cell-type pseudobulk (specificity)
fdr_cut       <- 0.05
lfc_cut       <- 1
min_cpm       <- 10      # detectability: median CPM in HPV- malignant pseudobulks
min_share     <- 0.5
min_fold      <- 2
cont_max      <- 1000    # per-cell CPM limit for PTPRC, COL1A2, PECAM1 and glandular markers
tcga_max_gb   <- 1.5
prior_primary <- c("BHN39", "BHN46", "BHN49", "BHN50")
nasopharynx   <- c("BHN77", "BHN34", "BHN61", "BHN73")
pos_up        <- c("MAGEA4", "PRAME", "MMP1", "LAMC2", "PDPN")
pos_down      <- c("KRT4", "KRT13", "CRNN")
tme_types     <- c("CD4T", "CD8T", "NK", "Bcell", "PlasmaCell", "Macrophage", "Monocyte",
                   "panDC", "Neutrophils", "Mast", "Fibroblast", "Endothelial")
mk <- list(squamous = c("KRT5", "KRT14", "KRT13", "KRT4", "KRT6A", "KRT15", "SFN"),
           glandular = c("PIGR", "MUC5B", "BPIFB1"), immune = "PTPRC", fibroblast = c("COL1A2", "DCN"),
           endothelial = c("PECAM1", "VWF"), muscle = c("ACTA2", "DES"))

## =========================================================================
## C. Cell-type pseudobulk with v2 labels (full atlas; run alone first)
## =========================================================================
if (!file.exists(paths$celltype_pb_v2)) {
  message("C. Cell-type pseudobulk from the full atlas")
  lab_c <- fread(paths$cell_labels_v2, select = c("cell", "dataset", "cell_state_primary"))
  am <- as.data.table(readRDS(paths$atlas_cells_meta))[, .(cell, Patient, Source, Cell_Type)]
  am <- am[study_of(Patient) != "Kurten2021" & !Patient %in% nasopharynx]   # Kurten: CD45-sorted
  mal <- lab_c[dataset == "HNSCatlas" & cell_state_primary == "malignant", cell]
  am[, group := fifelse(cell %in% mal, "malignant", fifelse(Cell_Type %in% tme_types, Cell_Type, "other"))]
  rm(lab_c); free()
  suppressPackageStartupMessages(library(SeuratObject))
  atlas <- readRDS(paths$atlas)
  counts <- LayerData(atlas, assay = "RNA", layer = "counts")
  rm(atlas); free()
  counts <- counts[, am$cell]; free()
  pbc <- pseudobulk(counts, am[, .(cell, key = paste(Patient, group, sep = "|"))])
  rm(counts); free()
  smp <- data.table(id = colnames(pbc$counts), n_cells = as.integer(pbc$n_cells))
  smp[, `:=`(patient = sub("\\|.*$", "", id), group = sub("^.*\\|", "", id))]
  setnames(smp, "id", "key")
  saveRDS(list(counts = pbc$counts, samples = smp,
               note = "malignant = cell_state_primary (step 04); 12 TME types = atlas Cell_Type; other = all remaining cells"),
          paths$celltype_pb_v2)
  rm(pbc, am); free()
}
if ("celltype" %in% commandArgs(TRUE)) quit(save = "no")

## =========================================================================
## B. Gene symbol harmonisation (HGNC complete set)
## =========================================================================
if (!file.exists(paths$hgnc)) {
  dir.create(dirname(paths$hgnc), recursive = TRUE, showWarnings = FALSE)
  download.file("https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",
                paths$hgnc, mode = "wb", quiet = TRUE)
}
hg <- read_hgnc(paths$hgnc)
ds_genes <- list(HNSCatlas = rownames(readRDS(paths$atlas_labelled)$counts))
free()
ds_genes$GSE181919_Choi <- rownames(readRDS(paths$choi_nl_lp_all)$counts)
ds_genes$GSE234933_Bill <- rownames(readRDS(paths$bill_epi_v2)$counts)
ds_genes$GSE172577_Peng <- rownames(readRDS(paths$peng_epi_v2)$counts)
free()
maps <- lapply(ds_genes, hgnc_map, hg = hg)
map_all <- rbindlist(lapply(names(maps), function(d) maps[[d]][, dataset := d]))
fwrite(map_all[how != "approved"], file.path(tdir, "B_symbol_changes_all_datasets.csv"))
sym_sum <- map_all[, .(rows = .N, approved = sum(how == "approved"), renamed_previous = sum(how == "previous"),
                       renamed_alias = sum(how == "alias"), ambiguous_dropped = sum(grepl("ambiguous", how)),
                       not_in_hgnc_kept = sum(how == "not_in_hgnc"),
                       symbols_after = uniqueN(mapped[!is.na(mapped)])), by = dataset]
sym_sum[, rows_merged := rows - ambiguous_dropped - symbols_after]
fwrite(sym_sum, file.path(tdir, "B_symbol_summary.csv"))
ex_genes <- c("PYCRL", "PYCR3", "HN1L", "JPT2", "ECE2", "EEF1AKMT4-ECE2", "C16orf74")
sym_ex <- map_all[original %in% ex_genes | mapped %in% ex_genes, .(dataset, original, how, mapped)][order(mapped, dataset)]
fwrite(sym_ex, file.path(tdir, "B_symbol_examples.csv"))
map_atlas <- maps$HNSCatlas
map_choi <- maps$GSE181919_Choi

## =========================================================================
## A. Re-annotation of Choi normal mucosa (NL) and leukoplakia (LP)
## =========================================================================
units <- readRDS(file.path(paths$cnv_cache, "cnv_units.rds"))$units
rq <- fread(paths$ref_quantiles)
pos <- fread(file.path(paths$cnv_cache, "cnv_genes.csv"))
setorder(pos, chr_n, start)
cpm_of <- function(x, genes, lib) colSums(x[intersect(genes, rownames(x)), , drop = FALSE]) / lib * 1e6
gate_cols <- function(x, lib) data.table(
  squamous_keratin_cpm = cpm_of(x, mk$squamous, lib), PTPRC_cpm = cpm_of(x, "PTPRC", lib),
  COL1A2_cpm = cpm_of(x, "COL1A2", lib), PECAM1_cpm = cpm_of(x, "PECAM1", lib),
  glandular_cpm = cpm_of(x, mk$glandular, lib))
call_at <- function(signal, corr, thr) fcase(
  is.na(signal), NA_character_,
  is.na(corr), fifelse(signal > thr, "uncertain", "non_malignant"),
  signal > thr & corr > cnv_corr_cut, "malignant",
  signal <= thr & corr <= cnv_corr_cut, "non_malignant",
  default = "uncertain")

if (!file.exists(paths$choi_reannotated)) {
  message("A. Choi NL/LP: clustering and marker annotation")
  suppressPackageStartupMessages(library(Seurat))
  ch <- readRDS(paths$choi_nl_lp_all)
  cm <- as.data.table(ch$meta)[, .(cell, patient, sample, tissue, hpv, choi_label = celltype)]
  seu <- CreateSeuratObject(ch$counts, min.cells = 0, min.features = 0)
  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(seu, nfeatures = 2000, verbose = FALSE)
  seu <- ScaleData(seu, verbose = FALSE)
  seu <- RunPCA(seu, npcs = 30, verbose = FALSE)
  seu <- FindNeighbors(seu, dims = 1:20, verbose = FALSE)
  seu <- FindClusters(seu, resolution = 0.5, random.seed = 1, verbose = FALSE)
  cm[, cluster := as.character(seu$seurat_clusters)]
  rm(seu); free()
  lib <- colSums(ch$counts)
  ## cluster label: marker class with the highest cluster-mean log1p(CP10k)
  ## (>= 0.5), as for Peng in step 02 (detection rates are not used)
  lx <- ch$counts[unlist(mk), ]
  lx <- t(t(lx) / lib) * 1e4
  lx@x <- log1p(lx@x)
  avg <- t(sapply(split(seq_len(ncol(lx)), cm$cluster), function(i) Matrix::rowMeans(lx[, i, drop = FALSE])))
  cl <- data.table(cluster = rownames(avg), cells = as.integer(table(cm$cluster)[rownames(avg)]),
                   squamous = rowMeans(avg[, mk$squamous]), glandular = rowMeans(avg[, mk$glandular]),
                   immune = avg[, "PTPRC"], fibroblast = pmax(avg[, "COL1A2"], avg[, "DCN"]),
                   endothelial = pmax(avg[, "PECAM1"], avg[, "VWF"]), muscle = pmax(avg[, "ACTA2"], avg[, "DES"]))
  sc <- as.matrix(cl[, .(squamous, glandular, immune, fibroblast, endothelial, muscle)])
  cl[, lineage := fifelse(apply(sc, 1, max) >= 0.5, colnames(sc)[apply(sc, 1, which.max)], "unassigned")]
  cl[lineage == "squamous", lineage := "squamous_epithelial"]
  cm[, cluster_lineage := cl$lineage[match(cluster, cl$cluster)]]
  fwrite(merge(cl, cm[, .(choi_top_label = names(which.max(table(choi_label))),
                          pct_top = round(100 * max(table(choi_label)) / .N)), by = cluster], by = "cluster"),
         file.path(tdir, "A_choi_cluster_annotation.csv"))
  ## per-cell squamous gate: keratin threshold = density minimum between the
  ## ambient and the squamous mode of log10(summed keratin CPM + 1)
  cm <- cbind(cm, gate_cols(ch$counts, lib))
  d <- density(log10(cm[squamous_keratin_cpm > 0, squamous_keratin_cpm] + 1), bw = 0.08)
  w <- d$x > 2.6 & d$x < 4.3
  krt_thr <- round(10^d$x[w][which.min(d$y[w])] - 1)
  cm[, squamous_gate := squamous_keratin_cpm >= krt_thr & PTPRC_cpm < cont_max & COL1A2_cpm < cont_max &
                        PECAM1_cpm < cont_max & glandular_cpm < cont_max]
  cm[, new_label := fcase(squamous_gate, "squamous_epithelial",
                          cluster_lineage == "squamous_epithelial", "squamous_cluster_gate_fail",
                          default = cluster_lineage)]
  cm[, keratin_threshold := krt_thr]
  ## CNV for every gated cell against the step 03 references of its patient
  message("A5. CNV scoring of gated cells (step 03 references replayed)")
  refs <- list()
  replay_step03_references(paths, units, function(u, cnt, lib) if (grepl("^P[0-9]", u)) refs[[u]] <<- cnt,
                           upto = "choi")
  sq_cells <- cm[squamous_gate == TRUE]
  scr <- rbindlist(lapply(split(sq_cells$cell, sq_cells$patient), function(cc) {
    p <- cm[match(cc[1], cell), patient]
    r <- score_new_cells(ch$counts[rownames(refs[[p]]), cc, drop = FALSE], refs[[p]], units[[p]], pos)
    r[, `:=`(unit = p, in_cnv_reference = cell %in% colnames(refs[[p]]))]
  }))
  ## previously scored cells must reproduce the cached step 03 scores
  old <- rbindlist(lapply(units, `[[`, "cells"))[, .(cell, signal_old = signal, corr_old = corr)]
  chk <- merge(scr, old, by = "cell")
  stopifnot(max(abs(chk$signal - chk$signal_old)) < 1e-12, max(abs(chk$corr - chk$corr_old), na.rm = TRUE) < 1e-12)
  scr[, `:=`(thr_p95 = rq$q95[match(unit, rq$unit)], thr_p99 = rq$q99[match(unit, rq$unit)])]
  scr[, `:=`(cnv_call_p95 = call_at(signal, corr, thr_p95), cnv_call_p99 = call_at(signal, corr, thr_p99))]
  cm <- merge(cm, scr[, .(cell, signal, corr, in_cnv_reference, cnv_call_p95, cnv_call_p99)], by = "cell", all.x = TRUE)
  state_4b <- function(call) with(cm, fcase(
    !squamous_gate, "not_squamous",
    tissue == "LP" & call %in% "malignant", "LP_CNVpos",
    tissue == "LP" & call %in% "non_malignant", "LP_CNVneg",
    tissue == "LP", "LP_uncertain",
    call %in% "non_malignant", "normal_epithelial",
    default = "uncertain"))
  cm[, `:=`(cell_state_4b = state_4b(cnv_call_p95), cell_state_strict_4b = state_4b(cnv_call_p99))]
  fwrite(cm, paths$choi_reannotated)
  rm(ch, refs, lx); free()
}
cm <- fread(paths$choi_reannotated)
krt_thr <- cm$keratin_threshold[1]

conf <- dcast(cm, sample + tissue + hpv + choi_label ~ new_label, value.var = "cell", fun.aggregate = length)
fwrite(conf, file.path(tdir, "A_confusion_choi_label_vs_new.csv"))
conf_tot <- dcast(cm, choi_label ~ new_label, value.var = "cell", fun.aggregate = length)
fwrite(conf_tot, file.path(tdir, "A_confusion_totals.csv"))
choi_states <- cm[, .N, by = .(sample, tissue, hpv, cell_state_4b)][order(sample)]
fwrite(dcast(choi_states, sample + tissue + hpv ~ cell_state_4b, value.var = "N", fill = 0),
       file.path(tdir, "A_choi_states_per_sample.csv"))

## A6. same gate on atlas normal epithelial cells; step 04b labels
lab <- fread(paths$cell_labels_v2)
s4b_cols <- c("squamous_keratin_cpm", "squamous_gate", "cell_state_4b", "cell_state_strict_4b", "primary_4b")
if (length(intersect(s4b_cols, names(lab)))) lab[, intersect(s4b_cols, names(lab)) := NULL]
atl <- readRDS(paths$atlas_labelled)
an <- lab[dataset == "HNSCatlas" & (cell_state_primary == "normal_epithelial" | cell_state_strict == "normal_epithelial"), cell]
ag <- cbind(data.table(cell = an), gate_cols(atl$counts[, an], colSums(atl$counts[, an])))
rm(atl); free()
ag[, squamous_gate := squamous_keratin_cpm >= krt_thr & PTPRC_cpm < cont_max & COL1A2_cpm < cont_max &
                      PECAM1_cpm < cont_max & glandular_cpm < cont_max]
gate_all <- rbind(ag[, .(cell, squamous_keratin_cpm, squamous_gate)],
                  cm[, .(cell, squamous_keratin_cpm, squamous_gate)])
lab <- merge(lab, gate_all, by = "cell", all.x = TRUE)
lab[, cell_state_4b := fcase(
  dataset == "GSE181919_Choi", cm$cell_state_4b[match(cell, cm$cell)],
  cell_state_primary == "normal_epithelial" & !squamous_gate %in% TRUE, "normal_failed_gate",
  default = cell_state_primary)]
lab[, cell_state_strict_4b := fcase(
  dataset == "GSE181919_Choi", cm$cell_state_strict_4b[match(cell, cm$cell)],
  cell_state_strict == "normal_epithelial" & !squamous_gate %in% TRUE, "normal_failed_gate",
  default = cell_state_strict)]
lab[, primary_4b := primary_step04 & !patient_id_unique %in% prior_primary]
fwrite(lab, paths$cell_labels_v2)
atlas_gate <- merge(ag, lab[, .(cell, patient_id_unique, hpv, sample_type, cell_state_primary)], by = "cell")
fwrite(atlas_gate[, .(cells = .N, pass = sum(squamous_gate), fail_keratin = sum(squamous_keratin_cpm < krt_thr),
                      fail_contaminant = sum(squamous_keratin_cpm >= krt_thr & !squamous_gate)),
                  by = .(hpv, patient_id_unique, sample_type)][order(hpv, patient_id_unique)],
       file.path(tdir, "A_atlas_normal_gate.csv"))

## A7. normal reference per patient (computed in Part D from the pseudobulks)

## =========================================================================
## D. Pseudobulk and DE (normal >= 10 cells, malignant >= 20)
## =========================================================================
excl_normal <- "P9"   # P9 normal mucosa: MAGEA4/MMP1 positive; sensitivity only
groups_atl <- list(malignant = quote(cell_state_primary == "malignant"),
                   malignant_strict = quote(cell_state_strict == "malignant"),
                   normal = quote(cell_state_4b == "normal_epithelial"),
                   normal_strict = quote(cell_state_strict_4b == "normal_epithelial"))
groups_choi <- list(normal = quote(cell_state_4b == "normal_epithelial"),
                    normal_strict = quote(cell_state_strict_4b == "normal_epithelial"),
                    LP_CNVpos = quote(cell_state_4b == "LP_CNVpos"), LP_CNVneg = quote(cell_state_4b == "LP_CNVneg"))
key_cells <- function(d, groups) rbindlist(lapply(names(groups), function(g)
  d[eval(groups[[g]]), .(cell, key = paste(patient_id_unique, g, sep = "|"))]))
if (!file.exists(paths$pseudobulk_v2b)) {
  message("D. Pseudobulk")
  atl <- readRDS(paths$atlas_labelled)
  pa <- pseudobulk(atl$counts, key_cells(lab[dataset == "HNSCatlas" & role == "discovery" & !flag_excluded], groups_atl))
  rm(atl); free()
  ch <- readRDS(paths$choi_nl_lp_all)
  pc <- pseudobulk(ch$counts, key_cells(copy(cm)[, patient_id_unique := patient], groups_choi))
  rm(ch); free()
  fa <- map_rows(pa$counts, map_atlas)
  fc <- map_rows(pc$counts, map_choi)
  g_sh <- intersect(rownames(fa), rownames(fc))
  keys <- union(colnames(fa), colnames(fc))
  sh <- Matrix(0, length(g_sh), length(keys), dimnames = list(g_sh, keys), sparse = TRUE)
  sh[, colnames(fa)] <- fa[g_sh, ]
  sh[, colnames(fc)] <- sh[, colnames(fc)] + fc[g_sh, ]
  n_sh <- setNames(numeric(length(keys)), keys)
  n_sh[names(pa$n_cells)] <- pa$n_cells
  n_sh[names(pc$n_cells)] <- n_sh[names(pc$n_cells)] + pc$n_cells
  saveRDS(list(full = list(counts = fa, n_cells = pa$n_cells), shared = list(counts = sh, n_cells = n_sh),
               source = list(atlas = colnames(fa), choi = colnames(fc))), paths$pseudobulk_v2b)
  rm(pa, pc, fa, fc, sh); free()
}
pb <- readRDS(paths$pseudobulk_v2b)
ptinfo <- lab[role == "discovery" & !flag_excluded,
              .(hpv = hpv[1], study = study[1], primary_4b = primary_4b[1], sensitivity = sensitivity_reason[1]),
              by = patient_id_unique]
ptinfo[patient_id_unique == "BHN50", sensitivity := "prior head-and-neck primary"]
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
sm_sh[, squamous_keratin_cpm := round(colSums(as.matrix(pb$shared$counts[intersect(mk$squamous, rownames(pb$shared$counts)), key, drop = FALSE])) / lib * 1e6)]
fwrite(sm_sh[order(hpv, group, patient)], file.path(tdir, "D_pseudobulk_samples.csv"))

## A7 / D1. normal reference per patient and patients gained with >= 10 cells
s04 <- fread(file.path(paths$step04, "tables", "A_discovery_patients.csv"))
nref <- sm_sh[group %in% c("normal", "normal_strict"), .(patient, hpv, group, source, n_cells, squamous_keratin_cpm, primary_4b)]
nref <- dcast(nref, patient + hpv + source + primary_4b ~ group, value.var = c("n_cells", "squamous_keratin_cpm"))
if (!"n_cells_normal" %in% names(nref)) nref[, n_cells_normal := 0L]
nref[is.na(n_cells_normal), n_cells_normal := 0L]
nref[, step04_normal_cells := s04$normal[match(patient, s04$patient_id_unique)]]
nref[, `:=`(pb_primary = n_cells_normal >= min_norm & primary_4b & patient != excl_normal,
            pb_at_20 = n_cells_normal >= 20,
            step04_pb = step04_normal_cells >= 20 & !is.na(step04_normal_cells))]
setorder(nref, hpv, -n_cells_normal)
fwrite(nref, file.path(tdir, "A_normal_reference_per_patient.csv"))

## DE sets
make_set <- function(hpv_arm, mal = "malignant", nor = "normal", pts_ok = NULL, paired = TRUE, nor_drop = excl_normal) {
  s <- sm_sh[hpv == hpv_arm & (primary_4b | patient %in% pts_ok)]
  s <- s[!(group %in% nor & patient %in% nor_drop)]
  agg <- function(gs, state, mn) {
    x <- s[group %in% gs, .(keys = list(key), n_cells = sum(n_cells)), by = .(patient, study)]
    x[n_cells >= mn][, state := state]
  }
  ss <- rbind(agg(mal, "malignant", min_mal), agg(nor, "normal", min_norm))
  if (paired) ss <- ss[patient %in% intersect(ss[state == "malignant", patient], ss[state == "normal", patient])]
  ss[, sample := paste(patient, state, sep = "|")]
  cnt <- do.call(cbind, lapply(ss$keys, function(k) Matrix::rowSums(pb$shared$counts[, k, drop = FALSE])))
  colnames(cnt) <- ss$sample
  ss[, `:=`(state = factor(state, levels = c("normal", "malignant")), patient = factor(patient), study = factor(study))]
  list(smp = ss[, !"keys"], counts = cnt)
}
run_c1 <- function(set) de_qlf(set$counts, set$smp, ~ patient + state, "statemalignant", fdr_cut, lfc_cut)
run_c2 <- function(set) {
  tab <- table(set$smp$study, set$smp$state)
  if (!any(tab[, "normal"] > 0 & tab[, "malignant"] > 0)) stop("C2: no study has both states")
  de_qlf(set$counts, set$smp, ~ study + state, "statemalignant", fdr_cut, lfc_cut)
}
save_de <- function(res, name) { fwrite(res$table, file.path(tdir, paste0(name, ".csv"))); res$table }

## study_zero: zero counts in every malignant pseudobulk (>= 20 cells) of a study
mz <- sm_full[group == "malignant" & n_cells >= min_mal]
st_sum <- sapply(split(mz$key, mz$study), function(k) Matrix::rowSums(pb$full$counts[, k, drop = FALSE]))
study_zero <- apply(st_sum == 0, 1, function(z) paste(colnames(st_sum)[z], collapse = " "))

## candidates: significant in C2 and same sign in C1
candidates <- function(c1, c2) {
  t <- merge(c2$table[, .(gene, C2_logFC = logFC, C2_logCPM = logCPM, C2_FDR = FDR, C2_dir = direction)],
             c1$table[, .(gene, C1_logFC = logFC, C1_FDR = FDR, C1_dir = direction)], by = "gene", all = TRUE)
  t[, study_zero := study_zero[gene]]
  t[is.na(study_zero), study_zero := ""]
  t[, candidate := fcase(C2_dir %in% "up" & C1_logFC > 0, "up", C2_dir %in% "down" & C1_logFC < 0, "down",
                         default = "")]
  t[, C2_sig_not_tested_C1 := C2_dir %in% c("up", "down") & is.na(C1_logFC)]
  t[, candidate_list := candidate != "" & study_zero == ""]
  t[]
}

message("D. DE")
res <- list()
for (arm in c("HPV-", "HPV+")) {
  s1 <- make_set(arm); s2 <- make_set(arm, paired = FALSE)
  r1 <- run_c1(s1); r2 <- run_c2(s2)
  tag <- if (arm == "HPV-") "HPVneg" else "HPVpos"
  save_de(r1, paste0("D_C1_", tag, "_paired")); save_de(r2, paste0("D_C2_", tag, "_allpatients"))
  cd <- candidates(r1, r2)
  fwrite(cd, file.path(tdir, paste0("D_candidates_", tag, ".csv")))
  res[[arm]] <- list(s1 = s1, s2 = s2, c1 = r1, c2 = r2, cand = cd)
}
de_counts <- rbindlist(lapply(names(res), function(a) with(res[[a]], data.table(
  hpv = a, C1_patients = uniqueN(s1$smp$patient), C1_genes = c1$n_genes_tested,
  C1_up = sum(c1$table$direction == "up"), C1_down = sum(c1$table$direction == "down"),
  C2_malignant_pb = sum(s2$smp$state == "malignant"), C2_normal_pb = sum(s2$smp$state == "normal"),
  C2_genes = c2$n_genes_tested, C2_up = sum(c2$table$direction == "up"), C2_down = sum(c2$table$direction == "down"),
  cand_up = sum(cand$candidate == "up"), cand_down = sum(cand$candidate == "down"),
  cand_up_list = sum(cand$candidate_list & cand$candidate == "up"),
  cand_down_list = sum(cand$candidate_list & cand$candidate == "down"),
  C2_sig_not_tested_C1 = sum(cand$C2_sig_not_tested_C1),
  study_zero_removed = sum(cand$candidate != "" & cand$study_zero != ""),
  normal_pb = paste(sort(as.character(s2$smp[state == "normal", patient])), collapse = " ")))))
fwrite(de_counts, file.path(tdir, "D_de_counts.csv"))
study_state <- rbindlist(lapply(names(res), function(a)
  dcast(res[[a]]$s2$smp, study ~ state, value.var = "sample", fun.aggregate = length)[, hpv := a]))
fwrite(study_state, file.path(tdir, "D_C2_study_by_state.csv"))
cand <- res[["HPV-"]]$cand

## =========================================================================
## E. Independent bulk replication: paired TCGA-HNSC (GDC STAR counts)
## =========================================================================
if (!file.exists(paths$tcga_counts)) {
  message("E. TCGA-HNSC: GDC file query")
  flt <- jsonlite::toJSON(list(op = "and", content = list(
    list(op = "in", content = list(field = "cases.project.project_id", value = list("TCGA-HNSC"))),
    list(op = "in", content = list(field = "data_type", value = list("Gene Expression Quantification"))),
    list(op = "in", content = list(field = "analysis.workflow_type", value = list("STAR - Counts"))),
    list(op = "in", content = list(field = "cases.samples.sample_type", value = list("Primary Tumor", "Solid Tissue Normal"))))),
    auto_unbox = TRUE)
  q <- httr::POST("https://api.gdc.cancer.gov/files", encode = "json",
                  body = list(filters = jsonlite::fromJSON(flt, simplifyVector = FALSE), format = "JSON", size = "5000",
                              fields = "file_id,file_name,file_size,cases.submitter_id,cases.samples.sample_type,cases.samples.submitter_id"))
  hits <- jsonlite::fromJSON(httr::content(q, as = "text", encoding = "UTF-8"))$data$hits
  man <- rbindlist(lapply(seq_len(nrow(hits)), function(i) data.table(
    file_id = hits$file_id[i], file_name = hits$file_name[i], file_size = hits$file_size[i],
    case = hits$cases[[i]]$submitter_id[1], sample_type = hits$cases[[i]]$samples[[1]]$sample_type[1],
    sample = hits$cases[[i]]$samples[[1]]$submitter_id[1])))
  paired <- intersect(man[sample_type == "Primary Tumor", case], man[sample_type == "Solid Tissue Normal", case])
  man <- man[case %in% paired]
  ## one file per case and sample type (first sample ID, vial A first)
  setorder(man, case, sample_type, sample, file_id)
  man <- man[, .SD[1], by = .(case, sample_type)]
  gb <- sum(man$file_size) / 1e9
  message(sprintf("E. %d paired cases, %d files, %.2f GB", length(paired), nrow(man), gb))
  if (gb > tcga_max_gb) stop(sprintf("TCGA download %.2f GB exceeds %.1f GB", gb, tcga_max_gb))
  dir.create(paths$tcga_dir, recursive = TRUE, showWarnings = FALSE)
  for (i in seq_len(nrow(man))) {
    f <- file.path(paths$tcga_dir, man$file_name[i])
    if (file.exists(f) && file.size(f) == man$file_size[i]) next
    for (try in 1:3) {
      ok <- tryCatch({ download.file(paste0("https://api.gdc.cancer.gov/data/", man$file_id[i]), f, mode = "wb", quiet = TRUE); TRUE },
                     error = function(e) FALSE)
      if (ok && file.size(f) == man$file_size[i]) break
    }
    if (file.size(f) != man$file_size[i]) stop("download failed: ", man$file_id[i])
  }
  rd <- lapply(man$file_name, function(fn) {
    x <- fread(file.path(paths$tcga_dir, fn), skip = 1)
    x[grepl("^ENSG", gene_id), .(gene_id, gene_name, unstranded)]
  })
  stopifnot(all(vapply(rd, function(x) identical(x$gene_id, rd[[1]]$gene_id), TRUE)))
  m <- do.call(cbind, lapply(rd, `[[`, "unstranded"))
  colnames(m) <- man$sample
  ens <- sub("\\..*$", "", rd[[1]]$gene_id)
  keep <- !grepl("_PAR_Y$", rd[[1]]$gene_id)
  ## Ensembl ID -> HGNC approved symbol; otherwise the GDC gene name, mapped
  sym <- hg$symbol[match(ens, hg$ensembl_gene_id)]
  gmap <- hgnc_map(unique(rd[[1]]$gene_name), hg)
  sym[is.na(sym)] <- gmap$mapped[match(rd[[1]]$gene_name[is.na(sym)], gmap$original)]
  keep <- keep & !is.na(sym)
  f <- factor(sym[keep])
  agg <- sparseMatrix(i = as.integer(f), j = which(keep), x = 1, dims = c(nlevels(f), nrow(m)), dimnames = list(levels(f), NULL))
  tc <- as.matrix(agg %*% m)
  fwrite(man, file.path(tdir, "E_tcga_manifest.csv"))
  saveRDS(list(counts = tc, samples = man, gb = gb), paths$tcga_counts)
  rm(rd, m, agg); free()
}
tc <- readRDS(paths$tcga_counts)
cb <- fread(paths$cbio_pancan, select = c("patientId", "SUBTYPE"))
ts <- copy(tc$samples)[, `:=`(hpv = fcase(cb$SUBTYPE[match(case, cb$patientId)] == "HNSC_HPV-", "HPV-",
                                          cb$SUBTYPE[match(case, cb$patientId)] == "HNSC_HPV+", "HPV+", default = "UNCLEAR"),
                              state = factor(fifelse(sample_type == "Primary Tumor", "tumour", "normal"), levels = c("normal", "tumour")))]
tcga_pairs <- ts[, .(cases = uniqueN(case)), by = hpv]
fwrite(tcga_pairs, file.path(tdir, "E_tcga_pairs_by_hpv.csv"))
tcga_de <- list()
for (arm in c("HPV-", "HPV+")) {
  s <- ts[hpv == arm][, patient := factor(case)]
  tcga_de[[arm]] <- de_qlf(tc$counts[, s$sample], s, ~ patient + state, "statetumour", fdr_cut, lfc_cut)
  save_de(tcga_de[[arm]], paste0("E_TCGA_", if (arm == "HPV-") "HPVneg" else "HPVpos", "_paired_tumour_vs_normal"))
}
replicate <- function(cd, tt) {
  x <- merge(cd, tt$table[, .(gene, TCGA_logFC = logFC, TCGA_FDR = FDR)], by = "gene", all.x = TRUE)
  x[, TCGA_replicated := !is.na(TCGA_FDR) & TCGA_FDR < fdr_cut &
                         sign(TCGA_logFC) == fifelse(candidate == "up" | C2_logFC > 0, 1, -1)]
  x[]
}
cand <- replicate(cand, tcga_de[["HPV-"]])
candp <- replicate(res[["HPV+"]]$cand, tcga_de[["HPV+"]])
rep_rate <- rbind(
  cand[candidate_list == TRUE, .(hpv = "HPV-", candidates = .N, tested_TCGA = sum(!is.na(TCGA_FDR)),
                                 replicated = sum(TCGA_replicated), pct = round(100 * mean(TCGA_replicated), 1)), by = candidate],
  candp[candidate_list == TRUE, .(hpv = "HPV+", candidates = .N, tested_TCGA = sum(!is.na(TCGA_FDR)),
                                  replicated = sum(TCGA_replicated), pct = round(100 * mean(TCGA_replicated), 1)), by = candidate])
fwrite(rep_rate, file.path(tdir, "E_replication_rate.csv"))
pc <- rbindlist(lapply(c(pos_up, pos_down), function(g) {
  r <- cand[gene == g]; t <- tcga_de[["HPV-"]]$table[gene == g]
  data.table(gene = g, expected = if (g %in% pos_up) "up" else "down",
             C1_logFC = if (nrow(r)) round(r$C1_logFC, 2) else NA, C1_FDR = if (nrow(r)) signif(r$C1_FDR, 2) else NA,
             C2_logFC = if (nrow(r)) round(r$C2_logFC, 2) else NA, C2_FDR = if (nrow(r)) signif(r$C2_FDR, 2) else NA,
             candidate = if (nrow(r)) r$candidate else "",
             TCGA_logFC = if (nrow(t)) round(t$logFC, 2) else NA, TCGA_FDR = if (nrow(t)) signif(t$FDR, 2) else NA)
}))
pc[, `:=`(sc_pass = candidate == expected,
          TCGA_pass = !is.na(TCGA_FDR) & TCGA_FDR < fdr_cut & sign(TCGA_logFC) == fifelse(expected == "up", 1, -1))]
fwrite(pc, file.path(tdir, "D_E_positive_controls.csv"))

## =========================================================================
## F. Tumour specificity (cell-type pseudobulk v2), detectability, proliferation
## =========================================================================
ct <- readRDS(paths$celltype_pb_v2)
ctc <- map_rows(ct$counts, map_atlas)
cs <- ct$samples[patient %in% ptinfo[hpv == "HPV-" & primary_4b == TRUE, patient_id_unique]]
lib_ct <- colSums(ctc)
genes_f <- intersect(unique(c(cand[candidate != "", gene], res[["HPV-"]]$c2$table[direction == "up", gene])), rownames(ctc))
cpm_ct <- as.matrix(ctc[genes_f, cs$key, drop = FALSE]) %*% diag(1e6 / lib_ct[cs$key], nrow = nrow(cs))
colnames(cpm_ct) <- cs$key
med_type <- function(g) {
  j <- cs[group == g & n_cells >= min_type, key]
  if (!length(j)) return(rep(NA_real_, length(genes_f)))
  apply(cpm_ct[, j, drop = FALSE], 1, median)
}
tum <- med_type("malignant")
tme <- sapply(tme_types, med_type)
max_tme <- apply(tme, 1, max, na.rm = TRUE)
max_tme_type <- tme_types[apply(tme, 1, which.max)]
share_pts <- cs[group == "malignant" & n_cells >= min_type, patient]
share <- sapply(share_pts, function(p) {
  tot <- Matrix::rowSums(ctc[genes_f, cs[patient == p, key], drop = FALSE])
  ifelse(tot > 0, ctc[genes_f, cs[patient == p & group == "malignant", key]] / tot, NA_real_)
})
tissue_share <- apply(share, 1, median, na.rm = TRUE)
mal_full <- sm_full[group == "malignant" & hpv == "HPV-" & primary_4b & n_cells >= min_mal]
lcpm_mal <- edgeR::cpm(DGEList(as.matrix(pb$full$counts[, mal_full$key])), log = TRUE, prior.count = 2)
cpm_mal <- as.matrix(pb$full$counts[genes_f, mal_full$key, drop = FALSE]) %*% diag(1e6 / mal_full$lib, nrow = nrow(mal_full))
med_mal <- apply(cpm_mal, 1, median)
## proliferation: Tirosh et al. 2016 G1/S and G2/M lists (Seurat
## cc.genes.updated.2019), mapped to HGNC; correlation with their mean
## log CPM across HPV- malignant pseudobulks
cc <- Seurat::cc.genes.updated.2019
cc_g <- unique(na.omit(hgnc_map(c(cc$s.genes, cc$g2m.genes), hg)$mapped))
cc_score <- colMeans(lcpm_mal[intersect(cc_g, rownames(lcpm_mal)), ])
cc_r <- apply(lcpm_mal[genes_f, , drop = FALSE], 1, function(v) suppressWarnings(cor(v, cc_score, method = "spearman")))
ann <- data.table(gene = genes_f, malignant_median_CPM = tum, max_TME_CPM = max_tme, max_TME_type = max_tme_type,
                  fold_vs_max_TME = tum / max_tme, tumour_share = tissue_share, malignant_pb_median_CPM = med_mal,
                  cell_cycle_list = genes_f %in% cc_g, cell_cycle_corr = round(cc_r, 2))
ann <- cbind(ann, setNames(as.data.table(round(tme, 2)), paste0("CPM_", tme_types)))
ann[, specificity := fcase(malignant_median_CPM >= min_fold * max_TME_CPM & tumour_share >= min_share, "tumour_specific",
                           malignant_median_CPM > max_TME_CPM, "tumour_enriched", default = "shared_with_TME")]
ann[, `:=`(detectable = malignant_pb_median_CPM >= min_cpm, proliferation = cell_cycle_list | cell_cycle_corr >= 0.7)]
final <- merge(cand[candidate != ""], ann, by = "gene", all.x = TRUE)
final[, tier := fcase(candidate == "up" & candidate_list & TCGA_replicated & specificity %in% "tumour_specific" & detectable %in% TRUE, "1_specific",
                      candidate == "up" & candidate_list & TCGA_replicated & specificity %in% "tumour_enriched" & detectable %in% TRUE, "2_enriched",
                      default = "")]
final[tier == "", tier := NA]
setorder(final, -candidate, tier, C2_FDR, na.last = TRUE)
fwrite(final, file.path(tdir, "F_final_candidates_HPVneg.csv"))
fwrite(ann, file.path(tdir, "F_specificity_annotation.csv"))
f_sum <- final[candidate == "up" & candidate_list == TRUE,
               .N, by = .(TCGA_replicated, specificity, detectable, proliferation)][order(-TCGA_replicated, specificity)]
fwrite(f_sum, file.path(tdir, "F_up_candidate_summary.csv"))
fwrite(data.table(patients_specificity = uniqueN(cs$patient), tumour_share_patients = length(share_pts),
                  malignant_pb_detectability = nrow(mal_full), cell_cycle_genes = length(cc_g)),
       file.path(tdir, "F_patients_used.csv"))

## =========================================================================
## G. Sensitivity (primary candidate definition re-applied)
## =========================================================================
cand_set <- function(s1, s2) candidates(run_c1(s1), run_c2(s2))
add_back <- c("BHN39", "BHN46", "BHN49", "OP8")
sens_runs <- list(
  "(a) strict CNV labels" = list(make_set("HPV-", mal = "malignant_strict", nor = "normal_strict"),
                                 make_set("HPV-", mal = "malignant_strict", nor = "normal_strict", paired = FALSE)),
  "(b) P9 normal included" = list(make_set("HPV-", nor_drop = NULL), make_set("HPV-", nor_drop = NULL, paired = FALSE)),
  "(c) BHN39/46/49 + OP8 added" = list(make_set("HPV-", pts_ok = add_back), make_set("HPV-", pts_ok = add_back, paired = FALSE)))
ref_up <- cand[candidate_list & candidate == "up", gene]
ref_dn <- cand[candidate_list & candidate == "down", gene]
ref_t1 <- final[tier %in% c("1_specific", "2_enriched"), gene]
sens <- rbindlist(lapply(names(sens_runs), function(n) {
  sr <- sens_runs[[n]]
  cd <- cand_set(sr[[1]], sr[[2]])
  fwrite(cd, file.path(tdir, paste0("G_candidates_", gsub("[^a-z0-9]+", "_", tolower(n)), ".csv")))
  up <- cd[candidate_list & candidate == "up", gene]; dn <- cd[candidate_list & candidate == "down", gene]
  ov <- function(a, b) c(length(intersect(a, b)), round(length(intersect(a, b)) / length(union(a, b)), 3),
                         round(100 * length(intersect(a, b)) / length(b), 1))
  data.table(analysis = n, C1_patients = uniqueN(sr[[1]]$smp$patient),
             C2_malignant = sum(sr[[2]]$smp$state == "malignant"), C2_normal = sum(sr[[2]]$smp$state == "normal"),
             up = length(up), up_shared = ov(up, ref_up)[1], up_jaccard = ov(up, ref_up)[2], up_pct_retained = ov(up, ref_up)[3],
             down = length(dn), down_jaccard = ov(dn, ref_dn)[2], down_pct_retained = ov(dn, ref_dn)[3],
             tier12_retained = sum(ref_t1 %in% up), tier12_total = length(ref_t1))
}))
fwrite(sens, file.path(tdir, "G_sensitivity_overlap.csv"))

## =========================================================================
## H. Figures
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
spec_cols <- c(tumour_specific = "#e34948", tumour_enriched = "#eda100", shared_with_TME = "#898781")

## Fig 17: squamous keratin gate
gd <- rbind(cm[, .(set = "Choi NL/LP, all cells", v = log10(squamous_keratin_cpm + 1), pass = squamous_gate)],
            ag[, .(set = "Atlas normal epithelial (step 04)", v = log10(squamous_keratin_cpm + 1), pass = squamous_gate)])
f17 <- ggplot(gd, aes(v, fill = pass)) + geom_histogram(bins = 70, linewidth = 0) +
  geom_vline(xintercept = log10(krt_thr + 1), linetype = "22", linewidth = 0.3) +
  facet_wrap(~ set, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c(`TRUE` = "#2a78d6", `FALSE` = "#c3c2b7"), name = "Squamous gate") +
  scale_y_continuous(trans = "log1p", breaks = c(0, 10, 100, 1000, 10000)) +
  labs(title = "Per-cell squamous keratin gate",
       subtitle = sprintf("summed CPM of %s\nthreshold %s CPM (density minimum); PTPRC, COL1A2, PECAM1 and glandular markers < %d CPM",
                          paste(mk$squamous, collapse = "/"), format(krt_thr, big.mark = ","), cont_max),
       x = "log10(squamous keratin CPM + 1)", y = "Cells (log1p scale)") + theme_fig()
save_fig(f17, "Fig17_squamous_keratin_gate", 160, 100)

## Fig 18: volcano C2 (HPV-)
v <- merge(res[["HPV-"]]$c2$table, cand[, .(gene, candidate, candidate_list)], by = "gene")
v[, cls := fcase(candidate_list & candidate == "up", "up candidate", candidate_list & candidate == "down", "down candidate",
                 direction != "ns", "C2 only / study_zero", default = "n.s.")]
lab_g <- unique(c(final[tier %in% c("1_specific", "2_enriched")][1:min(.N, 20), gene], intersect(c(pos_up, pos_down), v$gene)))
f18 <- ggplot(v, aes(logFC, -log10(PValue), colour = cls)) +
  ggrastr::geom_point_rast(size = 0.4, stroke = 0, shape = 16, raster.dpi = 300, dev = "ragg") +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "22", linewidth = 0.3, colour = ink2) +
  geom_text_repel(data = v[gene %in% lab_g], aes(label = gene, fontface = fifelse(gene %in% c(pos_up, pos_down), "bold", "plain")),
                  size = 2, colour = ink, max.overlaps = 40, segment.size = 0.2, min.segment.length = 0) +
  scale_colour_manual(values = c(`up candidate` = "#e34948", `down candidate` = "#2a78d6", `C2 only / study_zero` = "#eda100", n.s. = "#c3c2b7"), name = NULL) +
  labs(title = sprintf("C2: HPV- malignant vs normal epithelium (%d tumour, %d normal pseudobulks)",
                       de_counts[hpv == "HPV-", C2_malignant_pb], de_counts[hpv == "HPV-", C2_normal_pb]),
       subtitle = "edgeR QL, ~ study + state; candidate = C2 FDR < 0.05, |log2FC| >= 1, same sign in paired C1; bold = positive controls",
       x = "log2 fold change (malignant / normal)", y = "-log10 P") + theme_fig()
save_fig(f18, "Fig18_volcano_C2", 160, 110)

## Fig 19: specificity of up candidates
sp <- final[candidate == "up" & candidate_list == TRUE & !is.na(specificity)]
sp[, `:=`(x = log10(max_TME_CPM + 1), y = log10(malignant_median_CPM + 1),
          rep = factor(fifelse(TCGA_replicated, "TCGA replicated", "not replicated"), levels = c("TCGA replicated", "not replicated")))]
f19 <- ggplot(sp, aes(x, y, colour = specificity)) + geom_abline(intercept = 0, slope = 1, linewidth = 0.3, colour = ink2) +
  geom_point(aes(shape = detectable), size = 1, alpha = 0.85) +
  geom_text_repel(data = sp[tier %in% "1_specific"][1:min(.N, 25)], aes(label = gene), size = 2, colour = ink, segment.size = 0.2, max.overlaps = 40) +
  facet_wrap(~ rep) + scale_colour_manual(values = spec_cols, name = NULL) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), name = sprintf("median CPM >= %d", min_cpm)) +
  labs(title = "Tumour-cell specificity of HPV- up candidates (v2 CNV-malignant cells)",
       subtitle = sprintf("median CPM over HPV- patients of unsorted studies (>= %d cells per type); line: equal expression", min_type),
       x = "log10(CPM + 1), highest microenvironment cell type", y = "log10(CPM + 1), malignant cells") + theme_fig()
save_fig(f19, "Fig19_specificity_scatter", 180, 95)

## Fig 20: heatmap of final candidates (tier 1, then tier 2; up to 40)
hg_genes <- final[tier %in% c("1_specific", "2_enriched")][1:min(.N, 40), gene]
if (length(hg_genes)) {
  hs <- sm_sh[hpv == "HPV-" & primary_4b & ((group == "malignant" & n_cells >= min_mal) | (group == "normal" & n_cells >= min_norm))]
  lc <- edgeR::cpm(DGEList(as.matrix(pb$shared$counts[, hs$key])), log = TRUE, prior.count = 2)[hg_genes, , drop = FALSE]
  z <- t(scale(t(lc)))
  ord <- hs[order(-xtfrm(group), study, patient), key]
  hd <- data.table(gene = rep(rownames(z), ncol(z)), id = rep(colnames(z), each = nrow(z)), z = pmax(pmin(as.vector(z), 2.5), -2.5))
  setnames(hd, "id", "key")
  hd[, panel := factor(fifelse(grepl("[|]malignant$", key), "malignant", "normal epithelium"), levels = c("malignant", "normal epithelium"))]
  hd[, `:=`(gene = factor(gene, levels = rev(hg_genes)), key = factor(key, levels = ord))]
  tier_lab <- final$tier[match(hg_genes, final$gene)]
  f20 <- ggplot(hd, aes(key, gene, fill = z)) + geom_tile() +
    scale_fill_gradient2(low = "#2a78d6", mid = "#f0efec", high = "#e34948", name = "z (log CPM)") +
    scale_x_discrete(labels = function(k) paste0(sub("[|].*$", "", k), fifelse(grepl(paste0("^", excl_normal, "[|]normal"), k), "*", ""))) +
    scale_y_discrete(labels = function(g) paste0(g, fifelse(final$tier[match(g, final$gene)] == "1_specific", "", " (enriched)"))) +
    facet_grid(~ panel, scales = "free_x", space = "free_x") +
    labs(title = "Final HPV- candidates: up, TCGA-replicated, detectable (tumour-specific first)",
         subtitle = sprintf("HPV- primary pseudobulks (malignant >= %d, normal >= %d cells); %d HGNC-harmonised genes shared with Choi",
                            min_mal, min_norm, nrow(pb$shared$counts)), x = NULL, y = NULL) +
    theme_fig() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = base - 2),
                        axis.line = element_blank(), axis.ticks = element_blank())
  save_fig(f20, "Fig20_heatmap_final_candidates", 180, 30 + 3.2 * length(hg_genes))
}

## Fig 21: single-cell C2 vs TCGA log2FC
sv <- merge(res[["HPV-"]]$c2$table[, .(gene, sc = logFC)], tcga_de[["HPV-"]]$table[, .(gene, tcga = logFC)], by = "gene")
sv <- merge(sv, cand[, .(gene, candidate, candidate_list, TCGA_replicated)], by = "gene", all.x = TRUE)
sv[, cls := fcase(candidate_list & candidate == "up" & TCGA_replicated, "up candidate, replicated",
                  candidate_list & candidate == "down" & TCGA_replicated, "down candidate, replicated",
                  candidate_list %in% TRUE, "candidate, not replicated", default = "other genes")]
r_all <- cor(sv$sc, sv$tcga, method = "spearman"); r_c <- sv[candidate_list %in% TRUE, cor(sc, tcga, method = "spearman")]
f21 <- ggplot(sv[order(cls != "other genes")], aes(sc, tcga, colour = cls)) +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = rule) + geom_vline(xintercept = 0, linewidth = 0.2, colour = rule) +
  ggrastr::geom_point_rast(size = 0.5, stroke = 0, shape = 16, raster.dpi = 300, dev = "ragg") +
  geom_text_repel(data = sv[gene %in% lab_g], aes(label = gene), size = 2, colour = ink, segment.size = 0.2, max.overlaps = 40) +
  scale_colour_manual(values = c(`up candidate, replicated` = "#e34948", `down candidate, replicated` = "#2a78d6",
                                 `candidate, not replicated` = "#eda100", `other genes` = "#dcdbd5"), name = NULL) +
  labs(title = "Single-cell (C2) vs TCGA-HNSC paired bulk log2 fold change, HPV-",
       subtitle = sprintf("%d genes tested in both; Spearman r = %.2f (all), %.2f (candidates); TCGA: %d HPV- pairs, ~ patient + tissue",
                          nrow(sv), r_all, r_c, tcga_pairs[hpv == "HPV-", cases]),
       x = "log2FC malignant / normal epithelium (pseudobulk C2)", y = "log2FC tumour / normal (TCGA)") + theme_fig()
save_fig(f21, "Fig21_sc_vs_TCGA_log2FC", 150, 110)

## console summary
print(sym_sum); print(sym_ex); print(conf_tot); print(nref); print(de_counts); print(pc); print(tcga_pairs)
print(rep_rate); print(f_sum); print(sens)
message("Done.")
