## Step 02: data acquisition and assembly of the v2 epithelial object
##
## Purpose: gather raw-count epithelial cells from the HNSCatlas (discovery),
##   Choi GSE181919 normal mucosa / leukoplakia (reference), Bill GSE234933
##   recurrences / metastases (recurrence) and Peng GSE172577 (single-cell
##   validation); harmonise metadata and flags; no integration, no analysis.
## Inputs:  data/processed/HNSCatlas/atlas_epithelial.rds
##          data/raw/GSE181919/GSE181919_UMI_counts.txt.gz
##          data/raw/GSE234933/GSE234933_MGH_HNSCC_gex_raw_counts.tar.gz
##          data/raw/GSE172577/GSM*_SYSMH*_{matrix.mtx,features.tsv,barcodes.tsv}.gz
##          data/raw/step01_metadata/ (GEO cell and sample annotations)
##          results/step01_metadata/patient_clinical.csv
## Outputs: data/processed/v2/*.rds (one object per dataset + epithelial_v2.rds)
##          data/processed/v2/epithelial_v2_meta.csv.gz
##          results/step02_assembly/*.csv
## Date: 2026-10-12
## Run from the repo root:  Rscript scripts/02_assembly.R
## Each dataset is processed on its own and cached, so memory stays < 16 GB.

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
})
source("config/paths.R")
dir.create(paths$processed_v2, recursive = TRUE, showWarnings = FALSE)
out <- paths$step02
dir.create(out, recursive = TRUE, showWarnings = FALSE)
set.seed(1)
free <- function() invisible(gc(full = TRUE))
is_int <- function(m) if (inherits(m, "dgCMatrix")) all(m@x == round(m@x)) else all(m == round(m))

## canonical HNSCC markers (Puram et al. 2017)
mk <- list(epithelial  = c("KRT5", "KRT14", "KRT17", "EPCAM", "SFN"),
           immune      = "PTPRC",
           fibroblast  = c("COL1A2", "DCN"),
           endothelial = c("PECAM1", "VWF"))
## per-cell marker detection (counts > 0), the same for every dataset
marker_cols <- function(counts) {
  det <- function(g) {
    g <- intersect(g, rownames(counts))
    if (!length(g)) return(rep(NA_integer_, ncol(counts)))
    as.integer(colSums(counts[g, , drop = FALSE] > 0))
  }
  data.table(n_counts = as.integer(colSums(counts)), n_genes = as.integer(colSums(counts > 0)),
             epi_markers_n = det(mk$epithelial), ptprc_detected = det(mk$immune) > 0,
             fibro_detected = det(mk$fibroblast) > 0, endo_detected = det(mk$endothelial) > 0)
}

clin <- fread(paths$patient_clinical)
clin_cols <- clin[, .(patient = patient_id, site = site_group, subsite, stage_T = T, stage_N = N,
                      stage_M = M, hpv_test = hpv_method)]
treated   <- c("OP8", "OP13")
nasophar  <- "BHN77"

## =========================================================================
## A. Downloads: local size vs GEO size, raw-integer checks (filled below)
## =========================================================================
geo <- function(gse, f) sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%snnn/%s/suppl/%s",
                                substr(gse, 1, nchar(gse) - 3), gse, f)
remote_size <- function(u) {
  h <- tryCatch(curlGetHeaders(u), error = function(e) character(0))
  v <- as.numeric(sub(".*:\\s*", "", grep("^content-length", h, ignore.case = TRUE, value = TRUE)))
  if (length(v)) max(v) else NA_real_
}
dl <- data.table(
  file = c(paths$choi_counts, paths$bill_counts_tar, file.path(paths$peng_dir, "GSE172577_RAW.tar"),
           paths$bill_cells, paths$bill_samples),
  url = c(geo("GSE181919", "GSE181919_UMI_counts.txt.gz"),
          geo("GSE234933", "GSE234933_MGH_HNSCC_gex_raw_counts.tar.gz"),
          geo("GSE172577", "GSE172577_RAW.tar"),
          geo("GSE234933", "GSE234933_MGH_HNSCC_cell_annotation.txt.gz"),
          geo("GSE234933", "GSE234933_MGH_HNSCC_sample_annotation.txt.gz")))
for (i in seq_len(nrow(dl))) if (!file.exists(dl$file[i])) {
  dir.create(dirname(dl$file[i]), recursive = TRUE, showWarnings = FALSE)
  download.file(dl$url[i], dl$file[i], mode = "wb")
}
dl[, local_bytes := file.size(file)]
dl[, geo_bytes := vapply(url, remote_size, 0)]
dl[, size_ok := local_bytes == geo_bytes]
int_checks <- list()

## =========================================================================
## B. HNSCatlas epithelial cells (RNA counts only)
## =========================================================================
if (!file.exists(paths$atlas_epi_v2)) {
  message("B. Atlas: reading the epithelial object...")
  epi <- readRDS(paths$atlas_epithelial)
  counts <- LayerData(epi, assay = "RNA", layer = "counts")
  md <- epi@meta.data
  rm(epi); free()

  md <- data.table(cell = rownames(md), md)
  md[, sample := sub("_[ACGT]+[-.][0-9]+$", "", cell)]
  ## Choi cells carry no sample in their ID; all are the CA sample of the patient
  md[sample == "Choi", sample := paste0("Choi_", Patient)]
  md[, sample_type := fifelse(grepl("norm$", sample), "adjacent_normal",
                       fifelse(grepl("LN$", sample), "lymph_node", "primary_tumour"))]
  md[, compartment := fifelse(Cell_Type == "Tumor", "tumour_label", "nonmalignant_label")]
  md <- cbind(md, marker_cols(counts))
  saveRDS(list(counts = counts, meta = md), paths$atlas_epi_v2)
  rm(counts, md); free()
}
atlas <- readRDS(paths$atlas_epi_v2)
int_checks$atlas <- is_int(atlas$counts)
am <- atlas$meta
## original study label (46 values in the atlas): epithelial / unresolved /
## clearly non-epithelial; "suspect_label" = Tumor/Epithelial call on a cell
## the source study labelled as a non-epithelial type
epi_labels <- c("Epithelial", "Epithelial Cells", "Epithelial.cells", "Tumor", "Malignant.cells", "NormalEpith")
am[, original_label_class := fcase(Original_Cell_Type %in% epi_labels, "epithelial",
                                   Original_Cell_Type %in% c("Unresolved", "Unknown"), "unresolved",
                                   default = "non_epithelial")]
am[, suspect_label := original_label_class == "non_epithelial"]
am[, `:=`(dataset = "HNSCatlas", patient = Patient, hpv = HPV, original_label = Original_Cell_Type,
          role = "discovery")]

puram_split <- am[Patient %in% c("OP33", "OP34", "OP35"), .N, by = .(Patient, sample, sample_type, compartment)]
setorder(puram_split, Patient, sample, compartment)
fwrite(puram_split, file.path(out, "puram_OP33_35_split.csv"))
message("Puram split:"); print(puram_split)

suspect <- am[, .(tumour_label = sum(compartment == "tumour_label"),
                  tumour_suspect = sum(compartment == "tumour_label" & suspect_label),
                  tumour_unresolved = sum(compartment == "tumour_label" & original_label_class == "unresolved"),
                  nonmalignant_label = sum(compartment == "nonmalignant_label"),
                  nonmalignant_suspect = sum(compartment == "nonmalignant_label" & suspect_label)),
              by = .(patient, hpv)]
susp_labels <- am[suspect_label == TRUE, .N, by = .(compartment, original_label)][order(compartment, -N)]
## do suspect "Tumor" cells look epithelial by marker genes?
susp_markers <- am[compartment == "tumour_label", .(cells = .N,
                     pct_epi_markers_2plus = round(100 * mean(epi_markers_n >= 2), 1),
                     pct_PTPRC = round(100 * mean(ptprc_detected), 1),
                     pct_fibro = round(100 * mean(fibro_detected), 1),
                     pct_endo = round(100 * mean(endo_detected), 1)),
                   by = .(suspect_label)]
fwrite(suspect[order(-tumour_suspect)], file.path(out, "suspect_label_per_patient.csv"))
fwrite(susp_labels, file.path(out, "suspect_label_original_types.csv"))
fwrite(susp_markers, file.path(out, "suspect_label_marker_check.csv"))
print(susp_markers)

atlas_genes <- rownames(atlas$counts)
## CDKN2A reference: % of tumour-labelled cells expressing it, per atlas patient
cdkn2a_ref <- am[compartment == "tumour_label" & !suspect_label & !patient %in% c(nasophar, treated),
                 .(cells = .N, pct_CDKN2A = 100 * mean(atlas$counts["CDKN2A", cell] > 0)),
                 by = .(patient, hpv)][cells >= 50]
atlas_meta <- am
rm(atlas); free()

## =========================================================================
## C. Choi GSE181919 normal mucosa (NL) and leukoplakia (LP)
## =========================================================================
choi_bc <- fread(paths$choi_barcodes, skip = 1, header = FALSE,
                 col.names = c("barcode", "patient", "sample", "sex", "age", "tissue", "subsite", "hpv", "celltype"))
if (!file.exists(paths$choi_epi_v2)) {
  message("C. Choi: reading NL/LP columns from the dense UMI matrix...")
  hdr <- strsplit(readLines(gzfile(paths$choi_counts), n = 1), "\t")[[1]]
  sel <- choi_bc[tissue %in% c("NL", "LP")]
  idx <- match(gsub("-", ".", sel$barcode), hdr)
  stopifnot(!anyNA(idx))
  dt <- fread(cmd = sprintf("gzcat '%s'", paths$choi_counts), skip = 1, header = FALSE,
              select = c(1L, idx + 1L))
  genes <- dt[[1]]
  dt[, 1 := NULL]
  m <- as(as.matrix(dt), "CsparseMatrix")
  rm(dt); free()
  dimnames(m) <- list(genes, paste0("Choi_", sel$sample, "_", sel$barcode))
  sel[, cell := paste0("Choi_", sample, "_", barcode)]
  sel <- cbind(sel, marker_cols(m))
  saveRDS(list(counts = m, meta = sel), paths$choi_nl_lp_all)
  ep <- sel$celltype == "Epithelial.cells"
  saveRDS(list(counts = m[, ep], meta = sel[ep]), paths$choi_epi_v2)
  rm(m, sel); free()
}
choi <- readRDS(paths$choi_epi_v2)
int_checks$choi <- is_int(choi$counts)
choi_all_meta <- readRDS(paths$choi_nl_lp_all)$meta
free()
## barcodes already in the atlas? (atlas Choi IDs: Choi_<barcode with '.'>)
atlas_choi <- gsub("[.]", "-", sub("^Choi_", "", atlas_meta[dataset == "HNSCatlas" & grepl("^Choi_", cell), cell]))
choi_overlap <- sum(choi_all_meta$barcode %in% atlas_choi)
cm <- choi$meta
cm[, `:=`(dataset = "GSE181919_Choi", hpv = hpv, original_label = celltype,
          sample_type = fifelse(tissue == "NL", "normal_mucosa", "leukoplakia"),
          compartment = "nonmalignant_label", role = "discovery",
          site_choi = c(OC = "Oral cavity", OP = "Oropharynx", HP = "Hypopharynx")[subsite],
          suspect_label = FALSE)]
choi_tab <- merge(choi_all_meta[, .(total_cells = .N), by = sample],
                  cm[, .(epithelial_cells = .N), by = .(sample, patient, tissue, site = site_choi, hpv)],
                  by = "sample", all = TRUE)
choi_tab[, patient_in_atlas := patient %in% atlas_meta$patient]
fwrite(choi_tab, file.path(out, "choi_nl_lp_epithelial.csv"))
print(choi_tab)
choi_genes <- rownames(choi$counts)
choi_meta <- cm
rm(choi, choi_all_meta); free()

## =========================================================================
## D. Bill GSE234933 recurrences and metastases (samples not in the atlas)
## =========================================================================
bill_s <- fread(paths$bill_samples)
setnames(bill_s, c("sample", "sex", "age", "hpv", "smoking", "site", "specimen"))
bill_s[, in_atlas := paste0("B", sample) %in% atlas_meta$patient]
bill_cells <- fread(paths$bill_cells)
extra <- bill_s[in_atlas == FALSE]
if (!file.exists(paths$bill_epi_v2)) {
  message("D. Bill: extracting ", nrow(extra), " recurrence/metastasis samples...")
  cdir <- file.path(paths$bill_dir, "counts")
  dir.create(cdir, showWarnings = FALSE)
  want <- paste0(extra$sample, "_gex_raw_counts.rds")
  untar(paths$bill_counts_tar, files = want, exdir = cdir)
  mats <- list(); metas <- list(); stats <- list()
  for (s in extra$sample) {
    m <- readRDS(file.path(cdir, paste0(s, "_gex_raw_counts.rds")))
    a <- bill_cells[sample == s]
    keep <- a[subclustering == "Tumor cells", sample_barcode]
    stats[[s]] <- data.table(sample = s, barcodes_in_matrix = ncol(m), annotated_cells = nrow(a),
                             epithelial_cells = length(keep), integer_counts = is_int(m),
                             n_genes = nrow(m))
    mm <- m[, intersect(keep, colnames(m)), drop = FALSE]
    mats[[s]] <- mm
    metas[[s]] <- merge(data.table(cell = colnames(mm)), a, by.x = "cell", by.y = "sample_barcode")[match(colnames(mm), cell)]
    rm(m); free()
  }
  g <- Reduce(intersect, lapply(mats, rownames))
  counts <- do.call(cbind, lapply(mats, function(m) m[g, , drop = FALSE]))
  meta <- rbindlist(metas)
  meta <- cbind(meta, marker_cols(counts))
  colnames(counts) <- meta$cell <- paste0("Bill_", meta$cell)
  saveRDS(list(counts = counts, meta = meta, stats = rbindlist(stats)), paths$bill_epi_v2)
  rm(mats, metas, counts, meta); free()
}
bill <- readRDS(paths$bill_epi_v2)
int_checks$bill <- all(bill$stats$integer_counts)
bm <- bill$meta
bm <- merge(bm, extra[, .(sample, sex, age, hpv_b = hpv, smoking, site_b = site, specimen)], by = "sample")
bm[, `:=`(dataset = "GSE234933_Bill", patient = paste0("B", sample),
          hpv = fifelse(hpv_b == "Positive", "HPV+", "HPV-"),
          sample_type = fifelse(grepl("recurrence", specimen), "recurrence", "metastasis"),
          ## Bill epithelial subcluster: major_state Tumor / NormalEpith; cells
          ## without a major_state are normal epithelia of the metastatic site
          ## (hepatocyte, pneumocyte, exocrine, FAM3D), tumour programmes (TC*)
          ## without a malignancy call, or "unused"
          compartment = fcase(major_state == "Tumor", "tumour_label",
                              major_state == "NormalEpith", "nonmalignant_label",
                              grepl("^Hepat|^Pneumocyte|^Exocrine|^Epithelial_", minor_state), "nonmalignant_label",
                              default = "unannotated"),
          original_label = paste(subclustering, fcoalesce(major_state, "NA"), minor_state, sep = "/"),
          role = "recurrence",
          suspect_label = FALSE)]
bill_tab <- merge(extra[, .(sample, patient = paste0("B", sample), sex, age, hpv, smoking, site, specimen)],
                  bill$stats, by = "sample")
bill_tab <- merge(bill_tab, bm[, .(tumour_label = sum(compartment == "tumour_label"),
                                   nonmalignant_label = sum(compartment == "nonmalignant_label"),
                                   unannotated = sum(compartment == "unannotated")), by = sample],
                  by = "sample", all.x = TRUE)
bill_tab[, treatment := "UNCLEAR (not in GEO; paper Table S1 not openly downloadable)"]
## same patient as an atlas primary? GEO gives one sample per ID; demographic matches only
bill_s[, demo := paste(sex, age, hpv, smoking, site)]
pairs <- merge(bill_s[in_atlas == FALSE, .(rec_sample = sample, specimen, demo)],
               bill_s[in_atlas == TRUE, .(atlas_primary = paste0("B", sample), demo)], by = "demo")
bill_tab[, candidate_atlas_primary := pairs$atlas_primary[match(sample, pairs$rec_sample)]]
setorder(bill_tab, specimen, sample)
fwrite(bill_tab, file.path(out, "bill_recurrence_metastasis_samples.csv"))
fwrite(pairs, file.path(out, "bill_candidate_primary_recurrence_pairs.csv"))
print(bill_tab[, .(sample, hpv, site, specimen, annotated_cells, epithelial_cells, tumour_label, candidate_atlas_primary)])
bill_genes <- rownames(bill$counts)
bill_meta <- bm
rm(bill); free()

## =========================================================================
## E. Peng GSE172577: QC, cluster-level marker annotation, CDKN2A
## =========================================================================
if (!file.exists(paths$peng_epi_v2)) {
  message("E. Peng: reading 10x matrices...")
  tarf <- file.path(paths$peng_dir, "GSE172577_RAW.tar")
  if (!length(list.files(paths$peng_dir, "matrix.mtx.gz$"))) untar(tarf, exdir = paths$peng_dir)
  mfiles <- list.files(paths$peng_dir, "matrix.mtx.gz$", full.names = TRUE)
  mats <- list(); qc <- list()
  for (f in mfiles) {
    pre <- sub("matrix.mtx.gz$", "", f)
    s <- sub("^GSM[0-9]+_", "", basename(pre)); s <- sub("_$", "", s)
    m <- as(readMM(f), "CsparseMatrix")
    feat <- fread(paste0(pre, "features.tsv.gz"), header = FALSE)
    bc <- fread(paste0(pre, "barcodes.tsv.gz"), header = FALSE)[[1]]
    dimnames(m) <- list(make.unique(feat$V2), paste0("Peng_", s, "_", bc))
    nc <- colSums(m); ng <- colSums(m > 0)
    mt <- colSums(m[grep("^MT-", rownames(m)), , drop = FALSE]) / pmax(nc, 1) * 100
    ## QC (Luecken & Theis 2019): fixed floors + per-sample upper MAD outliers
    hi_nc <- 10^(median(log10(nc)) + 5 * mad(log10(nc)))
    hi_ng <- 10^(median(log10(ng)) + 5 * mad(log10(ng)))
    fail <- data.table(low_genes = ng < 200, low_counts = nc < 500, high_mt = mt > 20,
                       high_outlier = nc > hi_nc | ng > hi_ng)
    keep <- !Reduce(`|`, fail)
    qc[[s]] <- data.table(sample = s, barcodes = ncol(m), integer_counts = is_int(m),
                          removed_low_genes = sum(fail$low_genes), removed_low_counts = sum(fail$low_counts),
                          removed_high_mt = sum(fail$high_mt), removed_high_outlier = sum(fail$high_outlier),
                          upper_counts = round(hi_nc), upper_genes = round(hi_ng),
                          kept = sum(keep))
    mats[[s]] <- m[, keep]
    rm(m); free()
  }
  counts <- do.call(cbind, mats); rm(mats); free()
  qc <- rbindlist(qc)
  seu <- CreateSeuratObject(counts, min.cells = 0, min.features = 0)
  seu$sample <- sub("^Peng_([^_]+)_.*$", "\\1", colnames(seu))
  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(seu, nfeatures = 2000, verbose = FALSE)
  seu <- ScaleData(seu, verbose = FALSE)
  seu <- RunPCA(seu, npcs = 30, verbose = FALSE)
  seu <- FindNeighbors(seu, dims = 1:20, verbose = FALSE)
  seu <- FindClusters(seu, resolution = 0.5, random.seed = 1, verbose = FALSE)
  meta <- data.table(cell = colnames(seu), sample = seu$sample, cluster = as.character(seu$seurat_clusters))
  meta <- cbind(meta, marker_cols(counts))
  rm(seu); free()
  saveRDS(list(counts = counts, meta = meta, qc = qc), paths$peng_all_v2)
  rm(counts, meta); free()
}
## cluster label from mean log1p(CP10k) of the marker genes. Detection rates
## are not used: ambient keratin RNA is detected in most immune cells here.
## epithelial: mean(KRT5, KRT14, KRT17, SFN) >= 2.5 and PTPRC, COL1A2/DCN,
## PECAM1/VWF all < 0.5; clusters with high keratins AND PTPRC >= 0.5 are
## "mixed" (likely doublets) and not kept
if (!file.exists(paths$peng_epi_v2)) {
  pa <- readRDS(paths$peng_all_v2)
  x <- pa$counts[c(mk$epithelial, mk$immune, mk$fibroblast, mk$endothelial, "LYZ"), ]
  x <- t(t(x) / colSums(pa$counts)) * 1e4
  x@x <- log1p(x@x)
  avg <- t(sapply(split(seq_len(ncol(x)), pa$meta$cluster), function(i) Matrix::rowMeans(x[, i, drop = FALSE])))
  cl <- data.table(cluster = rownames(avg), cells = as.integer(table(pa$meta$cluster)[rownames(avg)]),
                   keratin = rowMeans(avg[, c("KRT5", "KRT14", "KRT17", "SFN")]), EPCAM = avg[, "EPCAM"],
                   PTPRC = avg[, "PTPRC"], LYZ = avg[, "LYZ"], fibro = pmax(avg[, "COL1A2"], avg[, "DCN"]),
                   endo = pmax(avg[, "PECAM1"], avg[, "VWF"]))
  cl[, lineage := fcase(keratin >= 2.5 & PTPRC < 0.5 & fibro < 0.5 & endo < 0.5, "epithelial",
                        keratin >= 2.5 & PTPRC >= 0.5, "mixed_epithelial_immune",
                        fibro >= 1, "fibroblast", endo >= 1, "endothelial",
                        PTPRC >= 0.5 | LYZ >= 1, "immune", default = "unassigned")]
  pa$meta[, lineage := cl$lineage[match(cluster, cl$cluster)]]
  ep <- pa$meta$lineage == "epithelial"
  saveRDS(list(counts = pa$counts[, ep], meta = pa$meta[ep], qc = pa$qc, clusters = cl), paths$peng_epi_v2)
  rm(pa, x); free()
}
peng <- readRDS(paths$peng_epi_v2)
int_checks$peng <- all(peng$qc$integer_counts)
fwrite(peng$qc, file.path(out, "peng_qc.csv"))
setorder(peng$clusters, -cells)
fwrite(peng$clusters, file.path(out, "peng_cluster_annotation.csv"))
print(peng$qc); print(peng$clusters)
pm <- peng$meta
pm[, CDKN2A_pos := as.numeric(peng$counts["CDKN2A", cell] > 0)]
peng_cdkn2a <- pm[, .(epithelial_cells = .N, pct_CDKN2A = round(100 * mean(CDKN2A_pos), 1)), by = sample]
ref <- cdkn2a_ref[, .(patients = .N, median_pct_CDKN2A = round(median(pct_CDKN2A), 1),
                      min = round(min(pct_CDKN2A), 1), max = round(max(pct_CDKN2A), 1)), by = hpv]
fwrite(peng_cdkn2a, file.path(out, "peng_cdkn2a.csv"))
fwrite(ref, file.path(out, "atlas_cdkn2a_reference.csv"))
print(peng_cdkn2a); print(ref)
pm[, `:=`(dataset = "GSE172577_Peng", patient = paste0("Peng_", sample), hpv = "UNCLEAR",
          sample_type = "primary_tumour", compartment = "unannotated", role = "sc_validation",
          original_label = paste0("cluster ", cluster, " (", lineage, ")"), suspect_label = FALSE)]
peng_genes <- rownames(peng$counts)
peng_meta <- pm
rm(peng); free()

## =========================================================================
## F. One v2 epithelial object: raw counts, shared genes, harmonised meta
## =========================================================================
genes <- Reduce(intersect, list(atlas_genes, choi_genes, bill_genes, peng_genes))
message(sprintf("F. Genes: atlas %d, Choi %d, Bill %d, Peng %d -> shared %d",
                length(atlas_genes), length(choi_genes), length(bill_genes), length(peng_genes), length(genes)))
gene_tab <- data.table(dataset = c("HNSCatlas", "GSE181919_Choi", "GSE234933_Bill", "GSE172577_Peng", "shared"),
                       genes = c(length(atlas_genes), length(choi_genes), length(bill_genes),
                                 length(peng_genes), length(genes)))
fwrite(gene_tab, file.path(out, "genes_per_dataset.csv"))

cols <- c("cell", "dataset", "patient", "sample", "sample_type", "compartment", "hpv", "role",
          "original_label", "suspect_label", "n_counts", "n_genes", "epi_markers_n",
          "ptprc_detected", "fibro_detected", "endo_detected")
bill_meta[, site := site_b]
meta <- rbindlist(list(atlas_meta[, ..cols], choi_meta[, c(cols, "site_choi"), with = FALSE],
                       bill_meta[, c(cols, "site"), with = FALSE], peng_meta[, ..cols]), fill = TRUE)
meta <- merge(meta, clin_cols, by = "patient", all.x = TRUE, suffixes = c("", "_clin"))
meta[dataset == "HNSCatlas", site := site_clin]
meta[dataset == "GSE181919_Choi", site := site_choi]
meta[dataset == "GSE234933_Bill", site := fifelse(is.na(site_clin), site, site_clin)]
meta[dataset == "GSE172577_Peng", site := "Oral cavity"]
meta[, c("site_choi", "site_clin") := NULL]
meta[, `:=`(flag_excluded = patient == nasophar | (dataset == "GSE234933_Bill" & site == "Nasopharynx"),
            flag_treated = patient %in% treated,
            flag_no_site_comparison = patient == "BHN17" | site %in% c("Unknown", NA))]
meta[, exclusion_reason := fcase(patient == nasophar, "nasopharynx primary (decision)",
                                 flag_excluded, "nasopharynx primary (same rule as BHN77)",
                                 flag_treated, "previously treated: primary-analysis exclusion, sensitivity only",
                                 default = "")]
meta[, primary_analysis := !flag_excluded & !flag_treated]
setkey(meta, cell)
fwrite(meta, paths$epithelial_v2_meta)

if (!file.exists(paths$epithelial_v2)) {
  message("F. Building the combined object...")
  parts <- list(
    readRDS(paths$atlas_epi_v2)$counts[genes, ],
    readRDS(paths$choi_epi_v2)$counts[genes, ],
    readRDS(paths$bill_epi_v2)$counts[genes, ],
    readRDS(paths$peng_epi_v2)$counts[genes, ])
  counts <- do.call(cbind, parts); rm(parts); free()
  storage.mode(counts@x) <- "double"
  md <- as.data.frame(meta[colnames(counts)])
  rownames(md) <- md$cell
  obj <- CreateSeuratObject(counts, meta.data = md, project = "epithelial_v2", min.cells = 0, min.features = 0)
  rm(counts); free()
  saveRDS(obj, paths$epithelial_v2)
  message("Saved ", paths$epithelial_v2, ": ", nrow(obj), " genes x ", ncol(obj), " cells")
  rm(obj); free()
}

## ---- summary tables ---------------------------------------------------------
dl[, raw_integer_counts := c(int_checks$choi, int_checks$bill, int_checks$peng, NA, NA)]
dl[, file := basename(file)]
fwrite(dl[, .(file, local_bytes, geo_bytes, size_ok, raw_integer_counts)], file.path(out, "downloads.csv"))
fwrite(data.table(dataset = names(int_checks), integer_counts = unlist(int_checks)), file.path(out, "integer_checks.csv"))

cells <- meta[, .N, by = .(dataset, sample_type, compartment, hpv, site)][order(dataset, sample_type, compartment, hpv, site)]
fwrite(cells, file.path(out, "cells_by_dataset_type_hpv_site.csv"))

hpvneg <- meta[hpv == "HPV-" & primary_analysis & compartment == "tumour_label",
               .(tumour_cells = .N, tumour_cells_not_suspect = sum(!suspect_label)),
               by = .(role, dataset, patient, sample_type, site)]
fwrite(hpvneg[order(role, -tumour_cells)], file.path(out, "hpvneg_tumour_cells_per_patient.csv"))
print(hpvneg[tumour_cells >= 50, .(patients = uniqueN(patient)), by = .(role, dataset)])
print(hpvneg[tumour_cells_not_suspect >= 50, .(patients_not_suspect = uniqueN(patient)), by = .(role, dataset)])
## Peng: malignancy not called yet; epithelial cells per patient
fwrite(meta[dataset == "GSE172577_Peng", .(epithelial_cells = .N), by = patient],
       file.path(out, "peng_epithelial_per_patient.csv"))
message("Done. Outputs in ", out, " and ", paths$processed_v2)
