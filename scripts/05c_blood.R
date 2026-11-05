## Step 05c (part A): blood background (Cillo 2020 PBMC scRNA-seq, GTEx v8)
##
## Purpose: label major immune cell types in the GSE139324 PBMC samples
##   (HNSCC HPV-, HNSCC HPV+, healthy donors) with standard markers (no author
##   annotation in GEO), pseudobulk per cell type x group and report CPM for
##   the 24 atlas genes; GTEx v8 median TPM in whole blood, esophagus mucosa
##   and minor salivary gland; the practical blood-low flag.
## Inputs:  data/raw/step05c_blood/{pbmc/*PBMC*, GSE139324_series_matrix.txt.gz, GTEx_v8_gene_median_tpm.gct.gz},
##          HGNC complete set
## Outputs: results/step05c_gene_atlas/tables/A_*.csv, data/processed/v2/step05c_cillo_pbmc_pseudobulk.rds
## Date: 2026-10-28

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
})
source("config/paths.R")
source("scripts/lib_de.R")

tdir <- file.path(paths$step05c, "tables")
dir.create(tdir, recursive = TRUE, showWarnings = FALSE)
hg <- read_hgnc(paths$hgnc)
genes24 <- c("LAMC2", "KRT13", "HOXD11", "HOXB7", "CASC9", "IGF2BP3", "LY6K", "SLC2A1", "SLC39A4", "E2F1", "MYBL2",
             "FOXM1", "HMGB3", "KIF23", "HMMR", "P4HA1", "NFIC", "CA9", "ARPC1A", "PSMC5", "TBP", "GUSB", "SDHA", "ALAS1")
min_genes <- 200
max_mito <- 15
min_cells_type <- 50
markers <- list(
  `CD4 T` = c("CD3D", "CD3E", "IL7R", "CD4"), `CD8 T` = c("CD3D", "CD8A", "CD8B"), NK = c("GNLY", "NKG7", "KLRD1", "KLRF1"),
  B = c("MS4A1", "CD79A", "CD79B"), plasma = c("MZB1", "JCHAIN", "TNFRSF17"), `CD14 monocyte` = c("CD14", "LYZ", "S100A8", "S100A9"),
  `CD16 monocyte` = c("FCGR3A", "MS4A7"), cDC = c("FCER1A", "CD1C", "CLEC10A"), pDC = c("LILRA4", "IL3RA", "CLEC4C"),
  platelet = c("PPBP", "PF4"))
free <- function() invisible(gc(full = TRUE))

## =========================================================================
## A1. GSE139324 PBMC samples and groups (GEO series matrix)
## =========================================================================
l <- readLines(gzfile(paths$gse139324_matrix))
getf <- function(key) gsub('"', "", strsplit(grep(key, l, value = TRUE)[1], "\t")[[1]][-1])
meta <- data.table(gsm = getf("^!Sample_geo_accession"), title = getf("^!Sample_title"), source = getf("^!Sample_source_name_ch1"))
hpv_line <- grep("hpv_status", grep("^!Sample_characteristics_ch1", l, value = TRUE), value = TRUE)
meta[, hpv := sub("hpv_status: ", "", gsub('"', "", strsplit(hpv_line, "\t")[[1]][-1]))]
meta <- meta[source == "peripheral blood"]
meta[, group := fcase(grepl("^HD", title), "healthy donor", hpv == "HPV_negative", "HNSCC HPV-", hpv == "HPV_positive", "HNSCC HPV+",
                      default = "UNCLEAR")]
files <- list.files(paths$gse139324_pbmc, pattern = "matrix.mtx.gz$", full.names = TRUE)
meta[, prefix := sapply(gsm, function(g) sub("_matrix.mtx.gz$", "", basename(grep(paste0("^", g, "_"), basename(files), value = TRUE)[1])))]
stopifnot(!anyNA(meta$prefix), all(meta$group != "UNCLEAR"))
fwrite(meta[, .N, by = group], file.path(tdir, "A_cillo_samples.csv"))
print(meta[, .N, by = group])

## =========================================================================
## A2. Read 10x matrices; Ensembl IDs -> HGNC symbols; QC
## =========================================================================
rd10x <- function(p) {
  f <- function(s) file.path(paths$gse139324_pbmc, paste0(p, s))
  m <- readMM(f("_matrix.mtx.gz"))
  g <- fread(f("_genes.tsv.gz"), header = FALSE)
  b <- fread(f("_barcodes.tsv.gz"), header = FALSE)$V1
  dimnames(m) <- list(g$V1, paste0(p, "_", b))
  list(m = as(m, "CsparseMatrix"), genes = g)
}
first <- rd10x(meta$prefix[1])
gtab <- first$genes[, .(ensembl = V1, symbol = V2)]
gtab[, hgnc := hg$symbol[match(ensembl, hg$ensembl_gene_id)]]
mp <- hgnc_map(gtab[is.na(hgnc), symbol], hg)
gtab[is.na(hgnc), hgnc := mp$mapped[match(symbol, mp$original)]]
fwrite(data.table(gene = genes24, in_reference = genes24 %in% gtab$hgnc,
                  ensembl = gtab$ensembl[match(genes24, gtab$hgnc)]), file.path(tdir, "A_cillo_gene_coverage.csv"))
mats <- list(); qc <- list()
for (i in seq_len(nrow(meta))) {
  x <- if (i == 1) first else rd10x(meta$prefix[i])
  stopifnot(identical(x$genes$V1, gtab$ensembl))
  m <- x$m
  ng <- colSums(m > 0); mito <- 100 * colSums(m[grepl("^MT-", gtab$symbol), , drop = FALSE]) / pmax(colSums(m), 1)
  keep <- ng >= min_genes & mito < max_mito
  qc[[i]] <- data.table(gsm = meta$gsm[i], title = meta$title[i], group = meta$group[i], cells_raw = ncol(m), cells_kept = sum(keep))
  mats[[i]] <- m[, keep]
}
rm(first); free()
counts <- do.call(cbind, mats); rm(mats); free()
cells <- data.table(cell = colnames(counts), title = rep(meta$title, sapply(qc, `[[`, "cells_kept")),
                    group = rep(meta$group, sapply(qc, `[[`, "cells_kept")))
fwrite(rbindlist(qc), file.path(tdir, "A_cillo_qc.csv"))
message(sprintf("A2. %d PBMC cells kept from %d samples", ncol(counts), nrow(meta)))

## =========================================================================
## A3. Marker-based labelling (Seurat, clusters labelled by marker sets)
## =========================================================================
sym_counts <- counts
rownames(sym_counts) <- make.unique(gtab$symbol)
so <- CreateSeuratObject(sym_counts, min.cells = 0, min.features = 0)
rm(sym_counts); free()
so <- NormalizeData(so, verbose = FALSE)
so <- FindVariableFeatures(so, nfeatures = 2000, verbose = FALSE)
so <- ScaleData(so, verbose = FALSE)
so <- RunPCA(so, npcs = 30, verbose = FALSE)
so <- FindNeighbors(so, dims = 1:30, verbose = FALSE)
so <- FindClusters(so, resolution = 0.8, verbose = FALSE)
cl <- as.character(so$seurat_clusters)
ln <- GetAssayData(so, layer = "data")
mk <- unique(unlist(markers))
mk <- mk[mk %in% rownames(ln)]
cm <- sapply(split(seq_along(cl), cl), function(j) Matrix::rowMeans(ln[mk, j, drop = FALSE]))
cz <- t(scale(t(cm)))
score <- sapply(markers, function(s) colMeans(cz[intersect(s, rownames(cz)), , drop = FALSE]))
lab <- colnames(score)[apply(score, 1, which.max)]
names(lab) <- rownames(score)
cl_tab <- data.table(cluster = rownames(score), cells = as.integer(table(cl)[rownames(score)]), label = lab,
                     top_score = apply(score, 1, max), second = apply(score, 1, function(v) colnames(score)[order(-v)[2]]),
                     round(score, 2))
fwrite(cl_tab[order(as.integer(cluster))], file.path(tdir, "A_cillo_cluster_labels.csv"))
cells[, cell_type := lab[cl]]
rm(so, ln); free()
print(cl_tab[order(as.integer(cluster)), .(cluster, cells, label, top_score = round(top_score, 2), second)])

## =========================================================================
## A4. Pseudobulk per cell type x group -> CPM for the 24 genes
## =========================================================================
key <- factor(paste(cells$group, cells$cell_type, sep = "|"))
ind <- sparseMatrix(i = seq_len(nrow(cells)), j = as.integer(key), x = 1, dims = c(nrow(cells), nlevels(key)),
                    dimnames = list(NULL, levels(key)))
pbm <- counts %*% ind
lib <- colSums(pbm)
rows <- match(genes24, gtab$hgnc)
cpm <- as.matrix(pbm[rows[!is.na(rows)], , drop = FALSE]) %*% diag(1e6 / lib)
dimnames(cpm) <- list(genes24[!is.na(rows)], levels(key))
ncell <- as.integer(table(key)[levels(key)])
pbl <- melt(as.data.table(cpm, keep.rownames = "gene"), id.vars = "gene", variable.name = "key", value.name = "CPM")
pbl[, `:=`(group = sub("\\|.*", "", key), cell_type = sub(".*\\|", "", key), n_cells = ncell[match(key, levels(key))])]
pbl[, `:=`(used_for_max = n_cells >= min_cells_type, library = lib[match(key, levels(key))])]
pbl[, UMI := CPM * library / 1e6]
fwrite(pbl[, .(gene, group, cell_type, n_cells, library, UMI, CPM, used_for_max)], file.path(tdir, "A_cillo_pseudobulk_cpm.csv"))
fwrite(cells[, .N, by = .(group, cell_type)][order(group, -N)], file.path(tdir, "A_cillo_cells_per_type.csv"))
saveRDS(list(counts = pbm, n_cells = ncell, genes = gtab), paths$step05c_pbmc)
rm(counts); free()

## =========================================================================
## A5. GTEx v8 median TPM; blood-low flag
## =========================================================================
gx <- fread(cmd = paste("gzip -dc", shQuote(paths$gtex_median)), skip = 2)
gx[, ensembl := sub("\\..*$", "", Name)]
gx <- gx[!grepl("_PAR_Y$", Name)]
gx[, hgnc := hg$symbol[match(ensembl, hg$ensembl_gene_id)]]
gt <- gx[hgnc %in% genes24, .(gene = hgnc, GTEx_whole_blood = `Whole Blood`, GTEx_esophagus_mucosa = `Esophagus - Mucosa`,
                              GTEx_minor_salivary_gland = `Minor Salivary Gland`)]
gt <- gt[, lapply(.SD, sum), by = gene]
fwrite(gt, file.path(tdir, "A_gtex_median_tpm.csv"))
mx <- pbl[used_for_max == TRUE, .(PBMC_max_CPM = max(CPM), PBMC_max_type = paste(group, cell_type, sep = ": ")[which.max(CPM)],
                                  PBMC_HPVneg_max = max(CPM[group == "HNSCC HPV-"]), PBMC_HPVpos_max = max(CPM[group == "HNSCC HPV+"]),
                                  PBMC_HD_max = max(CPM[group == "healthy donor"]), PBMC_max_UMI = UMI[which.max(CPM)],
                                  ## post hoc, descriptive: maximum over pseudobulks with >= 1 M counts (not the pre-specified flag)
                                  PBMC_max_CPM_lib1M = max(CPM[library >= 1e6])), by = gene]
bl <- merge(data.table(gene = genes24), merge(mx, gt, by = "gene", all = TRUE), by = "gene", all.x = TRUE)
bl[, blood_low := PBMC_max_CPM < 1 & GTEx_whole_blood < 1]
bl[, blood_low_lib1M_posthoc := PBMC_max_CPM_lib1M < 1 & GTEx_whole_blood < 1]
fwrite(bl, file.path(tdir, "A_blood_background.csv"))
print(bl[order(PBMC_max_CPM), .(gene, PBMC_max = round(PBMC_max_CPM, 2), PBMC_max_type, UMI = round(PBMC_max_UMI),
                                max_lib1M = round(PBMC_max_CPM_lib1M, 2), GTEx_blood = round(GTEx_whole_blood, 2), blood_low,
                                blood_low_lib1M_posthoc)])
