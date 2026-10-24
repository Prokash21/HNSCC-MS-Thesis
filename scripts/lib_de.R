## Helpers for step 04: reference-cell CNV signal quantiles, pseudobulk and
## edgeR quasi-likelihood DE.
##
## ref_signal_quantiles(): the step 03 cache keeps CNV scores of epithelial
##   cells and the 99th-percentile threshold per unit, but not the reference
##   cells' own signals. A reference cell's signal depends only on the
##   reference cells (gene means, smoothing, per-cell median and reference
##   centring are all computed from them), so the reference distribution is
##   rebuilt here from the reference cells alone with run_cnv() (lib_cnv.R),
##   without re-scoring any epithelial cell. Where step 03 sub-sampled
##   reference cells (Choi NL/LP, Bill, Peng), the random draws are replayed in
##   the order scripts/03_cnv.R made them. Each rebuilt unit is accepted only
##   if its 99th percentile equals the cached threshold.
## Date: 2026-10-24

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(edgeR)
})

## ---- reference-cell signal quantiles -------------------------------------
## number of random draws keep_heat() (scripts/03_cnv.R) made for one unit:
## sample(i, 80) for each group of > 80 cells, groups in sorted order
replay_heat <- function(group_sizes) {
  for (n in group_sizes[order(names(group_sizes))]) if (n > 80) sample(seq_len(n), 80)
  invisible(NULL)
}

ref_signal_quantiles <- function(paths, units, probs = c(0.95, 0.99), max_ref = 300,
                                 nasopharynx = c("BHN77", "BHN34", "BHN61", "BHN73")) {
  pos <- fread(file.path(paths$cnv_cache, "cnv_genes.csv"))
  setorder(pos, chr_n, start)
  meta2 <- fread(paths$epithelial_v2_meta)
  setkey(meta2, cell)
  dbl <- fread(file.path(paths$cnv_cache, "doublets.csv.gz"))
  out <- list()
  score <- function(u, cnt, lib) {
    r <- run_cnv(cnt, lib, rep(TRUE, ncol(cnt)), rep(FALSE, ncol(cnt)), rep(FALSE, ncol(cnt)), pos)
    q <- quantile(r$cells$signal, probs, names = FALSE)
    out[[u]] <<- data.table(unit = u, n_ref = ncol(cnt), q95 = q[1], q99 = q[2],
                            thr_cached = units[[u]]$thr, n_ref_cached = units[[u]]$n_ref)
  }
  group_sizes_epi <- function(u, prefix_fun) {
    cc <- units[[u]]$cells
    c(Reference = units[[u]]$n_ref, table(prefix_fun(cc)))
  }

  ## step 03 ran sections A (doublets) -> B (gene list) -> C (CNV) in one
  ## session after set.seed(1); only the last scDblFinder call matters
  ## because run_dbl() reset the seed before each sample
  set.seed(1)
  suppressPackageStartupMessages({ library(scDblFinder); library(SingleCellExperiment) })
  bc <- fread(paths$bill_cells)
  last_f <- tail(list.files(file.path(paths$bill_dir, "counts"), full.names = TRUE), 1)
  m <- readRDS(last_f)
  m <- m[, intersect(bc[sample == sub("_gex_raw_counts.rds$", "", basename(last_f)), sample_barcode], colnames(m))]
  set.seed(1)
  invisible(scDblFinder(SingleCellExperiment(list(counts = m)), verbose = FALSE))
  rm(m)
  n_atl <- ncol(readRDS(paths$atlas_epi_v2)$counts)
  invisible(sample(n_atl, min(20000, n_atl)))
  invisible(gc())

  ## C1 atlas (no sub-sampling; draws only for the BHN70 heatmap)
  ref <- readRDS(paths$cnv_input)
  ref_meta <- as.data.table(ref$meta)[group == "Reference"]
  am <- meta2[dataset == "HNSCatlas"]
  for (p in setdiff(unique(am$patient), nasopharynx)) {
    if (grepl("^P[0-9]", p) || !p %in% names(units)) next
    rc <- ref_meta[Patient == p, cell]
    score(p, ref$counts[, rc, drop = FALSE], colSums(ref$counts[, rc, drop = FALSE]))
    if (p == "BHN70") replay_heat(group_sizes_epi(p, function(cc) paste0("Epithelial: ", cc$cnv_call)))
  }
  ## C2 Choi patients
  ch <- readRDS(paths$choi_nl_lp_all)
  chm <- ch$meta
  g_choi <- intersect(rownames(ch$counts), rownames(ref$counts))
  choi_pts <- union(unique(am[grepl("^P[0-9]", patient), patient]), unique(chm$patient))
  for (p in choi_pts) {
    ca <- am[patient == p, cell]
    ca_ref <- ref_meta[Patient == p, cell]
    nl_ref <- chm[patient == p & celltype != "Epithelial.cells",
                  .(cell = if (.N > max_ref) sample(cell, max_ref) else cell), by = sample]$cell
    if (length(ca_ref) + length(nl_ref) < 30) next
    cnt <- cbind(ref$counts[g_choi, ca_ref, drop = FALSE], ch$counts[g_choi, nl_ref, drop = FALSE])
    score(p, cnt, colSums(cnt))
    if (p == "P15") replay_heat(group_sizes_epi(p, function(cc)
      paste0(fifelse(cc$cell %in% ca, "CA", sub("^Choi_([A-Z]+)[0-9]+_.*$", "\\1", cc$cell)), ": ", cc$cnv_call)))
  }
  rm(ref, ch, chm); invisible(gc())
  ## C3 Bill
  bm <- meta2[dataset == "GSE234933_Bill"]
  for (s in unique(sub("^Bill_([^_]+)_.*$", "\\1", bm$cell))) {
    u <- paste0("Bill_", s)
    if (paste0("B", s) %in% nasopharynx) next
    nonepi <- bc[sample == s & subclustering != "Tumor cells", sample_barcode]
    nonepi <- setdiff(nonepi, sub("^Bill_", "", dbl[dbl_class == "doublet", cell]))
    if (length(nonepi) > max_ref) nonepi <- sample(nonepi, max_ref)
    if (length(nonepi) < 30) next
    raw <- readRDS(file.path(paths$bill_dir, "counts", paste0(s, "_gex_raw_counts.rds")))
    score(u, raw[pos$gene[pos$gene %in% rownames(raw)], nonepi, drop = FALSE], colSums(raw[, nonepi, drop = FALSE]))
    rm(raw); invisible(gc())
    if (u == "Bill_HN23") replay_heat(group_sizes_epi(u, function(cc) paste0("Epithelial: ", cc$cnv_call)))
  }
  ## C4 Peng
  pa <- readRDS(paths$peng_all_v2)
  cl <- readRDS(paths$peng_epi_v2)$clusters
  pa$meta[, lineage := cl$lineage[match(cluster, cl$cluster)]]
  for (s in unique(pa$meta$sample)) {
    u <- paste0("Peng_", s)
    rc <- pa$meta[sample == s & lineage %in% c("immune", "fibroblast", "endothelial") &
                    !cell %in% dbl[dbl_class == "doublet", cell], cell]
    if (length(rc) > max_ref) rc <- sample(rc, max_ref)
    score(u, pa$counts[, rc, drop = FALSE], colSums(pa$counts[, rc, drop = FALSE]))
  }
  rm(pa); invisible(gc())
  res <- rbindlist(out)
  res[, match_99 := abs(q99 - thr_cached) <= 1e-9 * thr_cached & n_ref == n_ref_cached]
  res[]
}

## ---- pseudobulk -----------------------------------------------------------
## sum raw counts per key; `cells` is a data.table(cell, key) and a cell may
## appear under several keys (e.g. primary and strict labels)
pseudobulk <- function(counts, cells) {
  j <- match(cells$cell, colnames(counts))
  stopifnot(!anyNA(j))
  key <- factor(cells$key)
  ind <- sparseMatrix(i = j, j = as.integer(key), x = 1, dims = c(ncol(counts), nlevels(key)),
                      dimnames = list(NULL, levels(key)))
  m <- counts %*% ind
  list(counts = as(m, "CsparseMatrix"), n_cells = c(table(key))[levels(key)])
}

## study from the patient prefix (BHN Bill, OP Puram, HN Kurten, P Choi)
study_of <- function(p) fcase(grepl("^BHN", p), "Bill2023", grepl("^OP", p), "Puram2023",
                              grepl("^HN", p), "Kurten2021", grepl("^P[0-9]", p), "Choi2023",
                              default = NA_character_)

## ---- edgeR quasi-likelihood DE (Chen, Lun & Smyth 2016) -------------------
## filterByExpr on the state groups, TMM, robust QL fit, QL F-test of `coef`
de_qlf <- function(counts, smp, formula, coef, fdr = 0.05, lfc = 1) {
  d <- model.matrix(formula, data = smp)
  if (qr(d)$rank < ncol(d)) stop("design is not of full rank: ", deparse(formula))
  y <- DGEList(counts)
  keep <- filterByExpr(y, group = smp$state)
  y <- normLibSizes(y[keep, , keep.lib.sizes = FALSE])
  y <- estimateDisp(y, d)
  fit <- glmQLFit(y, d, robust = TRUE)
  tt <- topTags(glmQLFTest(fit, coef = coef), n = Inf, sort.by = "PValue")$table
  tt <- data.table(gene = rownames(tt), tt)
  tt[, sig := FDR < fdr & abs(logFC) >= lfc]
  tt[, direction := fifelse(!sig, "ns", fifelse(logFC > 0, "up", "down"))]
  list(table = tt, n_samples = ncol(y), n_genes_tested = nrow(y), design = colnames(d),
       df_prior = median(fit$df.prior))
}
