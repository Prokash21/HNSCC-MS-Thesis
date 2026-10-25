## Helpers for steps 04/04b: reference-cell CNV signal quantiles, scoring of
## new cells, pseudobulk and edgeR quasi-likelihood DE.
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
## Date: 2026-10-25

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

## replay of the step 03 reference-cell selection, in the order and with the
## random draws of scripts/03_cnv.R; on_unit(u, cnt, lib) is called with the
## reference counts (CNV genes or more) and library sizes of every unit.
## upto = "choi" stops after the Choi units.
replay_step03_references <- function(paths, units, on_unit, upto = c("all", "choi"), max_ref = 300,
                                     nasopharynx = c("BHN77", "BHN34", "BHN61", "BHN73")) {
  upto <- match.arg(upto)
  pos <- fread(file.path(paths$cnv_cache, "cnv_genes.csv"))
  meta2 <- fread(paths$epithelial_v2_meta)
  setkey(meta2, cell)
  dbl <- fread(file.path(paths$cnv_cache, "doublets.csv.gz"))
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
    on_unit(p, ref$counts[, rc, drop = FALSE], colSums(ref$counts[, rc, drop = FALSE]))
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
    on_unit(p, cnt, colSums(cnt))
    if (p == "P15") replay_heat(group_sizes_epi(p, function(cc)
      paste0(fifelse(cc$cell %in% ca, "CA", sub("^Choi_([A-Z]+)[0-9]+_.*$", "\\1", cc$cell)), ": ", cc$cnv_call)))
  }
  rm(ref, ch, chm); invisible(gc())
  if (upto == "choi") return(invisible(NULL))
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
    on_unit(u, raw[pos$gene[pos$gene %in% rownames(raw)], nonepi, drop = FALSE], colSums(raw[, nonepi, drop = FALSE]))
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
    on_unit(u, pa$counts[, rc, drop = FALSE], colSums(pa$counts[, rc, drop = FALSE]))
  }
  rm(pa); invisible(gc())
  invisible(NULL)
}

ref_signal_quantiles <- function(paths, units, probs = c(0.95, 0.99), ...) {
  pos <- fread(file.path(paths$cnv_cache, "cnv_genes.csv"))
  setorder(pos, chr_n, start)
  out <- list()
  replay_step03_references(paths, units, function(u, cnt, lib) {
    r <- run_cnv(cnt, lib, rep(TRUE, ncol(cnt)), rep(FALSE, ncol(cnt)), rep(FALSE, ncol(cnt)), pos)
    q <- quantile(r$cells$signal, probs, names = FALSE)
    out[[u]] <<- data.table(unit = u, n_ref = ncol(cnt), q95 = q[1], q99 = q[2],
                            thr_cached = units[[u]]$thr, n_ref_cached = units[[u]]$n_ref)
  }, ...)
  res <- rbindlist(out)
  res[, match_99 := abs(q99 - thr_cached) <= 1e-9 * thr_cached & n_ref == n_ref_cached]
  res[]
}

## CNV scores of new cells against a unit's step 03 reference cells and
## cached tumour-average profile (signal: run_cnv; correlation: with the
## cached profile, as step 03 did with its own profile)
score_new_cells <- function(new_counts, ref_counts, unit, pos) {
  cnt <- cbind(new_counts, ref_counts)
  is_ref <- rep(c(FALSE, TRUE), c(ncol(new_counts), ncol(ref_counts)))
  r <- run_cnv(cnt, colSums(cnt), is_ref, !is_ref, rep(FALSE, ncol(cnt)), pos, keep_s = TRUE)
  prof <- unit$profile
  corr <- if (is.null(prof)) rep(NA_real_, sum(!is_ref)) else
    as.vector(cor(r$s[, !is_ref, drop = FALSE], prof[r$genes]))
  data.table(cell = colnames(new_counts), signal = r$cells$signal[!is_ref], corr = corr, thr99_rebuilt = r$thr)
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

## ---- HGNC symbol harmonisation (HGNC; Seal et al. 2023) -------------------
## map symbols to current approved symbols: approved -> kept; else a unique
## previous symbol; else a unique alias; symbols that are previous/alias
## symbols of several approved genes are ambiguous and dropped; symbols HGNC
## does not know (e.g. clone-based lncRNA names) are kept unchanged
read_hgnc <- function(f) {
  hg <- fread(f, sep = "\t", quote = "", select = c("hgnc_id", "symbol", "status", "prev_symbol",
                                                   "alias_symbol", "ensembl_gene_id", "locus_group"))
  hg[, `:=`(prev_symbol = gsub('"', "", prev_symbol), alias_symbol = gsub('"', "", alias_symbol))]
  hg[status == "Approved"]
}
hgnc_map <- function(sym, hg) {
  long <- function(col) {
    x <- hg[get(col) != "", .(old = unlist(strsplit(get(col), "|", fixed = TRUE))), by = symbol]
    x[, .(n_targets = uniqueN(symbol), target = symbol[1]), by = old]
  }
  prev <- long("prev_symbol")
  alias <- long("alias_symbol")
  res <- data.table(original = sym)
  res[, how := fcase(original %in% hg$symbol, "approved",
                     original %in% prev[n_targets == 1, old], "previous",
                     original %in% prev[n_targets > 1, old], "ambiguous_previous",
                     original %in% alias[n_targets == 1, old], "alias",
                     original %in% alias[n_targets > 1, old], "ambiguous_alias",
                     default = "not_in_hgnc")]
  res[, mapped := fcase(how == "approved", original,
                        how == "previous", prev$target[match(original, prev$old)],
                        how == "alias", alias$target[match(original, alias$old)],
                        how == "not_in_hgnc", original,
                        default = NA_character_)]
  res[]
}
## sum rows of a genes x samples matrix by mapped symbol; unmapped rows dropped
map_rows <- function(m, map) {
  mp <- map$mapped[match(rownames(m), map$original)]
  keep <- !is.na(mp)
  f <- factor(mp[keep])
  agg <- sparseMatrix(i = as.integer(f), j = which(keep), x = 1, dims = c(nlevels(f), nrow(m)),
                      dimnames = list(levels(f), NULL))
  out <- agg %*% m
  if (!is(m, "sparseMatrix")) out <- as.matrix(out)
  out
}
