## Helpers for step 03: gene positions and per-sample CNV inference.
##
## Method (Tirosh et al. 2016; Puram et al. 2017), as in v1
## (archive_v1_HPVpos/scripts/05b_cnv_check.R):
##   log2(CPM/10 + 1) -> minus the mean of the sample's reference
##   (immune/stromal) cells -> clipped to +-3 -> 101-gene moving average along
##   each chromosome -> centred per cell (median) and on the reference.
##   CNV signal      = mean square of the smoothed profile
##   CNV correlation = Pearson correlation with the sample's tumour-average
##                     profile = mean of the top 5% of epithelial cells by CNV
##                     signal (at least 10 cells; Puram et al. 2017)
##   call: malignant      signal > 99th percentile of reference AND corr > 0.4
##         non_malignant  neither criterion met
##         uncertain      only one criterion met
## Date: 2026-10-12

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})

cnv_window   <- 101
cnv_corr_cut <- 0.4
cnv_sig_q    <- 0.99

## autosomal gene positions from Ensembl cDNA FASTA headers
gene_positions <- function(ensembl_dir) {
  fa <- list.files(ensembl_dir, pattern = "cdna.all.fa.gz$", full.names = TRUE)[1]
  hdr <- system(sprintf("gzcat '%s' | grep '^>'", fa), intern = TRUE)
  hdr <- hdr[grepl("chromosome:GRCh38:", hdr) & grepl("gene_symbol:", hdr)]
  loc <- regmatches(hdr, regexpr("chromosome:GRCh38:[^ ]+", hdr))
  sym <- sub("^gene_symbol:", "", regmatches(hdr, regexpr("gene_symbol:[^ ]+", hdr)))
  parts <- do.call(rbind, strsplit(sub("chromosome:GRCh38:", "", loc), ":"))
  pos <- data.table(gene = sym, chr = parts[, 1], start = as.numeric(parts[, 2]))
  pos <- pos[chr %in% as.character(1:22), .(start = min(start)), by = .(gene, chr)]
  pos <- pos[!duplicated(gene)]
  pos[, chr_n := as.integer(chr)]
  setorder(pos, chr_n, start)
  pos
}

## moving average within each chromosome (window truncated at the ends)
smooth_chr <- function(m, chr_runs, window = cnv_window) {
  half <- (window - 1) %/% 2
  out <- m
  for (idx in chr_runs) {
    if (length(idx) < 2) next
    cs <- rbind(0, apply(m[idx, , drop = FALSE], 2, cumsum))
    n <- length(idx)
    lo <- pmax(1, seq_len(n) - half)
    hi <- pmin(n, seq_len(n) + half)
    out[idx, ] <- (cs[hi + 1, , drop = FALSE] - cs[lo, , drop = FALSE]) / (hi - lo + 1)
  }
  out
}

## CNV for one sample.
##   counts   genes x cells (rows = the ordered CNV gene list, may be a subset)
##   lib      library size per cell (from the full gene set)
##   is_ref   reference (immune/stromal) cells
##   is_epi   epithelial cells to score
##   in_prof  epithelial cells allowed to define the tumour-average profile
##   pos      gene positions (ordered); keep_s: return smoothed matrix too
run_cnv <- function(counts, lib, is_ref, is_epi, in_prof = is_epi, pos, keep_s = FALSE) {
  genes <- pos$gene[pos$gene %in% rownames(counts)]
  p <- pos[match(genes, gene)]
  chr_runs <- split(seq_along(genes), factor(p$chr, levels = as.character(1:22)))
  x <- counts[genes, , drop = FALSE]
  e <- as.matrix(x %*% Diagonal(x = 1e5 / lib))
  e <- log2(e + 1)
  r <- e - rowMeans(e[, is_ref, drop = FALSE])
  rm(e)
  r[r > 3] <- 3
  r[r < -3] <- -3
  s <- smooth_chr(r, chr_runs)
  rm(r)
  s <- sweep(s, 2, apply(s, 2, median))
  s <- s - rowMeans(s[, is_ref, drop = FALSE])
  signal <- colMeans(s^2)
  thr <- unname(quantile(signal[is_ref], cnv_sig_q))
  cand <- which(in_prof & is_epi)
  n_top <- min(length(cand), max(10L, ceiling(0.05 * length(cand))))
  profile <- NULL
  corr <- rep(NA_real_, ncol(s))
  if (length(cand) >= 10) {
    top <- cand[order(signal[cand], decreasing = TRUE)][seq_len(n_top)]
    profile <- rowMeans(s[, top, drop = FALSE])
    corr <- as.vector(cor(s, profile))
  }
  call <- fifelse(is.na(corr), fifelse(signal > thr, "uncertain", "non_malignant"),
            fifelse(signal > thr & corr > cnv_corr_cut, "malignant",
            fifelse(signal <= thr & corr <= cnv_corr_cut, "non_malignant", "uncertain")))
  res <- data.table(cell = colnames(counts), is_ref = is_ref, is_epi = is_epi, signal = signal,
                    signal_rel = signal / thr, corr = corr, cnv_call = call)
  ## tumour-average profile of the cells called malignant (>= 10 cells)
  mal <- which(is_epi & call == "malignant")
  mal_profile <- if (length(mal) >= 10) rowMeans(s[, mal, drop = FALSE]) else NULL
  out <- list(cells = res, profile = profile, mal_profile = mal_profile, genes = genes, chr = p$chr,
              thr = thr, n_ref = sum(is_ref), n_profile = if (is.null(profile)) 0L else n_top)
  if (keep_s) out$s <- s
  out
}
