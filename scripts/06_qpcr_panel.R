## Step 06: qRT-PCR panel for HPV- tumour vs normal tissue
##
## Purpose: rank qPCR reference genes for tumour vs normal stability (geNorm,
##   NormFinder) in TCGA HPV- paired tissue and in single-cell pseudobulks;
##   rank target genes (tier 1/1b/2 candidates, positive controls, a normal-
##   mucosa control, CA9) by tumour vs normal effect, per-patient consistency,
##   expected expression and an anatomical-site check (oral cavity, larynx);
##   choose a core panel and a reserve list; assemble the primer ordering list.
## Inputs:  results/step04b_DE/tables/ (DE, TCGA paired DE), results/step05_TF/tables/
##          (A_candidates_final.csv, E_master_regulator_shortlist.csv)
##          data/processed/v2/{tcga_hnsc_paired_counts.rds, pseudobulk_step04b.rds,
##          celltype_pseudobulk_v2.rds, atlas_epithelial_labelled.rds, cell_labels_v2.csv.gz}
##          data/raw/tcga_hnsc_star/ (GDC STAR files: TPM), HGNC complete set,
##          cBioPortal hnsc_tcga primary site (downloaded)
##          results/step06_qPCR/tables/primer_*.csv (from scripts/06_primer_design.py)
## Outputs: results/step06_qPCR/ (tables, figures)
##          data/processed/v2/tcga_hnsc_paired_tpm.rds (cache)
## Date: 2026-10-27
## Run order (from the repo root):
##   1. Rscript scripts/06_qpcr_panel.R           (parts A-C; writes the gene list for primers)
##   2. envs/py310/bin/python scripts/06_primer_design.py
##   3. Rscript scripts/06_qpcr_panel.R           (again: part D, panel with primers, figures)

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(edgeR)
  library(ggplot2)
  library(patchwork)
})
source("config/paths.R")
source("scripts/lib_de.R")

out <- paths$step06
tdir <- file.path(out, "tables")
dir.create(tdir, recursive = TRUE, showWarnings = FALSE)
rd4 <- function(f) fread(file.path(paths$step04b, "tables", f))
rd5 <- function(f) fread(file.path(paths$step05, "tables", f))

fdr_cut  <- 0.05
min_cpm  <- 10       # detectable: median malignant pseudobulk CPM (as step 04b)
min_type <- 20
min_mal  <- 20
min_norm <- 10
excl_normal <- "P9"
tpm_min  <- 1        # expected expression for qPCR, TCGA median TPM in the tissue where the gene is high
tpm_ok   <- 10       #   < 1 very low, 1-10 low, >= 10 moderate/high
n_core   <- 10       # core panel targets (incl. positive and normal-mucosa controls)
n_reserve <- 8
n_ref_pilot <- 5     # reference genes to test on the first ~10 patient samples
ref_cand <- c("SDHA", "PSMC5", "TBP", "HPRT1", "ALAS1", "PPIA", "YWHAZ", "GUSB")
pos_up   <- c("MMP1", "PRAME", "MAGEA4", "LAMC2", "PDPN")
norm_ctrl <- c("KRT13", "KRT4")          # KRT13 preferred; KRT4 if KRT13 primers fail
hypoxia  <- "CA9"
down_tf  <- "NFIC"                       # step 05: lower activity in malignant cells (4 lines of evidence)
tf_pairs <- c(FOXM1 = "HMGB3")           # TF and consensus-regulon target, for co-expression in the cohort
priority <- c("NFIC", "ARPC1A", "LY6K")  # core if room, else top of the reserve (decision 2026-10-27)
tme_types <- c("CD4T", "CD8T", "NK", "Bcell", "PlasmaCell", "Macrophage", "Monocyte",
               "panDC", "Neutrophils", "Mast", "Fibroblast", "Endothelial")
oral_sites <- c("Oral Tongue", "Oral Cavity", "Floor of mouth", "Buccal Mucosa", "Alveolar Ridge",
                "Hard Palate", "Lip")

hg <- read_hgnc(paths$hgnc)

## =========================================================================
## A. Data: TCGA HPV- pairs (counts, TPM, site) and single-cell summaries
## =========================================================================
## A1. anatomical site (cBioPortal hnsc_tcga, PRIMARY_SITE_PATIENT)
if (!file.exists(paths$tcga_site)) {
  u <- paste0("https://www.cbioportal.org/api/studies/hnsc_tcga/clinical-data?clinicalDataType=PATIENT",
              "&attributeId=PRIMARY_SITE_PATIENT&projection=SUMMARY&pageSize=10000")
  st <- as.data.table(jsonlite::fromJSON(u))[, .(case = patientId, primary_site = value)]
  fwrite(st, paths$tcga_site)
}
site <- fread(paths$tcga_site)

tc <- readRDS(paths$tcga_counts)
cb <- fread(paths$cbio_pancan, select = c("patientId", "SUBTYPE"))
ts <- copy(tc$samples)[, `:=`(hpv = cb$SUBTYPE[match(case, cb$patientId)],
                              state = factor(fifelse(sample_type == "Primary Tumor", "tumour", "normal"),
                                             levels = c("normal", "tumour")),
                              primary_site = site$primary_site[match(case, site$case)])]
ts <- ts[hpv == "HNSC_HPV-"][, `:=`(patient = factor(case),
                                    site_group = fcase(primary_site %in% oral_sites, "oral_cavity",
                                                       primary_site == "Larynx", "larynx",
                                                       default = "other"))]
fwrite(ts[state == "tumour", .N, by = .(site_group, primary_site)][order(site_group, -N)],
       file.path(tdir, "A_tcga_pairs_by_site.csv"))

## A2. TPM from the same GDC STAR files (symbols mapped as in step 04b)
if (!file.exists(paths$tcga_tpm)) {
  man <- tc$samples
  rdl <- lapply(man$file_name, function(fn) {
    x <- fread(file.path(paths$tcga_dir, fn), skip = 1, select = c("gene_id", "gene_name", "tpm_unstranded"))
    x[grepl("^ENSG", gene_id)]
  })
  stopifnot(all(vapply(rdl, function(x) identical(x$gene_id, rdl[[1]]$gene_id), TRUE)))
  m <- do.call(cbind, lapply(rdl, `[[`, "tpm_unstranded"))
  colnames(m) <- man$sample
  ens <- sub("\\..*$", "", rdl[[1]]$gene_id)
  keep <- !grepl("_PAR_Y$", rdl[[1]]$gene_id)
  sym <- hg$symbol[match(ens, hg$ensembl_gene_id)]
  gmap <- hgnc_map(unique(rdl[[1]]$gene_name), hg)
  sym[is.na(sym)] <- gmap$mapped[match(rdl[[1]]$gene_name[is.na(sym)], gmap$original)]
  keep <- keep & !is.na(sym)
  f <- factor(sym[keep])
  agg <- sparseMatrix(i = as.integer(f), j = which(keep), x = 1, dims = c(nlevels(f), nrow(m)),
                      dimnames = list(levels(f), NULL))
  saveRDS(list(tpm = as.matrix(agg %*% m), samples = man), paths$tcga_tpm)
  rm(rdl, m, agg)
}
tpm <- readRDS(paths$tcga_tpm)$tpm[, ts$sample]
tpm_t <- tpm[, ts[state == "tumour"][order(case), sample]]
tpm_n <- tpm[, ts[state == "normal"][order(case), sample]]
pair_site <- ts[state == "tumour"][order(case), site_group]

## A3. paired tumour vs normal DE per site (edgeR QL, ~ patient + state, as step 04b)
tde_all <- rd4("E_TCGA_HPVneg_paired_tumour_vs_normal.csv")
tde_site <- lapply(c(oral_cavity = "oral_cavity", larynx = "larynx"), function(sg) {
  s <- ts[site_group == sg][, patient := droplevels(patient)]
  de_qlf(tc$counts[, s$sample], s, ~ patient + state, "statetumour", fdr_cut, 1)$table
})
fwrite(data.table(site = names(tde_site), pairs = ts[state == "tumour", .N, by = site_group][match(names(tde_site), site_group), N]),
       file.path(tdir, "A_site_de_design.csv"))

## A4. single-cell: DE (step 04b) and specificity / detectability (method of step 04b part F)
de2 <- rd4("D_C2_HPVneg_allpatients.csv"); de1 <- rd4("D_C1_HPVneg_paired.csv")
pb <- readRDS(paths$pseudobulk_v2b)
smp <- rd4("D_pseudobulk_samples.csv")
map_atlas <- hgnc_map(rownames(readRDS(paths$atlas_labelled)$counts), hg)
invisible(gc())
ct <- readRDS(paths$celltype_pb_v2)
ctc <- map_rows(ct$counts, map_atlas)
lab_pt <- fread(paths$cell_labels_v2, select = c("patient_id_unique", "hpv", "role", "flag_excluded", "primary_4b"))
neg_pts <- unique(lab_pt[role == "discovery" & !flag_excluded & hpv == "HPV-" & primary_4b == TRUE, patient_id_unique])
rm(lab_pt)
cs <- ct$samples[patient %in% neg_pts]
lib_ct <- colSums(ctc)
mal_k <- smp[group == "malignant" & hpv == "HPV-" & primary_4b == TRUE & n_cells >= min_mal]
nor_k <- smp[group == "normal" & hpv == "HPV-" & primary_4b == TRUE & n_cells >= min_norm & patient != excl_normal]
pb_cpm <- function(genes, k, m) {
  g <- intersect(genes, rownames(m))
  r <- setNames(rep(NA_real_, length(genes)), genes)
  r[g] <- apply(as.matrix(m[g, k$key, drop = FALSE]) %*% diag(1e6 / colSums(m[, k$key, drop = FALSE]), nrow = nrow(k)), 1, median)
  r
}
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
  d <- data.table(gene = genes, malignant_median_CPM = tum, max_TME_CPM = mx, max_TME_type = mxt,
                  tumour_share = apply(sp, 1, median, na.rm = TRUE),
                  malignant_pb_CPM = pb_cpm(genes, mal_k, pb$full$counts),
                  normal_pb_CPM = pb_cpm(genes, nor_k, pb$shared$counts))
  d[, specificity := fcase(malignant_median_CPM >= 2 * max_TME_CPM & tumour_share >= 0.5, "tumour_specific",
                           malignant_median_CPM > max_TME_CPM, "tumour_enriched", default = "shared_with_TME")]
  d[]
}

## =========================================================================
## B. Reference genes: tumour vs normal stability
## =========================================================================
## geNorm M (Vandesompele et al. 2002) with stepwise exclusion and pairwise
## variation V, and NormFinder rho (Andersen et al. 2004), transcribed from
## NormqPCR 1.58 (Perkins et al. 2012; LGPL-3): stabMeasureM, selectHKs
## (method = "geNorm"), stabMeasureRho. x: samples x genes, log2 scale.
genorm_M <- function(x) {
  n <- ncol(x)
  M <- vapply(seq_len(n), function(j) {
    A <- x[, j] - x[, -j, drop = FALSE]
    if (n > 2) mean(sqrt(rowSums((t(A) - colMeans(A))^2) / (nrow(A) - 1))) else sd(A)
  }, numeric(1))
  setNames(M, colnames(x))
}
genorm_steps <- function(x) {
  n <- ncol(x); rank <- setNames(rep(NA_integer_, n), colnames(x)); V <- c()
  while (ncol(x) > 2) {
    M <- genorm_M(x); worst <- which.max(M)
    rank[colnames(x)[worst]] <- ncol(x)
    V[paste0(ncol(x) - 1, "/", ncol(x))] <- sd(rowMeans(x[, -worst, drop = FALSE]) - rowMeans(x))
    x <- x[, -worst, drop = FALSE]
  }
  rank[colnames(x)] <- 1L
  list(rank = rank, V = V)
}
normfinder_rho <- function(x, group) {
  group <- factor(group); k <- ncol(x); m <- nlevels(group); Levels <- levels(group)
  ngr <- as.integer(table(group))
  mej <- rowMeans(x)
  meigr <- as.matrix(aggregate(x, by = list(group), mean)[, -1])
  megr <- rowMeans(meigr)
  x.split <- split(x, f = group)
  var.group.all <- matrix(0, ncol = k, nrow = m)
  for (i in seq_len(m)) {
    x.temp <- matrix(x.split[[i]], nrow = ngr[i])
    a <- colSums((t(t(x.temp) - meigr[i, ]) - mej[group == Levels[i]] + megr[i])^2) / (ngr[i] - 1)
    b <- sum(a)
    var.group.all[i, ] <- (a - b / (k * (k - 1))) / (1 - 2 / k)
  }
  var.group.all <- pmax(var.group.all, 0)
  m1i <- colMeans(meigr); m1j <- rowMeans(meigr); m1 <- mean(m1j)
  dif <- t(meigr - m1j) - m1i + m1
  va <- var.group.all / ngr
  tau <- max(sum(dif * dif) / ((m - 1) * (k - 1)) - mean(va), 0)
  dnew <- dif * tau / (tau + t(va))
  vanew <- t(va + tau * va / (tau + va))
  setNames(rowMeans(abs(dnew) + sqrt(vanew)), colnames(x))
}

## B1. TCGA HPV- tissue: log2 TMM-CPM of all 80 samples
y <- DGEList(tc$counts[, ts$sample])
y <- normLibSizes(y[filterByExpr(y, group = ts$state), , keep.lib.sizes = FALSE])
lc_t <- edgeR::cpm(y, log = TRUE, prior.count = 1)
stopifnot(all(ref_cand %in% rownames(lc_t)))
xr <- t(lc_t[ref_cand, ])
gs_t <- genorm_steps(xr)
os <- ts$site_group %in% c("oral_cavity", "larynx")
## B2. single-cell pseudobulks: 27 malignant + 7 normal (shared genes)
kk <- rbind(mal_k, nor_k)
lc_s <- edgeR::cpm(DGEList(as.matrix(pb$shared$counts[, kk$key])), log = TRUE, prior.count = 1)
stopifnot(all(ref_cand %in% rownames(lc_s)))
xs <- t(lc_s[ref_cand, ])
gs_s <- genorm_steps(xs)
rf <- data.table(gene = ref_cand,
                 TCGA_TPM_tumour = apply(tpm_t[ref_cand, ], 1, median), TCGA_TPM_normal = apply(tpm_n[ref_cand, ], 1, median),
                 TCGA_log2FC = tde_all$logFC[match(ref_cand, tde_all$gene)], TCGA_FDR = tde_all$FDR[match(ref_cand, tde_all$gene)],
                 oral_log2FC = tde_site$oral_cavity$logFC[match(ref_cand, tde_site$oral_cavity$gene)],
                 larynx_log2FC = tde_site$larynx$logFC[match(ref_cand, tde_site$larynx$gene)],
                 geNorm_M_TCGA = genorm_M(xr), geNorm_rank_TCGA = gs_t$rank[ref_cand],
                 NormFinder_rho_TCGA = normfinder_rho(xr, ts$state),
                 NormFinder_rho_TCGA_site = normfinder_rho(xr[os, ], paste(ts$site_group, ts$state)[os]),
                 sc_malignant_CPM = pb_cpm(ref_cand, mal_k, pb$shared$counts), sc_normal_CPM = pb_cpm(ref_cand, nor_k, pb$shared$counts),
                 sc_C2_log2FC = de2$logFC[match(ref_cand, de2$gene)], sc_C2_FDR = de2$FDR[match(ref_cand, de2$gene)],
                 geNorm_M_sc = genorm_M(xs), NormFinder_rho_sc = normfinder_rho(xs, kk$group))
## pre-specified combined rank: mean of the ranks of TCGA geNorm (stepwise), TCGA NormFinder,
## |TCGA tumour vs normal log2FC| and single-cell NormFinder
rf[, combined := (geNorm_rank_TCGA + frank(NormFinder_rho_TCGA) + frank(abs(TCGA_log2FC)) + frank(NormFinder_rho_sc)) / 4]
setorder(rf, combined)
rf[, `:=`(final_rank = seq_len(.N), pilot = seq_len(.N) <= n_ref_pilot)]
fwrite(rf, file.path(tdir, "B_reference_gene_stability.csv"))
fwrite(data.table(data = "TCGA", step = names(gs_t$V), V = gs_t$V), file.path(tdir, "B_geNorm_pairwise_variation.csv"))

## =========================================================================
## C. Target genes: ranking, site check, expected expression
## =========================================================================
cf <- rd5("A_candidates_final.csv")
tiers <- cf[tier %in% c("1", "1b", "2"), .(gene, tier, direction)]
tg <- rbind(tiers,
            data.table(gene = setdiff(pos_up, tiers$gene), tier = "positive_control", direction = "up"),
            data.table(gene = norm_ctrl, tier = "normal_control", direction = "down"),
            data.table(gene = setdiff(hypoxia, tiers$gene), tier = "hypoxia", direction = "up"),
            data.table(gene = setdiff(down_tf, tiers$gene), tier = "down_TF", direction = "down"))
genes <- tg$gene
sp <- specificity_of(genes)
tg <- merge(tg, sp, by = "gene", all.x = TRUE)
at <- function(tab, col, g = tg$gene) tab[[col]][match(g, tab$gene)]
tg[, `:=`(sc_C2_log2FC = at(de2, "logFC"), sc_C2_FDR = at(de2, "FDR"), sc_C1_log2FC = at(de1, "logFC"),
          TCGA_log2FC = at(tde_all, "logFC"), TCGA_FDR = at(tde_all, "FDR"),
          oral_log2FC = at(tde_site$oral_cavity, "logFC"), oral_FDR = at(tde_site$oral_cavity, "FDR"),
          larynx_log2FC = at(tde_site$larynx, "logFC"), larynx_FDR = at(tde_site$larynx, "FDR"))]
## per-pair consistency and expected expression from TPM
in_tpm <- tg$gene %in% rownames(tpm)
lr <- matrix(NA_real_, nrow(tg), ncol(tpm_t), dimnames = list(tg$gene, NULL))
lr[in_tpm, ] <- log2((tpm_t[tg$gene[in_tpm], ] + 0.5) / (tpm_n[tg$gene[in_tpm], ] + 0.5))
sgn <- fifelse(tg$direction == "up", 1, -1)
frac2 <- function(cols) rowMeans(sgn * lr[, cols, drop = FALSE] >= 1)
tg[, `:=`(pairs_2fold = frac2(TRUE), pairs_2fold_oral = frac2(pair_site == "oral_cavity"),
          pairs_2fold_larynx = frac2(pair_site == "larynx"),
          TCGA_TPM_tumour = NA_real_, TCGA_TPM_normal = NA_real_)]
tg[in_tpm, `:=`(TCGA_TPM_tumour = apply(tpm_t[gene, , drop = FALSE], 1, median),
                TCGA_TPM_normal = apply(tpm_n[gene, , drop = FALSE], 1, median))]
## expected expression where the gene is high (tumour for up, normal for down)
tg[, TPM_high := fifelse(direction == "up", TCGA_TPM_tumour, TCGA_TPM_normal)]
tg[, expression_class := fcase(is.na(TPM_high), "not in TCGA", TPM_high < tpm_min, "very low (<1 TPM)",
                               TPM_high < tpm_ok, "low (1-10 TPM)", default = "moderate/high (>=10 TPM)")]
## site check (paired TCGA). Oral cavity must pass (expected direction, FDR < 0.05): the study cohort is
## mostly oral cavity. Larynx (11 pairs, low power) excludes a gene only when it is significant in the
## opposite direction; same direction but not significant = "oral-cavity-supported, larynx underpowered"
site_class <- function(lfc, fdr, s) fcase(is.na(fdr), "not tested",
                                          fdr < fdr_cut & sign(lfc) == s, "pass",
                                          fdr < fdr_cut, "opposite direction (significant)",
                                          sign(lfc) == s, "same direction, not significant",
                                          default = "opposite sign, not significant")
tg[, `:=`(oral_class = site_class(oral_log2FC, oral_FDR, sgn), larynx_class = site_class(larynx_log2FC, larynx_FDR, sgn))]
tg[, `:=`(oral_pass = oral_class == "pass", larynx_pass = larynx_class == "pass",
          larynx_excludes = larynx_class == "opposite direction (significant)")]
tg[, site_note := fcase(!oral_pass, paste0("oral cavity: ", oral_class),
                        larynx_pass, "both sites pass",
                        larynx_excludes, "larynx: opposite direction, significant (excluded)",
                        larynx_class == "same direction, not significant", "oral-cavity-supported, larynx underpowered",
                        larynx_class == "opposite sign, not significant", "oral-cavity-supported, larynx no change (opposite sign, ns)",
                        default = paste0("oral-cavity-supported, larynx ", larynx_class))]
## eligibility (pre-specified)
detect_sc <- fifelse(tg$direction == "up", tg$malignant_pb_CPM, tg$normal_pb_CPM)
tg[, eligible := !is.na(TCGA_FDR) & TCGA_FDR < fdr_cut & sign(TCGA_log2FC) == sgn & oral_pass & !larynx_excludes &
                  !is.na(TPM_high) & TPM_high >= tpm_min & !is.na(detect_sc) & detect_sc >= min_cpm]
tg[, why_not := ""]
tg[!(TCGA_FDR < fdr_cut & sign(TCGA_log2FC) == sgn) %in% TRUE, why_not := paste0(why_not, "not up/down in TCGA; ")]
tg[!(oral_pass & !larynx_excludes), why_not := paste0(why_not, "site check: ", site_note, "; ")]
tg[!(TPM_high >= tpm_min) %in% TRUE, why_not := paste0(why_not, "TCGA TPM < 1; ")]
tg[!(detect_sc >= min_cpm) %in% TRUE, why_not := paste0(why_not, "single-cell CPM < 10; ")]
tg[direction == "down" & tier == "1b", why_not := paste0(why_not, "down candidate (the normal-control slot is KRT13/KRT4); ")]
## score among up targets: mean percentile rank of TCGA log2FC, share of pairs >= 2-fold up,
## single-cell C2 log2FC, tumour share of the tissue signal, log TPM in tumour
up <- tg$direction == "up"
pr <- function(v) { r <- frank(v, na.last = "keep"); (r - 1) / (sum(!is.na(v)) - 1) }
tg[, score := NA_real_]
tg[up, score := rowMeans(cbind(pr(TCGA_log2FC), pr(pairs_2fold), pr(sc_C2_log2FC), pr(tumour_share), pr(log2(TCGA_TPM_tumour + 1))),
                         na.rm = TRUE)]
tg[direction == "down" & tier == "1b", eligible := FALSE]

## =========================================================================
## D. Panel: core targets and reserve (needs primers in pass 2)
## =========================================================================
pf <- file.path(tdir, "primer_best.csv")
have_primers <- file.exists(pf)
prm <- if (have_primers) fread(pf) else NULL
primer_ok <- function(g) if (is.null(prm)) rep(TRUE, length(g)) else
  g %in% prm[off_target_products == "none" & detects_canonical == TRUE, gene]
master <- rd5("E_master_regulator_shortlist.csv")$TF
pick <- function(tg) {
  ok <- tg[eligible == TRUE & primer_ok(gene)]
  mtf <- intersect(master, ok$gene)
  core <- c(
    ok[tier == "positive_control"][order(-score)][1, gene],                         # positive control
    intersect(norm_ctrl, ok$gene)[1],                                                 # KRT13, else KRT4
    mtf,                                                                              # FOXM1 (step 05)
    intersect(tf_pairs[mtf], ok$gene),                                                # its regulon target HMGB3
    intersect(hypoxia, ok$gene))                                                      # CA9
  core <- unique(core[!is.na(core)])
  ## whole-tissue qPCR also measures immune and stromal cells: tumour_specific genes first
  fill <- c(ok[tier %in% c("1", "1b") & specificity == "tumour_specific" & !gene %in% core][order(-score), gene],
            ok[tier %in% c("1", "1b", "2") & specificity != "tumour_specific" & !gene %in% core][order(-score), gene],
            ok[tier == "2" & specificity == "tumour_specific" & !gene %in% core][order(-score), gene])
  fill <- unique(fill)
  core <- c(core, head(fill, n_core - length(core)))
  room <- n_core - length(core)                                                     # NFIC if room is left
  if (room > 0) core <- c(core, head(intersect(down_tf, ok$gene), room))
  reserve <- head(setdiff(c(intersect(priority, ok$gene), fill, ok[tier == "positive_control"][order(-score), gene]), core),
                  n_reserve)
  reserve <- c(reserve, setdiff(intersect(norm_ctrl, ok$gene), core))     # backup normal-mucosa control
  list(core = core, reserve = reserve)
}
sel <- pick(tg)
## optional Arm B add-on: CA9 when it fails only the single-cell CPM gate (hypoxia marks a subset of
## tumour cells) but is up in TCGA, passes both sites and is >= 10 TPM in tumour tissue
opt <- tg[gene %in% hypoxia & eligible == FALSE & TCGA_FDR < fdr_cut & TCGA_log2FC > 0 & oral_pass & larynx_pass &
            TPM_high >= tpm_ok & primer_ok(gene), gene]
tg[, panel := fcase(gene %in% sel$core, "core", gene %in% sel$reserve, "reserve", gene %in% opt, "optional (Arm B)",
                    default = "")]
## pilot reference genes: best combined rank among those with a clean primer pair
rf[, pilot := FALSE]
rf[final_rank %in% head(rf[primer_ok(gene), final_rank], n_ref_pilot), pilot := TRUE]
fwrite(rf, file.path(tdir, "B_reference_gene_stability.csv"))
tg[, role := fcase(tier == "positive_control", "positive control (up)", tier == "normal_control", "normal-mucosa control (down)",
                   tier == "hypoxia", "hypoxia marker (Arm B link)", gene %in% master, "target: step 05 master regulator",
                   tier == "down_TF", "down-regulated TF marker (step 05)", gene %in% tf_pairs, "target: tier 1, FOXM1 regulon target",
                   tier == "1", "target: tier 1", tier == "1b", "target: tier 1b", tier == "2", "target: tier 2", default = "")]
setorder(tg, -eligible, -score, na.last = TRUE)
tg[, rank := seq_len(.N)]
fwrite(tg, file.path(tdir, "C_target_ranking.csv"))

## TF-target pairs for co-expression in the cohort: which tier 1/1b genes are in the consensus
## pySCENIC regulon (step 05) or the CollecTRI regulon of a TF on the panel (core or reserve)
reg <- fread(file.path(paths$step05_scenic, "regulons_consensus.csv"))
ctri <- fread(paths$collectri, select = c("source_genesymbol", "target_genesymbol"))
panel_tfs <- intersect(c(sel$core, sel$reserve), unique(c(reg$TF, ctri$source_genesymbol, cf[is_TF == TRUE, gene])))
tt <- cf[tier %in% c("1", "1b"), .(gene, tier, direction)]
tt[, `:=`(panel = tg$panel[match(gene, tg$gene)],
          panel_TF_in_regulon = sapply(gene, function(g) paste(sort(reg[target == g & TF %in% panel_tfs, unique(TF)]), collapse = " ")),
          panel_TF_in_collectri = sapply(gene, function(g) paste(sort(ctri[target_genesymbol == g & source_genesymbol %in% panel_tfs,
                                                                            unique(source_genesymbol)]), collapse = " ")),
          any_TF_in_regulon = sapply(gene, function(g) paste(sort(reg[target == g, unique(TF)]), collapse = " ")))]
fwrite(tt, file.path(tdir, "C_tf_target_pairs.csv"))
fwrite(data.table(panel_TF = panel_tfs, consensus_regulon_targets = sapply(panel_tfs, function(t) reg[TF == t, uniqueN(target)])),
       file.path(tdir, "C_panel_tfs.csv"))

## gene list for primer design (pass 1): references, controls, CA9, all tier 1/1b up, eligible tier 2
if (!have_primers) {
  pl <- unique(rbind(
    data.table(gene = ref_cand, role = "reference gene"),
    tg[tier %in% c("positive_control", "normal_control", "hypoxia", "down_TF"), .(gene, role)],
    tg[tier %in% c("1", "1b") & direction == "up", .(gene, role)],
    tg[tier == "2" & eligible == TRUE][order(-score)][1:min(.N, 12), .(gene, role)]), by = "gene")
  fwrite(pl, file.path(tdir, "C_genes_for_primers.csv"))
  print(tg[panel != "", .(gene, tier, panel, score = round(score, 2), TCGA_log2FC = round(TCGA_log2FC, 2),
                          oral_pass, larynx_pass, TPM_high = round(TPM_high, 1), expression_class)])
  print(rf)
  message("Pass 1 done: run envs/py310/bin/python scripts/06_primer_design.py, then this script again")
  quit(save = "no")
}

## D2. final panel table and ordering list
core_tab <- tg[panel == "core"][order(match(gene, sel$core))]
fin <- rbind(core_tab[, .(set = "core", gene, role)],
             tg[panel == "reserve"][order(match(gene, sel$reserve)), .(set = "reserve", gene, role)],
             tg[panel == "optional (Arm B)", .(set = "optional (Arm B)", gene, role)],
             rf[pilot == TRUE, .(set = "reference (pilot)", gene, role = "reference gene")])
fin <- merge(fin, tg[, .(gene, tier, specificity, sc_C2_log2FC, TCGA_log2FC, TCGA_FDR, oral_log2FC, larynx_log2FC,
                         pairs_2fold, malignant_pb_CPM, normal_pb_CPM, TCGA_TPM_tumour, TCGA_TPM_normal,
                         expression_class, score)], by = "gene", all.x = TRUE, sort = FALSE)
fin <- merge(fin, rf[, .(gene, ref_TPM_tumour = TCGA_TPM_tumour, ref_TPM_normal = TCGA_TPM_normal, ref_final_rank = final_rank)],
             by = "gene", all.x = TRUE, sort = FALSE)
fin <- merge(fin, tg[, .(gene, site_note)], by = "gene", all.x = TRUE, sort = FALSE)
fin <- merge(fin, prm[, .(gene, forward, reverse, product_bp, design, main_isoforms_detected, isoform_specific,
                          off_target_products, pseudogene_products, readthrough_products, transcript, main_isoform_ids,
                          detected_transcript_ids)],
             by = "gene", all.x = TRUE, sort = FALSE)
fin[, set := factor(set, levels = c("core", "reference (pilot)", "optional (Arm B)", "reserve"))]
setorder(fin, set)
fwrite(fin, file.path(tdir, "D_panel_final.csv"))
ord <- fin[set %in% c("core", "reference (pilot)", "optional (Arm B)") & !is.na(forward)]
order_list <- rbind(ord[, .(oligo = paste0(gene, "_F"), sequence = forward, gene, set, product_bp)],
                    ord[, .(oligo = paste0(gene, "_R"), sequence = reverse, gene, set, product_bp)])
order_list[, `:=`(length = nchar(sequence), scale = "25 nmol", purification = "desalted")]
setorder(order_list, set, gene, oligo)
fwrite(order_list, file.path(tdir, "D_primer_order_list.csv"))

## Primer-BLAST-ready table (check manually on NCBI Primer-BLAST, RefSeq mRNA and genome databases);
## pseudogene flag: predicted product from an annotated pseudogene transcript, or HGNC pseudogenes of the gene
hg_pseudo <- hg[locus_group == "pseudogene", symbol]
pb_tab <- fin[!is.na(forward), .(gene, set, role, forward, reverse, product_bp, canonical_transcript = transcript,
                                 target_transcripts = detected_transcript_ids, main_isoform_ids, main_isoforms_detected,
                                 isoform_specific, off_target_products, pseudogene_products, readthrough_products)]
pb_tab[, hgnc_pseudogenes := sapply(gene, function(g) paste(grep(paste0("^", g, "P[0-9]+[A-Z]?$"), hg_pseudo, value = TRUE),
                                                           collapse = " "))]
pb_tab[, pseudogene_flag := fcase(pseudogene_products != "none", "PREDICTED pseudogene product: choose another pair",
                                  hgnc_pseudogenes != "", "gene has HGNC pseudogenes: check the genome database in Primer-BLAST; use DNase and no-RT controls",
                                  default = "none known")]
fwrite(pb_tab, file.path(tdir, "D_primer_blast_check.csv"))

## =========================================================================
## E. Figures
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

## Fig 26: reference gene stability
rl <- rf[, .(gene = factor(gene, levels = rev(gene)), pilot,
             `geNorm M (TCGA)` = geNorm_M_TCGA, `NormFinder rho (TCGA)` = NormFinder_rho_TCGA,
             `NormFinder rho (single cell)` = NormFinder_rho_sc, `|log2FC| tumour vs normal (TCGA)` = abs(TCGA_log2FC))]
rl <- melt(rl, id.vars = c("gene", "pilot"))
f26 <- ggplot(rl, aes(value, gene, fill = pilot)) + geom_col(width = 0.7) +
  facet_wrap(~ variable, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = c(`TRUE` = "#2a78d6", `FALSE` = "#c3c2b7"), labels = c(`TRUE` = "pilot set", `FALSE` = "not chosen"), name = NULL) +
  labs(title = "Reference genes: stability between tumour and normal tissue (HPV-)",
       subtitle = sprintf("TCGA: 40 paired HPV- cases; single cell: %d malignant and %d normal pseudobulks; lower = more stable; genes ordered by combined rank",
                          nrow(mal_k), nrow(nor_k)), x = NULL, y = NULL) +
  theme_fig() + theme(legend.position = "bottom")
save_fig(f26, "Fig26_reference_gene_stability", 180, 70)

## Fig 27: target ranking (core, reserve and controls)
show <- tg[panel != "" | tier %in% c("positive_control", "normal_control", "hypoxia", "down_TF") | tier == "1"]
show[, lab := factor(paste0(gene, fifelse(panel == "core", " *", "")), levels = rev(paste0(gene, fifelse(panel == "core", " *", ""))))]
sl <- melt(show[, .(lab, `single cell C2` = sc_C2_log2FC, `TCGA all` = TCGA_log2FC, `TCGA oral cavity` = oral_log2FC,
                    `TCGA larynx` = larynx_log2FC)], id.vars = "lab", value.name = "log2FC")
pass_l <- melt(show[, .(lab, `single cell C2` = "", `TCGA all` = fifelse(TCGA_FDR < fdr_cut, "", "x"),
                        `TCGA oral cavity` = fifelse(oral_pass, "", "x"),
                        `TCGA larynx` = fcase(larynx_pass, "", larynx_excludes, "x",
                                             larynx_class == "same direction, not significant", "ns",
                                             larynx_class == "opposite sign, not significant", "ns-", default = ""))],
               id.vars = "lab", value.name = "mark")
sl <- merge(sl, pass_l, by = c("lab", "variable"))
p27a <- ggplot(sl, aes(variable, lab, fill = pmax(pmin(log2FC, 6), -6))) + geom_tile(colour = "white") +
  geom_text(data = sl[!is.na(mark) & mark != ""], aes(label = mark), size = 2, colour = ink) +
  scale_fill_gradient2(low = "#2a78d6", mid = "#f0efec", high = "#e34948", name = "log2FC") +
  labs(x = NULL, y = NULL, subtitle = "tumour vs normal log2FC") + theme_fig() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.line = element_blank(), axis.ticks = element_blank())
p27b <- ggplot(show, aes(pmax(TPM_high, 0.1), lab, colour = expression_class)) + geom_point(size = 1.6) +
  geom_vline(xintercept = c(tpm_min, tpm_ok), linewidth = 0.3, linetype = 2, colour = rule) +
  scale_x_log10() + scale_colour_manual(values = c(`very low (<1 TPM)` = "#e34948", `low (1-10 TPM)` = "#eda100",
                                                   `moderate/high (>=10 TPM)` = "#2a78d6", `not in TCGA` = "#898781"), name = NULL) +
  labs(x = "median TPM (TCGA, where high)", y = NULL, subtitle = "expected expression") + theme_fig() +
  theme(axis.text.y = element_blank(), legend.position = "bottom", legend.direction = "vertical")
p27c <- ggplot(show, aes(pairs_2fold, lab)) + geom_col(width = 0.6, fill = "#898781") +
  labs(x = "share of TCGA pairs\n>= 2-fold in expected direction", y = NULL, subtitle = "per-patient consistency") +
  scale_x_continuous(limits = c(0, 1)) + theme_fig() + theme(axis.text.y = element_blank())
f27 <- (p27a | p27b | p27c) + plot_layout(widths = c(1.4, 1, 0.8)) +
  plot_annotation(title = "qPCR target ranking (HPV-): effect, site check, expected expression",
                  subtitle = paste0("* = core panel; tier 1 genes, controls, CA9 and the reserve shown.\n",
                                    "Heatmap: x = fails; ns = same direction, not significant; ",
                                    "ns- = opposite sign, not significant; grey = not tested"),
                  theme = theme(plot.title = element_text(size = base + 1, face = "bold"),
                                plot.subtitle = element_text(size = base - 0.5, colour = ink2)))
save_fig(f27, "Fig27_target_ranking", 180, 30 + 4 * nrow(show))

## Fig 28: core panel and pilot references in TCGA pairs (log2 TPM, tumour vs normal, by site)
fg <- c(sel$core, rf[pilot == TRUE, gene])
pt <- rbindlist(lapply(fg, function(g) data.table(gene = g, case = rep(seq_len(ncol(tpm_t)), 2),
                                                  state = rep(c("normal", "tumour"), each = ncol(tpm_t)),
                                                  site = rep(pair_site, 2), log2TPM = log2(c(tpm_n[g, ], tpm_t[g, ]) + 1))))
pt[, gene := factor(gene, levels = fg)]
f28 <- ggplot(pt, aes(state, log2TPM)) + geom_line(aes(group = case, colour = site), linewidth = 0.2, alpha = 0.5) +
  geom_boxplot(width = 0.35, outlier.shape = NA, fill = NA, linewidth = 0.3) +
  facet_wrap(~ gene, nrow = 2, scales = "free_y") +
  scale_colour_manual(values = c(oral_cavity = "#2a78d6", larynx = "#e34948", other = "#898781"), name = NULL) +
  labs(title = "Core panel and pilot reference genes in paired TCGA HPV- tissue",
       subtitle = "lines join tumour and normal of one case; coloured by anatomical site", x = NULL, y = "log2(TPM + 1)") +
  theme_fig() + theme(legend.position = "bottom")
save_fig(f28, "Fig28_panel_tcga_pairs", 180, 95)

## Fig 29: primer map (amplicon zoom on the canonical transcript)
mp <- fread(file.path(tdir, "primer_map.csv"))
pg <- intersect(fin[set %in% c("core", "reference (pilot)") & !is.na(forward), gene], mp$gene)
flank <- 120
rows <- lapply(pg, function(g) {
  m <- mp[gene == g]; f <- m[kind == "forward"]; r <- m[kind == "reverse"]; tx <- m[kind == "transcript"]
  lo <- max(0, f$start - flank); hi <- min(tx$end, r$end + flank)
  ex <- m[kind == "exon" & end > lo & start < hi][, `:=`(start = pmax(start, lo), end = pmin(end, hi))]
  list(ex = ex[, .(gene, label, x0 = start - f$start, x1 = end - f$start)],
       pr = data.table(gene = g, kind = c("forward", "reverse"), x0 = c(0, r$start - f$start), x1 = c(f$end - f$start, r$end - f$start)),
       amp = data.table(gene = g, x0 = 0, x1 = r$end - f$start))
})
ex <- rbindlist(lapply(rows, `[[`, "ex")); prr <- rbindlist(lapply(rows, `[[`, "pr")); amp <- rbindlist(lapply(rows, `[[`, "amp"))
strip <- setNames(sprintf("%s (%d bp, %s)", pg, fin$product_bp[match(pg, fin$gene)], fin$main_isoforms_detected[match(pg, fin$gene)]), pg)
ex[, gene := factor(strip[gene], levels = strip)]
prr[, gene := factor(strip[gene], levels = strip)]
amp[, gene := factor(strip[gene], levels = strip)]
ex[, shade := factor(as.integer(label) %% 2)]
f29 <- ggplot() +
  geom_rect(data = ex, aes(xmin = x0, xmax = x1, ymin = -0.25, ymax = 0.25, fill = shade), colour = "white", linewidth = 0.4) +
  geom_text(data = ex, aes(x = (x0 + x1) / 2, y = 0, label = paste0("exon ", label)), size = (base - 1.5) / .pt, colour = ink2) +
  geom_segment(data = amp, aes(x = x0, xend = x1, y = -0.42, yend = -0.42), colour = "#e34948", linewidth = 1.2) +
  geom_segment(data = prr[kind == "forward"], aes(x = x0, xend = x1, y = 0.4, yend = 0.4), colour = ink, linewidth = 0.5,
               arrow = arrow(length = unit(1.2, "mm"), type = "closed")) +
  geom_segment(data = prr[kind == "reverse"], aes(x = x1, xend = x0, y = 0.4, yend = 0.4), colour = ink, linewidth = 0.5,
               arrow = arrow(length = unit(1.2, "mm"), type = "closed")) +
  facet_wrap(~ gene, ncol = 3, scales = "free_x") +
  scale_fill_manual(values = c("0" = "#d6d5cf", "1" = "#c3c2b7"), guide = "none") +
  scale_y_continuous(limits = c(-0.55, 0.55)) +
  labs(title = "qPCR primer positions on the canonical transcript (core panel and pilot references)",
       subtitle = "arrows = primers; red = amplicon; boxes = exons; strip: product size, main isoforms detected",
       x = "position relative to the forward primer (bp)", y = NULL) +
  theme_fig() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
                      strip.text = element_text(hjust = 0))
save_fig(f29, "Fig29_primer_map", 180, 25 + 22 * ceiling(length(pg) / 3))

print(fin[, .(set, gene, role, TCGA_log2FC = round(TCGA_log2FC, 2), expression_class, product_bp, main_isoforms_detected,
              off_target_products)])
message("Done.")
