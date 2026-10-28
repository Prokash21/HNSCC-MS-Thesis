## Step 05c (part B-C): gene atlas matrix, annotations, correlations
##
## Purpose: assemble the 24-gene atlas matrix (log2FC and FDR from step 04b and
##   recomputed step 06 site DE; single-cell expression by cell type with the
##   step 04b part F method; bulk TPM in TCGA and GTEx; blood background;
##   outcome log HR from steps 07a-07c), the right annotations (specificity,
##   network membership, blood-low), LFC concordance statistics and gene-gene
##   Spearman co-expression in TCGA (all, oral) and single-cell malignant
##   pseudobulks.
## Inputs:  results/step04b_DE/tables/{D_C2_*, D_C1_*, E_TCGA_*, D_pseudobulk_samples}.csv,
##          data/processed/v2/{tcga_hnsc_paired_counts, tcga_hnsc_paired_tpm, pseudobulk_step04b, celltype_pseudobulk_v2,
##          atlas_epithelial_labelled, step07a_tcga_hpvneg}.rds, data/processed/v2/cell_labels_v2.csv.gz,
##          data/raw/step01_metadata/{cbioportal_hnsc_tcga_pan_can_atlas_2018_patient, cbioportal_hnsc_tcga_primary_site}.csv,
##          results/step06_qPCR/tables/D_panel_final.csv, results/step07a_recurrence/tables/S_gene_cox.csv,
##          results/step07b_validation/tables/V_gene_cox.csv, results/step07c_oral/tables/O_gene_cox.csv,
##          data/processed/v2/step05_scenic/regulons_consensus.csv, data/reference/collectri/collectri_omnipath.tsv,
##          results/step05b_functional/tables/E_gene_target_membership.csv, results/step05c_gene_atlas/tables/A_*.csv
## Outputs: results/step05c_gene_atlas/tables/{B_*, C_*}.csv
## Date: 2026-10-28

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(edgeR)
})
source("config/paths.R")
source("scripts/lib_de.R")

tdir <- file.path(paths$step05c, "tables")
rd4 <- function(f) fread(file.path(paths$step04b, "tables", f))
rd <- function(f) fread(file.path(tdir, f))
fdr_cut <- 0.05
min_type <- 20; min_mal <- 20; min_norm <- 10; excl_normal <- "P9"
tme_types <- c("CD4T", "CD8T", "NK", "Bcell", "PlasmaCell", "Macrophage", "Monocyte",
               "panDC", "Neutrophils", "Mast", "Fibroblast", "Endothelial")
oral_sites <- c("Oral Tongue", "Oral Cavity", "Floor of mouth", "Buccal Mucosa", "Alveolar Ridge",
                "Hard Palate", "Lip")
groups <- list(identity = c("LAMC2", "HOXD11", "HOXB7", "CASC9", "IGF2BP3", "LY6K", "SLC2A1", "SLC39A4", "ARPC1A"),
               network = c("E2F1", "MYBL2", "FOXM1", "HMGB3"), outcome = c("KIF23", "HMMR", "P4HA1"),
               `down marker` = c("KRT13", "NFIC"), optional = "CA9", reference = c("PSMC5", "TBP", "GUSB", "SDHA", "ALAS1"))
gl <- data.table(gene = unlist(groups), group = rep(names(groups), lengths(groups)))
genes <- gl$gene
targets <- gl[group != "reference", gene]
net_block <- c("E2F1", "MYBL2", "FOXM1", "KIF23", "HMGB3", "HMMR")
hg <- read_hgnc(paths$hgnc)
free <- function() invisible(gc(full = TRUE))

## =========================================================================
## B1. log2FC and FDR: step 04b tables + TCGA site DE recomputed (step 06 code)
## =========================================================================
lfc_of <- function(d, nm) d[gene %in% genes, setNames(.(gene, logFC, FDR), c("gene", paste0("lfc_", nm), paste0("fdr_", nm)))]
tc <- readRDS(paths$tcga_counts)
cb <- fread(paths$cbio_pancan, select = c("patientId", "SUBTYPE"))
site <- fread(paths$tcga_site)
ts <- copy(tc$samples)[, `:=`(hpv = cb$SUBTYPE[match(case, cb$patientId)],
                              state = factor(fifelse(sample_type == "Primary Tumor", "tumour", "normal"), levels = c("normal", "tumour")),
                              primary_site = site$primary_site[match(case, site$case)])]
ts <- ts[hpv == "HNSC_HPV-"][, `:=`(patient = factor(case),
                                    site_group = fcase(primary_site %in% oral_sites, "oral_cavity", primary_site == "Larynx", "larynx",
                                                       default = "other"))]
site_de <- lapply(c(oral_cavity = "oral_cavity", larynx = "larynx"), function(sg) {
  s <- ts[site_group == sg][, patient := droplevels(patient)]
  de_qlf(tc$counts[, s$sample], s, ~ patient + state, "statetumour", fdr_cut, 1)$table
})
lfc <- Reduce(function(a, b) merge(a, b, by = "gene", all = TRUE), list(
  data.table(gene = genes),
  lfc_of(rd4("D_C2_HPVneg_allpatients.csv"), "C2"), lfc_of(rd4("D_C1_HPVneg_paired.csv"), "C1"),
  lfc_of(rd4("E_TCGA_HPVneg_paired_tumour_vs_normal.csv"), "TCGA"), lfc_of(site_de$oral_cavity, "TCGA_oral"),
  lfc_of(site_de$larynx, "TCGA_larynx"), lfc_of(rd4("D_C2_HPVpos_allpatients.csv"), "C2_HPVpos")))
## check against the step 06 panel values
pan <- fread(file.path(paths$step06, "tables", "D_panel_final.csv"))[!is.na(oral_log2FC)]
chk <- merge(pan[, .(gene, oral_log2FC, larynx_log2FC)], lfc[, .(gene, lfc_TCGA_oral, lfc_TCGA_larynx)], by = "gene")
chk[, `:=`(max_abs_diff = pmax(abs(oral_log2FC - lfc_TCGA_oral), abs(larynx_log2FC - lfc_TCGA_larynx)))]
fwrite(chk, file.path(tdir, "B_site_de_check.csv"))
stopifnot(max(chk$max_abs_diff) < 1e-6)
fwrite(data.table(site = c("oral_cavity", "larynx"), pairs = ts[state == "tumour", .N, by = site_group][match(c("oral_cavity", "larynx"), site_group), N],
                  genes_tested = sapply(site_de, nrow)), file.path(tdir, "B_site_de_design.csv"))
rm(tc); free()

## =========================================================================
## B2. Single-cell expression and specificity (step 04b part F method, as step 06)
## =========================================================================
pb <- readRDS(paths$pseudobulk_v2b)
smp <- rd4("D_pseudobulk_samples.csv")
map_atlas <- hgnc_map(rownames(readRDS(paths$atlas_labelled)$counts), hg)
free()
ct <- readRDS(paths$celltype_pb_v2)
ctc <- map_rows(ct$counts, map_atlas)
lab_pt <- fread(paths$cell_labels_v2, select = c("patient_id_unique", "hpv", "role", "flag_excluded", "primary_4b"))
neg_pts <- unique(lab_pt[role == "discovery" & !flag_excluded & hpv == "HPV-" & primary_4b == TRUE, patient_id_unique])
rm(lab_pt)
cs <- ct$samples[patient %in% neg_pts]
lib_ct <- colSums(ctc)
mal_k <- smp[group == "malignant" & hpv == "HPV-" & primary_4b == TRUE & n_cells >= min_mal]
nor_k <- smp[group == "normal" & hpv == "HPV-" & primary_4b == TRUE & n_cells >= min_norm & patient != excl_normal]
pb_cpm <- function(g, k, m) {
  r <- setNames(rep(NA_real_, length(g)), g); gg <- intersect(g, rownames(m))
  r[gg] <- apply(as.matrix(m[gg, k$key, drop = FALSE]) %*% diag(1e6 / colSums(m[, k$key, drop = FALSE]), nrow = nrow(k)), 1, median)
  r
}
gg <- intersect(genes, rownames(ctc))
cp <- as.matrix(ctc[gg, cs$key, drop = FALSE]) %*% diag(1e6 / lib_ct[cs$key], nrow = nrow(cs))
colnames(cp) <- cs$key
med <- function(t) { j <- cs[group == t & n_cells >= min_type, key]; apply(cp[, j, drop = FALSE], 1, median) }
tum <- med("malignant"); tme <- sapply(tme_types, med)
spx <- sapply(cs[group == "malignant" & n_cells >= min_type, patient], function(p) {
  tot <- Matrix::rowSums(ctc[gg, cs[patient == p, key], drop = FALSE])
  ifelse(tot > 0, ctc[gg, cs[patient == p & group == "malignant", key]] / tot, NA_real_)
})
sc <- data.table(gene = gg, malignant_median_CPM = tum, max_TME_CPM = apply(tme, 1, max),
                 max_TME_type = tme_types[apply(tme, 1, which.max)], tumour_share = apply(spx, 1, median, na.rm = TRUE))
sc[, specificity := fcase(malignant_median_CPM >= 2 * max_TME_CPM & tumour_share >= 0.5, "tumour_specific",
                          malignant_median_CPM > max_TME_CPM, "tumour_enriched", default = "shared_with_TME")]
sc <- cbind(sc, setNames(as.data.table(tme), paste0("sc_", tme_types)))
sc[, `:=`(sc_malignant = pb_cpm(gene, mal_k, pb$full$counts), sc_normal = pb_cpm(gene, nor_k, pb$shared$counts))]
## check against step 06 panel values (same method)
chk2 <- merge(pan[, .(gene, mal6 = malignant_pb_CPM, nor6 = normal_pb_CPM, spec6 = specificity)],
              sc[, .(gene, sc_malignant, sc_normal, specificity)], by = "gene")
fwrite(chk2, file.path(tdir, "B_single_cell_check.csv"))
stopifnot(max(abs(chk2$mal6 - chk2$sc_malignant)) < 1e-6, all(chk2$spec6 == chk2$specificity))

## pseudobulk log CPM of malignant pseudobulks for co-expression (TMM, prior 2)
ym <- DGEList(as.matrix(pb$full$counts[, mal_k$key]))
ym <- normLibSizes(ym)
lc_sc <- edgeR::cpm(ym, log = TRUE, prior.count = 2)[intersect(targets, rownames(ym)), ]
rm(pb, ct, ctc, cp); free()

## =========================================================================
## B3. Bulk tissue (TCGA, GTEx), blood background, outcome, network
## =========================================================================
t7 <- readRDS(paths$step07a_tcga)
tpm <- readRDS(paths$tcga_tpm)$tpm
nor_s <- ts[state == "normal", sample]
bulk <- data.table(gene = genes, TCGA_tumour_median_TPM = t7$tpm_median[genes],
                   TCGA_normal_median_TPM = apply(tpm[genes, nor_s, drop = FALSE], 1, median))
gt <- rd("A_gtex_median_tpm.csv"); bl <- rd("A_blood_background.csv")
s7a <- fread(file.path(paths$step07a, "tables", "S_gene_cox.csv"))[endpoint == "PFI" & model == "multivariable", .(gene, hr_TCGA_PFI = HR, p_TCGA_PFI = p)]
s7c <- fread(file.path(paths$step07c, "tables", "O_gene_cox.csv"))[endpoint == "PFI" & model == "main", .(gene, hr_TCGA_oral_PFI = HR, p_TCGA_oral_PFI = p)]
s7b <- fread(file.path(paths$step07b, "tables", "V_gene_cox.csv"))[cohort == "GSE41613" & endpoint == "DSS" & model == "multivariable",
                                                                   .(gene, hr_GSE41613_DSS = HR, p_GSE41613_DSS = p)]
## network membership (step 05b networks): TF / target of E2F1, MYBL2 or FOXM1 / none
scen <- unique(fread(file.path(paths$step05_scenic, "regulons_consensus.csv"))[, .(TF, target)])
raw <- fread(paths$collectri)
cl <- raw[!grepl("_", source_genesymbol), .(TF = source_genesymbol, target = target_genesymbol)]
nm <- hgnc_map(unique(c(cl$TF, cl$target)), hg)
cl[, `:=`(TF = nm$mapped[match(TF, nm$original)], target = nm$mapped[match(target, nm$original)])]
nets <- unique(rbind(scen, cl[!is.na(TF) & !is.na(target)]))
core3 <- c("E2F1", "MYBL2", "FOXM1")
netm <- data.table(gene = genes)
netm[, regulated_by := sapply(gene, function(g) paste(sort(unique(nets[TF %in% core3 & target == g & TF != g, TF])), collapse = ";"))]
netm[, network := fcase(gene %in% core3, "TF (E2F1/MYBL2/FOXM1)", regulated_by != "", "target", default = "none")]
m5b <- fread(file.path(paths$step05b, "tables", "E_gene_target_membership.csv"))
m5b[, core3_5b := sapply(targeted_by, function(x) paste(sort(unique(setdiff(sub("_.*", "", trimws(strsplit(x, ";")[[1]])), "NFIC"))), collapse = ";"))]
chk3 <- merge(m5b[, .(gene, core3_5b)], netm[, .(gene, regulated_by)], by = "gene")
fwrite(chk3, file.path(tdir, "B_network_check.csv"))

## =========================================================================
## B4. Atlas matrix
## =========================================================================
atlas <- Reduce(function(a, b) merge(a, b, by = "gene", all.x = TRUE),
                list(gl, lfc, sc, bulk, gt, bl[, .(gene, PBMC_HPVneg_max, PBMC_HPVpos_max, PBMC_HD_max, PBMC_max_CPM, PBMC_max_type,
                                                    PBMC_max_UMI, blood_low, blood_low_lib1M_posthoc)],
                     s7a, s7c, s7b, netm))
atlas[, gene := factor(gene, genes)]
setorder(atlas, gene)
fwrite(atlas, file.path(tdir, "B_gene_atlas_matrix.csv"))
print(atlas[, .(gene, group, C2 = round(lfc_C2, 2), TCGA = round(lfc_TCGA, 2), oral = round(lfc_TCGA_oral, 2),
                lar = round(lfc_TCGA_larynx, 2), specificity, network, blood_low)])

## =========================================================================
## C1. LFC concordance (all genes tested in both)
## =========================================================================
c2 <- rd4("D_C2_HPVneg_allpatients.csv"); c1 <- rd4("D_C1_HPVneg_paired.csv"); tt <- rd4("E_TCGA_HPVneg_paired_tumour_vs_normal.csv")
conc <- rbindlist(lapply(list(C2 = c2, C1 = c1), function(d) {
  x <- merge(d[, .(gene, sc = logFC)], tt[, .(gene, tcga = logFC)], by = "gene")
  data.table(n = nrow(x), pearson = cor(x$sc, x$tcga), spearman = cor(x$sc, x$tcga, method = "spearman"),
             pearson_24 = cor(x[gene %in% genes]$sc, x[gene %in% genes]$tcga),
             spearman_24 = cor(x[gene %in% genes]$sc, x[gene %in% genes]$tcga, method = "spearman"), n_24 = sum(x$gene %in% genes))
}), idcol = "single_cell")
fwrite(conc, file.path(tdir, "C_lfc_concordance.csv"))
print(conc)

## =========================================================================
## C2. Gene-gene co-expression (Spearman) among the 19 target genes
## =========================================================================
oral_s <- t7$samples[case %in% site[primary_site %in% oral_sites, case], sample]
mats <- list(TCGA_all = t7$lcpm[targets, ], TCGA_oral = t7$lcpm[targets, oral_s], sc_malignant = lc_sc)
co <- rbindlist(lapply(names(mats), function(k) {
  m <- mats[[k]]
  pr <- t(combn(rownames(m), 2))
  r <- rbindlist(lapply(seq_len(nrow(pr)), function(i) {
    ct <- suppressWarnings(cor.test(m[pr[i, 1], ], m[pr[i, 2], ], method = "spearman", exact = FALSE))
    data.table(gene1 = pr[i, 1], gene2 = pr[i, 2], rho = unname(ct$estimate), p = ct$p.value)
  }))
  r[, `:=`(dataset = k, n = ncol(m), FDR = p.adjust(p, "BH"))]
}))
co[, network_pair := gene1 %in% net_block & gene2 %in% net_block]
fwrite(co, file.path(tdir, "C_coexpression.csv"))
print(co[, .(n = n[1], pairs = .N, p05 = sum(p < 0.05), fdr05 = sum(FDR < 0.05), net_pairs = sum(network_pair),
             net_median_rho = median(rho[network_pair]), other_median_rho = median(rho[!network_pair])), by = dataset])
