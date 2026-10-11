## Pseudobulk differential expression on the HNSCatlas epithelial compartment.
##
##   Plan A (paired, within patient): tumour vs normal squamous epithelium
##     A1  all paired patients
##     A2  HPV+ patients
##     A3  HPV- patients
##     A4  A2 vs A3 (does the tumour change differ by HPV status?)
##   Plan B (between patients): HPV+ vs HPV- tumour cells, adjusted for study
##
## The normal reference is the clean squamous subset of the atlas'
## non-malignant epithelial cells (scripts/05): gland, ciliated and
## likely-malignant clusters and immune/stromal doublets are excluded, and
## cells with a tumour-like copy-number profile are removed (scripts/05b).
## Raw counts are summed per patient x group (pseudobulk) and tested with
## edgeR quasi-likelihood. A patient needs >= 20 tumour cells, and for Plan A
## also >= 10 normal squamous cells (muscat's default minimum).
## Significant = FDR < 0.05 and |log2FC| >= 1.
##
## Inputs:  data/processed/HNSCatlas/atlas_epithelial.rds      (scripts/02)
##          HNSCC-HPV/results/01_atlas/tables/patient_metadata.csv               (scripts/04)
##          data/processed/HNSCatlas/normal_reference_cells.rds (scripts/05b)
## Outputs: HNSCC-HPV/results/03_DGE/tables/*.csv
##          data/processed/HNSCatlas/pseudobulk.rds (cache, rebuilt when the
##          grouping changes), dge_results.rds
## Run from the repo root:  Rscript scripts/06_dge_pseudobulk.R

suppressPackageStartupMessages({
  library(Matrix)
  library(edgeR)
})

in_dir     <- "data/processed/HNSCatlas"
out_dir    <- "HNSCC-HPV/results/03_DGE/tables"
min_tumour <- 20
min_normal <- 10
fdr_cut    <- 0.05
lfc_cut    <- 1
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ---- pseudobulk: sum raw counts per patient x group -----------------------
pb_file    <- file.path(in_dir, "pseudobulk.rds")
pb_version <- "Tumor + normal squamous, CNV-positive cells removed (scripts/05, 05b)"
pb <- if (file.exists(pb_file)) readRDS(pb_file)
if (is.null(pb) || !identical(pb$version, pb_version)) {
  message("Building pseudobulk from the epithelial object...")
  suppressPackageStartupMessages(library(Seurat))
  epi <- readRDS(file.path(in_dir, "atlas_epithelial.rds"))
  normal_ref <- readRDS(file.path(in_dir, "normal_reference_cells.rds"))
  group <- ifelse(epi$Compartment == "Tumor", "Tumor",
                  ifelse(colnames(epi) %in% normal_ref, "Normal squamous", NA))
  use <- !is.na(group)
  counts <- LayerData(epi, assay = "RNA", layer = "counts")[, use]
  grp <- factor(paste(epi$Patient[use], group[use], sep = "|"))
  ind <- sparseMatrix(i = seq_along(grp), j = as.integer(grp), x = 1,
                      dims = c(length(grp), nlevels(grp)),
                      dimnames = list(colnames(counts), levels(grp)))
  pb <- list(version = pb_version,
             counts  = as.matrix(counts %*% ind),
             samples = data.frame(sample      = levels(grp),
                                  Patient     = sub("\\|.*$", "", levels(grp)),
                                  Compartment = sub("^.*\\|", "", levels(grp)),
                                  n_cells     = as.integer(table(grp))))
  saveRDS(pb, pb_file)
  rm(epi, counts)
  invisible(gc())
}

meta <- read.csv("HNSCC-HPV/results/01_atlas/tables/patient_metadata.csv")
smp <- merge(pb$samples, meta[, c("Patient", "Study", "HPV", "Sex")], by = "Patient")
smp$Compartment <- factor(smp$Compartment, levels = c("Normal squamous", "Tumor"))

## one tidy table per test, written to CSV
de_table <- function(res, name) {
  tab <- topTags(res, n = Inf, sort.by = "PValue")$table
  tab <- data.frame(gene = rownames(tab), tab, row.names = NULL)
  tab$significant <- tab$FDR < fdr_cut & abs(tab$logFC) >= lfc_cut
  write.csv(tab, file.path(out_dir, paste0(name, ".csv")), row.names = FALSE)
  message(sprintf("%-38s up %5d | down %5d | genes tested %d", name,
                  sum(tab$significant & tab$logFC > 0),
                  sum(tab$significant & tab$logFC < 0), nrow(tab)))
  tab
}

## ---- Plan A: tumour vs normal squamous, paired within patient -------------
n_of <- function(comp) {
  tapply(smp$n_cells * (smp$Compartment == comp), smp$Patient, sum)
}
paired <- names(which(n_of("Tumor") >= min_tumour & n_of("Normal squamous") >= min_normal))
sA <- smp[smp$Patient %in% paired, ]
sA <- sA[order(sA$Patient, sA$Compartment), ]
smp$in_planA <- smp$Patient %in% paired
write.csv(smp, file.path(out_dir, "pseudobulk_samples.csv"), row.names = FALSE)

yA <- DGEList(pb$counts[, sA$sample])
keep <- filterByExpr(yA, group = interaction(sA$HPV, sA$Compartment))
yA <- normLibSizes(yA[keep, , keep.lib.sizes = FALSE])
message(sprintf("\nPlan A: %d patients (HPV+ %d, HPV- %d)", length(paired),
                sum(sA$HPV == "HPV+") / 2, sum(sA$HPV == "HPV-") / 2))

patient <- factor(sA$Patient)
tumor   <- as.numeric(sA$Compartment == "Tumor")
hpv_pos <- as.numeric(sA$HPV == "HPV+")

## A1: one tumour effect across all paired patients
d1   <- model.matrix(~ patient + tumor)
fit1 <- glmQLFit(estimateDisp(yA, d1), d1, robust = TRUE)
A1   <- de_table(glmQLFTest(fit1, coef = "tumor"), "A1_tumour_vs_normal_all")

## A2-A4: tumour effect estimated separately within HPV+ and HPV- patients
d2 <- cbind(model.matrix(~ patient),
            tumor_HPVpos = tumor * hpv_pos, tumor_HPVneg = tumor * (1 - hpv_pos))
fit2 <- glmQLFit(estimateDisp(yA, d2), d2, robust = TRUE)
A2 <- de_table(glmQLFTest(fit2, coef = "tumor_HPVpos"), "A2_tumour_vs_normal_HPVpos")
A3 <- de_table(glmQLFTest(fit2, coef = "tumor_HPVneg"), "A3_tumour_vs_normal_HPVneg")
con <- (colnames(d2) == "tumor_HPVpos") - (colnames(d2) == "tumor_HPVneg")
A4 <- de_table(glmQLFTest(fit2, contrast = con), "A4_tumour_change_HPVpos_vs_HPVneg")

## ---- Plan B: HPV+ vs HPV- tumour cells, between patients ------------------
sB <- smp[smp$Compartment == "Tumor" & smp$n_cells >= min_tumour, ]
hpv   <- factor(sB$HPV, levels = c("HPV-", "HPV+"))
study <- factor(sB$Study)
yB <- DGEList(pb$counts[, sB$sample])
keep <- filterByExpr(yB, group = hpv)
yB <- normLibSizes(yB[keep, , keep.lib.sizes = FALSE])
message(sprintf("\nPlan B: %d patients (HPV+ %d, HPV- %d), adjusted for study",
                nrow(sB), sum(hpv == "HPV+"), sum(hpv == "HPV-")))
dB   <- model.matrix(~ study + hpv)
fitB <- glmQLFit(estimateDisp(yB, dB), dB, robust = TRUE)
B    <- de_table(glmQLFTest(fitB, coef = "hpvHPV+"), "B_HPVpos_vs_HPVneg_tumour")

## ---- candidates: cancer-specific AND HPV-linked ---------------------------
wide <- function(genes) {
  out <- data.frame(gene = genes)
  for (nm in c("A1", "A2", "A3", "A4", "B")) {
    t <- get(nm)[, c("gene", "logFC", "FDR")]
    names(t)[2:3] <- paste0(nm, c("_log2FC", "_FDR"))
    out <- merge(out, t, by = "gene", all.x = TRUE)
  }
  out[order(out$B_FDR), ]
}
up <- function(t) t$gene[t$significant & t$logFC > 0]
dn <- function(t) t$gene[t$significant & t$logFC < 0]
cand_pos <- wide(intersect(up(A2), up(B)))   # tumour-up in HPV+, higher in HPV+ tumours
cand_neg <- wide(intersect(up(A3), dn(B)))   # tumour-up in HPV-, higher in HPV- tumours
write.csv(cand_pos, file.path(out_dir, "candidates_HPVpos.csv"), row.names = FALSE)
write.csv(cand_neg, file.path(out_dir, "candidates_HPVneg.csv"), row.names = FALSE)
message(sprintf("\nCandidates: HPV+ %d | HPV- %d", nrow(cand_pos), nrow(cand_neg)))
message("HPV+ candidates: ", paste(cand_pos$gene, collapse = ", "))
message("HPV- candidates: ", paste(cand_neg$gene, collapse = ", "))

## ---- built-in checks ------------------------------------------------------
message("\nCheck: CDKN2A (p16) should be strongly up in HPV+ tumours (B)")
print(B[B$gene == "CDKN2A", c("gene", "logFC", "FDR")], row.names = FALSE)
message("Check: cancer-testis MAGEA4 and viral HPV16-E6 should be up in tumour (A1, A2)")
print(A1[A1$gene == "MAGEA4", c("gene", "logFC", "FDR")], row.names = FALSE)
if ("HPV16-E6" %in% A2$gene) {
  print(A2[A2$gene == "HPV16-E6", c("gene", "logFC", "FDR")], row.names = FALSE)
} else {
  message("HPV16-E6 not tested in A2 (filterByExpr: too few samples express it)")
}

saveRDS(list(A1 = A1, A2 = A2, A3 = A3, A4 = A4, B = B, yA = yA, sA = sA,
             yB = yB, sB = sB, cand_pos = cand_pos, cand_neg = cand_neg),
        file.path(in_dir, "dge_results.rds"))
message("\nTables in ", out_dir)
