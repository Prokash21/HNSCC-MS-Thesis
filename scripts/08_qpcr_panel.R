## qRT-PCR panel screen for the HPV+ candidates from scripts/06.
##
## For each candidate: is it expressed mainly by tumour cells (or also by
## immune / stromal cells, which whole-tissue qPCR would pick up too), how much
## of a whole-tissue signal would come from tumour cells, is it expressed highly
## enough for qPCR, and does it still differ by HPV status in whole tissue?
## Also pre-screens qPCR reference genes for stability.
##
## The panel follows the strategy note in HNSCC-HPV/docs/: one established
## HPV+ gene as a positive control, ZFR2, one lncRNA and one little-studied
## gene, plus the planned immune / stromal markers and reference genes.
## "Novelty" classes are copied from that note, not verified here.
##
## Inputs:  HNSCC-HPV/atlas/HNSCC_Atlas.rds (only to build the cache)
##          data/processed/HNSCatlas/normal_reference_cells.rds    (scripts/05b)
##          HNSCC-HPV/results/01_atlas/tables/patient_metadata.csv (scripts/04)
##          HNSCC-HPV/results/03_DGE/tables/*.csv                  (scripts/06)
## Outputs: HNSCC-HPV/results/04_qPCR_panel/ (figures + tables/)
##          data/processed/HNSCatlas/celltype_pseudobulk.rds (cache)
## Run from the repo root:  Rscript scripts/08_qpcr_panel.R

suppressPackageStartupMessages({
  library(Matrix)
  library(edgeR)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

in_dir  <- "data/processed/HNSCatlas"
dge_dir <- "HNSCC-HPV/results/03_DGE/tables"
out_dir <- "HNSCC-HPV/results/04_qPCR_panel"
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

min_cells  <- 20     # cells per patient x group to use that group
min_normal <- 10     # as in scripts/06
min_cpm    <- 10     # "detectable by qPCR": median CPM in HPV+ tumour cells
min_share  <- 0.5    # tumour cells give >= half of the whole-tissue signal
min_spec   <- 1      # >= 2x higher in tumour than in any immune/stromal type
sorted_study <- "Kurten 2021 (GSE164690)"   # CD45-sorted: no whole-tissue view

## ---- cache: counts per patient x cell group, expressing cells per group ---
cache_file    <- file.path(in_dir, "celltype_pseudobulk.rds")
cache_version <- "normal squamous = normal_reference_cells.rds (scripts/05b)"
pbc <- if (file.exists(cache_file)) readRDS(cache_file)
if (is.null(pbc) || !identical(pbc$version, cache_version)) {
  message("Building cell-type pseudobulk from the full atlas (slow)...")
  suppressPackageStartupMessages(library(Seurat))
  atlas  <- readRDS("HNSCC-HPV/atlas/HNSCC_Atlas.rds")
  counts <- LayerData(atlas, assay = "RNA", layer = "counts")
  md     <- atlas@meta.data[, c("Patient", "Source", "Cell_Type")]
  rm(atlas)
  invisible(gc())
  normal_ref <- readRDS(file.path(in_dir, "normal_reference_cells.rds"))
  ct <- ifelse(is.na(md$Cell_Type), "Unassigned", md$Cell_Type)
  group <- ct
  group[ct == "Tumor"] <- "Tumour"
  group[ct == "Tumor" & md$Source == "CD45p"] <- "Other"   # doublets in CD45+ sorts
  group[ct == "Epithelial"] <- "Epithelial, other"
  group[rownames(md) %in% normal_ref] <- "Normal squamous"
  group[ct %in% c("Multi", "Unassigned")] <- "Other"
  key <- factor(paste(md$Patient, group, sep = "|"))
  ind <- sparseMatrix(i = seq_along(key), j = as.integer(key), x = 1,
                      dims = c(length(key), nlevels(key)),
                      dimnames = list(colnames(counts), levels(key)))
  ## stored entries of a counts matrix are > 0, so row hits = expressing cells
  detected <- vapply(split(seq_len(ncol(counts)), group), function(j) {
    tabulate(counts[, j, drop = FALSE]@i + 1L, nbins = nrow(counts))
  }, integer(nrow(counts)))
  rownames(detected) <- rownames(counts)
  saveRDS(list(version  = cache_version,
               counts   = counts %*% ind,
               samples  = data.frame(sample  = levels(key),
                                     Patient = sub("\\|.*$", "", levels(key)),
                                     Group   = sub("^.*\\|", "", levels(key)),
                                     n_cells = as.integer(table(key))),
               detected = detected,
               n_group  = c(table(group))),
          cache_file)
  rm(counts, ind)
  invisible(gc())
  pbc <- readRDS(cache_file)
}

meta <- read.csv("HNSCC-HPV/results/01_atlas/tables/patient_metadata.csv")
smp  <- merge(pbc$samples, meta[, c("Patient", "Study", "HPV")], by = "Patient")
smp  <- smp[match(colnames(pbc$counts), smp$sample), ]
hpv_pos <- meta$Patient[meta$HPV == "HPV+"]

dge  <- function(f) read.csv(file.path(dge_dir, f))
A1   <- dge("A1_tumour_vs_normal_all.csv")
A2   <- dge("A2_tumour_vs_normal_HPVpos.csv")
B    <- dge("B_HPVpos_vs_HPVneg_tumour.csv")
cand <- dge("candidates_HPVpos.csv")$gene

markers   <- c("GZMK", "CD8A", "COL1A2", "SFRP2")       # planned TME markers
ref_genes <- c("RPLP0", "HPRT1", "PSMC5", "RPL30", "RPL27", "OAZ1", "MTCH1", "GAPDH",
               "ACTB", "TBP", "B2M", "PPIA", "ALAS1", "RPS18", "SDHA", "YWHAZ",
               "GUSB", "PGK1", "HMBS", "UBC", "TFRC")
## novelty verdicts as given in the strategy note (HNSCC-HPV/docs/)
novelty <- c(SYCP2 = "Established", STAG3 = "Established", SMC1B = "Established",
             APOBEC3B = "Established", CDKN2A = "Established", E2F1 = "Established",
             MCM2 = "Established", MCM5 = "Established", TYMS = "Established",
             CCNE2 = "Established", SYNGR3 = "Partial", "MIR9-3HG" = "Partial",
             LINC00885 = "Partial", ZFR2 = "Partial", "FOXD3-AS1" = "Partial",
             KLHL35 = "Partial", NTS = "Open", WDR76 = "Open", NEURL1 = "Open",
             PODXL2 = "Open", RIBC2 = "Open")

genes <- intersect(unique(c(cand, markers, ref_genes)), rownames(pbc$counts))
lib   <- colSums(pbc$counts)
cpm_g <- sweep(as.matrix(pbc$counts[genes, , drop = FALSE]), 2, lib, "/") * 1e6

## median CPM of a cell group across patients with enough cells of that group
med_cpm <- function(group, patients = NULL, min_n = min_cells) {
  s <- smp$Group == group & smp$n_cells >= min_n &
    (is.null(patients) | smp$Patient %in% patients)
  apply(cpm_g[, s, drop = FALSE], 1, median)
}
other_groups <- c("CD4T", "CD8T", "NK", "Bcell", "PlasmaCell", "Macrophage", "Monocyte",
                  "panDC", "Neutrophils", "Mast", "Fibroblast", "Endothelial")
tum_pos <- med_cpm("Tumour", hpv_pos)
tum_neg <- med_cpm("Tumour", setdiff(meta$Patient, hpv_pos))
norm_sq <- med_cpm("Normal squamous", min_n = min_normal)
other   <- sapply(other_groups, med_cpm)
max_other      <- apply(other, 1, max)
max_other_type <- other_groups[apply(other, 1, which.max)]

## share of each gene's whole-sample counts coming from tumour cells, per HPV+
## patient of the unsorted studies, then the median across patients
share_pts <- smp$Patient[smp$Group == "Tumour" & smp$n_cells >= min_cells &
                         smp$HPV == "HPV+" & smp$Study != sorted_study]
share <- sapply(share_pts, function(p) {
  tot <- rowSums(cpm_g[, smp$Patient == p, drop = FALSE] *
                 rep(lib[smp$Patient == p], each = length(genes)) / 1e6)
  tum <- cpm_g[, smp$Patient == p & smp$Group == "Tumour"] *
         lib[smp$Patient == p & smp$Group == "Tumour"] / 1e6
  ifelse(tot > 0, tum / tot, NA)
})
tissue_share <- apply(share, 1, median, na.rm = TRUE)

## ---- whole tissue: all cells of a patient summed (unsorted studies) -------
wt_pts <- unique(smp$Patient[smp$Study != sorted_study])
wt_counts <- sapply(wt_pts, function(p) rowSums(pbc$counts[, smp$Patient == p, drop = FALSE]))
wt_meta <- meta[match(wt_pts, meta$Patient), ]
hpv   <- factor(wt_meta$HPV, levels = c("HPV-", "HPV+"))
study <- factor(wt_meta$Study)
y <- DGEList(wt_counts)
y <- normLibSizes(y[filterByExpr(y, group = hpv), , keep.lib.sizes = FALSE])
d <- model.matrix(~ study + hpv)
fit <- glmQLFit(estimateDisp(y, d), d, robust = TRUE)
wt <- topTags(glmQLFTest(fit, coef = "hpvHPV+"), n = Inf)$table
wt <- data.frame(gene = rownames(wt), wt, row.names = NULL)
write.csv(wt, file.path(out_dir, "tables", "wholetissue_HPVpos_vs_HPVneg.csv"), row.names = FALSE)
lc_wt <- cpm(DGEList(wt_counts), log = TRUE, prior.count = 2)
message(sprintf("Whole tissue: %d patients (HPV+ %d, HPV- %d), study-adjusted",
                length(wt_pts), sum(hpv == "HPV+"), sum(hpv == "HPV-")))

## ---- candidate screen table -----------------------------------------------
at <- function(tab, col) tab[[col]][match(cand, tab$gene)]
scr <- data.frame(gene = cand,
                  novelty = ifelse(cand %in% names(novelty), novelty[cand], "Not assessed"),
                  A2_log2FC = at(A2, "logFC"), A2_FDR = at(A2, "FDR"),
                  B_log2FC = at(B, "logFC"), B_FDR = at(B, "FDR"),
                  tumour_HPVpos_CPM = tum_pos[cand], tumour_HPVneg_CPM = tum_neg[cand],
                  normal_squamous_CPM = norm_sq[cand],
                  max_immune_stromal_CPM = max_other[cand],
                  max_immune_stromal_type = max_other_type[match(cand, genes)],
                  specificity_log2 = log2((tum_pos[cand] + 1) / (max_other[cand] + 1)),
                  tissue_tumour_share = tissue_share[cand],
                  pct_tumour_cells_expressing = 100 * pbc$detected[cand, "Tumour"] / pbc$n_group[["Tumour"]],
                  wholetissue_log2FC = at(wt, "logFC"), wholetissue_FDR = at(wt, "FDR"),
                  row.names = NULL)
scr$detectable      <- scr$tumour_HPVpos_CPM >= min_cpm
scr$tumour_dominant <- scr$tissue_tumour_share >= min_share & scr$specificity_log2 >= min_spec
scr$wholetissue_HPV <- !is.na(scr$wholetissue_FDR) & scr$wholetissue_FDR < 0.05 &
                       scr$wholetissue_log2FC > 0
scr$passes_all <- with(scr, detectable & tumour_dominant & wholetissue_HPV)
scr <- scr[order(!scr$passes_all, scr$wholetissue_FDR), ]
write.csv(scr, file.path(out_dir, "tables", "candidate_screen.csv"), row.names = FALSE)
message(sprintf("Candidates: %d | detectable %d | tumour-dominant %d | HPV+ higher in whole tissue %d | all three %d",
                nrow(scr), sum(scr$detectable), sum(scr$tumour_dominant, na.rm = TRUE),
                sum(scr$wholetissue_HPV), sum(scr$passes_all, na.rm = TRUE)))

## ---- panel: the note's preferred gene per role if it passes, else the best
## passing gene of that kind (smallest whole-tissue FDR) --------------------
is_lnc <- grepl("^LINC|-AS[0-9]+$|HG$", scr$gene)
ok <- scr$passes_all %in% TRUE
## which filters a gene fails, for the panel notes
why <- function(g) {
  r <- scr[scr$gene == g, ]
  paste(c(if (!r$detectable) sprintf("expression %.1f CPM < %d", r$tumour_HPVpos_CPM, min_cpm),
          if (!isTRUE(r$tumour_dominant))
            sprintf("tumour share %.0f%%, specificity %.1f", 100 * r$tissue_tumour_share, r$specificity_log2),
          if (!r$wholetissue_HPV) sprintf("whole-tissue FDR %.2g", r$wholetissue_FDR)),
        collapse = "; ")
}
chosen <- character(0)
choose <- function(role, preferred, rows) {
  pass <- setdiff(intersect(preferred, scr$gene[ok]), chosen)
  if (length(pass)) {
    g <- pass[1]
    note <- "preferred in strategy note"
  } else {
    r <- scr[rows & ok & !scr$gene %in% chosen, ]
    g <- if (nrow(r)) r$gene[1] else NA_character_
    note <- sprintf("%s failed (%s); best passing alternative", preferred[1], why(preferred[1]))
  }
  chosen <<- c(chosen, g)
  data.frame(role = role, gene = g, note = note)
}
panel <- rbind(
  choose("Positive control (established HPV+ gene)", "SYCP2", scr$novelty == "Established"),
  choose("Novelty: under-studied, seen in HPV+ lists", "ZFR2", scr$novelty == "Partial" & !is_lnc),
  choose("Novelty: lncRNA", c("MIR9-3HG", "LINC00885"), is_lnc),
  choose("Novelty: little studied (highest risk)", "RIBC2", scr$novelty != "Established"))
panel <- cbind(panel, scr[match(panel$gene, scr$gene),
                          c("novelty", "tumour_HPVpos_CPM", "tissue_tumour_share",
                            "specificity_log2", "B_log2FC", "wholetissue_log2FC", "wholetissue_FDR")])
write.csv(panel, file.path(out_dir, "tables", "qPCR_panel.csv"), row.names = FALSE)
message("\nSuggested panel:")
print(panel[, c("role", "gene", "note")], row.names = FALSE)

## ---- reference-gene pre-screen (confirm with geNorm/NormFinder on tissue) -
ref <- intersect(ref_genes, rownames(lc_wt))
rs <- data.frame(gene = ref,
                 mean_logCPM = rowMeans(lc_wt[ref, ]),
                 sd_wholetissue = apply(lc_wt[ref, ], 1, sd),
                 HPV_log2FC = wt$logFC[match(ref, wt$gene)],
                 tumour_vs_normal_log2FC = A1$logFC[match(ref, A1$gene)])
rs$stability_rank <- rank(rank(rs$sd_wholetissue) + rank(abs(rs$HPV_log2FC)) +
                          rank(abs(rs$tumour_vs_normal_log2FC)), ties.method = "first")
rs <- rs[order(rs$stability_rank), ]
write.csv(rs, file.path(out_dir, "tables", "reference_gene_screen.csv"), row.names = FALSE)
message("\nMost stable reference genes: ", paste(head(rs$gene, 5), collapse = ", "))
print(rs, row.names = FALSE, digits = 2)

message("\nPlanned TME markers in whole tissue (HPV+ vs HPV-):")
print(wt[wt$gene %in% markers, c("gene", "logFC", "FDR")], row.names = FALSE, digits = 2)

saveRDS(list(scr = scr, panel = panel, rs = rs, lc_wt = lc_wt, wt_meta = wt_meta, wt = wt,
             cpm_g = cpm_g, smp = smp, detected = pbc$detected, n_group = pbc$n_group,
             novelty = novelty, markers = markers,
             thresholds = list(min_cells = min_cells, min_normal = min_normal,
                               min_cpm = min_cpm, min_share = min_share, min_spec = min_spec)),
        file.path(in_dir, "qpcr_screen.rds"))
