## Extract the epithelial compartment (malignant "Tumor" + non-malignant
## "Epithelial" cells) from the HNSCatlas Seurat object, for tumour-vs-normal
## comparison / DGE and downstream substate discovery.
##
## Atlas facts this relies on (checked on the object, not the README):
##   - metadata: Patient, Source, Sex, HPV, Original_Cell_Type, Cell_Type
##     (no "Dataset" column)
##   - Source has NO normal-tissue samples (CA, CD45n, CD45p, Primary,
##     Unknown primary), so the normal comparator is the non-malignant
##     "Epithelial" cells from within the tumours
##   - assays: RNA (genes) + integrated (1000 features, default); pca/umap
##     were computed on the integrated assay
##
## Also writes a light per-cell table (metadata + atlas UMAP coordinates) so
## figures can be redrawn without reloading the 4 GB atlas.
##
## Run from the repo root:  Rscript scripts/02_atlas_extract_epithelial.R

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

atlas_path <- "HNSCC-HPV/atlas/HNSCC_Atlas.rds"
out_dir    <- "data/processed/HNSCatlas"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

t0 <- Sys.time()
message("Reading atlas (xz-compressed, this is slow)...")
atlas <- readRDS(atlas_path)
message("Loaded in ", format(round(Sys.time() - t0, 1)))
md <- atlas@meta.data

## one HPV status per patient. OP16 carries a single stray HPV+ cell among its
## 3,329 HPV- cells (atlas metadata glitch); majority vote per patient restores
## the paper's cohort (46 HPV-: 32 M / 14 F; 32 HPV+: 30 M / 2 F)
hpv_major <- tapply(md$HPV, md$Patient, function(x) names(which.max(table(x))))
fix <- md$HPV != hpv_major[md$Patient]
message(sprintf("HPV harmonised for %d cell(s) in: %s", sum(fix),
                paste(unique(md$Patient[fix]), collapse = ", ")))
md$HPV <- unname(hpv_major[md$Patient])
atlas$HPV <- md$HPV
message("\n== Patients: HPV x Sex ==")
print(table(unique(md[, c("Patient", "HPV", "Sex")])[, c("HPV", "Sex")]))

message("\n== seq_folder ==")
print(table(md$seq_folder, useNA = "ifany"))
message("\n== Patients per Source x HPV ==")
print(tapply(md$Patient, list(md$Source, md$HPV), function(p) length(unique(p))))

## epithelial lineage. Tumor/Epithelial calls inside CD45+ sorted fractions
## (immune-only by design) are most likely doublets -> excluded.
epi_md <- md[md$Cell_Type %in% c("Tumor", "Epithelial") & md$Source != "CD45p", ]
message(sprintf("\nEpithelial lineage kept: %d cells (dropped %d Tumor/Epithelial cells from CD45p sorts)",
                nrow(epi_md), sum(md$Cell_Type %in% c("Tumor", "Epithelial") & md$Source == "CD45p")))
message("\n== Epithelial lineage: Cell_Type x Source x HPV ==")
print(ftable(table(epi_md$Cell_Type, epi_md$Source, epi_md$HPV)))

## per-patient availability -- decides what a tumour-vs-normal DGE can be
per_patient <- do.call(rbind, lapply(split(epi_md, epi_md$Patient), function(d) data.frame(
  Patient = d$Patient[1], HPV = d$HPV[1], Sex = d$Sex[1],
  Source = paste(sort(unique(d$Source)), collapse = "+"),
  n_Tumor = sum(d$Cell_Type == "Tumor"), n_Epithelial = sum(d$Cell_Type == "Epithelial"))))
per_patient <- per_patient[order(per_patient$HPV, -per_patient$n_Tumor), ]
message("\n== Epithelial lineage per patient ==")
print(per_patient, row.names = FALSE)
message(sprintf("\nPatients with >=20 Tumor AND >=20 Epithelial cells (paired pseudobulk): %d (HPV+ %d, HPV- %d)",
                sum(per_patient$n_Tumor >= 20 & per_patient$n_Epithelial >= 20),
                sum(per_patient$n_Tumor >= 20 & per_patient$n_Epithelial >= 20 & per_patient$HPV == "HPV+"),
                sum(per_patient$n_Tumor >= 20 & per_patient$n_Epithelial >= 20 & per_patient$HPV == "HPV-")))
write.csv(per_patient, file.path(out_dir, "epithelial_cells_per_patient.csv"), row.names = FALSE)

## scVI / pseudobulk DGE need raw integer counts -- check them per Source
lay <- Layers(atlas[["RNA"]])
message("\nRNA layers: ", paste(lay, collapse = ", "))
for (l in grep("^counts", lay, value = TRUE)) {
  m <- LayerData(atlas, assay = "RNA", layer = l)
  if (!inherits(m, "dgCMatrix")) { message(l, ": class ", class(m)[1]); next }
  message(sprintf("%s: %d genes x %d cells", l, nrow(m), ncol(m)))
  if (ncol(m) == nrow(md)) {
    for (s in sort(unique(md$Source))) {
      idx <- which(md$Source == s)
      x <- m[, idx[seq_len(min(3000, length(idx)))]]@x
      message(sprintf("  %-38s non-integer fraction %.4f | max %.1f", s, mean(x %% 1 != 0), max(x)))
    }
  }
}

## light per-cell table for plotting
emb <- Embeddings(atlas, "umap")
cells <- data.frame(cell = rownames(md), md,
                    UMAP_1 = emb[rownames(md), 1], UMAP_2 = emb[rownames(md), 2],
                    check.names = FALSE)
saveRDS(cells, file.path(out_dir, "atlas_cells_meta_umap.rds"))

## epithelial compartment object (keeps both assays, reductions and metadata)
epi <- subset(atlas, cells = rownames(epi_md))
rm(atlas); invisible(gc())
epi$Compartment <- ifelse(epi$Cell_Type == "Tumor", "Tumor", "Non-malignant epithelial")
print(epi)
saveRDS(epi, file.path(out_dir, "atlas_epithelial.rds"))

message("\nDone in ", format(round(Sys.time() - t0, 1)), ". Outputs in ", out_dir, ":")
message(" - atlas_cells_meta_umap.rds        (all cells: metadata + atlas UMAP)")
message(" - atlas_epithelial.rds             (Seurat object: Tumor + non-malignant Epithelial)")
message(" - epithelial_cells_per_patient.csv (cells per patient, HPV, Source)")
