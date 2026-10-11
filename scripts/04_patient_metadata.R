## Per-patient metadata for the HNSCatlas: study, sample type, HPV, sex and how
## many tumour / non-malignant epithelial cells each patient contributes.
##
## Input:  data/processed/HNSCatlas/atlas_cells_meta_umap.rds (scripts/02)
## Output: HNSCC-HPV/results/01_atlas/tables/patient_metadata.csv
## Run from the repo root:  Rscript scripts/04_patient_metadata.R

cells   <- readRDS("data/processed/HNSCatlas/atlas_cells_meta_umap.rds")
out_dir <- "HNSCC-HPV/results/01_atlas/tables"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## the atlas has no Dataset column; patient-ID prefixes identify the study
## (HN = CD45-/CD45+ sorted fractions, i.e. GSE164690; P = "CA" samples)
study <- c(BHN = "Bill 2023 (GSE234933)", OP = "Puram 2023 (GSE182227)",
           HN = "Kurten 2021 (GSE164690)", P = "Choi 2023 (GSE181919)")
prefix <- sub("^(BHN|OP|HN|P).*$", "\\1", unique(cells$Patient))
stopifnot(all(prefix %in% names(study)))

per <- do.call(rbind, lapply(split(cells, cells$Patient), function(d) {
  ## tumour/epithelial calls inside CD45+ sorts are excluded, as in scripts/02
  e <- d[d$Source != "CD45p", ]
  data.frame(
    Patient     = d$Patient[1],
    Study       = study[[sub("^(BHN|OP|HN|P).*$", "\\1", d$Patient[1])]],
    Sample_type = paste(sort(unique(d$Source)), collapse = " + "),
    HPV         = d$HPV[1],
    Sex         = d$Sex[1],
    Total_cells = nrow(d),
    Tumor_cells = sum(e$Cell_Type == "Tumor", na.rm = TRUE),
    NonMalignant_epithelial_cells = sum(e$Cell_Type == "Epithelial", na.rm = TRUE))
}))
## usable for a paired (within-patient) tumour vs non-malignant comparison
per$Paired_20plus <- per$Tumor_cells >= 20 & per$NonMalignant_epithelial_cells >= 20
per <- per[order(per$HPV, per$Study, -per$Tumor_cells), ]
write.csv(per, file.path(out_dir, "patient_metadata.csv"), row.names = FALSE)

message("Patients: ", nrow(per))
message("\n== Patients: HPV x Sample type ==")
print(table(per$Sample_type, per$HPV))
message("\n== Patients: HPV x Study ==")
print(table(per$Study, per$HPV))
message("\n== Patients with >=20 tumour AND >=20 non-malignant epithelial cells ==")
print(table(Paired = per$Paired_20plus, HPV = per$HPV))
message("\n== Cells (CD45+ sorts excluded): HPV x compartment ==")
print(aggregate(cbind(Tumor_cells, NonMalignant_epithelial_cells) ~ HPV, per, sum))
message("\nWrote ", file.path(out_dir, "patient_metadata.csv"))
