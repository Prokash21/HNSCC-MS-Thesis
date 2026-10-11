## Project paths, kept in one place. Source from the repo root:
##   source("config/paths.R")
## Date: 2026-10-11

paths <- list(
  ## ---- data (git-ignored) ----------------------------------------------
  atlas              = "data/atlas/HNSCC_Atlas.rds",
  raw_dir            = "data/raw",
  raw_GSE103322      = "data/raw/GSE103322/GSE103322_HNSCC_all_data.txt.gz",
  external_dir       = "data/external",
  reference_dir      = "data/reference",
  collectri          = "data/reference/collectri/collectri_omnipath.tsv",
  ensembl_dir        = "data/reference/ensembl",

  ## processed objects from v1, reused in v2
  processed_atlas    = "data/processed/HNSCatlas",
  atlas_cells_meta   = "data/processed/HNSCatlas/atlas_cells_meta_umap.rds",
  atlas_epithelial   = "data/processed/HNSCatlas/atlas_epithelial.rds",
  epithelial_umap    = "data/processed/HNSCatlas/epithelial_umap.rds",
  nonmalignant_subtypes = "data/processed/HNSCatlas/nonmalignant_subtypes.rds",
  normal_reference   = "data/processed/HNSCatlas/normal_reference_cells.rds",
  cnv_input          = "data/processed/HNSCatlas/cnv_input.rds",
  pseudobulk_v1      = "data/processed/HNSCatlas/pseudobulk.rds",
  celltype_pseudobulk = "data/processed/HNSCatlas/celltype_pseudobulk.rds",
  processed_GSE103322 = "data/processed/GSE103322",

  ## ---- v1 results (archived, tracked) ----------------------------------
  v1_results         = "archive_v1_HPVpos/results",
  patient_metadata   = "archive_v1_HPVpos/results/01_atlas/tables/patient_metadata.csv",
  cnv_scores_v1      = "archive_v1_HPVpos/results/02_normal_cleaning/tables/cnv_scores_per_cell.csv",
  dge_v1             = "archive_v1_HPVpos/results/03_DGE/tables",

  ## ---- v2 results -------------------------------------------------------
  results            = "results",
  step01             = "results/step01_metadata",
  step01_raw         = "data/raw/step01_metadata"     # downloaded GEO/TCGA metadata
)
