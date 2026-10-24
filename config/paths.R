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
  step01_raw         = "data/raw/step01_metadata",    # downloaded GEO/TCGA metadata
  patient_clinical   = "results/step01_metadata/patient_clinical.csv",

  ## ---- step 02: raw downloads and assembled v2 objects -------------------
  choi_counts        = "data/raw/GSE181919/GSE181919_UMI_counts.txt.gz",
  choi_barcodes      = "data/raw/step01_metadata/GSE181919_Barcode_metadata.txt.gz",
  bill_counts_tar    = "data/raw/GSE234933/GSE234933_MGH_HNSCC_gex_raw_counts.tar.gz",
  bill_dir           = "data/raw/GSE234933",
  bill_cells         = "data/raw/step01_metadata/GSE234933_MGH_HNSCC_cell_annotation.txt.gz",
  bill_samples       = "data/raw/step01_metadata/GSE234933_MGH_HNSCC_sample_annotation.txt.gz",
  peng_dir           = "data/raw/GSE172577",
  processed_v2       = "data/processed/v2",
  atlas_epi_v2       = "data/processed/v2/atlas_epithelial_counts.rds",
  choi_nl_lp_all     = "data/processed/v2/choi_nl_lp_allcells.rds",
  choi_epi_v2        = "data/processed/v2/choi_nl_lp_epithelial.rds",
  bill_epi_v2        = "data/processed/v2/bill_recurrence_epithelial.rds",
  peng_all_v2        = "data/processed/v2/peng_allcells_qc.rds",
  peng_epi_v2        = "data/processed/v2/peng_epithelial.rds",
  epithelial_v2      = "data/processed/v2/epithelial_v2.rds",
  epithelial_v2_meta = "data/processed/v2/epithelial_v2_meta.csv.gz",
  step02             = "results/step02_assembly",

  ## ---- step 03: CNV calls and final labels -------------------------------
  cnv_cache          = "data/processed/v2/cnv",
  atlas_labelled     = "data/processed/v2/atlas_epithelial_labelled.rds",
  choi_labelled      = "data/processed/v2/choi_nl_lp_epithelial_labelled.rds",
  bill_labelled      = "data/processed/v2/bill_recurrence_epithelial_labelled.rds",
  peng_labelled      = "data/processed/v2/peng_epithelial_labelled.rds",
  cell_labels_v2     = "data/processed/v2/cell_labels_v2.csv.gz",
  step03             = "results/step03_cnv",

  ## ---- step 04: relabelling, pseudobulk, differential expression ---------
  ref_quantiles      = "data/processed/v2/cnv/ref_signal_quantiles.csv",
  pseudobulk_v2      = "data/processed/v2/pseudobulk_step04.rds",
  step04             = "results/step04_DE"
)
