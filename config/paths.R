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
  step04             = "results/step04_DE",

  ## ---- step 04b: normal reference fix, symbols, bulk replication ---------
  hgnc               = "data/reference/hgnc/hgnc_complete_set.txt",
  cbio_pancan        = "data/raw/step01_metadata/cbioportal_hnsc_tcga_pan_can_atlas_2018_patient.csv",
  celltype_pb_v2     = "data/processed/v2/celltype_pseudobulk_v2.rds",
  choi_reannotated   = "data/processed/v2/choi_nl_lp_reannotated.csv.gz",
  pseudobulk_v2b     = "data/processed/v2/pseudobulk_step04b.rds",
  tcga_dir           = "data/raw/tcga_hnsc_star",
  tcga_counts        = "data/processed/v2/tcga_hnsc_paired_counts.rds",
  step04b            = "results/step04b_DE",

  ## ---- step 05: TF activity and regulatory networks ----------------------
  lambert_tfs        = "data/reference/tf_lists/Lambert2018_TF_names_v_1.01.txt",
  cistarget_dir      = "data/reference/cistarget",
  step05_export      = "data/processed/v2/step05_export",
  step05_scenic      = "data/processed/v2/step05_scenic",
  step05             = "results/step05_TF",

  ## ---- step 06: qRT-PCR panel and primers ---------------------------------
  tcga_site          = "data/raw/step01_metadata/cbioportal_hnsc_tcga_primary_site.csv",
  tcga_tpm           = "data/processed/v2/tcga_hnsc_paired_tpm.rds",
  step06             = "results/step06_qPCR",

  ## ---- step 07a: recurrence arm and survival -----------------------------
  papers_dir         = "data/raw/step01_metadata/papers",
  gene_sets_dir      = "config/gene_sets",
  tcga_all_counts    = "data/raw/step07_bulk/TCGA-HNSC.star_counts.tsv.gz",
  tcga_cdr           = "data/raw/tcga_cdr/TCGA-CDR-SupplementalTableS1.xlsx",
  gpl570             = "data/raw/step07_bulk/GPL570.annot.gz",
  step07a_tcga       = "data/processed/v2/step07a_tcga_hpvneg.rds",
  step07a_bill       = "data/processed/v2/step07a_bill_cells.rds",
  step07a            = "results/step07a_recurrence",

  ## ---- step 07b: independent validation (GEO cohorts) --------------------
  gse65858           = "data/raw/step07_bulk/GSE65858_series_matrix.txt.gz",
  gse41613           = "data/raw/step07_bulk/GSE41613_series_matrix.txt.gz",
  gpl10558           = "data/raw/step07_bulk/GPL10558.annot.gz",
  step07b_geo        = "data/processed/v2/step07b_geo.rds",
  step07b            = "results/step07b_validation",

  ## ---- step 07c: TCGA oral cavity subgroup (exploratory) -----------------
  step07c            = "results/step07c_oral",

  ## ---- step 05b: functional analysis ------------------------------------
  msigdb_h           = "data/reference/msigdb/h.all.v2024.1.Hs.symbols.gmt",
  msigdb_gobp        = "data/reference/msigdb/c5.go.bp.v2024.1.Hs.symbols.gmt",
  msigdb_gobp_json   = "data/reference/msigdb/c5.go.bp.v2024.1.Hs.json",
  msigdb_reactome    = "data/reference/msigdb/c2.cp.reactome.v2024.1.Hs.symbols.gmt",
  step05b_sets       = "data/processed/v2/step05b_gene_sets.rds",
  step05b            = "results/step05b_functional",

  ## ---- step 05c: gene atlas, blood background -----------------------------
  gse139324_pbmc     = "data/raw/step05c_blood/pbmc",
  gse139324_matrix   = "data/raw/step05c_blood/GSE139324_series_matrix.txt.gz",
  gtex_median        = "data/raw/step05c_blood/GTEx_v8_gene_median_tpm.gct.gz",
  step05c_pbmc       = "data/processed/v2/step05c_cillo_pbmc_pseudobulk.rds",
  step05c            = "results/step05c_gene_atlas",

  ## ---- step 08: final qPCR list, primers, ordering, wet-lab sheets -------
  step08             = "results/step08_primers"
)
