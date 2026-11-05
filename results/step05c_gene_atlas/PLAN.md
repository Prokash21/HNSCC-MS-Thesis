# Step 05c plan: gene atlas heatmap, blood background check, correlation plots

Written 2026-10-28, before any figure was made and before the blood data were processed. This step is **descriptive only**: no primer design, no gene list or selection change, and no new outcome models. Later changes are listed as post hoc in `STEP05c_SUMMARY.md`.

## 0. Genes and groups (24)

Group assignment is mine; the user gave the six group names but not a per-gene mapping.

| Group | Genes |
|---|---|
| identity (tumour-cell markers, up) | LAMC2, HOXD11, HOXB7, CASC9, IGF2BP3, LY6K, SLC2A1, SLC39A4, ARPC1A |
| network | E2F1, MYBL2, FOXM1, HMGB3 |
| outcome (step 07a eligible) | KIF23, HMMR, P4HA1 |
| down marker | KRT13, NFIC |
| optional | CA9 |
| reference (QC: log2FC should be near 0) | PSMC5, TBP, GUSB, SDHA, ALAS1 |

- **Target genes:** the 19 non-reference genes.
- **Co-expression network block:** E2F1, MYBL2, FOXM1, KIF23, HMGB3, HMMR, as specified by the user.

## 1. Data already in the repository (reused, not re-derived)

- **log2FC and FDR:**
  - step 04b `D_C2_HPVneg_allpatients`, `D_C1_HPVneg_paired`, `E_TCGA_HPVneg_paired_tumour_vs_normal`, `D_C2_HPVpos_allpatients`.
  - **TCGA oral cavity and larynx:** step 06 computed these but kept only the panel genes' log2FC. They are recomputed with the step 06 code exactly: HPV- pairs, cBioPortal site, edgeR QL `~ patient + state`, oral 28 pairs and larynx 11 pairs. They are checked against `D_panel_final.csv` (`oral_log2FC`, `larynx_log2FC`), where the values must agree.
- **Single-cell expression:** the step 04b part F method, as used in step 06 (`specificity_of`).
  - malignant and normal epithelium: median CPM of the step 04b HPV- pseudobulks (malignant ≥ 20 cells; normal ≥ 10 cells, P9 excluded);
  - 12 TME types: median CPM over HPV- patients' cell-type pseudobulks with ≥ 20 cells (`celltype_pseudobulk_v2.rds`);
  - tumour-specificity class: the same rule (tumour_specific / tumour_enriched / shared_with_TME).
- **Bulk tissue:**
  - TCGA HPV- tumour median TPM: step 07a, 415 tumours;
  - TCGA normal median TPM: the 40 HPV- paired normals (step 06 TPM file).
- **Outcome:** log HR per SD, as tested in earlier steps.
  - TCGA PFI all-site: step 07a `S_gene_cox`, multivariable;
  - TCGA oral PFI: step 07c `O_gene_cox`, main model;
  - GSE41613 DSS: step 07b `V_gene_cox`, multivariable.

  Genes not tested in a step are left empty (grey): E2F1 in 07b and 07c; CASC9 in 07b; the reference genes everywhere. No new models are fitted.
- **Network membership:** step 05b networks (pySCENIC consensus regulons plus HGNC-mapped CollecTRI, as in `05b_network.R`). Classes:
  - **TF:** E2F1, MYBL2, FOXM1;
  - **target:** a target of any of the three in either source;
  - **none.**

  Checked against `E_gene_target_membership.csv` for the overlapping genes.

## 2. New external data (A)

1. **Cillo 2020 (GSE139324), CD45+ scRNA-seq.**
   - `GSE139324_RAW.tar`, 543 MB; series matrix for metadata. Only the **PBMC** samples are used: HNSCC HPV- and HPV+ (`hpv_status`) and healthy donors (HD). TIL and tonsil samples are not used.
   - **Annotation:** GEO provides no author cell annotation (per-sample 10x matrices only). Major immune types are labelled with standard markers, and this is stated in the summary.
   - **QC:** cells with ≥ 200 detected genes and < 15% mitochondrial counts. Doublets are not removed, because only pseudobulks are used.
   - **Labelling:** Seurat v5, all PBMC cells together (no integration):
     - LogNormalize, 2,000 variable genes, 30 PCs, Louvain clusters (resolution 0.8);
     - each cluster gets the label of the marker set with the highest mean scaled expression across clusters.
   - **Marker sets** (standard PBMC markers, e.g. Hao 2021):

     | Type | Markers |
     |---|---|
     | CD4 T | CD3D, CD3E, IL7R, CD4 |
     | CD8 T | CD3D, CD8A, CD8B |
     | NK | GNLY, NKG7, KLRD1, KLRF1 |
     | B | MS4A1, CD79A, CD79B |
     | plasma | MZB1, JCHAIN, TNFRSF17 |
     | CD14 monocyte | CD14, LYZ, S100A8, S100A9 |
     | CD16 monocyte | FCGR3A, MS4A7 |
     | cDC | FCER1A, CD1C, CLEC10A |
     | pDC | LILRA4, IL3RA, CLEC4C |
     | platelet | PPBP, PF4 |

     The cluster × marker table is saved for review.
   - **Pseudobulk:** raw counts summed per cell type × group (HPV-, HPV+, HD), then CPM for the 24 genes, with cell numbers reported. Gene IDs are mapped from Ensembl IDs to HGNC; genes absent from the reference are reported.
2. **GTEx v8 gene median TPM** (`GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz`): Whole Blood, Esophagus - Mucosa and Minor Salivary Gland. Ensembl IDs are mapped to HGNC.
3. **"blood-low" flag (practical, pre-specified, not literature-based):**
   - **Rule:** max PBMC cell-type CPM < 1 **and** GTEx whole blood TPM < 1. The maximum is over all cell types with ≥ 50 cells in a group and over the three groups.
   - **Scope:** CD45+ PBMC data show the **background** expression of blood immune cells only. They cannot show whether tumour RNA or circulating tumour cells are detectable in blood.

## 3. Gene atlas heatmap (B; ComplexHeatmap)

- **Rows:** 24 genes, split by group.
- **Column blocks**, each with its own colour scale and legend:
  1. **log2FC** (diverging): single-cell C2, C1, TCGA paired all-site, TCGA oral, TCGA larynx, HPV+ C2. A dot marks FDR < 0.05; empty = not tested.
  2. **Single-cell log2(CPM + 1):** malignant, normal epithelium, CD4 T, CD8 T, NK, B, plasma, macrophage, monocyte, DC, neutrophil, mast, fibroblast, endothelial.
  3. **Bulk log2(TPM + 1):** TCGA HPV- tumour median, TCGA normal median, GTEx esophagus mucosa, GTEx minor salivary gland.
  4. **Blood log2(CPM or TPM + 1):** Cillo PBMC HPV- patients (max over cell types), healthy donors (max over cell types), GTEx whole blood.
  5. **Outcome log HR** (diverging): TCGA PFI all-site, TCGA oral PFI, GSE41613 DSS. A dot marks p < 0.05.
- **Right annotations:** specificity class, network membership, blood-low flag.
- **Matrix file:** the full matrix is saved as `tables/B_gene_atlas_matrix.csv`.

## 4. Correlation plots (C)

1. **LFC concordance:**
   - C2 log2FC vs TCGA paired log2FC over all genes tested in both (grey), with the 24 genes labelled; Pearson and Spearman r.
   - The same for C1 vs TCGA.
2. **Gene–gene co-expression**, Spearman rho among the 19 target genes, in:
   - TCGA HPV- tumours (step 07a TMM log CPM, n = 415);
   - the TCGA oral subgroup (n = 277; site as in step 07a);
   - single-cell HPV- malignant pseudobulks (n = 27; `pb$full` counts, TMM, log CPM prior 2).

   Display and reporting:
   - heatmaps clustered by hierarchical clustering on 1 − rho (average linkage);
   - cells with nominal p < 0.05 are marked, and BH FDR is reported in the table;
   - the network block is outlined.

   n = 27 is small and is stated.

## 5. Outputs

- **Scripts:** `scripts/05c_blood.R` (Cillo + GTEx), `scripts/05c_atlas.R` (matrix, site DE, annotations), `scripts/05c_figures.R`.
- **Tables:** `results/step05c_gene_atlas/tables/`.
- **`STEP05c_SUMMARY.md`** (local):
  - plain-language reading of the heatmap;
  - blood-low genes;
  - network co-expression;
  - limitations and UNCLEAR items;
  - references.
- **Figures:**
  - Fig47 gene atlas heatmap;
  - Fig48 blood panel (Cillo per cell type × group, and GTEx);
  - Fig49 LFC concordance scatters;
  - Fig50 co-expression heatmaps.
