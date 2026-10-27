# PROJECT STATUS: HNSCC single-cell analysis (HNSCatlas)

Status as of 2026-10-11. Everything below was read from files in this repo or from lightweight inspection of the saved objects. No analysis was re-run. **UNCLEAR** marks anything that could not be verified.

**Important framing note.** The existing analyses were designed with **HPV+ as the primary focus**: the candidate genes, qPCR panel and primers are all HPV+ genes. The new plan makes **HPV-negative tumour cells primary**. For that plan, nearly every downstream step (candidates, qPCR screen, primers, TF interpretation) has **not been done yet**. Only the HPV- tumour-vs-normal contrast (A3) exists, and it uses just 5 patients.

---

## 1. FOLDER STRUCTURE

Project root: `/Users/prokash/DWCT/Thesis-Cancer/Thesis` (about 6.1 GB in total; `data/` is about 2.2 GB). `data/`, `*.rds` and the atlas are git-ignored.

```text
.
├── README.md, LICENSE, .gitignore, PROJECT_STATUS.md
├── PLAN/
│   ├── HNSCC_Research_Plan.docx          (new HPV-negative plan, v1.0, 11 Oct 2026)
│   └── Sample_Collection_Protocol.docx   (untracked, not in git)
├── scripts/                              (01–12, see section 4)
├── HNSCC-HPV/
│   ├── atlas/
│   │   ├── HNSCC_Atlas.rds               4,034,095,316 B (4.0 GB)  ← source atlas
│   │   ├── HNSCatlas.R                   (exploratory: DimPlot + tables)
│   │   └── .Rhistory, .RDataTmp
│   ├── docs/                             (pipeline .docx/.pdf, strategy-note PDF, 2 Elsevier PDFs, umap.png)
│   └── results/
│       ├── 01_atlas/                     Fig1, Fig2, tables/patient_metadata.csv (5.8 KB)
│       ├── 02_normal_cleaning/           Fig3, Fig3b, tables/{cnv_scores_per_cell.csv 5.0 MB, cnv_summary.csv 180 B, nonmalignant_clusters.csv 3.5 KB}
│       ├── 03_DGE/                       Fig4–6, tables/ (A1–A4, B ≈1.1 MB each; candidates_HPVpos 18 KB; candidates_HPVneg 0.3 KB; pseudobulk_samples 10 KB)
│       │   └── first_pass_all_nonmalignant/  older run: Fig4–6 + tables/ (≈1.1–1.3 MB each)
│       ├── 04_qPCR_panel/                Fig7, Fig8, tables/{candidate_screen 27 KB, qPCR_panel 0.9 KB, reference_gene_screen 1.8 KB, wholetissue_HPVpos_vs_HPVneg 1.6 MB}
│       ├── 05_primers/                   Fig9, tables/{primer_candidates 22 KB, primer_map 6 KB, primer_order_sheet 2 KB}
│       └── 06_TF_activity/               Fig10, tables/{TF_activity_HPVpos_vs_HPVneg 67 KB, _all_contrasts 121 KB, _per_patient 762 KB}
└── data/
    ├── raw/GSE103322/GSE103322_HNSCC_all_data.txt.gz          90.2 MB
    ├── processed/GSE103322/
    │   ├── GSE103322_seurat.rds          256.6 MB
    │   ├── expr_log2tpm_matrix.rds       128.3 MB
    │   ├── expr_log2tpm_matrix.csv.gz     90.1 MB
    │   └── cell_metadata.csv               0.37 MB
    ├── processed/HNSCatlas/
    │   ├── atlas_epithelial.rds          1,163 MB   (Seurat: Tumor + non-malignant Epithelial)
    │   ├── cnv_input.rds                   207.8 MB (cache for 05b)
    │   ├── nonmalignant_clustered.rds       80.8 MB
    │   ├── nonmalignant_epithelial.rds      65.2 MB
    │   ├── celltype_pseudobulk.rds          46.0 MB (cache for 08)
    │   ├── atlas_cells_meta_umap.rds         7.1 MB (all 274,911 cells: metadata + UMAP)
    │   ├── qpcr_screen.rds                   6.9 MB
    │   ├── dge_results.rds                   4.3 MB
    │   ├── pseudobulk.rds                    3.6 MB
    │   ├── epithelial_umap.rds               1.0 MB
    │   ├── nonmalignant_subtypes.rds         0.30 MB
    │   ├── normal_reference_cells.rds        8 KB
    │   └── epithelial_cells_per_patient.csv  2.8 KB
    └── reference/
        ├── collectri/collectri_omnipath.tsv              4.1 MB
        └── ensembl/Homo_sapiens.GRCh38.{cdna 183.9 MB, ncrna 41.0 MB}.fa.gz
```

There are no `.h5ad`, `.parquet`, `.xlsx` or `.ipynb` files. All work is in R (Seurat), plus one Python script.

---

## 2. DATA OBJECTS

### 2.1 `HNSCC-HPV/atlas/HNSCC_Atlas.rds`: the HNSCatlas (Conde-Lopez et al., iScience 2026)
- **Format:** Seurat v5 object, xz-compressed RDS, 4.0 GB. The full object was **not reloaded** for this report. Structure was taken from its subset `atlas_epithelial.rds` (same assays and commands) and from the full-metadata export `atlas_cells_meta_umap.rds`.
- **Size:** **274,911 cells**, 78 patients. Genes: RNA assay **45,075** features. The integrated assay has **1,000** features.
- **Raw counts:** **yes**, in Seurat assay `RNA`, layer `counts`. Checked on the epithelial subset: 0 % non-integer values in every Source (CD45n, Primary, CA, Unknown primary); max 19,671. The `data` layer holds log-normalised values.
- **Integrated assay:** `integrated` (layers data, scale.data). It is the **default assay**, and `pca`/`umap` were computed on it. **Do not use it for DE.**
- **How the embedding was made** (from `@commands` stored in the object, dated 2024-04-07, i.e. by the atlas authors and not in this project):
  - `FindIntegrationAnchors`: LogNormalize, k.anchor 5, k.score 30, max.features 200, reduction CCA (Seurat default). There were **4 objects in the sample tree**, so integration was across 4 groups, presumably the 4 studies (batch key **UNCLEAR**, as no dataset column exists).
  - `IntegrateData`: dims 1:20, k.weight 100. This produced 1,000 integrated features.
  - `ScaleData` → `RunPCA` with npcs 20 → `RunUMAP` on pca, 30 neighbours, cosine metric, min.dist 0.3, seed 42.
  - **No scVI** embedding exists anywhere.
- **Metadata columns:** `orig.ident` (UNCLEAR whether present in the full object; it is absent in the exports), `seq_folder`, `nCount_RNA`, `nFeature_RNA`, `Patient`, `Source`, `Sex`, `HPV`, `Original_Cell_Type`, `Cell_Type`.
  - There is **no** study/dataset column. Study is derived from the patient-ID prefix in `scripts/04`: BHN = Bill 2023 GSE234933, OP = Puram 2023 GSE182227, HN = Kurten 2021 GSE164690, P = Choi 2023 GSE181919.
  - There is **no** anatomical-site, tissue, stage, treatment or recurrence column.
  - `nCount_RNA` ranges 500–241,583 (median 3,405). `nFeature_RNA` ranges 10–11,264 (median 1,255).

| Column | Values (cells) |
|---|---|
| `Source` (sample/sort type) | Primary 172,987 · CD45p 49,437 · CD45n 23,643 · CA 23,088 · Unknown primary of the head and neck 5,756 |
| `HPV` | HPV- 147,744 · HPV+ 127,167 (after per-patient majority harmonisation; see §3.5) |
| `Sex` | Male 226,784 · Female 48,127 |
| `Cell_Type` | CD4T 50,435 · **Tumor 48,718** · CD8T 43,467 · Fibroblast 19,957 · Bcell 16,747 · Endothelial 16,149 · Macrophage 10,640 · PlasmaCell 8,388 · NK 8,137 · Monocyte 8,128 · panDC 8,023 · **Epithelial 4,818** · Neutrophils 2,868 · Multi 2,417 · Mast 1,524 · **NA 24,495** |
| `Original_Cell_Type` | 46 unique study-specific labels, e.g. T_cells 46,416, Epithelial 23,770, Tumor 18,275, Malignant.cells 4,829, Unresolved 4,149, NormalEpith 418 … |
| `seq_folder` | "SeuratProject" 73,080 (= all CD45n + CD45p, i.e. Kurten); NA otherwise |

- **Malignant label:** `Cell_Type == "Tumor"`. The non-malignant epithelial label is `Cell_Type == "Epithelial"`. The atlas authors decided how cells were called malignant. This project did not re-call malignant cells. It only used a CNV check on the normal side (§5.2).

### 2.2 `data/processed/HNSCatlas/atlas_epithelial.rds` (made by scripts/02)
- Seurat object with **52,909 cells** × 45,075 RNA genes + 1,000 integrated features. Default assay is `integrated`. Reductions: `pca` (20 dims, integrated) and `umap` (atlas UMAP).
- RNA layers: `counts` (raw integers) and `data`. The integrated layers are `data` and `scale.data`.
- Metadata is the atlas columns plus `Compartment` (Tumor / Non-malignant epithelial). HPV is harmonised per patient.
- Selection: `Cell_Type ∈ {Tumor, Epithelial}` **and** `Source != "CD45p"`. This drops 624 Tumor + 3 Epithelial cells found in CD45+ sorts.

| | HPV- | HPV+ |
|---|---|---|
| Tumor | 19,629 | 28,465 |
| Non-malignant epithelial | 2,554 | 2,261 |

- HPV16 viral genes are present as rows: `HPVC1, HPV16-E1, -E2, -E5, -L2, -L1, -E6, -E7`.

### 2.3 Other HNSCatlas-derived objects
| File | Content |
|---|---|
| `atlas_cells_meta_umap.rds` | data.frame, 274,911 rows: all metadata + atlas UMAP_1/UMAP_2 |
| `epithelial_umap.rds` | sub-UMAP of the 52,909 epithelial cells (re-run on the integrated assay: ScaleData → PCA 30 → UMAP dims 1:30, seed 42) |
| `nonmalignant_epithelial.rds` | Seurat object of the 4,815 non-malignant epithelial cells, RNA only (counts + data) |
| `nonmalignant_clustered.rds` | the same, after re-clustering (scripts/05) |
| `nonmalignant_subtypes.rds` | per cell: cluster, subtype, squamous_clean, module scores, UMAP |
| `normal_reference_cells.rds` | character vector of **1,136** cell IDs, the final "normal squamous" reference (HPV- 377, HPV+ 759) |
| `cnv_input.rds` | raw counts + meta for the CNV check (cache) |
| `pseudobulk.rds` | raw-count pseudobulk, patient × {Tumor, Normal squamous}; version string "Tumor + normal squamous, CNV-positive cells removed"; integer counts verified |
| `celltype_pseudobulk.rds` | raw-count pseudobulk, patient × every cell type (for the qPCR screen) |
| `dge_results.rds` | list: A1–A4, B tables, yA/sA, yB/sB (edgeR DGELists + sample tables), cand_pos, cand_neg |
| `qpcr_screen.rds` | list: screen, panel, reference genes, whole-tissue DE, CPMs, thresholds |

### 2.4 `data/processed/GSE103322/` (Puram 2017; **not** part of the atlas)
- `GSE103322_seurat.rds`: Seurat object with **5,902 cells, 21 patients**. The `counts` **and** `data` layers both hold **log2(TPM/10+1)**. **There are no raw counts** because GEO supplies only TPM.
- Metadata: `cell_id, patient_id, lymph_node, classified_cancer, classified_normal, noncancer_type, cell_type`.
- `cell_type`: malignant 2,215 · fibroblast 1,440 · T cell 1,237 · unresolved 324 · endothelial 260 · B cell 138 · mast 120 · macrophage 98 · dendritic 51 · myocyte 19. `lymph_node`: 0 = 4,541, 1 = 1,361 cells.
- HPV status per patient is **not** in this object (UNCLEAR, would need the paper's supplementary table). This dataset has not been used downstream yet. The README calls it the independent validation set.

---

## 3. COHORT NUMBERS (from `results/01_atlas/tables/patient_metadata.csv` and the pseudobulk sample tables)

Definitions: "Tumor" means `Cell_Type == Tumor`. "NM epi" means `Cell_Type == Epithelial`. Both exclude CD45p-sorted cells. "Clean normal" means cells left in the final normal-squamous reference after cleaning and CNV removal (§5.2).

### 3.1 Patients per study × HPV
| Study (GEO) | HPV- | HPV+ | Total |
|---|---|---|---|
| Bill 2023 (GSE234933), prefix BHN | 16 | 8 | 24 |
| Choi 2023 (GSE181919), prefix P | 13 | 7 | 20 |
| Kurten 2021 (GSE164690), prefix HN | 12 | 6 | 18 |
| Puram 2023 (GSE182227), prefix OP | 5 | 11 | 16 |
| **Total** | **46** | **32** | **78** |

Sex: HPV- 32 M / 14 F; HPV+ 30 M / 2 F.

### 3.2 Sample type (Source) per patient × HPV
| Source | HPV- | HPV+ |
|---|---|---|
| Primary (Bill, Puram) | 20 | 19 |
| CA (Choi) | 13 | 7 |
| CD45n + CD45p (Kurten) | 9 | 6 |
| CD45p only (Kurten: HN02, HN03, HN04, which have no epithelial cells) | 3 | 0 |
| Unknown primary of the head and neck (BHN17) | 1 | 0 |

**Anatomical site (oral cavity / oropharynx / larynx …): UNCLEAR.** The atlas has no site column, and no other file in the repo records site. It would have to come from each study's GEO/supplementary clinical tables. This matters a lot for the HPV- arm, where oral cavity and larynx tumours differ.

### 3.3 HPV-negative patients (n = 46)
| Patient | Study | Source | Tumor | NM epi | Clean normal | ≥50 T | ≥20 T & ≥20 NM |
|---|---|---|---|---|---|---|---|
| BHN70 | Bill | Primary | 2,644 | 38 | 32 | ✓ | ✓ |
| BHN64 | Bill | Primary | 1,620 | 64 | 5 | ✓ | ✓ |
| BHN60 | Bill | Primary | 1,391 | 37 | 34 | ✓ | ✓ |
| BHN7 | Bill | Primary | 518 | 33 | 4 | ✓ | ✓ |
| BHN1 | Bill | Primary | 505 | 22 | 0 | ✓ | ✓ |
| BHN39 | Bill | Primary | 316 | 1 | 1 | ✓ | |
| BHN17 | Bill | Unknown primary | 236 | 18 | 0 | ✓ | |
| BHN46 | Bill | Primary | 177 | 8 | 8 | ✓ | |
| BHN67 | Bill | Primary | 120 | 0 | 0 | ✓ | |
| BHN75 | Bill | Primary | 105 | 2 | 0 | ✓ | |
| BHN72 | Bill | Primary | 98 | 2 | 0 | ✓ | |
| BHN31TS | Bill | Primary | 39 | 7 | 7 | | |
| BHN49 | Bill | Primary | 26 | 0 | 0 | | |
| BHN30 | Bill | Primary | 20 | 6 | 6 | | |
| BHN40-3-05 | Bill | Primary | 9 | 7 | 5 | | |
| BHN74 | Bill | Primary | 0 | 1 | 1 | | |
| P8 | Choi | CA | 731 | 146 | 0 | ✓ | ✓ |
| P4 | Choi | CA | 427 | 78 | 0 | ✓ | ✓ |
| P51 | Choi | CA | 274 | 3 | 0 | ✓ | |
| P21 | Choi | CA | 272 | 17 | 0 | ✓ | |
| P9 | Choi | CA | 178 | 2 | 0 | ✓ | |
| P30 | Choi | CA | 94 | 10 | 1 | ✓ | |
| P26 | Choi | CA | 83 | 4 | 0 | ✓ | |
| P38 | Choi | CA | 81 | 2 | 0 | ✓ | |
| P15 | Choi | CA | 59 | 0 | 0 | ✓ | |
| P6 | Choi | CA | 52 | 0 | 0 | ✓ | |
| P31 | Choi | CA | 45 | 1 | 1 | | |
| P7 | Choi | CA | 43 | 0 | 0 | | |
| P60 | Choi | CA | 27 | 0 | 0 | | |
| HN15 | Kurten | CD45n+CD45p | 478 | 0 | 0 | ✓ | |
| HN01 | Kurten | CD45n+CD45p | 445 | 1 | 0 | ✓ | |
| HN08 | Kurten | CD45n+CD45p | 173 | 0 | 0 | ✓ | |
| HN07 | Kurten | CD45n+CD45p | 135 | 1 | 1 | ✓ | |
| HN11 | Kurten | CD45n+CD45p | 67 | 1 | 1 | ✓ | |
| HN10 | Kurten | CD45n+CD45p | 65 | 0 | 0 | ✓ | |
| HN09 | Kurten | CD45n+CD45p | 48 | 0 | 0 | | |
| HN06 | Kurten | CD45n+CD45p | 18 | 0 | 0 | | |
| HN05 | Kurten | CD45n+CD45p | 7 | 0 | 0 | | |
| HN02 | Kurten | CD45p | 0 | 0 | 0 | | |
| HN03 | Kurten | CD45p | 0 | 0 | 0 | | |
| HN04 | Kurten | CD45p | 0 | 0 | 0 | | |
| OP8 | Puram | Primary | 5,548 | 94 | 9 | ✓ | ✓ |
| OP12 | Puram | Primary | 1,187 | 123 | 104 | ✓ | ✓ |
| OP10 | Puram | Primary | 649 | 61 | 34 | ✓ | ✓ |
| OP19 | Puram | Primary | 359 | 491 | 123 | ✓ | ✓ |
| OP16 | Puram | Primary | 260 | 1,273 | 0 | ✓ | ✓ |

**HPV- summary:**
- 46 patients; 19,629 tumour cells; 2,554 NM epithelial cells.
- **≥ 50 tumour cells: 32 patients** (Bill 11, Choi 10, Kurten 6, Puram 5). ≥ 20 tumour cells: 39. 0 tumour cells: 4 (HN02, HN03, HN04, BHN74).
- **≥ 20 tumour AND ≥ 20 NM epi: 12 patients** (BHN70, BHN64, BHN60, BHN7, BHN1, P8, P4, OP8, OP12, OP10, OP19, OP16).
- After normal-reference cleaning, only **5 HPV- patients** keep ≥ 20 tumour and ≥ 10 clean normal cells (the Plan A rule): BHN60, BHN70, OP10, OP12, OP19. OP16's 1,273 "NM epi" cells were all removed (glandular or likely malignant).

### 3.4 HPV-positive patients (n = 32)
| Patient | Study | Source | Tumor | NM epi | Clean normal | ≥50 T | ≥20 T & ≥20 NM |
|---|---|---|---|---|---|---|---|
| BHN63 | Bill | Primary | 5,523 | 622 | 90 | ✓ | ✓ |
| BHN68 | Bill | Primary | 2,072 | 5 | 1 | ✓ | |
| BHN76 | Bill | Primary | 1,295 | 68 | 66 | ✓ | ✓ |
| BHN77 | Bill | Primary | 823 | 176 | 99 | ✓ | ✓ |
| BHN58 | Bill | Primary | 322 | 89 | 79 | ✓ | ✓ |
| BHN59 | Bill | Primary | 263 | 12 | 3 | ✓ | |
| BHN50 | Bill | Primary | 129 | 19 | 5 | ✓ | |
| BHN52 | Bill | Primary | 24 | 13 | 12 | | |
| P43 | Choi | CA | 1,367 | 197 | 1 | ✓ | ✓ |
| P57 | Choi | CA | 619 | 12 | 0 | ✓ | |
| P46 | Choi | CA | 297 | 4 | 1 | ✓ | |
| P22 | Choi | CA | 107 | 5 | 0 | ✓ | |
| P59 | Choi | CA | 41 | 0 | 0 | | |
| P86 | Choi | CA | 26 | 0 | 0 | | |
| P84 | Choi | CA | 6 | 0 | 0 | | |
| HN13 | Kurten | CD45n+CD45p | 1,193 | 27 | 20 | ✓ | ✓ |
| HN12 | Kurten | CD45n+CD45p | 729 | 2 | 2 | ✓ | |
| HN17 | Kurten | CD45n+CD45p | 542 | 2 | 2 | ✓ | |
| HN16 | Kurten | CD45n+CD45p | 416 | 1 | 0 | ✓ | |
| HN14 | Kurten | CD45n+CD45p | 273 | 26 | 23 | ✓ | ✓ |
| HN18 | Kurten | CD45n+CD45p | 78 | 15 | 14 | ✓ | |
| OP33 | Puram | Primary | 3,163 | 207 | 91 | ✓ | ✓ |
| OP20 | Puram | Primary | 2,439 | 170 | 43 | ✓ | ✓ |
| OP9 | Puram | Primary | 1,925 | 38 | 14 | ✓ | ✓ |
| OP4 | Puram | Primary | 1,899 | 85 | 26 | ✓ | ✓ |
| OP34 | Puram | Primary | 1,579 | 90 | 22 | ✓ | ✓ |
| OP6 | Puram | Primary | 410 | 62 | 45 | ✓ | ✓ |
| OP5 | Puram | Primary | 301 | 44 | 1 | ✓ | ✓ |
| OP14 | Puram | Primary | 157 | 80 | 64 | ✓ | ✓ |
| OP17 | Puram | Primary | 152 | 23 | 11 | ✓ | ✓ |
| OP35 | Puram | Primary | 150 | 155 | 18 | ✓ | ✓ |
| OP13 | Puram | Primary | 145 | 12 | 6 | ✓ | |

**HPV+ summary:**
- 32 patients; 28,465 tumour cells; 2,261 NM epithelial cells.
- **≥ 50 tumour cells: 28** (Bill 7, Choi 4, Kurten 6, Puram 11). ≥ 20 tumour: 31. 0 tumour: 0.
- **≥ 20 tumour AND ≥ 20 NM epi: 17.**
- Plan A (≥ 10 clean normal) used **17** HPV+ patients: Bill 5, Kurten 3, Puram 9.

### 3.5 HPV status: missing or ambiguous
- No patient has missing HPV.
- **OP16** had 1 HPV+ cell among 3,329 HPV- cells (an atlas metadata glitch). scripts/02 set HPV per patient by majority vote, giving HPV-. The result matches the paper (46 HPV- / 32 HPV+).
- **HPV status cannot be checked from viral reads outside Puram 2023.** HPV16 transcripts are non-zero only in Puram tumours: all 11 HPV+ Puram patients have viral CPM > 10 (median 330 CPM) and all 5 HPV- Puram patients have 0. Every Bill, Choi and Kurten sample has **0** viral counts, whether labelled HPV+ or HPV-. The README states those studies did not map HPV reads. For 21 HPV+ and 41 HPV- patients, HPV status therefore rests only on the clinical metadata in the atlas (method UNCLEAR: p16 IHC versus HPV DNA/RNA).

---

## 4. SCRIPTS AND NOTEBOOKS (run order = numbering; all are run from the repo root)

| # | File (lang) | What it does | Key inputs → outputs | Key parameters |
|---|---|---|---|---|
| 0 | `HNSCC-HPV/atlas/HNSCatlas.R` (R) | Exploratory: atlas DimPlot by Cell_Type, tables of Source and cell types | atlas → screen only (`docs/umap.png` is probably from this; UNCLEAR) | none |
| 01 | `01_parse_GSE103322.R` (R) | Parses the Puram 2017 GEO text file into metadata, a log2TPM matrix and a Seurat object | raw txt.gz → `data/processed/GSE103322/*` | drops duplicate gene symbols; stores log2(TPM/10+1) in both counts and data |
| 02 | `02_atlas_extract_epithelial.R` (R) | Loads the atlas, harmonises HPV per patient, extracts Tumor + Epithelial (no CD45p), checks raw counts are integers, exports per-cell metadata/UMAP | atlas → `atlas_epithelial.rds`, `atlas_cells_meta_umap.rds`, `epithelial_cells_per_patient.csv` | HPV = majority per patient; exclude `Source == CD45p` |
| 03 | `03_atlas_figures.R` (R) | Fig1 (all cells), Fig2 (epithelial compartment, re-embedded) | meta/UMAP, atlas_epithelial → Fig1, Fig2, `epithelial_umap.rds` | sub-UMAP on the **integrated** assay, PCA 30, dims 1:30, seed 42 |
| 04 | `04_patient_metadata.R` (R) | One row per patient: study (from ID prefix), source, HPV, sex, cell counts | meta → `patient_metadata.csv` | Paired_20plus = ≥ 20 T & ≥ 20 NM epi |
| 05 | `05_nonmalignant_subtypes.R` (R) | Re-clusters the 4,815 NM epithelial cells; labels clusters squamous / glandular / ciliated / likely malignant; flags immune/stromal doublets | atlas_epithelial → `nonmalignant_*.rds`, `nonmalignant_clusters.csv`, Fig3 | RNA assay; NormalizeData; 2,000 HVGs; PCA 30; neighbours/UMAP dims 1:20; Louvain **res 0.5**; rules: glandular if ≥ 30 % express PIGR/MUC5B/BPIFB1; likely malignant if one patient ≥ 80 % of the cluster, or ≥ 25 % HPV16+ or ≥ 25 % cancer-testis-antigen+; squamous if score > 0.25; doublet if immune or stromal score > 0.5. **No batch correction** |
| 05b | `05b_cnv_check.R` (R) | inferCNV-style CNV scoring, implemented by hand, to remove CNV-positive cells from the normal reference | atlas counts, subtypes, Ensembl gene positions → `cnv_scores_per_cell.csv`, `cnv_summary.csv`, `normal_reference_cells.rds`, Fig3b | per patient; reference = immune + stromal cells (≥ 30 needed; ≤ 300 sampled); ≤ 300 tumour cells; top 7,000 expressed autosomal genes; log2(CPM/10+1); clip ±3; 101-gene window; tumour-like = signal > 99th pct of reference AND correlation > 0.4; patient informative if ≥ 50 % of its tumour cells are tumour-like |
| 06 | `06_dge_pseudobulk.R` (R) | Pseudobulk edgeR QL: tumour vs normal squamous (paired: A1–A4) and HPV+ vs HPV- tumour (B); candidate lists | atlas_epithelial (RNA counts), normal_ref, patient_metadata → `03_DGE/tables/*`, `pseudobulk.rds`, `dge_results.rds` | min 20 tumour cells; min 10 normal cells (Plan A); filterByExpr; TMM; robust QL; FDR < 0.05 & \|log2FC\| ≥ 1 |
| 07 | `07_dge_figures.R` (R) | Fig4 volcanoes, Fig5 UpSet/Venn, Fig6 checks | `dge_results.rds` → Fig4–6 | none |
| 08 | `08_qpcr_panel.R` (R) | Screens the HPV+ candidates for qPCR: tumour specificity versus immune/stroma, tumour share of the whole-tissue signal, expression level, whole-tissue HPV effect; picks a 4-gene panel; reference-gene pre-screen | atlas counts → `celltype_pseudobulk.rds`, `04_qPCR_panel/tables/*`, `qpcr_screen.rds` | min_cells 20; detectable = median CPM ≥ 10 in HPV+ tumour; tumour share ≥ 0.5; specificity log2 ≥ 1; Kurten excluded from whole-tissue analysis (sorted); whole-tissue edgeR `~ study + hpv` |
| 09 | `09_qpcr_figures.R` (R) | Fig7, Fig8 | `qpcr_screen.rds` → Fig7, Fig8 | none |
| 10 | `10_primer_design.py` (Python via `uv run --with primer3-py`) | Primer3 SYBR primers on the Ensembl canonical transcript (exon-junction spanning), checked against the whole transcriptome | Ensembl REST + cDNA/ncRNA FASTA → `05_primers/tables/*` | product 70–150 bp; Tm 58–62 (opt 60); GC 40–60 %; 18–24 nt; off-target check: 12-nt 3′ seed, ≤ 4 mismatches, product ≤ 2 kb |
| 11 | `11_primer_map.R` (R) | Fig9 primer map | primer tables → Fig9 | flank 120 bp |
| 12 | `12_tf_activity.R` (R) | TF activity with decoupleR ULM + CollecTRI at two levels: (1) on edgeR signed statistics per contrast; (2) per patient on tumour pseudobulk, then `lm(activity ~ study + hpv)` | `dge_results.rds`, CollecTRI → `06_TF_activity/tables/*`, Fig10 | minsize 5; robust = both levels FDR < 0.05 with the same sign |

The **script that produced `03_DGE/first_pass_all_nonmalignant/`** is not in the repo. It was an earlier version of 06 that used all NM epithelial cells as "normal" (UNCLEAR exact parameters).

---

## 5. ANALYSES COMPLETED

### 5.1 QC
- This project did **no cell-level QC**. It relies on the atlas authors' QC (nCount_RNA ≥ 500 in the data; other thresholds UNCLEAR).
- Raw counts were confirmed to be integers per Source (scripts/02 and this report).

### 5.2 Normal-reference cleaning (scripts/05 + 05b)
Sub-clustering of the 4,815 NM epithelial cells gave 21 clusters:

| Subtype | Cells |
|---|---|
| Squamous | 1,406 → 1,309 after removing doublets |
| Glandular | 1,680 |
| Likely malignant | 1,692 |
| Ciliated | 37 |

The CNV check ran on 74 patients, of which **48 were informative** (HPV- 23, HPV+ 25). Cells with a tumour-like CNV profile, by group:

| Group (informative patients) | Cells | Patients | % tumour-like |
|---|---|---|---|
| Reference (immune/stromal) | 13,721 | 48 | 0.1 % |
| Normal squamous | 1,098 | 38 | **15.8 %** |
| Glandular | 1,680 | 13 | 0.2 % |
| Likely malignant | 1,500 | 32 | 63.9 % |
| Tumour | 10,261 | 48 | 81.4 % |

**Final normal reference: 1,136 cells** (HPV- 377, HPV+ 759). The CNV-positive squamous cells were removed.

### 5.3 Pseudobulk DGE (scripts/06): edgeR quasi-likelihood on raw-count pseudobulk
Thresholds: FDR < 0.05 and |log2FC| ≥ 1. Genes tested: 10,627 (Plan A) and 10,857 (Plan B).

| Contrast | Design | Patients (HPV- / HPV+) | Up | Down |
|---|---|---|---|---|
| A1 tumour vs normal, all | `~ patient + tumor` | 22 (5 / 17) | 518 | 237 |
| A2 tumour vs normal, HPV+ | `~ patient + tumor:HPV+ + tumor:HPV-`, coef tumor_HPVpos | 17 HPV+ | 439 | 146 |
| **A3 tumour vs normal, HPV-** | same model, coef tumor_HPVneg | **5 HPV-** (BHN60, BHN70, OP10, OP12, OP19) | **31** | **29** |
| A4 (A2 − A3) | contrast | 22 | 1 | 0 |
| B HPV+ vs HPV- tumour | `~ study + hpv` (tumour pseudobulk only) | 70 (39 / 31) | 400 | 452 |

Plan B patients by study: Bill 14 / 8, Choi 13 / 6, Kurten 7 / 6, Puram 5 / 11 (HPV- / HPV+).

The **first pass** (normal = all NM epithelial cells; not the current analysis) gave A1 64 up / 287 down, A2 30 / 171, A3 16 / 40, A4 1 / 2, B unchanged.

**Top 20 by p-value per contrast (gene, log2FC, FDR)**

| A1 all | A2 HPV+ | A3 HPV- | A4 HPV+ vs HPV- change | B HPV+ vs HPV- tumour |
|---|---|---|---|---|
| NTS 3.54 1.0e-5 | HPV16-E1 6.15 1.4e-5 | PAX9 −3.45 7.6e-4 | PAX9 3.21 0.025 | SMC1B 5.34 3.3e-21 |
| HPV16-E1 6.11 1.0e-5 | HSD17B2 −4.31 1.4e-5 | SYT8 −4.90 2.0e-3 | DHRS3 2.44 0.054 | SYNGR3 4.24 1.3e-18 |
| HMSD 3.74 1.0e-5 | HPV16-E7 5.72 1.6e-5 | HOXA10 5.54 4.2e-3 | RASIP1 3.93 0.054 | ZFR2 6.39 6.9e-18 |
| HOXA10 4.07 1.3e-5 | NTS 3.45 3.3e-5 | MAGEA4 4.88 4.2e-3 | CDKN2A 5.26 0.061 | MIR9-3HG 3.43 4.9e-16 |
| HPV16-E7 5.68 1.5e-5 | HMSD 3.81 4.7e-5 | VAMP5 −2.76 7.5e-3 | SELENOP 5.57 0.061 | KRTDAP −6.32 5.2e-14 |
| FADS2 3.21 2.8e-5 | SMARCD3 −3.71 4.7e-5 | MAP1B 2.91 0.010 | CTSL −3.06 0.072 | SYCP2 4.32 3.9e-13 |
| RFC4 1.53 3.2e-5 | ZFR2 3.35 4.7e-5 | ENO2 3.89 0.021 | BCL2A1 4.63 0.072 | ACTN1 −2.26 5.1e-12 |
| RSRC1 1.50 3.2e-5 | MCM2 2.21 8.7e-5 | FOXA1 −3.09 0.021 | KRT15 3.17 0.072 | STAG3 2.70 3.1e-11 |
| ADAM23 3.34 3.7e-5 | RFC4 1.69 1.1e-4 | SH3YL1 −2.60 0.023 | ZBED2 −4.19 0.083 | MEI1 4.10 3.1e-11 |
| EPCAM 2.42 3.7e-5 | OVOS2 4.27 1.1e-4 | HOXC10 5.52 0.023 | SSBP2 2.16 0.099 | RPA2 1.40 1.4e-10 |
| ACTL6A 1.31 3.7e-5 | EFEMP2 −3.37 1.1e-4 | PRAME 4.08 0.023 | CDA −5.77 0.10 | CDKN2C 2.72 2.0e-10 |
| OVOS2 4.05 3.7e-5 | MLF1 2.01 1.2e-4 | DEPTOR −2.97 0.024 | HTRA1 −3.05 0.10 | ARHGAP4 2.77 3.9e-10 |
| RP11-357H14.17 3.77 3.9e-5 | HPV16-E6 6.28 1.2e-4 | SELENOP −4.89 0.024 | MYH10 −3.49 0.10 | MMP2 −4.58 4.5e-10 |
| ALG1L 2.30 5.2e-5 | RIBC2 3.86 1.2e-4 | DHRS3 −2.01 0.024 | FOXA1 2.93 0.10 | EREG −5.16 7.6e-10 |
| EFEMP2 −2.88 5.2e-5 | ALG1L 2.62 1.2e-4 | THSD4 −3.79 0.024 | SIM2 3.33 0.10 | FOXRED2 1.96 1.5e-9 |
| HOXB7 3.12 5.6e-5 | CDKN2A 2.39 2.0e-4 | PHYHD1 −5.35 0.024 | KRT75 −5.55 0.10 | KLHL35 3.18 1.9e-9 |
| PRAME 3.91 7.0e-5 | ABCC5 1.60 2.0e-4 | HAGLROS 3.48 0.024 | ASS1 2.73 0.12 | VCAM1 4.57 2.5e-9 |
| ASRGL1 3.10 8.3e-5 | ZBED2 −3.63 2.3e-4 | PNCK 5.10 0.024 | IGF2BP2 −3.71 0.12 | SIM2 3.40 2.7e-9 |
| GPC3 4.56 8.3e-5 | RSRC1 1.54 2.4e-4 | NDUFA4L2 3.76 0.024 | PARVB −2.95 0.13 | SMTN −1.99 2.7e-9 |
| ZIC2 3.56 9.7e-5 | NUSAP1 1.90 2.4e-4 | MGST1 −5.75 0.024 | IGFL2-AS1 −10.37 0.14 | MMP1 −5.98 2.7e-9 |

(In A1–A3 and B, positive log2FC = higher in tumour or HPV+. In A4, positive = larger tumour-vs-normal change in HPV+.)

**Candidate lists:**
- **HPV+ candidates** (up in A2 AND up in B): **95** genes, e.g. SMC1B, SYNGR3, ZFR2, MIR9-3HG, SYCP2, STAG3, KLHL35, RIBC2, FOXD3-AS1, CDKN2A, E2F1 …
- **HPV- candidates** (up in A3 AND down in B): **1** gene, **HOXD11** (A3 log2FC 4.41, FDR 0.030; B log2FC −2.91, FDR 4e-6).

### 5.4 Built-in QC checks
| Check | Result | Outcome |
|---|---|---|
| CDKN2A (p16) higher in HPV+ tumours (B) | log2FC 1.87, FDR 0.0069 | ✓ expected direction |
| CDKN2A tumour vs normal | A2 +2.39 (FDR 2e-4); A3 −2.87 (FDR 0.19, n.s.) | consistent with p16 loss in HPV- (not significant) |
| MAGEA4 up in tumour | A1 +4.69 (FDR 1.2e-4); A3 +4.88 (FDR 0.004); A2 +3.94 (FDR 0.11) | ✓ |
| HPV16-E6 / E7 up in HPV+ tumour (A2) | E6 +6.28 (FDR 1.2e-4); E7 +5.72 (FDR 1.6e-5); 0 in A3 | ✓, but this signal comes from Puram only |
| HPV16 genes in B | not tested (removed by filterByExpr; zero in 3 of 4 studies) | see §3.5 |
| CNV check of the normal reference | 15.8 % of "clean" squamous cells were tumour-like; removed | done for informative patients only (48 / 74) |

### 5.5 qRT-PCR candidate screen (scripts/08): **HPV+ only**
- Of the 95 HPV+ candidates: 78 detectable, 38 tumour-dominant, 62 HPV+-higher in whole tissue, **25 pass all three filters**.
- Whole-tissue pseudobulk HPV+ vs HPV- (`~ study + hpv`, unsorted studies only): 16,407 genes tested, 3,811 with FDR < 0.05.
- **Suggested panel:**

| Role | Gene | Note |
|---|---|---|
| Positive control | SYCP2 | B log2FC 4.32; whole-tissue log2FC 2.98 |
| Novelty | ZFR2 | B log2FC 6.39; whole-tissue log2FC 4.86 |
| Novelty (lncRNA) | MIR9-3HG | B log2FC 3.43; whole-tissue log2FC 2.09 |
| Novelty (highest risk) | KLHL35 | replaces RIBC2, which failed the expression filter (8.3 CPM < 10) |

- Novelty labels are copied from the strategy-note PDF and were not verified.
- **Reference genes** (most stable first): SDHA, RPS18, RPLP0, PSMC5, RPL27, PPIA, YWHAZ, GUSB. To be confirmed with geNorm/NormFinder.
- TME markers carried along: GZMK, CD8A, COL1A2, SFRP2.

### 5.6 Primer design (scripts/10–11)
- Best primer pair per gene for 13 genes: SYCP2, ZFR2, MIR9-3HG, KLHL35, GZMK, COL1A2, SFRP2, RPS18, SDHA, PSMC5, RPLP0, FOXD3-AS1, RIBC2.
- Products are 78–150 bp; Tm is 59–61 °C.
- All pairs span an exon junction except FOXD3-AS1, whose primers sit in different exons.
- **Possible off-target products:**
  - RPS18: pseudogenes and LINC01641, with 1 mismatch.
  - RPLP0: RPLP0P6 (0 mismatches) and RPLP0P10.
  - FOXD3-AS1: ITGB3BP (0 mismatches).
- None of the pairs has been checked in Primer-BLAST, and none has been ordered (UNCLEAR whether ordered).

### 5.7 TF activity (scripts/12): decoupleR ULM + CollecTRI (599 TFs scored)
- **HPV+ vs HPV- tumour:** 45 robust TFs.
  - 11 higher in HPV+: E2F4, ZMYND8, E2F1, E2F2, E2F3, SIM2, KLF1, MECOM, PPARG, MAML1, MEIS2.
  - 34 higher in HPV-, top: SMAD3, JUN, RUNX2, SMAD4, AP1, JUNB, JUND, ESR1, ETS1, HHEX, HOXD3, SNAI1, ETV4, SMAD2, SRF.
- **Positive controls:**
  - E2F1, E2F2 and E2F3 are robust and higher in HPV+ ✓.
  - TFDP1 is significant only at patient level.
  - FOXM1 and MYBL2 are not significant.
- **Tumour vs normal (A1, contrast level only):** 71 TFs with FDR < 0.05.
  - Higher in tumour, top: E2F4, E2F1, E2F2, E2F3, HCFC1, SRSF2, ARID3A, SALL4, E2F5, BARX2.
  - Higher in normal, top: JUNB, JUN, JUND, AP1, RELA, SP1, NFKB1, RFXANK, TP53, RFXAP.
- A2: 69 TFs with FDR < 0.05. **A3 (HPV-) was not run** through TF activity.
- **No GRN/regulon inference** has been done (no pySCENIC; no TF-target network beyond the CollecTRI prior).

### 5.8 Not done at all
- scVI
- Re-clustering of tumour cells / malignant substates
- Malignant re-calling (e.g. inferCNV on the tumour side)
- Pathway / gene-set enrichment
- Recurrence analysis

---

## 6. FIGURES (PNG at 600 dpi + PDF; `first_pass_all_nonmalignant/` holds older Fig4–6)

| File | Description |
|---|---|
| `01_atlas/Fig1_atlas_all_cells` | Atlas UMAP of all 274,911 cells by annotated cell type |
| `01_atlas/Fig2_epithelial_compartment` | Tumour + NM epithelial cells: position in the atlas, re-embedded UMAP by compartment and by HPV, cells per patient (20-cell line) |
| `02_normal_cleaning/Fig3_nonmalignant_cleaning` | Sub-clustered NM epithelial cells by subtype + marker dot plot per cluster |
| `02_normal_cleaning/Fig3b_CNV_check` | CNV heatmaps (BHN63, P8), signal-vs-correlation scatter per group, % tumour-like per group (PDF is git-ignored, 5.2 MB) |
| `03_DGE/Fig4_DGE_volcano` | Volcano plots A1–A3 and B + number of DE genes per contrast |
| `03_DGE/Fig5_DGE_overlaps` | UpSet/Venn of up-regulated genes across contrasts |
| `03_DGE/Fig6_DGE_checks` | CDKN2A by HPV; HPV16-E6 and MAGEA4 tumour vs normal, paired per patient |
| `04_qPCR_panel/Fig7_qPCR_candidate_screen` | Expression of the named candidates across atlas cell types; screen of the 95 HPV+ candidates (tumour share vs expression) |
| `04_qPCR_panel/Fig8_qPCR_panel_expectations` | Expected whole-tissue ("virtual bulk") HPV+ vs HPV- expression of the panel + TME markers; reference-gene ranking |
| `05_primers/Fig9_primer_map` | Primer positions and exons on each canonical transcript |
| `06_TF_activity/Fig10_TF_activity` | Top robust HPV-associated TFs; positive-control TFs per patient; top tumour-vs-normal TFs (A1) |
| `docs/umap.png` | Atlas UMAP screenshot (origin UNCLEAR, probably HNSCatlas.R) |

---

## 7. ENVIRONMENT

| Component | Version |
|---|---|
| R | 4.6.1 (2026-06-24) |
| Seurat / SeuratObject | 5.5.1 / 5.4.0 |
| edgeR / limma | 4.10.3 / 3.68.5 |
| decoupleR (R) | 2.17.0 |
| DESeq2 | 1.52.0 (installed, not used) |
| SingleCellExperiment | 1.34.0 |
| Matrix / ggplot2 | 1.7.6 / 4.0.3 |
| muscat, infercnv, zellkonverter | not installed |
| Python (system) | 3.9.6; scanpy, anndata, scvi-tools, decoupler, pydeseq2 and pySCENIC are **not installed** |
| Python for primers | run via `uv` 0.11.15 (`uv run --with primer3-py`), a temporary environment (primer3-py version UNCLEAR) |
| Machine | macOS (Darwin 25.6), 16 GB RAM. Loading the 4 GB atlas takes about 3 min and most of the RAM; the 1.1 GB epithelial object took 31 s and 2 GB |

---

## 8. PROBLEMS, CAVEATS AND OPEN ISSUES

1. **The HPV- tumour-vs-normal analysis is very underpowered.** A3 has only 5 HPV- patients (from 2 studies), against 17 HPV+. The single "HPV- candidate" (HOXD11) rests on this. For the new HPV-negative primary aim this is the main bottleneck. The normal reference is thin for HPV- (377 cells), and most Choi and Kurten HPV- patients have no clean normal cells.
2. **The analysis to date targets HPV+.** The candidate definition, qPCR screen, panel and primers are all HPV+-directed and cannot be reused as-is for an HPV- primary aim.
3. **HPV status is unverifiable outside Puram 2023.** There are no viral reads in Bill, Choi or Kurten. All HPV16-based signals (A1/A2 top hits HPV16-E1/E6/E7, and the "≥ 25 % viral" rule for likely-malignant clusters in scripts/05) only work for Puram patients.
4. **The atlas has no adjacent-normal tissue.** The "normal" group is NM epithelial cells from inside tumour samples. Before cleaning it mixed glandular (1,680), likely malignant (1,692) and ciliated cells. After cleaning, 15.8 % of the squamous cells were still CNV-positive in informative patients. In 26 of 74 patients CNV could not be checked, so residual tumour cells may remain in the reference for those patients.
5. **The atlas's Tumor / Epithelial labels are noisy.** "Tumor" cells include cells originally labelled B-cells (188), Endothelial (757), Erythroid (625), Fibroblast (419), Macrophages (206) and others. "Epithelial" includes 424 cells originally labelled T.cells and 920 "Unresolved". This project never re-called malignancy on the tumour side.
6. **CD45-sorted samples (Kurten):**
   - CD45p Tumor/Epithelial calls (627 cells) were excluded as likely doublets. That is reasonable.
   - The CD45n fractions remain in the tumour pseudobulk. They are excluded only from the whole-tissue (virtual-bulk) analysis.
   - 3 Kurten HPV- patients (HN02–HN04) are CD45p-only and have no epithelial cells.
7. **Integrated values:** the atlas UMAP and the Fig2 sub-UMAP use the Seurat-integrated assay. This is fine for visualisation. **DE used raw RNA counts (correct).** The scripts/05 NM sub-clustering used the uncorrected RNA assay, so clusters can be driven by patient or study. The "one patient ≥ 80 %" rule then labels such clusters "likely malignant", which could discard genuine normal cells from single-patient-dominated clusters (e.g. OP16, P8, P4).
8. **Unequal groups and study confounding in Plan B:** HPV+ is 11 / 31 from Puram, HPV- is 14 / 39 from Bill. The model adjusts for study, but within-study HPV contrasts are thin (e.g. Puram 5 HPV- vs 11 HPV+).
9. **Patients with few cells:** Plan B includes tumour pseudobulks with only 20–50 cells (e.g. BHN30 20, BHN49 26, P60 27, BHN52 24, P86 26), which makes them noisier.
10. **Anatomical site is unknown** (§3.2). HPV- oral cavity, larynx and hypopharynx cannot be separated, and the Bangladeshi target (gingivobuccal oral cancer) cannot be matched.
11. **There is no recurrence or outcome metadata** in the atlas, so Aim 2 needs external bulk cohorts.
12. **GSE103322 (Puram 2017)** has only log2TPM (no counts), so it is suitable for marker/score validation, not count-based DE. HPV labels for it are not yet attached.
13. **Reproducibility:**
    - The first-pass DGE script is not in the repo.
    - The HPV+ "novelty" classes in scripts/08 are hard-coded from a strategy-note PDF.
    - Package versions are not pinned (no renv lockfile).
    - The CollecTRI file was downloaded live from OmniPath (date UNCLEAR).
14. **Primer off-targets:** RPLP0 (exact match to the RPLP0P6 pseudogene) and FOXD3-AS1 (exact match to ITGB3BP) need Primer-BLAST review before ordering.

---

## 9. NOT YET DONE (no code or output found)

- **HPV-negative primary pipeline:** malignant vs normal re-definition for HPV- with enough patients; an HPV- candidate list; an HPV- qPCR screen and primers. Only A3 and HOXD11 exist.
- **Malignant cell calling by CNV on the tumour side** (inferCNV/copyKAT), as the new plan specifies. Only the normal side was checked.
- **scVI integration, Leiden clustering and malignant substate discovery** (named in `docs/HNSCC_Research_Pipeline_E.docx`). Not started; no Python single-cell stack is installed.
- **Gene regulatory networks** (pySCENIC regulons; ChIP-Atlas/Cistrome binding support). Only TF activity with a CollecTRI prior exists, and **not for A3/HPV-**.
- **Pathway / gene-set enrichment** (GSEA, ORA, PROGENy).
- **Recurrence arm (Aim 2):** cell-state programmes (p-EMT, hypoxia, stem-like), scoring in bulk cohorts with outcome data (ssGSEA/GSVA), survival models.
- **External validation:** GSE103322 single-cell; TCGA-HNSC bulk (cBioPortal/UALCAN/GEPIA).
- **Anatomical-site and clinical metadata** pulled from the 4 source GEO series.
- **Wet lab:** primer verification (Primer-BLAST), ordering, qRT-PCR. Nothing is recorded in the repo.
