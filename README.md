# HPV-associated tumour-cell genes in HNSCC

Thesis analysis: single-cell RNA-seq of head and neck squamous cell carcinoma
(the HNSCatlas, Conde-Lopez et al., iScience 2026) to find tumour-cell genes
linked to HPV status, then shortlist them for qRT-PCR validation in a
Bangladeshi HNSCC cohort.

## Folder layout

```text
scripts/                    analysis pipeline, run in order from the repo root
data/                       large inputs and intermediate files (not in git)
  raw/GSE103322/            Puram 2017 raw data (independent validation set)
  processed/GSE103322/      parsed Puram 2017 data
  processed/HNSCatlas/      atlas subsets, pseudobulk caches, DGE objects
HNSCC-HPV/
  atlas/                    HNSCC_Atlas.rds (4 GB download, not in git)
  docs/                     research plans and the strategy note
  results/
    01_atlas/               Fig1-2, per-patient metadata
    02_normal_cleaning/     Fig3 sub-clusters of the non-malignant cells, Fig3b CNV check
    03_DGE/                 Fig4-6, DGE tables (+ first-pass run for comparison)
    04_qPCR_panel/          Fig7-8, candidate screen, suggested panel, reference genes
    05_primers/             Fig9 primer map, primer candidates and order sheet
    06_TF_activity/         Fig10, TF activity (HPV+ vs HPV-, tumour vs normal)
```

Each results folder holds its figures (PNG at 600 dpi and PDF) and a
`tables/` subfolder with the CSV files.

## Pipeline

| Script | What it does | Main outputs |
| --- | --- | --- |
| `01_parse_GSE103322.R` | Parses Puram 2017 (GSE103322) into a Seurat object | `data/processed/GSE103322/` |
| `02_atlas_extract_epithelial.R` | Loads the atlas, harmonises HPV per patient, extracts tumour and epithelial cells | `data/processed/HNSCatlas/` |
| `03_atlas_figures.R` | Fig1 (all cells), Fig2 (epithelial compartment) | `results/01_atlas/` |
| `04_patient_metadata.R` | One row per patient: study, sample type, HPV, sex, cell counts | `results/01_atlas/tables/` |
| `05_nonmalignant_subtypes.R` | Sub-clusters the non-malignant epithelial cells and keeps clean squamous cells as the normal reference (Fig3) | `results/02_normal_cleaning/` |
| `05b_cnv_check.R` | inferCNV-style copy-number check of that cleaning, independent of marker genes; removes CNV-positive cells from the normal reference (Fig3b) | `results/02_normal_cleaning/` |
| `06_dge_pseudobulk.R` | Pseudobulk DGE with edgeR: tumour vs normal (paired, A1-A4), HPV+ vs HPV- tumour (B) | `results/03_DGE/tables/` |
| `07_dge_figures.R` | Fig4 volcano plots, Fig5 UpSet/Venn, Fig6 checks | `results/03_DGE/` |
| `08_qpcr_panel.R` | qPCR screen of the HPV+ candidates (cell-type specificity, tumour share of the whole-tissue signal, expression level, whole-tissue HPV effect), suggested panel, reference-gene pre-screen | `results/04_qPCR_panel/tables/` |
| `09_qpcr_figures.R` | Fig7 candidate screen, Fig8 expected tissue qPCR results and reference genes | `results/04_qPCR_panel/` |
| `10_primer_design.py` | SYBR Green qPCR primers (Primer3, exon-junction spanning) for the panel, screened against the human transcriptome (Ensembl) | `results/05_primers/tables/` |
| `11_primer_map.R` | Fig9: primer positions on each transcript | `results/05_primers/` |
| `12_tf_activity.R` | TF activity with decoupleR (ULM) + CollecTRI: HPV+ vs HPV- and tumour vs normal (Fig10) | `results/06_TF_activity/` |

Run a step with `Rscript scripts/<script>.R` from the repo root. Scripts 02, 05b
and 08 (first run) load the 4 GB atlas: about 3 minutes and most of 16 GB RAM.
Later steps reuse cached files in `data/processed/HNSCatlas/`.

Requirements: R 4.6 with Seurat 5, edgeR, decoupleR, ggplot2, ggrastr,
ggrepel, patchwork and ragg. The primer script is Python, run through uv
(`uv run --with primer3-py python scripts/10_primer_design.py`), which uses a
temporary environment and leaves the system Python untouched. Scripts 10 and
12 download reference data once (Ensembl transcriptome, CollecTRI from
OmniPath) into `data/reference/`.

## Data

- **HNSCatlas**: `HNSCC_Atlas.rds` from the link in the paper
  (<https://hifis-storage.desy.de:2880/Helmholtz/E220-Radioonc_biol-DKFZ/HNSCC_Atlas.rds>),
  placed in `HNSCC-HPV/atlas/`. It combines GSE234933, GSE182227, GSE164690
  and GSE181919 (78 patients, 274,911 cells).
- **GSE103322** (Puram 2017): not part of the atlas, so usable as an
  independent dataset. Raw file goes in `data/raw/GSE103322/`.

## Key decisions

- The patient is the unit of replication: counts are summed per patient
  (pseudobulk) and tested with edgeR quasi-likelihood. Significant means
  FDR < 0.05 and |log2FC| >= 1.
- The atlas has no normal-tissue samples. "Normal" is the clean squamous
  subset of the atlas' non-malignant epithelial cells; gland, ciliated and
  likely-malignant clusters and immune/stromal doublets are removed, and so
  are cells with a tumour-like copy-number profile (`results/02_normal_cleaning/`).
- Patient OP16 had one stray HPV+ cell among 3,329 HPV- cells; HPV is set per
  patient, which matches the paper (46 HPV-, 32 HPV+).
- Tumour/epithelial calls inside CD45+ sorted fractions are dropped as likely
  doublets.
- Viral HPV16 reads exist only for Puram 2023 (GSE182227); the other studies
  did not map them.
