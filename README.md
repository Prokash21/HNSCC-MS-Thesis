# Tumour-cell transcription factor networks and recurrence genes in HNSCC

**Single-cell discovery of tumor cell-specific transcription factor networks and
recurrence-associated genes in head and neck squamous cell carcinoma, with
validation in a Bangladeshi cohort.**

- **Primary analysis:** HPV-negative tumour cells.
- **Secondary arm:** HPV-positive tumour cells, run with the same pipeline.
- **Arm A:** tumour cell-specific genes, transcription factors and gene
  regulatory networks.
- **Arm B:** recurrence-associated programmes and genes, validated in bulk
  cohorts with outcome data.
- **Validation:** qRT-PCR on tissue from Bangladeshi HNSCC patients (tumour,
  adjacent normal, margins, contralateral and healthy mucosa).
- **Covariate:** anatomical site (oral cavity, pharynx, larynx) is treated as a
  covariate, never as a filter.

Discovery data: the HNSCatlas (Conde-Lopez et al., iScience 2026; 78 patients,
274,911 cells from GSE234933, GSE181919, GSE164690 and GSE182227). GSE103322
(Puram 2017) is kept apart as an independent dataset.

## Folder guide

```text
config/paths.R        every data and result path, sourced by all scripts
scripts/              v2 analysis scripts, numbered in run order
results/stepNN_*/     v2 outputs, one folder per step
archive_v1_HPVpos/    finished v1 HPV+ analysis (scripts + results), see its README
docs/
  PROJECT_STATUS_v1.md   status report of the v1 analysis
  pipeline/              earlier study plans
  references/            reference papers (PDFs not in git)
PLAN/                 research plan and sample collection protocol
data/                 not in git
  atlas/              HNSCC_Atlas.rds (4 GB)
  raw/                raw downloads (GEO, metadata)
  processed/          HNSCatlas/ and GSE103322/ intermediate objects
  external/           TCGA and other bulk cohorts
  reference/          Ensembl transcriptome, CollecTRI network
```

## Running

Run every script from the repo root, e.g. `Rscript scripts/01_metadata.R`.
Each script sources `config/paths.R`, so no path is hard-coded elsewhere.

| Step | Script | Output |
|---|---|---|
| 01 | `scripts/01_metadata.R`: clinical metadata and data availability | `results/step01_metadata/` |

## Methods principles

- Differential expression uses patient-level pseudobulk on raw counts (edgeR or
  DESeq2); integrated or batch-corrected values are never used for DE.
- The patient is the unit of replication.
- Discovery and validation datasets stay separate.
- FDR control throughout; simple pre-specified scores rather than fitted signatures.
- Every step follows a published reference method, cited in its summary.

## Data sources

- HNSCatlas: <https://hifis-storage.desy.de:2880/Helmholtz/E220-Radioonc_biol-DKFZ/HNSCC_Atlas.rds>
  → `data/atlas/`
- GEO series are downloaded by the scripts into `data/raw/`.

Requirements: R 4.6 (Seurat 5, edgeR, decoupleR, ggplot2, patchwork, ragg).
