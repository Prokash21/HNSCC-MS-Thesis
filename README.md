<h1 align="center">Recurrence genes and TF networks in HNSCC</h1>

<p align="center">
  <b>Single-cell RNA sequencing reveals recurrence-associated genes and their
  transcription factor networks in head and neck squamous cell carcinoma (HNSCC),
  with validation in a Bangladeshi cohort.</b>
</p>

<p align="center">
  <img src="https://img.shields.io/badge/R-4.6-276DC3?logo=r&logoColor=white" alt="R 4.6">
  <img src="https://img.shields.io/badge/Seurat-v5-6A5ACD" alt="Seurat v5">
  <img src="https://img.shields.io/badge/scRNA--seq-78%20patients%20%7C%20~275k%20cells-2E8B57" alt="78 patients, ~275k cells">
  <img src="https://img.shields.io/badge/focus-HPV--negative-C0392B" alt="HPV-negative">
</p>

---

## Overview

| Part | What |
| --- | --- |
| 🧬 **Discovery** | HNSCatlas single-cell data (GSE234933, GSE181919, GSE164690, GSE182227). Primary analysis on HPV-negative tumour cells. |
| 🎯 **Tumour-cell genes and TFs** | Malignant cells identified by inferred CNV, compared with normal epithelium by patient-level pseudobulk (edgeR), and scored for TF activity (decoupleR, pySCENIC). |
| 📈 **Recurrence** | Programme and gene scores tested against outcome in bulk cohorts (TCGA-HNSC, GEO). |
| 🧪 **Validation** | qRT-PCR panel on tumour and matched normal tissue from Bangladeshi HNSCC patients. |

```mermaid
flowchart LR
    A[Single-cell atlas] --> B[CNV: malignant cells]
    B --> C[Pseudobulk DE<br/>tumour vs normal]
    C --> D[TF activity<br/>and regulons]
    C --> E[Recurrence<br/>TCGA / GEO]
    D --> F[qRT-PCR panel]
    E --> F
    F --> G[Bangladeshi cohort]
```

## Highlights

<table>
  <tr>
    <td width="50%" align="center">
      <b>Malignant cells by inferred copy number</b><br>
      <img src="results/step03_cnv/Fig11_cnv_heatmaps.png" alt="CNV heatmaps">
    </td>
    <td width="50%" align="center">
      <b>HPV-negative tumour vs normal epithelium</b><br>
      <img src="results/step04b_DE/Fig18_volcano_C2.png" alt="Volcano plot">
    </td>
  </tr>
</table>

<p align="center"><b>Pathways → transcription factors → candidate genes</b></p>
<p align="center"><img src="results/step05b_functional/Fig46_story_diagram.png" alt="Pathway TF gene diagram" width="90%"></p>

<p align="center"><b>Tumour programme scores and outcome in TCGA HPV-negative HNSCC</b></p>
<p align="center"><img src="results/step07a_recurrence/Fig31_tcga_cox_forest.png" alt="Cox forest plot" width="90%"></p>

## Repository

```text
config/paths.R      all data and result paths
scripts/            analysis scripts, numbered in run order
results/stepNN_*/   outputs and figures, one folder per step
PLAN/               research plan and sample collection protocol
archive_v1_HPVpos/  earlier HPV-positive analysis
```

Run scripts from the repo root, e.g. `Rscript scripts/01_metadata.R`.
Requires R 4.6 (Seurat 5, edgeR, DESeq2, decoupleR, survival, ggplot2).
Data files are not in git; the atlas is available from
<https://hifis-storage.desy.de:2880/Helmholtz/E220-Radioonc_biol-DKFZ/HNSCC_Atlas.rds>.
