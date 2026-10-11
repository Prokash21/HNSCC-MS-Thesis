# v1: HPV-positive tumour-cell analysis (archived)

Completed analysis from the first round of the thesis (Sept–Oct 2026), kept
unchanged for reference. Tagged in git as `v1-HPVpos-complete`.

## What v1 did

Single-cell analysis of the HNSCatlas (78 patients, 274,911 cells; four GEO
studies) focused on tumour-cell genes linked to HPV status:

1. Extracted the epithelial compartment (Tumor + non-malignant Epithelial,
   CD45+ sorts excluded) and harmonised HPV status per patient.
2. Cleaned the non-malignant epithelial cells into a normal squamous reference
   (sub-clustering, marker rules, inferCNV-style CNV check) → 1,136 cells.
3. Patient-level pseudobulk DE with edgeR QL (raw counts):
   - A1–A4 tumour vs normal squamous, paired within patient (22 patients: 5 HPV-, 17 HPV+)
   - B HPV+ vs HPV- tumour, adjusted for study (70 patients)
4. qRT-PCR screen of the 95 HPV+ candidates, a 4-gene panel and reference genes.
5. SYBR Green primer design (Primer3, exon-junction spanning).
6. TF activity with decoupleR ULM + CollecTRI.

## Key results

| Analysis | Result |
|---|---|
| A1 tumour vs normal (all) | 518 up / 237 down (FDR < 0.05, \|log2FC\| ≥ 1) |
| A2 tumour vs normal (HPV+) | 439 up / 146 down |
| A3 tumour vs normal (HPV-) | 31 up / 29 down (only 5 patients) |
| B HPV+ vs HPV- tumour | 400 up / 452 down |
| HPV+ candidates (A2 up ∩ B up) | 95 genes; 25 pass all qPCR filters |
| qPCR panel | SYCP2 (positive control), ZFR2, MIR9-3HG, KLHL35 |
| Reference genes | SDHA, RPS18, RPLP0, PSMC5 |
| TF activity HPV+ vs HPV- | 45 robust TFs (E2F1–4 higher in HPV+; SMAD3, AP-1, RUNX2 higher in HPV-) |

Full details: `docs/PROJECT_STATUS_v1.md`.

## Contents

- `scripts/00_HNSCatlas_explore.R`: exploratory atlas plot (was `HNSCC-HPV/atlas/HNSCatlas.R`)
- `scripts/01`–`12`: pipeline, run in numeric order from the repo root
- `results/01_atlas` … `results/06_TF_activity`: figures (PNG + PDF) and `tables/`

## Rerunning v1

The scripts were **not edited** after the move, so they still use the old
paths. To rerun one, first replace the old paths in it (or check out the
`v1-HPVpos-complete` tag, where the old layout still exists):

| Old path | New path |
|---|---|
| `scripts/NN_*.R`, `scripts/10_primer_design.py` | `archive_v1_HPVpos/scripts/` |
| `HNSCC-HPV/atlas/HNSCatlas.R` | `archive_v1_HPVpos/scripts/00_HNSCatlas_explore.R` |
| `HNSCC-HPV/atlas/HNSCC_Atlas.rds` | `data/atlas/HNSCC_Atlas.rds` |
| `HNSCC-HPV/results/<step>/` | `archive_v1_HPVpos/results/<step>/` |
| `HNSCC-HPV/docs/*` | `docs/pipeline/*` (reference papers → `docs/references/`) |
| `PROJECT_STATUS.md` | `docs/PROJECT_STATUS_v1.md` |
| `data/processed/`, `data/raw/`, `data/reference/` | unchanged |

`HNSCatlas.R` also calls `setwd("~/DWCT/Thesis-Cancer/Thesis/HNSCC-HPV/atlas")`,
which no longer exists.
