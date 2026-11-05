# Step 07a plan: recurrence arm and survival (HPV- only)

Written 2026-10-28, before any step 07a outcome result was computed. Every later change is listed as post hoc in `STEP07a_SUMMARY.md`. There is no primer design in this step.

## 0. Seen before planning (transparency)

An earlier exploratory run was made on 2026-10-27 and removed on 2026-10-28 at the user's request (commits 3246bb9 and 82bc904, reset to 095a265). It had already looked at these outcome results:

- **TCGA HPV- PFI** (Cox with age + stage + site):
  - the precomputed Buffa hypoxia score from cBioPortal: HR 1.25 per SD, p = 0.011;
  - the MSigDB Hallmark hypoxia score: not significant;
  - de novo NMF meta-programmes from our HPV- malignant cells, including a G2/M-like and a p-EMT/hypoxia-like programme: none at FDR < 0.1;
  - every step 06 panel, reserve and CA9 gene: none at p < 0.05.
- **Bill relapse vs primary:** malignant-cell pseudobulk DE (no gene at FDR < 0.05) and de novo programme scores (none at FDR < 0.05; the p-EMT/hypoxia-like programme was nominally higher in locoregional recurrences).
- **Two GEO cohorts** (GSE65858, GSE41613): not part of this plan.

These results are therefore **seen before planning**. The gene sets, models and thresholds below come from the user's specification of 2026-10-28.

## 1. Data and exclusions (same as steps 03–06)

- **Single cell:**
  - HPV- malignant cells: `cell_state_4b == "malignant"`, `primary_4b == TRUE`.
  - Bill primaries: BHN patients in the step 05 export. This excludes the second/third primaries BHN39, BHN46 and BHN49, and the nasopharynx samples.
  - Bill relapses: role `recurrence`, HPV-, `primary_4b == TRUE`.
  - HN23 and HN45 are two recurrences of one patient. They are pooled as patient BHN23, since the patient is the unit.
  - A patient enters the Bill tests with ≥ 20 malignant cells (as in steps 05/06).
  - The full atlas is not loaded.
- **TCGA:**
  - HPV- primary tumours only (cBioPortal `SUBTYPE` = HNSC_HPV-, as in steps 01/04b/06).
  - One primary tumour (`-01`) RNA-seq sample per patient (first by sample ID).
  - Counts: GDC STAR unstranded counts from the UCSC Xena GDC hub (`data/raw/step07_bulk/TCGA-HNSC.star_counts.tsv.gz`). These were verified equal to the raw GDC STAR files for 5 samples.
- **Endpoints:** TCGA-CDR Table S1 (Liu 2018), PFI (primary) and OS (secondary), in days. Patients with redaction flags or a missing or zero time are excluded and counted in the REMARK flow.
- **Covariates:**
  - Stage: CDR `ajcc_pathologic_tumor_stage`, grouped I–II / III / IV.
  - Site: cBioPortal primary site, grouped oral cavity / larynx / other (oropharynx and hypopharynx), as in step 06.
  - Age: CDR `age_at_initial_pathologic_diagnosis`.
  - Missing stage is coded as its own level ("unknown") and the number is reported. Site is a covariate, never a filter.

## 2. Gene sets (taken from the papers' own tables)

| Programme | Role | Set | Source |
|---|---|---|---|
| 1 Hypoxia | primary | Toustrup 2011, 15 genes | Table S2, genes marked "**" ("final 15-Gene Hypoxia Classifier"); transcribed to `config/gene_sets/toustrup2011_tableS2.txt` |
| 1 Hypoxia | primary | Buffa 2010 metagene, 51 genes | Table S5 "Common Hypoxia Signature Genes", file `41416_2010_BF6605450_MOESM10_ESM.xls` |
| 1 Hypoxia | concordance | Gavish 2023 MP6 "Hypoxia", 50 genes | `41586_2023_6130_MOESM6_ESM.xlsx`, sheet "Cancer MPs" |
| 2 p-EMT | primary | Puram 2017 p-EMT, 100 genes | Table S7, column "p-EMT", genes 1–100 (`mmc1.pdf`, p. 6); transcribed to `config/gene_sets/puram2017_tableS7_pEMT.txt` |
| 2 p-EMT | concordance | Gavish 2023 MP13 "EMT-II", 50 genes | same sheet; the Gavish supplement (`MOESM3`, p. 14) describes MP13 as "highly consistent with the EMT program described previously by Puram et al." |
| 3 Proliferation | primary | Gavish 2023 MP1 "Cell Cycle - G2/M" and MP2 "Cell Cycle - G1/S", 50 genes each | same sheet |
| 3 Proliferation | concordance | Tirosh 2016 G1/S (43 genes) and G2/M (55 genes) | Table S5 "Cell cycle gene-sets", phase-specific genes (`aad0501_table_s5.xlsx`) |

- **Symbol mapping:** all symbols are mapped to current HGNC symbols with the step 04b function `hgnc_map`. Each change is listed.
- **Scoring rule:** genes absent from a dataset are dropped from that score, and the coverage is reported.
- **Exploratory set dropped: tumour-budding signature** (Ourailidis 2026). The 28 genes are given only as figure labels (main Fig. 3a; Additional File 1 Fig. S8), not as a table. The paper reports that its score correlates with the Puram p-EMT score (Spearman rho 0.94 spatial, 0.90 scRNA-seq) and shares 7 of the 15 Puram p-EMT genes. Programme 2 therefore covers this biology. Ourailidis 2026 is cited only in the rationale and discussion.
- **Panel overlap:** the overlap of every set with the step 06 core panel, reserve and CA9 is reported.

## 3. Scores

- **Single cells:** Tirosh 2016 method.
  - Relative expression = log2(CPM/10 + 1) minus the gene mean over the scored cells.
  - Score = mean relative expression of the set minus that of control genes. Controls are 100 per set gene, drawn from the same of 30 expression bins.
  - Patient score = mean over that patient's malignant cells.
- **Pseudobulk and TCGA:** mean z-score of log expression.
  - log2 CPM with TMM (edgeR, prior count 2 for pseudobulks, 1 for TCGA), after `filterByExpr`.
  - z per gene across the samples of that analysis; score = mean z of the set genes present.
  - In Cox models the score is standardised, so HR is per 1 SD.

## 4. Bill GSE234933: relapse vs primary (B)

- **Groups:** relapse = locoregional recurrences + distant metastases; primary = Bill HPV- primaries. No patient has both a primary and a relapse sample in the eligible data (HN13 is not BHN46, step 03), so all tests are unpaired. n per group is stated.
- **Programme scores:** two-sided Wilcoxon rank-sum, relapse vs primary, for every set (primary and concordance), using the patient mean single-cell score. BH FDR within the primary family and within the concordance family. The pseudobulk mean-z score is reported alongside.
- **DE:** malignant-cell pseudobulk on raw counts, edgeR QL as in step 04b: `filterByExpr`, TMM, robust QL fit, design `~ site + group`. Site is oral cavity vs other/unknown. The log2FC relapse vs primary per gene is used in D1.4 and D2.
- **Sensitivity:** locoregional recurrences only.
- **Limitation:** treatment history and time to relapse are unknown, so differences may reflect treatment effects.

## 5. TCGA survival (C)

1. **Programme Cox models**, for every set, for PFI and OS:
   - univariable: `Surv ~ score`;
   - multivariable (primary): `Surv ~ score + stage + site`;
   - pre-specified sensitivity: `+ age`.

   BH FDR on the multivariable p values, per endpoint, within the 5 primary sets (Toustrup, Buffa, Puram p-EMT, Gavish G1/S, Gavish G2/M), and separately within the 4 concordance sets.
2. **Per-gene Cox** (z-scored log CPM; univariable and `~ gene + stage + site`), for PFI and OS. The genes are all primary-set genes plus every step 06 core, reserve and CA9 gene (reference genes excluded). Report HR, 95% CI, p and BH FDR, with FDR across all tested genes per endpoint.
3. **Proportional hazards:** `cox.zph` (score term and global) for every model. Violations are reported, not corrected.
4. **Kaplan–Meier:** median split of each primary-set score (fixed cut-off, no optimal cut-point search), PFI and OS, with log-rank p. Descriptive only.
5. **REMARK profile table:**
   - patient flow (CDR HNSC → HPV- → RNA-seq → valid endpoint);
   - n and events per endpoint;
   - median follow-up (reverse Kaplan–Meier);
   - cut-off rule, missing-covariate handling;
   - univariable and multivariable HR with 95% CI, and the PH check.

## 6. Selection (D), programmes 1–3 only

**D1. Eligibility.** A gene is eligible if ALL hold:
1. A primary set containing the gene is associated with PFI in TCGA: multivariable HR > 1 and FDR < 0.10 (primary family). A gene in several primary sets passes if any of them passes. For hypoxia, Toustrup and Buffa are judged separately.
2. Per-gene multivariable Cox for PFI: HR > 1 (same direction) and p < 0.05.
3. Malignant pseudobulk CPM ≥ 10: median over the step 04b HPV- malignant pseudobulks with ≥ 20 cells, as in step 06. The tumour-cell share is reported as in step 04b part F: the median over HPV- patients of malignant counts divided by counts in all cell-type pseudobulks (`celltype_pseudobulk_v2.rds`), recomputed for every set gene with the step 04b code.
4. Bill relapse vs primary log2FC > 0 (significance not required). A gene removed by `filterByExpr` in the Bill DE counts as "not evaluable" and fails.
5. TCGA median TPM ≥ 10 over the HPV- primary tumours.
   - TPM is computed from the counts with GDC gene lengths. Lengths are derived as counts / `tpm_unstranded` from the 86 GDC STAR files already downloaded (step 04b).
   - The computed TPM is checked against `tpm_unstranded` for those samples.

Output `D_recurrence_eligible_all.csv`: every set gene with all metrics, its programme(s), a flag for the step 06 panel/reserve/CA9, the pass/fail of each criterion and an `eligible` column, without a cap. If a programme fails criterion 1, its genes are still listed, marked "programme not prognostic".

**D2. Ranking.**
- Eligible genes are ranked by the mean of 4 percentile ranks, each with higher = better:
  - log HR for PFI (multivariable);
  - Bill relapse log2FC;
  - tumour-cell share;
  - log TCGA median TPM.
- Output `D_recurrence_ranked.csv` in rank order.
- Suggested shortlist: walk down the ranking and take genes in strict order, at most 2 per programme, up to a suggested size of 4. A gene in two programmes counts toward both. Genes already in the panel count toward the shortlist.
- There is no minimum. If nothing is eligible this is stated, and thresholds are not relaxed.

**D3. Regulators.**
- For every primary set that passes D1.1, test enrichment of each TF's targets among the set genes (one-sided hypergeometric). TF targets come from two sources:
  - consensus pySCENIC regulons (step 05, `regulons_consensus.csv`);
  - CollecTRI (`collectri_omnipath.tsv`; complexes such as AP1/NFKB are excluded because they are not single genes).
- Universe: the 19,965 genes of the step 05 single-cell export; targets and set genes are restricted to it.
- BH FDR across all TF × set tests per source.
- Eligible: FDR < 0.05 and the TF expressed in malignant cells (malignant pseudobulk CPM ≥ 10).
- Ranked by enrichment FDR (best over sources); suggested shortlist at most 2.
- If no programme passes D1.1, D3 has no input and this is stated.

**D4.** Write `candidates_for_step08.csv` with columns gene, role (programme gene / regulator), source (set or TF source), rank, shortlist (yes/no). This file is the input to the separate final primer step.

## 7. Outputs

- `STEP07a_SUMMARY.md` covering:
  - results with n per group;
  - decisions, post hoc changes and UNCLEAR items;
  - primary vs concordance agreement per programme (direction and FDR of both, gene overlap);
  - the full eligible list size and the shortlist;
  - a reference for every method.
- Figures:
  - Bill programme scores, relapse vs primary;
  - TCGA Cox forest plots, PFI and OS;
  - Kaplan–Meier curves;
  - the REMARK profile table;
  - the eligible-gene ranking plot.
