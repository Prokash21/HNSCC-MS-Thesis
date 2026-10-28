# Step 07c plan: TCGA HPV- oral cavity subgroup (EXPLORATORY)

Written 2026-10-28, before any step 07c outcome result was computed. Every later change is listed as post hoc in `STEP07c_SUMMARY.md`. There is no primer design in this step, and this step does not change the step 07a shortlist (`results/step07a_recurrence/candidates_for_step08.csv`).

## 0. Status of this analysis

- **This analysis is EXPLORATORY and post hoc.** It was not part of the step 07a or 07b plans.
- **Why it is done:**
  - Step 07b found the programme signal strong in the oral-only cohort GSE41613 and absent in the mixed-site cohort GSE65858.
  - Our own Bangladeshi cohort is oral cavity only.
- **Results already seen before planning:**
  - TCGA HPV- all-site results (step 07a);
  - GSE41613 and GSE65858 results (step 07b).
- **Results not yet seen:** to my knowledge, no TCGA site-specific (oral or larynx) outcome result has been computed. The removed 2026-10-27 run is not kept, so whether it looked at a site subgroup is UNCLEAR.
- **Deviation from the project rule.** The rule is "anatomical site is a covariate, never a filter". This step filters on site at the user's request, for exploration only.
  - The primary analysis stays the all-site step 07a analysis with site as a covariate.
  - Part C (site interaction in all-site data) is the formal check of whether the effect differs by site.
  - No subgroup result here is used to select or drop genes.

## 1. Population and data (reused from step 07a, nothing recomputed)

- **Patients:** TCGA HPV- primary tumours, 415 with RNA-seq (`data/processed/v2/step07a_tcga_hpvneg.rds`). Endpoints and exclusions are as in step 07a: TCGA-CDR PFI (primary) and OS (secondary), with redacted, missing or zero-time patients excluded.
- **Site:** cBioPortal primary site, grouped exactly as in step 07a.
  - **Oral cavity** (expected about 277): Oral Tongue, Oral Cavity, Floor of mouth, Buccal Mucosa, Alveolar Ridge, Hard Palate, Lip.
  - **Larynx** (about 106).
  - **Other** (about 32): base of tongue, tonsil, oropharynx, hypopharynx. This group is not analysed alone.
- **Covariates:**
  - stage, the CDR `ajcc_pathologic_tumor_stage` grouped I–II / III / IV, with unknown as its own level, as in step 07a;
  - age, the CDR `age_at_initial_pathologic_diagnosis`.
- **Scores:** the step 07a programme scores (`B_tcga_programme_scores.csv`). These are the mean z of TMM log CPM, z computed across all 415 HPV- tumours, and are not recomputed in the subgroup.
- **Gene values:** step 07a log CPM.
- **HR scale:** in each model the score or gene is standardised within the analysed patients, so HR is per 1 SD within that subgroup.
- **REMARK profile** (Sauerbrei 2018), for oral and larynx: flow, n, PFI and OS events, median follow-up (reverse Kaplan–Meier), stage and age distribution, missing data (complete case per model), cut-off rule.

## 2. Programmes (B), oral cavity; larynx for contrast (descriptive)

- **Sets:** all 9 step 07a sets: 5 primary (Toustrup, Buffa, Puram p-EMT, Gavish G1/S, Gavish G2/M) and 4 concordance (Gavish MP6, Gavish MP13, Tirosh G1/S, Tirosh G2/M).
- **Cox models**, for PFI (primary) and OS (secondary):
  - univariable `Surv ~ score` (reported for REMARK);
  - **main model** `Surv ~ score + stage`;
  - sensitivity `Surv ~ score + stage + age`.
- **Reporting:** HR per SD, 95% CI, p. BH FDR on the main-model p values per subgroup and endpoint, within the 5 primary sets and within the 4 concordance sets.
- **Proportional hazards:** `cox.zph` (score term and global) for every model.
- **Kaplan–Meier:** median split (median of the score in the analysed oral patients; no cut-point search), PFI and OS, log-rank p. Descriptive only.
- **Larynx:** the same models, descriptive only, with the number of events stated. There is no FDR-based claim for larynx.
- **No validation verdict is given.** This step is exploratory.

## 3. Site interaction (C), all-site TCGA HPV-

- **Variable:** `oral` = oral cavity vs not (larynx + other).
- **Model**, for every set, PFI and OS: `Surv ~ score * oral + stage`, with the score standardised over all analysed patients.
- **Reporting:**
  - the interaction term (ratio of the HR per SD in oral to the HR per SD in non-oral) with 95% CI and p;
  - the HR per SD within oral and within non-oral derived from the same model.
- **Power:** interaction tests have low power, so a non-significant interaction is not evidence that the effect is the same across sites, and a nominal p < 0.05 across 9 sets × 2 endpoints is not a finding by itself. No FDR-based claim is made.

## 4. Genes (D), oral cavity, descriptive

- **Genes:**
  - the 7 eligible step 07a genes: P4HA1, HMMR, KIF23, SEC61G, CCNA2, ECT2, CKAP2;
  - MYBL2, E2F4, HIF1A, EPAS1;
  - every step 06 core, reserve and CA9 gene (`D_panel_final.csv`; reference genes excluded).
- **Model:** per-gene Cox on standardised log CPM, main model `~ gene + stage` and sensitivity `+ age`, PFI and OS.
- **Reporting:** HR, 95% CI, p, and BH FDR across the listed genes per endpoint (main model). Genes not in the step 07a log CPM matrix are listed.
- **Network table:** one small table with NFIC, FOXM1, HMGB3 and MYBL2 together.

## 5. Oral side by side (E): TCGA oral vs GSE41613

- **Rows:**
  - programmes: TCGA oral PFI (main and +age models) next to GSE41613 DSS (step 07b multivariable, age + stage);
  - secondary: TCGA oral OS next to GSE41613 OS;
  - genes: the part D genes, where present in GSE41613 (step 07b `V_gene_cox.csv` covers the 7 eligible genes, MYBL2, E2F4 and the panel genes).
- **Not available from step 07b:** HIF1A and EPAS1 in GSE41613. They are marked "not run in step 07b" and not computed here.
- **No pooling or meta-analysis.**

## 6. Outputs

- **Scripts:** `scripts/07c_oral.R` (models, tables) and `scripts/07c_figures.R`.
- **Tables:** `results/step07c_oral/tables/`.
- **`STEP07c_SUMMARY.md`** (local, git-ignored like the other summaries), clearly labelled EXPLORATORY:
  - n and events;
  - results;
  - limitations and UNCLEAR items;
  - a reference for every method.
- **Figures:**
  - forest plot, oral vs larynx vs all-site TCGA;
  - Kaplan–Meier curves for the 3 TCGA-prognostic sets in oral;
  - oral side by side, TCGA vs GSE41613.
