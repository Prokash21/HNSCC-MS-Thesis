# Step 07b plan: independent validation of the step 07a programmes (HPV- only)

Written 2026-10-28, before any step 07b outcome result was computed. Every later change is listed as post hoc in `STEP07b_SUMMARY.md`. There is no primer design in this step, and this step does not change the step 07a shortlist (`results/step07a_recurrence/candidates_for_step08.csv`).

## 0. Seen before planning (transparency)

- **Removed exploratory run.** An earlier exploratory run was made on 2026-10-27 and removed on 2026-10-28 at the user's request (commits 3246bb9 and 82bc904). It had already tested GSE65858 and GSE41613 for:
  - the MSigDB Hallmark hypoxia score;
  - de novo NMF programmes from our HPV- malignant cells;
  - 61 single genes, including every step 06 panel, reserve and CA9 gene.

  Its numbers are not kept in the repository. The step 07a sets (Toustrup, Buffa, Gavish MP1/MP2/MP6/MP13, Puram, Tirosh) were **not** tested in these two cohorts there. The panel, reserve and CA9 gene results in part D below are therefore **seen before planning**, and part D is descriptive only.
- **Step 07a TCGA results** are known (`STEP07a_SUMMARY.md`): Buffa, Gavish G1/S and Gavish G2/M passed the PFI criterion; Toustrup and Puram p-EMT did not. This plan tests those results; it does not re-select sets.
- **Metadata read for this plan.** To name the GEO fields, I tabulated the sample characteristics of both series matrices, which include the overall event counts. No expression value was read, and no association between expression and outcome was computed. Counts seen, over all samples in each series:
  - GSE65858 (270 samples): PFS events 133 and OS events 94;
  - GSE41613 (97 samples): `vital` = Alive 46, Dead-oral ca 30, Dead-non OC 14, Dead -unk cause 7.

## 1. Cohorts, fields and exclusions

Fields are the `!Sample_characteristics_ch1` keys of the GEO series matrices in `data/raw/step07_bulk/`.

### GSE65858 (Wichmann 2015; Illumina HumanHT-12 v4, GPL10558)

- **Patients:**
  - `hpv16_dna_rna` = "DNA-" (HPV16 DNA-negative; n = 196). These are the same samples as `hpv_dna` = "Negative".
  - Samples with "Other HPV", HPV16 DNA+ or NA are excluded.
  - `tumor_type` = "Primary" only. "Relapse" and "Secondary" tumours are excluded, as TCGA used primary tumours only.
  - The series is already one array per patient: the submitters removed duplicate pairs.
- **Primary endpoint PFS:** time `pfs`, event `pfs_event` (TRUE/FALSE). **Secondary endpoint OS:** `os`, `os_event`.
  - Times are taken as days. The values reach about 2,400, which is not plausible in months. The unit is not stated in GEO, so this is checked against the paper and reported.
  - Patients with time ≤ 0 or a missing time or event are excluded and counted in the flow.
- **Covariates:**
  - age: `age`, continuous years;
  - stage: `uicc_stage`, grouped I–II / III / IV (IVA, IVB and IVC → IV);
  - site: `tumor_site`, grouped oral cavity ("Cavum Oris") / larynx / other (oropharynx, hypopharynx).
- **Expression:** series matrix values. These are log2, RSN-normalised and ComBat-corrected, with 31,330 features (`!Sample_data_processing`). They are used as given.

### GSE41613 (Lohavanichbutr 2013; Affymetrix HG-U133 Plus 2.0, GPL570)

- **Patients:** all 97 samples. Every one is "HPV-negative oral squamous cell carcinoma (OSCC)" (`tissue type`), one per patient.
- **Primary endpoint, disease-specific survival (DSS):** time `fu time`; event = `vital` "Dead-oral ca".
  - "Alive", "Dead-non OC" and "Dead -unk cause" are censored at `fu time`.
  - Times are taken as months: the maximum is about 85, which fits the study's follow-up. This is checked against the paper and reported.
- **Secondary endpoint OS:** event = any `vital` starting with "Dead".
- **Pre-specified sensitivity (descriptive, not part of the validation rule):** DSS with "Dead -unk cause" counted as events.
- **Covariates:**
  - age: `age`, given only as 4 groups (19–39, 40–49, 50–59, 60–88). Entered as an ordinal score 1–4 (1 df), to keep the number of model terms low for about 30 events.
  - stage: `tumor stage`, I/II vs III/IV.
  - site: none, since all tumours are oral cavity.
- **Expression:** series matrix values (gcRMA, `!Sample_data_processing`), used as given. If the values are not on a log2 scale, they are log2-transformed and this is reported.

### Both cohorts

- **Missing covariates:** complete-case analysis in each model. Numbers dropped are reported in the REMARK profile. Site is a covariate, never a filter.
- **Discovery and validation stay separate:** no parameter (gene weight, cut-off, gene list) is estimated in these cohorts.

## 2. Probes to genes

- **Probe annotation:** the `Gene symbol` column of `GPL10558.annot.gz` / `GPL570.annot.gz` (GEO platform annotation files).
  - Probes with an empty symbol, or with several symbols ("///"), are dropped.
  - Symbols are mapped to current HGNC symbols with the step 04b function `hgnc_map`; probes that do not map are dropped.
- **Collapse:** one probe per gene, the probe with the highest mean expression over all samples in the series (`maxMean`, Miller 2011). This is computed before the patient filter, so the probe choice does not depend on the analysed subset.
- **Gene sets:** the step 07a lists in `results/step07a_recurrence/tables/A_gene_sets.csv`, which are built from `config/gene_sets/` and the paper tables. They are not changed.
  - Coverage (genes present / genes in set) is reported per set and platform.
  - Set genes absent from a platform are dropped from that score.

## 3. Scores and models (as in step 07a)

- **Score:** for each gene, z-score the log expression across the analysed patients of that cohort and endpoint data set. The score is the mean z of the set genes present. In Cox models the score is standardised again, so HR is per 1 SD.
- **Cox models**, for every set (5 primary + 4 concordance), for the primary and secondary endpoint of each cohort:
  - univariable: `Surv ~ score`;
  - multivariable (primary model):
    - GSE65858: `Surv ~ score + age + stage + site`;
    - GSE41613: `Surv ~ score + age + stage`.
- **Proportional hazards:** `cox.zph`, for the score term and global, on every model. Violations are reported, not corrected.
- **Kaplan–Meier:** median split of each primary-set score (fixed cut-off, no cut-point search), for both endpoints, with log-rank p. Descriptive only.
- **REMARK profile per cohort** (Sauerbrei 2018):
  - patient flow (series → HPV filter → primary tumour → valid endpoint);
  - n and events per endpoint;
  - median follow-up (reverse Kaplan–Meier);
  - covariate distribution;
  - missing-data handling;
  - cut-off rule;
  - univariable and multivariable HR with 95% CI, and the PH check.
- **Multiple testing:** BH FDR on the multivariable p values is reported per cohort and endpoint, within the 5 primary sets and within the 4 concordance sets. The validation rule (section 4) uses the raw p values, as specified.

## 4. Pre-specified validation rule (primary sets only)

- **Rule:** a programme set that was prognostic in TCGA (Buffa, Gavish G1/S, Gavish G2/M) is **validated** if, on the primary endpoint (GSE65858 PFS; GSE41613 DSS) in the multivariable model:
  - HR > 1 in **both** cohorts; **and**
  - p < 0.05 in **at least one**.
- **Otherwise** it is "not validated". The verdict also reports which part failed (direction or significance).
- **Other sets:**
  - Toustrup and Puram p-EMT (not prognostic in TCGA) are reported with the same models, without a verdict.
  - The concordance sets (Gavish MP6 hypoxia, Gavish MP13 EMT-II, Tirosh G1/S, Tirosh G2/M) are reported alongside, descriptive only.
- **Side-by-side table:** TCGA results are copied from step 07a (`S_programme_cox.csv`), not recomputed.
  - Both 07a models are shown: the multivariable model (stage + site, the 07a primary) and the +age sensitivity model, whose covariates match the GEO models.
  - The endpoints differ (TCGA PFI, GSE65858 PFS, GSE41613 DSS). This is stated as a limitation.

## 5. Genes (descriptive only)

- **Genes:**
  - the 7 eligible step 07a genes: P4HA1, HMMR, KIF23, SEC61G, CCNA2, ECT2, CKAP2;
  - the shortlisted regulators MYBL2 and E2F4;
  - every step 06 core, reserve and CA9 gene (`D_panel_final.csv`, sets core / reserve / optional (Arm B); reference genes excluded).
- **Model:** per-gene Cox (z-scored log expression; univariable and the multivariable model of section 3), both cohorts, primary and secondary endpoints.
- **Reporting:** HR, 95% CI, p, and BH FDR across the listed genes per cohort and endpoint. There is no selection or verdict. Genes absent from a platform are listed.

## 6. Outputs

- **Scripts:** `scripts/07b_cohorts.R` (parse, filter, collapse, scores, coverage), `scripts/07b_survival.R` (Cox, KM, REMARK, genes), `scripts/07b_figures.R`.
- **Tables:** `results/step07b_validation/tables/`.
- **`STEP07b_SUMMARY.md`** (kept local, git-ignored like the other step summaries):
  - n and events;
  - results table with TCGA and both cohorts side by side;
  - validation verdict per programme;
  - limitations, post hoc changes and UNCLEAR items;
  - a reference for every method.
- **Figures:**
  - forest plot of programme HRs across TCGA and both cohorts;
  - Kaplan–Meier curves per cohort for the 3 TCGA-prognostic sets (primary endpoint);
  - REMARK profiles.
