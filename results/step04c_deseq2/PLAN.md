# Step 04c plan: DESeq2 sensitivity analysis of the step 04b pseudobulk DE

Written 2026-10-28, before DESeq2 was run. This step is **descriptive only**: no candidate list, tier, panel or result of steps 04b–08 is changed. Later changes are listed as post hoc in `STEP04c_SUMMARY.md`.

## 1. Inputs (exactly as step 04b)

- **Single cell:**
  - `pseudobulk_step04b.rds` (`shared` count matrix) and the step 04b sample table (`D_pseudobulk_samples.csv`);
  - the step 04b `make_set` logic: primary patients only; malignant ≥ 20 cells; normal ≥ 10 cells; P9 normal excluded;
  - C1 = patients with both states (paired); C2 = all patients.
  - `study_zero` is taken from the `full` matrix as in step 04b.
- **TCGA:** `tcga_hnsc_paired_counts.rds`, HPV- pairs (cBioPortal subtype), 40 cases.
- **Arms:** HPV- (main) and HPV+ (secondary) for C1 and C2. TCGA is HPV- only (HPV+ has 3 pairs).
- **Reconstruction check:** before DESeq2, edgeR QL is re-run on the rebuilt sets with the step 04b function `de_qlf`. The log2FC and FDR must equal the saved step 04b tables, which proves the inputs are identical.

## 2. DESeq2 (A)

- **Genes:** the genes kept by step 04b `filterByExpr` (the genes in each saved edgeR table), so both methods test the same genes.
- **Fit:** DESeq2 defaults otherwise: median-of-ratios size factors on these genes, parametric dispersion fit, Wald test, BH (`padj`).
- **Independent filtering and Cook's distance:** kept at the defaults (independent filtering on; Cook's outlier handling and replacement for groups with ≥ 7 replicates). The report states how many genes get `padj = NA` from independent filtering and how many get `pvalue = NA` as Cook's outliers, and how many were replaced. For comparisons, `padj = NA` counts as not significant.
- **Designs** (state reference = normal):
  - C2: `~ study + state`;
  - C1: `~ patient + state`;
  - TCGA: `~ patient + state`.
- **Reported:**
  - unshrunken log2FC (MLE), for comparison with edgeR;
  - apeglm-shrunken log2FC (`lfcShrink(type = "apeglm")`), used for plots only.

## 3. Pre-specified comparison with edgeR QL (B)

1. **log2FC agreement:** Pearson and Spearman r between edgeR logFC and DESeq2 unshrunken log2FC, over all tested genes and over the step 04b candidate genes, for C2, C1 and TCGA (HPV-). HPV+ C2 and C1 are reported as secondary. Scatter plots (Fig51).
2. **Significance agreement:** genes with FDR < 0.05 (edgeR `FDR`, DESeq2 `padj`), separately up and down (by the sign of each method's own log2FC). Counts in edgeR only, both and DESeq2 only; Jaccard index; Venn diagrams (Fig52).
   - This uses FDR < 0.05 alone, without a fold-change threshold.
3. **Candidate rule re-applied** (HPV-):
   - **rule:** C2 significant (FDR < 0.05 and |log2FC| ≥ 1, unshrunken) and C1 log2FC of the same sign, excluding `study_zero` genes. These are the step 04b thresholds.
   - **compared with** the step 04b `candidate_list` (142 up, 73 down);
   - **reported:** numbers recovered, and genes gained and lost by name, each with its edgeR and DESeq2 numbers.
4. **The 24 atlas genes** (step 08 list: 17 targets, NFIC, ARPC1A, 5 references):
   - edgeR vs DESeq2 log2FC and FDR side by side for C2, C1 and TCGA (Fig53);
   - flag any gene whose direction or significance (FDR < 0.05) changes between methods.
5. **Reference genes** (PSMC5, TBP, GUSB, SDHA, ALAS1): DESeq2 log2FC is expected near 0; values are reported against the edgeR values.

For genes that differ, the likely reason is checked and reported: low counts (mean normalised count), Cook's outlier or replacement, and independent filtering.

## 4. Outputs

- **Scripts:** `scripts/04c_deseq2.R` (fits and tables), `scripts/04c_figures.R`.
- **Tables:** `results/step04c_deseq2/tables/`.
- **`STEP04c_SUMMARY.md`** (local):
  - verdict;
  - B1–B5 numbers;
  - differing genes and reasons;
  - limitations;
  - references.
- **Figures:**
  - Fig51 log2FC scatter, edgeR vs DESeq2 (C2, C1, TCGA);
  - Fig52 Venn of significant genes;
  - Fig53 the 24 genes side by side.

  Fig51 is already used by step 08 (`Fig51_primer_map_new`). The numbers requested for this step are kept, in this step's own folder.
