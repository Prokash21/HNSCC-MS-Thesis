# Step 05b plan: functional analysis (pathways, network check, story)

Written 2026-10-28, before any enrichment was run. Every later change is listed as post hoc in `STEP05b_SUMMARY.md`. There is no primer design in this step. No candidate list, panel or shortlist is changed.

## 0. Context and what was seen before planning

- **Already known:**
  - the step 04b DE results, step 05 TF activity and regulons, step 06 panel and step 07a–07c outcome results;
  - the step 05 TF findings (e.g. FOXM1/E2F/MYBL2 activity up and TP53 activity down in tumour), which motivate parts E and F.
- **Not yet run:** no gene-set enrichment (GSEA, ORA) and no PROGENy analysis has been run on these data in v2.
- **"23 planned qPCR genes."** No list of 23 genes exists in the repository. The step 06 panel (`D_panel_final.csv`) has 25 genes: 20 target genes (10 core, 9 reserve, CA9) and 5 reference genes (PSMC5, TBP, GUSB, SDHA, ALAS1).
  - This step uses the **20 target genes** as "qPCR genes". The reference genes were chosen for stability, not biology, and are listed separately in the story table.
  - The mismatch with "23" is UNCLEAR and is reported.

## 1. Gene sets (one release, read from local files)

- **Collections, all MSigDB v2024.1 human symbol files** in `data/reference/msigdb/`:
  - **Hallmark (primary):** `h.all.v2024.1.Hs.symbols.gmt`, 50 sets (already present).
  - **GO BP (secondary):** `c5.go.bp.v2024.1.Hs.symbols.gmt`, 7,608 sets, downloaded 2026-10-28 from data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/. GO IDs come from the matching `c5.go.bp.v2024.1.Hs.json` (`exactSource`).
  - **Reactome (secondary):** `c2.cp.reactome.v2024.1.Hs.symbols.gmt`, 1,736 sets, same source.
- **Symbols:** gene-set symbols are mapped to current HGNC symbols with the step 04b function `hgnc_map`, the same mapping as the DE tables. Unmapped or ambiguous symbols are dropped and counted.
- **Set size:** after restricting each set to the genes tested in a contrast, only sets with 15–500 genes are used, for all collections.
- **Recurrence programmes:** the 9 step 07a sets (`A_gene_sets.csv`) are tested as a separate custom collection, with their own BH FDR. They are not mixed into Hallmark.
- **Fixed in advance:** no collection, size limit or threshold is changed after results are seen.

## 2. Ranking statistics (A)

- **Statistic:** for each contrast and gene, `sign(logFC) * sqrt(F)`, from the edgeR QL F-test (Squair 2021).
- **HPV- (main):**
  - C2 single-cell pseudobulk, malignant vs normal, `~ study + state` (`D_C2_HPVneg_allpatients.csv`, 11,617 genes);
  - C1, paired patients (`D_C1_HPVneg_paired.csv`, 9,915 genes);
  - TCGA paired tumour vs normal, 40 pairs (`E_TCGA_HPVneg_paired_tumour_vs_normal.csv`).
- **HPV+ (secondary arm, same pipeline):** C2 and C1 HPV+ (`D_C2_HPVpos_allpatients.csv`, `D_C1_HPVpos_paired.csv`).
  - TCGA HPV+ has only 3 pairs (step 04b), so it is not used.
- **Universe:** the genes tested in that contrast.

## 3. GSEA (B)

- **Method:** `fgsea::fgseaMultilevel` (Korotkevich), `eps = 0`, `set.seed(1)`, minSize 15, maxSize 500, on each ranking.
- **Significance:** BH FDR (`padj`) < 0.05 within each collection and contrast.
- **Replicated Hallmark set:** same NES sign and FDR < 0.05 in **both** C2 HPV- and TCGA HPV-. C1 is reported alongside (descriptive; only 7 pairs).
- **Also reported:** GO BP and Reactome, with the same rule for "replicated".
- **Leading edge:**
  - leading-edge genes are reported for every significant set;
  - the 20 qPCR target genes and the 7 step 07a eligible genes are marked inside them.

## 4. ORA (C)

- **Gene lists** (step 05 `A_candidates_final.csv`, HPV-):
  - **up:** tiers 1, 1b and 2 with direction up;
  - **down:** tier 1b down plus `down_replicated`.
- **Universe:** the genes tested in C2 HPV- (11,617).
- **Method:** `clusterProfiler::enricher` (TERM2GENE from the GMT files), GO BP and Hallmark, minGSSize 15, maxGSSize 500, BH. Significant: `p.adjust` < 0.05.
- **GO redundancy:** `simplifyEnrichment` (Gu 2022). GO semantic similarity of the significant GO BP terms (`GO_similarity`, ont BP, default measure), then `simplifyGO` (default binary cut). Done separately for up and down terms, and for the replicated GSEA GO BP terms (up and down).
  - If fewer than 5 significant terms are available, clustering is skipped and this is stated.
  - rrvgo is not used.

## 5. PROGENy pathway activity (D)

- **Model:** `progeny::model_human_full`, top 500 genes per pathway by p-value, weight = coefficient. This is the same network that `decoupleR::get_progeny(top = 500)` builds.
- **Scoring:** activity with `decoupleR::run_mlm` on gene z-scores.
- **Single cell, HPV- (and HPV+ secondary):**
  - pseudobulks as in step 05 C2: malignant with ≥ 20 cells, normal with ≥ 10 cells, P9 normal excluded;
  - counts from the `shared` pseudobulk matrix, `filterByExpr`, TMM, log CPM (prior 2), z per gene;
  - test: `lm(activity ~ study + state)` over all pseudobulks (main), plus a paired t-test on the patients with both a malignant and a normal pseudobulk.
- **TCGA, HPV- 40 pairs:**
  - step 04b counts (`tcga_hnsc_paired_counts.rds`), `filterByExpr`, TMM, log CPM (prior 2), z per gene;
  - test: paired t-test (tumour minus normal) per pathway.
- **Multiple testing:** BH FDR over the 14 pathways per test. Significant: FDR < 0.05.
- **Replicated pathway:** same direction and FDR < 0.05 in the single-cell `lm` test and in the TCGA paired test.

## 6. Network check (E), pre-specified, descriptive

- **Sources:**
  - step 05 consensus pySCENIC regulons (`regulons_consensus.csv`; activating regulons only);
  - CollecTRI (`collectri_omnipath.tsv`, HGNC-mapped as in step 05, complexes excluded, sign kept).
- **E1. Edges among E2F1, MYBL2, FOXM1, NFIC:** every directed edge TF → TF, with its source (pySCENIC, CollecTRI or both) and the CollecTRI sign.
- **E2. Target overlaps:** targets of E2F1, MYBL2 and FOXM1, per source, compared with:
  - Hallmark G2M checkpoint and E2F targets;
  - the step 07a Gavish G1/S and G2/M sets (Tirosh G1/S and G2/M as concordance).

  The test is a one-sided hypergeometric test. The universe is the genes tested in C2 HPV- (11,617); targets and sets are restricted to it. BH FDR is applied over all TF × set × source tests.
- **E3. Membership:**
  - which of the 7 eligible step 07a genes (P4HA1, HMMR, KIF23, SEC61G, CCNA2, ECT2, CKAP2) and the 20 qPCR target genes are targets of E2F1, MYBL2, FOXM1 or NFIC;
  - all regulators of P4HA1 and HMMR in either source.
- **Caveat:** regulons are inferred (co-expression plus motif enrichment, pySCENIC) or curated from literature and text mining (CollecTRI). They are not experimentally proven in these samples.

## 7. Link to earlier steps (F)

1. **TF activity vs pathways** (step 05 C2 and TCGA activity vs NES):
   - FOXM1, E2F1, E2F4 and MYBL2 vs Hallmark G2M checkpoint and E2F targets;
   - TP53 activity vs Hallmark P53 pathway;
   - HOX TFs (HOXD11, HOXB7) vs GO BP sets whose name contains "ANTERIOR_POSTERIOR_PATTERN";
   - circadian TFs (BMAL1 and any CLOCK/NR1D/PER/CRY/DBP/TEF/HLF TF in step 05) vs GO BP and Reactome sets whose name contains "CIRCADIAN";
   - NFIC: no matching set is expected; any set with "NFI" in its name is reported.
2. **Recurrence programmes:** the 9 step 07a sets in tumour vs normal (GSEA on C2, C1, TCGA), together with the Hallmark hypoxia, EMT and G2M/E2F sets.
3. **Story table:** one row per qPCR target gene (20; reference genes listed separately), with:
   - the significant Hallmark sets where the gene is in the leading edge (C2 and TCGA);
   - regulating TFs: all pySCENIC consensus regulators, plus CollecTRI regulators restricted to TFs with step 05 C2 activity FDR < 0.05;
   - outcome HR and p from step 07a (TCGA all-site PFI, multivariable), step 07c (TCGA oral PFI, stage) and step 07b (GSE41613 DSS, multivariable).

   The 7 eligible step 07a genes are added in a second block.

## 8. Outputs

- **Scripts:** `scripts/05b_gsea.R`, `05b_ora.R`, `05b_progeny.R`, `05b_network.R`, `05b_figures.R`.
- **Tables:** `results/step05b_functional/tables/`.
- **`STEP05b_SUMMARY.md`** (local, git-ignored):
  - main findings in plain words;
  - replicated pathways;
  - network check;
  - story table;
  - limitations: bulk TCGA includes the microenvironment; 7 normal patients in single cell; inferred regulons;
  - UNCLEAR items;
  - references.
- **Figures:**
  - Hallmark dot plot (C2, C1, TCGA side by side);
  - Venn of HPV- up and down genes (single-cell C2 candidates vs TCGA paired FDR < 0.05, same direction);
  - GO cluster heatmap;
  - PROGENy heatmap;
  - network diagram (E2F1/MYBL2/FOXM1/NFIC → targets, eligible and qPCR genes highlighted);
  - story diagram (pathways → TFs → panel genes).
