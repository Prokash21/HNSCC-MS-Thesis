## Step 08 (part A): final qPCR gene list (22 assays + 2 optional)
##
## Purpose: the final assay list with group, role, source step, reason and the
##   step 05c blood-candidate flag; the input list for the primer step (which
##   genes reuse a step 06 pair and which need a new design).
## Inputs:  results/step06_qPCR/tables/{D_panel_final.csv, primer_best.csv}
## Outputs: results/step08_primers/tables/{A_final_gene_list.csv, A_genes_for_primers.csv}
## Date: 2026-10-28

suppressPackageStartupMessages(library(data.table))
source("config/paths.R")

tdir <- file.path(paths$step08, "tables")
dir.create(tdir, recursive = TRUE, showWarnings = FALSE)

gl <- rbindlist(list(
  data.table(group = "identity", gene = c("LAMC2", "HOXD11", "HOXB7", "CASC9", "IGF2BP3", "LY6K", "SLC2A1", "SLC39A4")),
  data.table(group = "proliferation network", gene = c("E2F1", "MYBL2", "FOXM1", "KIF23", "HMMR", "HMGB3")),
  data.table(group = "hypoxia", gene = c("P4HA1", "CA9")),
  data.table(group = "differentiation / down", gene = "KRT13"),
  data.table(group = "pilot reference", gene = c("PSMC5", "TBP", "GUSB", "SDHA", "ALAS1")),
  data.table(group = "optional", gene = c("NFIC", "ARPC1A"))))
gl[, order_set := fifelse(group == "optional", "optional (order if budget allows)", "core (22 assays)")]

info <- rbindlist(list(
  list("LAMC2", "tumour-cell identity, up; literature positive control", "04b (tier), 06 core",
       "replicated tumour-specific up gene; positive control for the assay set"),
  list("HOXD11", "tumour-cell identity, up", "04b tier 1, 06 core", "tumour-specific, absent in normal epithelium and TME; replicated in TCGA"),
  list("HOXB7", "tumour-cell identity, up", "04b tier 1, 06 core", "tumour-specific, replicated in TCGA"),
  list("CASC9", "tumour-cell identity, up (lncRNA)", "04b tier 1, 06 core", "tumour-specific lncRNA, absent in normal epithelium and TME"),
  list("IGF2BP3", "tumour-cell identity, up", "04b tier 1, 06 core", "tumour-specific, replicated in TCGA"),
  list("LY6K", "tumour-cell identity, up", "04b, 06 core", "tumour-specific, replicated in TCGA"),
  list("SLC2A1", "tumour-cell identity, up (glycolysis / hypoxia link)", "04b, 06 core", "tumour-specific, high expression; mTORC1/PI3K leading edge (05b)"),
  list("SLC39A4", "tumour-cell identity, up", "04b, 06 reserve -> core", "tumour-specific; moved from the step 06 reserve to the final list (user list)"),
  list("E2F1", "proliferation network TF (hub)", "05b, 05c (deviation)", "E2F1 -> MYBL2 in pySCENIC and CollecTRI; co-expressed with MYBL2/FOXM1/KIF23/HMMR"),
  list("MYBL2", "proliferation network TF", "07a regulator shortlist (D3 rank 1)", "targets enriched in the prognostic proliferation sets"),
  list("FOXM1", "proliferation network TF", "05 master regulator, 06 core", "TF activity and expression up in tumour; G2/M target enrichment"),
  list("KIF23", "proliferation, recurrence programme gene", "07a eligible + shortlist", "prognostic G2/M programme gene; E2F1/MYBL2/FOXM1 target"),
  list("HMMR", "proliferation, recurrence programme gene", "07a eligible + shortlist", "prognostic G2/M programme gene; associated with outcome in 07a-07c"),
  list("HMGB3", "proliferation network target", "05 (FOXM1 regulon), 06 core", "tumour-specific up; FOXM1 pySCENIC target (weak co-expression in 05c)"),
  list("P4HA1", "hypoxia, recurrence programme gene", "07a eligible + shortlist (rank 1)", "prognostic hypoxia gene (Buffa); mostly fibroblast-expressed (tumour share 0.16)"),
  list("CA9", "hypoxia marker", "06 optional -> core (deviation)", "largest tumour-vs-normal change; absent in normal, immune cells and blood; co-expressed with P4HA1 (05c)"),
  list("KRT13", "normal-mucosa differentiation marker, down", "06 core (normal control)", "strongly down in tumour; high in normal epithelium"),
  list("PSMC5", "pilot reference gene", "06 pilot reference", "stable tumour vs normal (geNorm/NormFinder, step 06)"),
  list("TBP", "pilot reference gene", "06 pilot reference", "stable tumour vs normal (step 06)"),
  list("GUSB", "pilot reference gene", "06 pilot reference", "stable tumour vs normal (step 06)"),
  list("SDHA", "pilot reference gene", "06 pilot reference", "stable in step 06; low TCGA TPM, UNCLEAR (05c)"),
  list("ALAS1", "pilot reference gene", "06 pilot reference", "least stable in bulk TCGA (log2FC -0.54, 05c); keep in pilot"),
  list("NFIC", "down TF (optional)", "05, 06 reserve -> optional (deviation)", "stroma-confounded; no network role after the 05b correction"),
  list("ARPC1A", "tumour up (optional)", "06 reserve -> optional (deviation)", "small effect; blood and stroma background (05c)")))
setnames(info, c("gene", "role", "source_step", "reason"))
gl <- merge(gl, info, by = "gene", sort = FALSE)
gl[, blood_candidate := gene %in% c("LAMC2", "LY6K", "CASC9", "HOXD11", "CA9")]
stopifnot(nrow(gl) == 24, gl[order_set == "core (22 assays)", .N] == 22)

## primer source: step 06 pair if one exists for the gene, else new design
p6 <- fread(file.path(paths$step06, "tables", "primer_best.csv"))
fin6 <- fread(file.path(paths$step06, "tables", "D_panel_final.csv"))
gl[, primer_source := fifelse(gene %in% fin6[!is.na(forward), gene], "step 06 pair (reuse)", "new design (step 06 settings)")]
stopifnot(all(gl[primer_source == "step 06 pair (reuse)", gene] %in% p6$gene))
gl[, idx := .I]
fwrite(gl[, .(gene, group, role, source_step, reason, blood_candidate, order_set, primer_source)], file.path(tdir, "A_final_gene_list.csv"))
fwrite(gl[, .(gene, role = group, primer_source)], file.path(tdir, "A_genes_for_primers.csv"))
print(gl[, .(gene, group, primer_source, blood_candidate)])
print(gl[, .N, by = .(order_set, primer_source)])
