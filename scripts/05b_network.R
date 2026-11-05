## Step 05b (parts E-F): network check and links to earlier steps
##
## Purpose: edges among E2F1, MYBL2, FOXM1 and NFIC (pySCENIC consensus
##   regulons, CollecTRI); hypergeometric overlap of their targets with
##   proliferation sets; target membership of eligible and qPCR genes;
##   regulators of P4HA1 and HMMR; TF activity vs pathway NES; recurrence
##   programmes in tumour vs normal; the gene story table (pathway -> TF ->
##   outcome).
## Inputs:  data/processed/v2/step05_scenic/regulons_consensus.csv, data/reference/collectri/collectri_omnipath.tsv,
##          data/processed/v2/step05b_gene_sets.rds, results/step04b_DE/tables/D_C2_HPVneg_allpatients.csv,
##          results/step05_TF/tables/{C_contrast_activity_HPVneg.csv, C_TCGA_activity_HPVneg.csv},
##          results/step05b_functional/tables/{B_hpvneg_wide.csv, B_gsea_all.csv},
##          results/step06_qPCR/tables/D_panel_final.csv, results/step07a_recurrence/tables/{D_recurrence_ranked.csv, S_gene_cox.csv},
##          results/step07b_validation/tables/V_gene_cox.csv, results/step07c_oral/tables/O_gene_cox.csv, HGNC complete set
## Outputs: results/step05b_functional/tables/{E_*, F_*}.csv
## Date: 2026-10-28

suppressPackageStartupMessages({
  library(data.table)
})
source("config/paths.R")
source("scripts/lib_de.R")

tdir <- file.path(paths$step05b, "tables")
rd <- function(f) fread(file.path(tdir, f))
fdr_cut <- 0.05
hg <- read_hgnc(paths$hgnc)
core_tfs <- c("E2F1", "MYBL2", "FOXM1", "NFIC")
target_tfs <- c("E2F1", "MYBL2", "FOXM1")
sets <- readRDS(paths$step05b_sets)
universe <- fread(file.path(paths$step04b, "tables", "D_C2_HPVneg_allpatients.csv"))$gene
panel <- fread(file.path(paths$step06, "tables", "D_panel_final.csv"))
qpcr <- panel[set %in% c("core", "reserve", "optional (Arm B)"), .(gene, panel_set = set)]
refs <- panel[set == "reference (pilot)", gene]
eligible <- fread(file.path(paths$step07a, "tables", "D_recurrence_ranked.csv"))$gene

## =========================================================================
## E0. Networks: pySCENIC consensus regulons and CollecTRI (as step 05)
## =========================================================================
scen <- unique(fread(file.path(paths$step05_scenic, "regulons_consensus.csv"))[, .(TF, target, sign = 1, source = "pySCENIC")])
raw <- fread(paths$collectri)
cl <- raw[!grepl("_", source_genesymbol), .(TF = source_genesymbol, target = target_genesymbol,
                                             sign = fifelse(as.logical(is_inhibition) & !as.logical(is_stimulation), -1, 1))]
nm <- hgnc_map(unique(c(cl$TF, cl$target)), hg)
cl[, `:=`(TF = nm$mapped[match(TF, nm$original)], target = nm$mapped[match(target, nm$original)], source = "CollecTRI")]
cl <- unique(cl[!is.na(TF) & !is.na(target)], by = c("TF", "target"))
nets <- rbind(scen, cl)
fwrite(nets[, .(TFs = uniqueN(TF), edges = .N), by = source], file.path(tdir, "E_network_size.csv"))

## =========================================================================
## E1. Edges among E2F1, MYBL2, FOXM1, NFIC
## =========================================================================
e1 <- nets[TF %in% core_tfs & target %in% core_tfs & TF != target]
e1 <- e1[, .(sources = paste(sort(unique(source)), collapse = "+"),
             collectri_sign = if (any(source == "CollecTRI")) sign[source == "CollecTRI"][1] else NA_real_), by = .(TF, target)]
fwrite(e1[order(TF, target)], file.path(tdir, "E_core_tf_edges.csv"))
print(e1)

## =========================================================================
## E2. Target overlap with proliferation sets (hypergeometric, universe = C2 tested genes)
## =========================================================================
ref_sets <- list(HALLMARK_G2M_CHECKPOINT = sets$H$HALLMARK_G2M_CHECKPOINT, HALLMARK_E2F_TARGETS = sets$H$HALLMARK_E2F_TARGETS,
                 Gavish2023_MP2_G1S = sets$RECUR$Gavish2023_MP2_G1S, Gavish2023_MP1_G2M = sets$RECUR$Gavish2023_MP1_G2M,
                 Tirosh2016_G1S = sets$RECUR$Tirosh2016_G1S, Tirosh2016_G2M = sets$RECUR$Tirosh2016_G2M)
N <- length(universe)
e2 <- rbindlist(lapply(c("pySCENIC", "CollecTRI"), function(src) rbindlist(lapply(target_tfs, function(tf) {
  tg <- intersect(nets[source == src & TF == tf, target], universe)
  rbindlist(lapply(names(ref_sets), function(s) {
    sg <- intersect(ref_sets[[s]], universe); k <- intersect(tg, sg)
    data.table(source = src, TF = tf, set = s, targets = length(tg), set_genes = length(sg), overlap = length(k),
               expected = length(tg) * length(sg) / N,
               p = if (length(tg)) phyper(length(k) - 1, length(sg), N - length(sg), length(tg), lower.tail = FALSE) else NA_real_,
               overlap_genes = paste(sort(k), collapse = ";"))
  }))
}))))
e2[, `:=`(fold = overlap / expected, FDR = p.adjust(p, "BH"))]
fwrite(e2, file.path(tdir, "E_target_set_overlap.csv"))
print(e2[, .(source, TF, set, targets, overlap, fold = round(fold, 1), FDR = signif(FDR, 2))])

## =========================================================================
## E3. Eligible and qPCR genes among targets; regulators of P4HA1 and HMMR
## =========================================================================
genes_e3 <- rbind(data.table(gene = eligible, group = "07a eligible"), data.table(gene = qpcr$gene, group = paste("qPCR", qpcr$panel_set)))
e3 <- copy(genes_e3)
for (tf in core_tfs) for (src in c("pySCENIC", "CollecTRI"))
  e3[, (paste(tf, src, sep = "_")) := gene %in% nets[source == src & TF == tf, target]]
e3[, targeted_by := apply(.SD, 1, function(r) paste(names(r)[r], collapse = "; ")), .SDcols = patterns("_(pySCENIC|CollecTRI)$")]
fwrite(e3, file.path(tdir, "E_gene_target_membership.csv"))
act <- fread(file.path(paths$step05, "tables", "C_contrast_activity_HPVneg.csv"))[, .(TF, act_C2 = score_C2, act_C2_FDR = FDR_C2)]
act_t <- fread(file.path(paths$step05, "tables", "C_TCGA_activity_HPVneg.csv"))[, .(TF, act_TCGA, act_TCGA_FDR)]
regs <- nets[target %in% c("P4HA1", "HMMR"), .(sources = paste(sort(unique(source)), collapse = "+"),
                                              collectri_sign = if (any(source == "CollecTRI")) sign[source == "CollecTRI"][1] else NA_real_),
             by = .(gene = target, TF)]
regs <- merge(merge(regs, act, by = "TF", all.x = TRUE), act_t, by = "TF", all.x = TRUE)
setorder(regs, gene, act_C2_FDR)
fwrite(regs, file.path(tdir, "E_regulators_P4HA1_HMMR.csv"))
print(e3[, .(gene, group, targeted_by)])
print(regs[, .(gene, TF, sources, act_C2 = round(act_C2, 2), act_C2_FDR = signif(act_C2_FDR, 2), act_TCGA_FDR = signif(act_TCGA_FDR, 2))])

## =========================================================================
## F1. TF activity (step 05) vs matching pathways (GSEA NES)
## =========================================================================
w <- rd("B_hpvneg_wide.csv")
circ <- act[grepl("^(BMAL1|CLOCK|NR1D[12]|PER[1-3]|CRY[12]|DBP|TEF|HLF)$", TF), TF]
pairs <- rbind(data.table(TF = c("FOXM1", "E2F1", "E2F4", "MYBL2"), match = "HALLMARK_G2M_CHECKPOINT|HALLMARK_E2F_TARGETS"),
               data.table(TF = "TP53", match = "HALLMARK_P53_PATHWAY"),
               data.table(TF = c("HOXD11", "HOXB7"), match = "ANTERIOR_POSTERIOR_PATTERN"),
               data.table(TF = circ, match = "CIRCADIAN"),
               data.table(TF = "NFIC", match = "NFI"))
f1 <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
  m <- w[grepl(pairs$match[i], set) & collection %in% c("H", "GOBP", "REACTOME")]
  if (!nrow(m)) return(data.table(TF = pairs$TF[i], set = "(no matching set tested)"))
  m[, .(TF = pairs$TF[i], collection, set, NES_C2 = NES_C2_HPVneg, FDR_C2 = FDR_C2_HPVneg, NES_TCGA = NES_TCGA_HPVneg,
        FDR_TCGA = FDR_TCGA_HPVneg, NES_C1 = NES_C1_HPVneg, replicated)]
}), fill = TRUE)
f1 <- merge(f1, merge(act, act_t, by = "TF", all = TRUE), by = "TF", all.x = TRUE)
fwrite(f1[order(TF, FDR_C2)], file.path(tdir, "F_tf_vs_pathway.csv"))

## =========================================================================
## F2. Recurrence programmes and matching Hallmark sets, tumour vs normal
## =========================================================================
f2 <- w[collection == "RECUR" | set %in% c("HALLMARK_HYPOXIA", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                                          "HALLMARK_G2M_CHECKPOINT", "HALLMARK_E2F_TARGETS"),
        .(collection, set, NES_C2 = NES_C2_HPVneg, FDR_C2 = FDR_C2_HPVneg, NES_C1 = NES_C1_HPVneg, FDR_C1 = FDR_C1_HPVneg,
          NES_TCGA = NES_TCGA_HPVneg, FDR_TCGA = FDR_TCGA_HPVneg, replicated)]
fwrite(f2, file.path(tdir, "F_recurrence_programmes_tumour_vs_normal.csv"))

## =========================================================================
## F3. Story table: gene -> pathways (leading edge) -> regulating TF -> outcome
## =========================================================================
g_all <- fread(file.path(tdir, "B_gsea_all.csv"))[collection == "H" & FDR < fdr_cut & contrast %in% c("C2_HPVneg", "TCGA_HPVneg")]
le_sets <- function(g, k) {
  x <- g_all[contrast == k & sapply(strsplit(leadingEdge_genes, ";"), function(v) g %in% v)]
  paste(sprintf("%s (%s)", sub("^HALLMARK_", "", x$set), fifelse(x$NES > 0, "+", "-")), collapse = "; ")
}
rep_h <- w[collection == "H" & replicated == TRUE]
sig_tfs <- act[act_C2_FDR < fdr_cut, TF]
s07a <- fread(file.path(paths$step07a, "tables", "S_gene_cox.csv"))[endpoint == "PFI" & model == "multivariable"]
s07c <- fread(file.path(paths$step07c, "tables", "O_gene_cox.csv"))[endpoint == "PFI" & model == "main"]
s07b <- fread(file.path(paths$step07b, "tables", "V_gene_cox.csv"))[cohort == "GSE41613" & endpoint == "DSS" & model == "multivariable"]
fmt <- function(d, g) if (g %in% d$gene) d[gene == g, sprintf("%.2f (p %.2g)", HR[1], p[1])] else "not tested"
story <- rbind(qpcr[, .(gene, block = paste("qPCR", panel_set))], data.table(gene = eligible, block = "07a eligible"),
               data.table(gene = refs, block = "qPCR reference"))
story <- story[!duplicated(gene)]
story[, `:=`(
  hallmark_LE_C2 = sapply(gene, le_sets, k = "C2_HPVneg"),
  hallmark_LE_TCGA = sapply(gene, le_sets, k = "TCGA_HPVneg"),
  replicated_hallmark_LE = sapply(gene, function(g) paste(sub("^HALLMARK_", "", rep_h[sapply(strsplit(LE_shared, ";"), function(v) g %in% v), set]),
                                                          collapse = "; ")),
  regulators_pySCENIC = sapply(gene, function(g) paste(sort(nets[source == "pySCENIC" & target == g, TF]), collapse = ";")),
  regulators_CollecTRI_activeTF = sapply(gene, function(g) paste(sort(intersect(nets[source == "CollecTRI" & target == g, TF], sig_tfs)),
                                                                 collapse = ";")),
  core_tf_targeting = e3$targeted_by[match(gene, e3$gene)],
  TCGA_allsite_PFI_HR = sapply(gene, fmt, d = s07a),
  TCGA_oral_PFI_HR = sapply(gene, fmt, d = s07c),
  GSE41613_DSS_HR = sapply(gene, fmt, d = s07b))]
story[block == "qPCR reference", c("TCGA_allsite_PFI_HR", "TCGA_oral_PFI_HR", "GSE41613_DSS_HR") := "not applicable (reference gene)"]
fwrite(story, file.path(tdir, "F_story_table.csv"))
print(story[, .(gene, block, replicated_hallmark_LE, regulators_pySCENIC, core_tf_targeting, TCGA_allsite_PFI_HR,
                TCGA_oral_PFI_HR, GSE41613_DSS_HR)])
