## Step 07a (part E): recurrence candidate selection (pre-specified, PLAN.md section 6)
##
## Purpose: D1 eligibility of every primary-set gene (programmes 1-3), D2
##   ranking of eligible genes and a suggested shortlist, D3 regulators whose
##   targets are enriched in a prognostic programme, D4 the candidate file for
##   the final primer step.
## Inputs:  results/step07a_recurrence/tables/{A_gene_sets.csv, B_gene_annotation.csv, C_bill_DE.csv,
##          S_programme_cox.csv, S_gene_cox.csv}
##          data/processed/v2/step05_scenic/regulons_consensus.csv, data/processed/v2/step05_export/genes.csv,
##          data/reference/collectri/collectri_omnipath.tsv, data/processed/v2/pseudobulk_step04b.rds,
##          results/step04b_DE/tables/D_pseudobulk_samples.csv, HGNC complete set
## Outputs: results/step07a_recurrence/tables/D_*.csv, results/step07a_recurrence/candidates_for_step08.csv
## Date: 2026-10-28

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})
source("config/paths.R")
source("scripts/lib_de.R")

out <- paths$step07a
tdir <- file.path(out, "tables")
rd <- function(f) fread(file.path(tdir, f))
fdr_prog <- 0.10
p_gene   <- 0.05
min_cpm  <- 10
min_tpm  <- 10
fdr_tf   <- 0.05
max_per_programme <- 2
shortlist_size    <- 4
max_tf_short      <- 2
min_mal  <- 20
hg <- read_hgnc(paths$hgnc)

gs <- rd("A_gene_sets.csv")[role == "primary"]
ann <- rd("B_gene_annotation.csv")
bill <- rd("C_bill_DE.csv")
pc <- rd("S_programme_cox.csv")[model == "multivariable" & endpoint == "PFI" & role == "primary"]
gc_ <- rd("S_gene_cox.csv")[model == "multivariable" & endpoint == "PFI"]

## =========================================================================
## D1. Eligibility
## =========================================================================
pc[, prognostic := HR > 1 & FDR < fdr_prog]
pass_sets <- pc[prognostic == TRUE, set]
message("D1. prognostic primary sets (PFI, multivariable HR > 1, FDR < 0.10): ", paste(pass_sets, collapse = ", "))
el <- gs[, .(sets = paste(set, collapse = ";"), programmes = paste(sort(unique(programme)), collapse = ";"),
             prognostic_sets = paste(set[set %in% pass_sets], collapse = ";"),
             prognostic_programmes = paste(sort(unique(programme[set %in% pass_sets])), collapse = ";")), by = gene]
el <- merge(el, gc_[, .(gene, PFI_HR = HR, PFI_lo = lo, PFI_hi = hi, PFI_logHR = logHR, PFI_p = p, PFI_FDR = FDR,
                        PFI_ph_p = ph_p_term)], by = "gene", all.x = TRUE)
os_ <- rd("S_gene_cox.csv")[model == "multivariable" & endpoint == "OS", .(gene, OS_HR = HR, OS_p = p, OS_FDR = FDR)]
el <- merge(el, os_, by = "gene", all.x = TRUE)
el <- merge(el, ann[, .(gene, malignant_pb_CPM, tumour_share, max_TME_type, tcga_median_TPM, panel_set)], by = "gene", all.x = TRUE)
el <- merge(el, bill[, .(gene, bill_log2FC = logFC, bill_FDR = FDR)], by = "gene", all.x = TRUE)
el[, `:=`(c1_programme_prognostic = prognostic_sets != "",
          c2_gene_PFI = (PFI_HR > 1 & PFI_p < p_gene) %in% TRUE,
          c3_malignant_CPM = (malignant_pb_CPM >= min_cpm) %in% TRUE,
          c4_bill_same_direction = (bill_log2FC > 0) %in% TRUE,
          c5_tcga_TPM = (tcga_median_TPM >= min_tpm) %in% TRUE)]
el[, eligible := c1_programme_prognostic & c2_gene_PFI & c3_malignant_CPM & c4_bill_same_direction & c5_tcga_TPM]
fails <- function(c2, c3, c4, c5, fc) paste(c(if (!c2) "gene PFI", if (!c3) "malignant CPM",
                                               if (!c4) (if (is.na(fc)) "Bill not evaluable" else "Bill direction"),
                                               if (!c5) "TCGA TPM"), collapse = "; ")
el[, status := fifelse(!c1_programme_prognostic, "programme not prognostic",
                       fifelse(eligible, "eligible", paste("fails:", fails(c2_gene_PFI, c3_malignant_CPM, c4_bill_same_direction,
                                                                           c5_tcga_TPM, bill_log2FC)))), by = gene]
el[, in_panel := !is.na(panel_set)]
setorder(el, -eligible, -c1_programme_prognostic, PFI_p)
fwrite(el, file.path(tdir, "D_recurrence_eligible_all.csv"))
crit <- el[c1_programme_prognostic == TRUE, .(genes_in_prognostic_sets = .N, c2 = sum(c2_gene_PFI), c3 = sum(c3_malignant_CPM),
                                              c4 = sum(c4_bill_same_direction), c5 = sum(c5_tcga_TPM),
                                              c2_c3 = sum(c2_gene_PFI & c3_malignant_CPM),
                                              c2_c3_c5 = sum(c2_gene_PFI & c3_malignant_CPM & c5_tcga_TPM),
                                              eligible = sum(eligible))]
fwrite(crit, file.path(tdir, "D_criteria_counts.csv"))
print(crit)

## =========================================================================
## D2. Ranking: mean percentile rank (higher = better) and suggested shortlist
## =========================================================================
rk <- el[eligible == TRUE]
pct <- function(v) (frank(v, ties.method = "average") - 1) / max(1, length(v) - 1)
if (nrow(rk)) {
  rk[, `:=`(pr_logHR = pct(PFI_logHR), pr_bill = pct(bill_log2FC), pr_share = pct(tumour_share),
            pr_TPM = pct(log(tcga_median_TPM)))]
  rk[, rank_score := (pr_logHR + pr_bill + pr_share + pr_TPM) / 4]
  setorder(rk, -rank_score)
  rk[, rank := seq_len(.N)]
  used <- c(); rk[, shortlist := FALSE]
  for (i in seq_len(nrow(rk))) {
    if (sum(rk$shortlist) >= shortlist_size) break
    pr <- strsplit(rk$prognostic_programmes[i], ";")[[1]]
    if (all(vapply(pr, function(p) sum(used == p) < max_per_programme, TRUE))) {
      rk$shortlist[i] <- TRUE; used <- c(used, pr)
    }
  }
} else rk[, `:=`(rank_score = numeric(), rank = integer(), shortlist = logical())]
fwrite(rk, file.path(tdir, "D_recurrence_ranked.csv"))
message(sprintf("D2. eligible genes: %d; shortlist: %s", nrow(rk), paste(rk[shortlist == TRUE, gene], collapse = ", ")))

## =========================================================================
## D3. Regulators: TF target enrichment in prognostic primary sets
## =========================================================================
universe <- fread(file.path(paths$step05_export, "genes.csv"))$gene
reg <- fread(file.path(paths$step05_scenic, "regulons_consensus.csv"))[, .(TF, target)]
raw <- fread(paths$collectri)
cl <- raw[!grepl("_", source_genesymbol), .(TF = source_genesymbol, target = target_genesymbol)]
nm <- hgnc_map(unique(c(cl$TF, cl$target)), hg)
cl[, `:=`(TF = nm$mapped[match(TF, nm$original)], target = nm$mapped[match(target, nm$original)])]
nets <- list(pySCENIC_consensus = unique(reg), CollecTRI = unique(cl[!is.na(TF) & !is.na(target)]))
N <- length(universe)
enr <- rbindlist(lapply(names(nets), function(src) {
  net <- nets[[src]][target %in% universe]
  rbindlist(lapply(pass_sets, function(s) {
    sg <- intersect(gs[set == s, gene], universe)
    net[, {
      tg <- unique(target); k <- length(intersect(tg, sg))
      .(source = src, set = s, n_targets = length(tg), set_genes = length(sg), overlap = k,
        overlap_genes = paste(sort(intersect(tg, sg)), collapse = ";"),
        p = phyper(k - 1, length(tg), N - length(tg), length(sg), lower.tail = FALSE))
    }, by = TF]
  }))
}))
if (nrow(enr)) enr[, FDR := p.adjust(p, "BH"), by = source]
## TF expression in malignant cells: step 04b HPV- malignant pseudobulks (>= 20 cells), median CPM
pb <- readRDS(paths$pseudobulk_v2b)
smp4 <- fread(file.path(paths$step04b, "tables", "D_pseudobulk_samples.csv"))
mal_k <- smp4[group == "malignant" & hpv == "HPV-" & primary_4b == TRUE & n_cells >= min_mal]
m <- pb$full$counts
tfs <- intersect(unique(enr$TF), rownames(m))
tf_cpm <- apply(as.matrix(m[tfs, mal_k$key, drop = FALSE]) %*% diag(1e6 / colSums(m[, mal_k$key, drop = FALSE]), nrow = nrow(mal_k)),
                1, median)
enr[, TF_malignant_pb_CPM := tf_cpm[TF]]
enr[, eligible := FDR < fdr_tf & (TF_malignant_pb_CPM >= min_cpm) %in% TRUE]
setorder(enr, FDR, p)
fwrite(enr[p < 0.05 | eligible], file.path(tdir, "D_regulator_enrichment.csv"))
tf_rank <- enr[eligible == TRUE, .(best_FDR = min(FDR), sources = paste(sort(unique(source)), collapse = ";"),
                                   sets = paste(sort(unique(set)), collapse = ";"), TF_malignant_pb_CPM = TF_malignant_pb_CPM[1],
                                   best_overlap = overlap[which.min(FDR)]), by = TF][order(best_FDR)]
tf_rank[, `:=`(rank = seq_len(.N), shortlist = seq_len(.N) <= max_tf_short)]
fwrite(tf_rank, file.path(tdir, "D_regulators_ranked.csv"))
message(sprintf("D3. eligible regulators: %d; shortlist: %s", nrow(tf_rank), paste(tf_rank[shortlist == TRUE, TF], collapse = ", ")))
print(tf_rank)

## =========================================================================
## D4. Candidate file for the final primer step
## =========================================================================
cand <- rbind(rk[, .(gene, role = "programme gene", source = prognostic_sets, rank, shortlist = fifelse(shortlist, "yes", "no"))],
              tf_rank[, .(gene = TF, role = "regulator", source = paste(sources, sets, sep = ": "), rank,
                          shortlist = fifelse(shortlist, "yes", "no"))])
fwrite(cand, file.path(out, "candidates_for_step08.csv"))
print(cand)
