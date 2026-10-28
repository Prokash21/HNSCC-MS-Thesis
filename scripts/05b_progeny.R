## Step 05b (part D): PROGENy pathway activity
##
## Purpose: PROGENy pathway activity (top 500 genes per pathway, decoupleR
##   multivariate linear model) per single-cell pseudobulk (malignant vs
##   normal; lm ~ study + state and paired t-test, HPV- and HPV+) and per TCGA
##   HPV- tumour/normal sample (paired t-test over 40 pairs).
## Inputs:  data/processed/v2/pseudobulk_step04b.rds, results/step04b_DE/tables/D_pseudobulk_samples.csv,
##          data/processed/v2/tcga_hnsc_paired_counts.rds, data/raw/step01_metadata/cbioportal_hnsc_tcga_pan_can_atlas_2018_patient.csv,
##          progeny::model_human_full, HGNC complete set
## Outputs: results/step05b_functional/tables/D_*.csv
## Date: 2026-10-28

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(edgeR)
  library(decoupleR)
})
source("config/paths.R")
source("scripts/lib_de.R")

tdir <- file.path(paths$step05b, "tables")
fdr_cut <- 0.05
top_n <- 500
min_mal <- 20
min_norm <- 10
excl_normal <- "P9"
hg <- read_hgnc(paths$hgnc)

## =========================================================================
## D1. PROGENy network: top 500 genes per pathway by p-value (as get_progeny)
## =========================================================================
data(model_human_full, package = "progeny")
net <- as.data.table(model_human_full)[, `:=`(gene = as.character(gene), pathway = as.character(pathway))]
net <- net[order(p.value), head(.SD, top_n), by = pathway]
nm <- hgnc_map(unique(net$gene), hg)
net[, target := nm$mapped[match(gene, nm$original)]]
net <- unique(net[!is.na(target), .(source = pathway, target, weight)], by = c("source", "target"))
fwrite(net[, .N, by = source], file.path(tdir, "D_progeny_network_size.csv"))
mlm <- function(z) {
  a <- as.data.table(run_mlm(mat = z, net = net, .source = "source", .target = "target", .mor = "weight", minsize = 5))
  m <- dcast(a, source ~ condition, value.var = "score")
  out <- as.matrix(m[, -1]); rownames(out) <- m$source
  out[, colnames(z), drop = FALSE]
}
zlog <- function(counts, group) {
  y <- DGEList(as.matrix(counts))
  y <- normLibSizes(y[filterByExpr(y, group = group), , keep.lib.sizes = FALSE])
  t(scale(t(edgeR::cpm(y, log = TRUE, prior.count = 2))))
}

## =========================================================================
## D2. Single cell: pseudobulks as step 05 C2 (patient level)
## =========================================================================
pb <- readRDS(paths$pseudobulk_v2b)
smp <- fread(file.path(paths$step04b, "tables", "D_pseudobulk_samples.csv"))
sc_res <- list(); sc_act <- list(); sc_design <- list()
for (arm in c("HPV-", "HPV+")) {
  ps <- smp[hpv == arm & primary_4b == TRUE &
              ((group == "malignant" & n_cells >= min_mal) | (group == "normal" & n_cells >= min_norm & patient != excl_normal))]
  ps[, state := factor(fifelse(group == "malignant", "malignant", "normal"), levels = c("normal", "malignant"))]
  z <- zlog(pb$shared$counts[, ps$key], ps$state)
  act <- mlm(z)
  pairs <- intersect(ps[state == "malignant", patient], ps[state == "normal", patient])
  pk_m <- ps[state == "malignant"][match(pairs, patient), key]
  pk_n <- ps[state == "normal"][match(pairs, patient), key]
  study <- factor(ps$study)
  r <- rbindlist(lapply(rownames(act), function(p) {
    a <- act[p, ]
    d <- a[pk_m] - a[pk_n]
    co <- summary(lm(a ~ study + ps$state))$coefficients
    data.table(pathway = p, all_effect = co["ps$statemalignant", "Estimate"], all_p = co["ps$statemalignant", "Pr(>|t|)"],
               paired_diff = mean(d), paired_p = t.test(d)$p.value)
  }))
  r[, `:=`(arm = arm, all_FDR = p.adjust(all_p, "BH"), paired_FDR = p.adjust(paired_p, "BH"))]
  sc_res[[arm]] <- r
  sc_act[[arm]] <- melt(as.data.table(act, keep.rownames = "pathway"), id.vars = "pathway", variable.name = "key",
                        value.name = "activity")[, `:=`(arm = arm, state = ps$state[match(key, ps$key)],
                                                        patient = ps$patient[match(key, ps$key)])]
  sc_design[[arm]] <- data.table(arm = arm, malignant = sum(ps$state == "malignant"), normal = sum(ps$state == "normal"),
                                 pairs = length(pairs), pair_patients = paste(pairs, collapse = " "), genes = nrow(z))
}
rm(pb); invisible(gc())

## =========================================================================
## D3. TCGA HPV- 40 pairs: paired t-test tumour - normal
## =========================================================================
tc <- readRDS(paths$tcga_counts)
cb <- fread(paths$cbio_pancan, select = c("patientId", "SUBTYPE"))
ts <- copy(tc$samples)[, `:=`(hpv = cb$SUBTYPE[match(case, cb$patientId)],
                              state = factor(fifelse(sample_type == "Primary Tumor", "tumour", "normal"), c("normal", "tumour")))]
ts <- ts[hpv == "HNSC_HPV-"]
stopifnot(uniqueN(ts$case) == 40, all(ts[, .N, by = case]$N == 2))
zt <- zlog(tc$counts[, ts$sample], ts$state)
at <- mlm(zt)
tum <- ts[state == "tumour"][order(case), sample]; nor <- ts[state == "normal"][order(case), sample]
tr <- rbindlist(lapply(rownames(at), function(p) {
  d <- at[p, tum] - at[p, nor]
  data.table(pathway = p, tcga_paired_diff = mean(d), tcga_paired_p = t.test(d)$p.value)
}))
tr[, tcga_paired_FDR := p.adjust(tcga_paired_p, "BH")]
tcga_act <- melt(as.data.table(at, keep.rownames = "pathway"), id.vars = "pathway", variable.name = "sample",
                 value.name = "activity")[, `:=`(state = ts$state[match(sample, ts$sample)], case = ts$case[match(sample, ts$sample)])]

## =========================================================================
## D4. Tables
## =========================================================================
sc_all <- rbindlist(sc_res)
fwrite(sc_all, file.path(tdir, "D_progeny_single_cell.csv"))
fwrite(rbindlist(sc_design), file.path(tdir, "D_progeny_design.csv"))
fwrite(rbindlist(sc_act), file.path(tdir, "D_progeny_activity_pseudobulk.csv"))
fwrite(tcga_act, file.path(tdir, "D_progeny_activity_tcga.csv"))
neg <- merge(sc_all[arm == "HPV-", .(pathway, sc_effect = all_effect, sc_p = all_p, sc_FDR = all_FDR,
                                     sc_paired_diff = paired_diff, sc_paired_FDR = paired_FDR)], tr, by = "pathway")
neg[, replicated := sc_FDR < fdr_cut & tcga_paired_FDR < fdr_cut & sign(sc_effect) == sign(tcga_paired_diff)]
neg[, direction := fifelse(replicated, fifelse(sc_effect > 0, "up in tumour", "down in tumour"), "")]
setorder(neg, -replicated, sc_p)
fwrite(neg, file.path(tdir, "D_progeny_hpvneg_summary.csv"))
print(rbindlist(sc_design))
print(neg[, .(pathway, sc = round(sc_effect, 2), sc_FDR = signif(sc_FDR, 2), sc_pair = round(sc_paired_diff, 2),
              sc_pair_FDR = signif(sc_paired_FDR, 2), tcga = round(tcga_paired_diff, 2), tcga_FDR = signif(tcga_paired_FDR, 2),
              replicated)])
print(sc_all[arm == "HPV+", .(pathway, eff = round(all_effect, 2), FDR = signif(all_FDR, 2))][order(FDR)])
