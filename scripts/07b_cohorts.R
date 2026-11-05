## Step 07b (part A): GEO validation cohorts, probes to genes, programme scores
##
## Purpose: parse GSE65858 and GSE41613 series matrices, apply the pre-specified
##   HPV- and primary-tumour filters, define endpoints and covariates, collapse
##   probes to genes (maxMean), report gene-set coverage and compute programme
##   scores (mean z of log expression) per cohort.
## Inputs:  data/raw/step07_bulk/{GSE65858,GSE41613}_series_matrix.txt.gz, {GPL10558,GPL570}.annot.gz,
##          HGNC complete set, results/step07a_recurrence/tables/A_gene_sets.csv
## Outputs: results/step07b_validation/tables/A_*.csv
##          data/processed/v2/step07b_geo.rds (gene-level expression, clinical tables, scores)
## Date: 2026-10-28

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})
source("config/paths.R")
source("scripts/lib_de.R")

tdir <- file.path(paths$step07b, "tables")
dir.create(tdir, recursive = TRUE, showWarnings = FALSE)
hg <- read_hgnc(paths$hgnc)
gs <- fread(file.path(paths$step07a, "tables", "A_gene_sets.csv"))

## =========================================================================
## A1. Series matrix: characteristics and expression
## =========================================================================
read_series <- function(f) {
  stopifnot(file.exists(f))
  l <- readLines(gzfile(f))
  ch <- grep("^!Sample_characteristics_ch1", l, value = TRUE)
  v <- lapply(strsplit(ch, "\t"), function(x) gsub('"', "", x[-1]))
  key <- sapply(v, function(x) sub(":.*", "", x[which(x != "")[1]]))
  cl <- as.data.table(setNames(lapply(v, function(x) sub("^[^:]*: ?", "", x)), make.unique(key)))
  gsm <- gsub('"', "", strsplit(grep("^!Sample_geo_accession", l, value = TRUE), "\t")[[1]][-1])
  cl[, gsm := gsm]
  b <- grep("^!series_matrix_table_(begin|end)", l)
  ex <- fread(text = l[(b[1] + 1):(b[2] - 1)], quote = '"')
  m <- as.matrix(ex[, -1]); rownames(m) <- as.character(ex[[1]])
  stopifnot(identical(colnames(m), gsm))
  list(clin = cl, expr = m)
}

## =========================================================================
## A2. Probes to genes: GEO annotation symbol, HGNC, maxMean (Miller 2011)
## =========================================================================
collapse_max_mean <- function(m, annot_file, platform) {
  stopifnot(file.exists(annot_file))
  an <- fread(cmd = paste("gzip -dc", shQuote(annot_file)), skip = "ID\tGene title", select = c("ID", "Gene symbol"),
              colClasses = "character")
  setnames(an, c("probe", "symbol"))
  an <- an[probe %in% rownames(m)]
  st <- data.table(platform = platform, probes_in_matrix = nrow(m), probes_annotated = nrow(an),
                   probes_empty_symbol = an[symbol == "", .N], probes_multi_symbol = an[grepl("///", symbol), .N])
  an <- an[symbol != "" & !grepl("///", symbol)]
  mp <- hgnc_map(unique(an$symbol), hg)
  an[, `:=`(gene = mp$mapped[match(symbol, mp$original)], how = mp$how[match(symbol, mp$original)])]
  st[, `:=`(probes_symbol_not_mapped = an[is.na(gene) | how == "not_in_hgnc", .N])]
  an <- an[!is.na(gene) & how != "not_in_hgnc"]
  an[, mean_expr := rowMeans(m[probe, , drop = FALSE])]
  pick <- an[order(-mean_expr), .SD[1], by = gene]
  st[, `:=`(probes_used = nrow(an), genes = nrow(pick))]
  g <- m[pick$probe, , drop = FALSE]; rownames(g) <- pick$gene
  list(expr = g, probes = pick[, .(gene, probe, symbol, how, mean_expr)], stats = st)
}

## =========================================================================
## A3. Cohorts: filters, endpoints, covariates
## =========================================================================
stage_grp <- function(s) factor(fcase(s %in% c("I", "II"), "I-II", s == "III", "III", grepl("^IV", s), "IV"),
                                c("I-II", "III", "IV"))
num <- function(v) suppressWarnings(as.numeric(v))

## GSE65858 (Wichmann 2015)
a <- read_series(paths$gse65858)
c1 <- a$clin
f1 <- data.table(step = c("GSE65858 series samples", "HPV16 DNA-negative (hpv16_dna_rna = DNA-)",
                          "  and primary tumour (tumor_type = Primary)"),
                 n = c(nrow(c1), c1[hpv16_dna_rna == "DNA-", .N], c1[hpv16_dna_rna == "DNA-" & tumor_type == "Primary", .N]))
c1 <- c1[hpv16_dna_rna == "DNA-" & tumor_type == "Primary"]
stopifnot(c1[hpv_dna != "Negative", .N] == 0)
c1 <- c1[, .(cohort = "GSE65858", gsm, patient = ID,
             PFS.time = num(pfs), PFS = as.integer(pfs_event == "TRUE"),
             OS.time = num(os), OS = as.integer(os_event == "TRUE"),
             age = num(age), stage = stage_grp(uicc_stage), stage_raw = uicc_stage,
             site = factor(fcase(tumor_site == "Cavum Oris", "oral cavity", tumor_site == "Larynx", "larynx",
                                 tumor_site %in% c("Oropharynx", "Hypopharynx"), "other"),
                           c("oral cavity", "larynx", "other")),
             site_raw = tumor_site, sex = gender, treatment)]
stopifnot(all(c1$PFS %in% 0:1), all(c1$OS %in% 0:1))
x1 <- collapse_max_mean(a$expr, paths$gpl10558, "GPL10558")

## GSE41613 (Lohavanichbutr 2013)
b <- read_series(paths$gse41613)
c2 <- b$clin
f2 <- data.table(step = c("GSE41613 series samples", "HPV-negative OSCC (tissue type)"),
                 n = c(nrow(c2), c2[grepl("^HPV-negative", `tissue type`), .N]))
c2 <- c2[grepl("^HPV-negative", `tissue type`)]
age_lv <- c("19-39", "40-49", "50-59", "60-88")
stopifnot(all(c2$age %in% age_lv), all(c2$vital %in% c("Alive", "Dead-oral ca", "Dead-non OC", "Dead -unk cause")))
c2 <- c2[, .(cohort = "GSE41613", gsm, patient = gsm,
             DSS.time = num(`fu time`), DSS = as.integer(vital == "Dead-oral ca"),
             OS.time = num(`fu time`), OS = as.integer(grepl("^Dead", vital)),
             DSS_unk.time = num(`fu time`), DSS_unk = as.integer(vital %in% c("Dead-oral ca", "Dead -unk cause")),
             age = match(age, age_lv), age_group = age,
             stage = factor(`tumor stage`, c("I/II", "III/IV")), stage_raw = `tumor stage`,
             site = factor("oral cavity"), site_raw = "oral cavity", sex = Sex, treatment, vital)]
## log2 check: gcRMA output is log2; transform only if the scale is clearly linear
ex2 <- b$expr
log2_applied <- max(ex2, na.rm = TRUE) > 100
if (log2_applied) ex2 <- log2(pmax(ex2, 0) + 1)
x2 <- collapse_max_mean(ex2, paths$gpl570, "GPL570")

scale_info <- rbind(data.table(cohort = "GSE65858", min = min(a$expr), median = median(a$expr), max = max(a$expr), log2_applied = FALSE),
                    data.table(cohort = "GSE41613", min = min(b$expr), median = median(b$expr), max = max(b$expr),
                               log2_applied = log2_applied))
fwrite(scale_info, file.path(tdir, "A_expression_scale.csv"))
fwrite(rbind(x1$stats, x2$stats), file.path(tdir, "A_probe_collapse.csv"))
fwrite(rbind(cbind(platform = "GPL10558", x1$probes), cbind(platform = "GPL570", x2$probes)),
       file.path(tdir, "A_probe_per_gene.csv"))
print(scale_info); print(rbind(x1$stats, x2$stats))

## =========================================================================
## A4. Coverage and scores (mean z over set genes present)
## =========================================================================
sets <- split(gs$gene, gs$set)
cohorts <- list(GSE65858 = list(clin = c1, expr = x1$expr[, c1$gsm], platform = "GPL10558", flow = f1),
                GSE41613 = list(clin = c2, expr = x2$expr[, c2$gsm], platform = "GPL570", flow = f2))
cov <- rbindlist(lapply(names(cohorts), function(k) {
  g <- rownames(cohorts[[k]]$expr)
  data.table(cohort = k, platform = cohorts[[k]]$platform, set = names(sets), genes_in_set = lengths(sets),
             genes_present = sapply(sets, function(s) length(intersect(s, g))),
             genes_missing = sapply(sets, function(s) paste(sort(setdiff(s, g)), collapse = ";")))
}))
cov <- merge(unique(gs[, .(set, programme, role)]), cov, by = "set")
cov[, coverage := genes_present / genes_in_set]
setorder(cov, cohort, role, programme, set)
fwrite(cov, file.path(tdir, "A_set_coverage.csv"))
print(cov[, .(cohort, set, genes_in_set, genes_present, coverage = round(coverage, 2))])

mean_z <- function(m) {
  z <- (m - rowMeans(m)) / matrixStats::rowSds(m)
  z[is.finite(rowSums(z)), , drop = FALSE]
}
for (k in names(cohorts)) {
  z <- mean_z(cohorts[[k]]$expr)
  cohorts[[k]]$scores <- data.table(gsm = colnames(z), sapply(sets, function(s) colMeans(z[intersect(s, rownames(z)), , drop = FALSE])))
}
fwrite(rbindlist(lapply(names(cohorts), function(k) cbind(cohort = k, cohorts[[k]]$scores))), file.path(tdir, "A_programme_scores.csv"))
fwrite(rbindlist(lapply(cohorts, `[[`, "flow")), file.path(tdir, "A_cohort_flow.csv"))
saveRDS(cohorts, paths$step07b_geo)
message("A. saved ", paths$step07b_geo)
