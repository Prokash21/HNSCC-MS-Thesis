## Step 07a (part D): TCGA HPV- survival, PFI (primary) and OS (secondary)
##
## Purpose: Cox models of programme scores and of single genes (set genes and
##   step 06 panel, reserve and CA9 genes) in TCGA HPV- primary tumours;
##   univariable, multivariable (stage + site) and sensitivity (+ age);
##   proportional hazards checks; median-split Kaplan-Meier data; REMARK
##   profile (patient flow, events, follow-up, missing covariates).
## Inputs:  data/raw/tcga_cdr/TCGA-CDR-SupplementalTableS1.xlsx (Liu et al. 2018; downloaded from GDC)
##          data/raw/step01_metadata/cbioportal_hnsc_tcga_primary_site.csv
##          data/processed/v2/step07a_tcga_hpvneg.rds, results/step07a_recurrence/tables/{A_gene_sets.csv,
##          B_tcga_programme_scores.csv}, results/step06_qPCR/tables/D_panel_final.csv
## Outputs: results/step07a_recurrence/tables/S_*.csv
## Date: 2026-10-28

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(readxl)
})
source("config/paths.R")

tdir <- file.path(paths$step07a, "tables")
oral_sites <- c("Oral Tongue", "Oral Cavity", "Floor of mouth", "Buccal Mucosa", "Alveolar Ridge",
                "Hard Palate", "Lip")
gs <- fread(file.path(tdir, "A_gene_sets.csv"))
tc <- readRDS(paths$step07a_tcga)
sc <- fread(file.path(tdir, "B_tcga_programme_scores.csv"))
panel <- fread(file.path(paths$step06, "tables", "D_panel_final.csv"))[set %in% c("core", "reserve", "optional (Arm B)")]

## =========================================================================
## S1. Endpoints and covariates (TCGA-CDR), patient flow
## =========================================================================
cdr <- as.data.table(read_excel(paths$tcga_cdr, sheet = 1, guess_max = 12000, na = c("", "#N/A", "NA")))
cdr <- cdr[type == "HNSC"]
num <- function(v) suppressWarnings(as.numeric(v))
site <- fread(paths$tcga_site)
cdr[, `:=`(case = bcr_patient_barcode, PFI = num(PFI), PFI.time = num(PFI.time), OS = num(OS), OS.time = num(OS.time),
           age = num(age_at_initial_pathologic_diagnosis),
           stage = factor(fcase(ajcc_pathologic_tumor_stage %in% c("Stage I", "Stage II"), "I-II",
                                ajcc_pathologic_tumor_stage == "Stage III", "III",
                                grepl("^Stage IV", ajcc_pathologic_tumor_stage), "IV", default = "unknown"),
                          c("I-II", "III", "IV", "unknown")))]
cdr[, primary_site := site$primary_site[match(case, site$case)]]
cdr[, site := factor(fcase(primary_site %in% oral_sites, "oral cavity", primary_site == "Larynx", "larynx",
                           !is.na(primary_site), "other"), c("oral cavity", "larynx", "other"))]
cdr[, redacted := !is.na(Redaction) & Redaction != ""]

samples <- tc$samples
flow <- rbind(tc$flow,
              data.table(step = "HPV- with a TCGA-CDR HNSC record", n = sum(samples$case %in% cdr$case)),
              data.table(step = "  of these, not redacted", n = sum(samples$case %in% cdr[redacted == FALSE, case])))
base <- merge(samples[, .(case, sample)], cdr[redacted == FALSE], by = "case")
base <- merge(base, sc[, -"case"], by = "sample")
ep_data <- function(ep) {
  d <- copy(base)[, `:=`(time = get(paste0(ep, ".time")), event = get(ep))]
  d[!is.na(time) & time > 0 & !is.na(event)]
}
eps <- c("PFI", "OS")
dd <- setNames(lapply(eps, ep_data), eps)
for (ep in eps) flow <- rbind(flow, data.table(step = sprintf("analysed for %s (time > 0, event known)", ep), n = nrow(dd[[ep]])))
fwrite(flow, file.path(tdir, "S_patient_flow.csv"))

## REMARK profile: n, events, median follow-up (reverse Kaplan-Meier), covariates
remark <- rbindlist(lapply(eps, function(ep) {
  d <- dd[[ep]]
  rk <- summary(survfit(Surv(time, 1 - event) ~ 1, data = d))$table
  data.table(endpoint = ep, patients = nrow(d), events = sum(d$event),
             median_followup_days = rk[["median"]], median_followup_months = rk[["median"]] / 30.44,
             stage_I_II = sum(d$stage == "I-II"), stage_III = sum(d$stage == "III"), stage_IV = sum(d$stage == "IV"),
             stage_unknown = sum(d$stage == "unknown"),
             site_oral = sum(d$site %in% "oral cavity"), site_larynx = sum(d$site %in% "larynx"),
             site_other = sum(d$site %in% "other"), site_missing = sum(is.na(d$site)),
             age_missing = sum(is.na(d$age)), age_median = median(d$age, na.rm = TRUE))
}))
fwrite(remark, file.path(tdir, "S_remark_profile.csv"))
print(flow); print(remark)

## =========================================================================
## S2. Cox models
## =========================================================================
models <- list(univariable = character(), multivariable = c("stage", "site"), sensitivity_age = c("stage", "site", "age"))
cox_fit <- function(d, adj) {
  d <- d[complete.cases(d[, c("time", "event", "x", adj), with = FALSE])]
  f <- coxph(as.formula(paste("Surv(time, event) ~", paste(c("x", adj), collapse = " + "))), data = d)
  s <- summary(f)$coefficients["x", ]
  zph <- tryCatch(cox.zph(f)$table, error = function(e) NULL)
  data.table(n = f$n, events = f$nevent, HR = exp(s[["coef"]]), lo = exp(s[["coef"]] - qnorm(0.975) * s[["se(coef)"]]),
             hi = exp(s[["coef"]] + qnorm(0.975) * s[["se(coef)"]]), logHR = s[["coef"]], p = s[["Pr(>|z|)"]],
             ph_p_term = if (is.null(zph)) NA_real_ else zph["x", "p"],
             ph_p_global = if (is.null(zph)) NA_real_ else zph["GLOBAL", "p"])
}
run_all <- function(values, label_dt) rbindlist(lapply(eps, function(ep) {
  d <- copy(dd[[ep]])
  rbindlist(lapply(names(values), function(v) {
    xv <- values[[v]][match(d$sample, names(values[[v]]))]
    stopifnot(!anyNA(xv))
    d[, x := as.numeric(scale(xv))]
    rbindlist(lapply(names(models), function(m) cbind(data.table(feature = v, endpoint = ep, model = m),
                                                     cox_fit(d, models[[m]]))))
  }))
}))

## programme scores (HR per SD of score)
set_info <- unique(gs[, .(set, programme, role)])
prog_vals <- setNames(lapply(set_info$set, function(s) setNames(sc[[s]], sc$sample)), set_info$set)
prog <- run_all(prog_vals)
prog <- merge(set_info, prog, by.x = "set", by.y = "feature")
prog[model == "multivariable", FDR := p.adjust(p, "BH"), by = .(endpoint, role)]
setorder(prog, endpoint, role, programme, set, model)
fwrite(prog, file.path(tdir, "S_programme_cox.csv"))
print(prog[model == "multivariable", .(endpoint, set, role, n, events, HR = round(HR, 2), lo = round(lo, 2), hi = round(hi, 2),
                                       p = signif(p, 2), FDR = signif(FDR, 2), ph = signif(ph_p_term, 2))])

## single genes (HR per SD of log CPM): primary-set genes + panel, reserve, CA9
gene_list <- unique(c(gs[role == "primary", gene], panel$gene))
tested <- intersect(gene_list, rownames(tc$lcpm))
gene_vals <- setNames(lapply(tested, function(g) tc$lcpm[g, ]), tested)
gcox <- run_all(gene_vals)
gcox[model == "multivariable", FDR := p.adjust(p, "BH"), by = endpoint]
gcox[, `:=`(primary_sets = sapply(feature, function(g) paste(gs[role == "primary" & gene == g, set], collapse = ";")),
            panel_set = panel$set[match(feature, panel$gene)])]
setnames(gcox, "feature", "gene")
fwrite(gcox, file.path(tdir, "S_gene_cox.csv"))
fwrite(data.table(gene = setdiff(gene_list, tested), reason = "not in TCGA log CPM after filterByExpr"),
       file.path(tdir, "S_gene_not_tested.csv"))
message(sprintf("S2. genes: %d listed, %d tested; multivariable PFI p < 0.05: %d (FDR < 0.10: %d)", length(gene_list),
                length(tested), gcox[model == "multivariable" & endpoint == "PFI", sum(p < 0.05)],
                gcox[model == "multivariable" & endpoint == "PFI", sum(FDR < 0.1)]))

## =========================================================================
## S3. Kaplan-Meier, median split of each primary-set score (descriptive)
## =========================================================================
km <- list(); lr <- list()
for (ep in eps) for (s in set_info[role == "primary", set]) {
  d <- copy(dd[[ep]])[, grp := factor(fifelse(get(s) > median(get(s)), "high", "low"), c("low", "high"))]
  f <- survfit(Surv(time, event) ~ grp, data = d)
  km[[paste(ep, s)]] <- rbind(data.table(endpoint = ep, set = s, group = c("low", "high"), time = 0, surv = 1, n_risk = NA_real_),
                              data.table(endpoint = ep, set = s, group = rep(sub("grp=", "", names(f$strata)), f$strata),
                                         time = f$time, surv = f$surv, n_risk = f$n.risk))
  t <- survdiff(Surv(time, event) ~ grp, data = d)
  lr[[paste(ep, s)]] <- data.table(endpoint = ep, set = s, n_low = sum(d$grp == "low"), n_high = sum(d$grp == "high"),
                                   events_low = sum(d[grp == "low", event]), events_high = sum(d[grp == "high", event]),
                                   cutoff = "median of the score in the analysed patients",
                                   logrank_p = pchisq(t$chisq, 1, lower.tail = FALSE))
}
fwrite(rbindlist(km), file.path(tdir, "S_km_curves.csv"))
fwrite(rbindlist(lr), file.path(tdir, "S_km_logrank.csv"))
