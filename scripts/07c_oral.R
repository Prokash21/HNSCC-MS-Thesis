## Step 07c (part A): TCGA HPV- oral cavity subgroup (exploratory)
##
## Purpose: Cox models of the step 07a programme scores in the oral cavity
##   subgroup (larynx for contrast), site interaction in all-site data,
##   descriptive per-gene models in oral, REMARK profile, median-split
##   Kaplan-Meier data, and an oral side-by-side table with GSE41613 (step 07b).
## Inputs:  data/raw/tcga_cdr/TCGA-CDR-SupplementalTableS1.xlsx, data/raw/step01_metadata/cbioportal_hnsc_tcga_primary_site.csv
##          data/processed/v2/step07a_tcga_hpvneg.rds, results/step07a_recurrence/tables/{A_gene_sets.csv,
##          B_tcga_programme_scores.csv, S_programme_cox.csv, D_recurrence_ranked.csv}
##          results/step07b_validation/tables/{V_programme_cox.csv, V_gene_cox.csv}, results/step06_qPCR/tables/D_panel_final.csv
## Outputs: results/step07c_oral/tables/O_*.csv
## Date: 2026-10-28

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(readxl)
})
source("config/paths.R")

tdir <- file.path(paths$step07c, "tables")
dir.create(tdir, recursive = TRUE, showWarnings = FALSE)
t07a <- file.path(paths$step07a, "tables")
t07b <- file.path(paths$step07b, "tables")
oral_sites <- c("Oral Tongue", "Oral Cavity", "Floor of mouth", "Buccal Mucosa", "Alveolar Ridge",
                "Hard Palate", "Lip")
gs <- fread(file.path(t07a, "A_gene_sets.csv"))
tc <- readRDS(paths$step07a_tcga)
sc <- fread(file.path(t07a, "B_tcga_programme_scores.csv"))
set_info <- unique(gs[, .(set, programme, role)])

## =========================================================================
## O1. Endpoints and covariates (TCGA-CDR, as in step 07a), subgroups, REMARK
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
base <- merge(tc$samples[, .(case, sample)], cdr[redacted == FALSE], by = "case")
base <- merge(base, sc[, -"case"], by = "sample")
base[, oral := factor(fifelse(site == "oral cavity", "oral", "non-oral"), c("non-oral", "oral"))]
eps <- c("PFI", "OS")
ep_data <- function(ep, sub = NULL) {
  d <- copy(base)[, `:=`(time = get(paste0(ep, ".time")), event = get(ep))]
  d <- d[!is.na(time) & time > 0 & !is.na(event)]
  if (!is.null(sub)) d <- d[site == sub]
  droplevels(d)
}
groups <- c(oral = "oral cavity", larynx = "larynx")
dd <- list()
for (g in names(groups)) for (ep in eps) dd[[g]][[ep]] <- ep_data(ep, groups[[g]])
for (ep in eps) dd[["all"]][[ep]] <- ep_data(ep)

flow <- rbind(data.table(step = "TCGA HPV- primary tumours with RNA-seq (step 07a)", n = nrow(tc$samples)),
              data.table(step = "  with a non-redacted TCGA-CDR record", n = nrow(base)),
              data.table(step = sprintf("  site %s", levels(base$site)), n = as.integer(table(base$site))),
              data.table(step = "  site missing", n = sum(is.na(base$site))),
              rbindlist(lapply(names(groups), function(g) data.table(
                step = sprintf("%s: analysed for %s (time > 0, event known)", g, eps),
                n = sapply(eps, function(ep) nrow(dd[[g]][[ep]]))))))
fwrite(flow, file.path(tdir, "O_patient_flow.csv"))

tab <- function(v) paste(sprintf("%s %d", names(table(v, useNA = "ifany")), table(v, useNA = "ifany")), collapse = "; ")
remark <- rbindlist(lapply(names(groups), function(g) rbindlist(lapply(eps, function(ep) {
  d <- dd[[g]][[ep]]
  rk <- summary(survfit(Surv(time, 1 - event) ~ 1, data = d))$table
  data.table(subgroup = g, endpoint = ep, patients = nrow(d), events = sum(d$event),
             median_followup_months = rk[["median"]] / 30.44,
             stage = tab(d$stage), primary_site = tab(d$primary_site),
             age = sprintf("median %.0f (range %.0f-%.0f)", median(d$age, na.rm = TRUE), min(d$age, na.rm = TRUE),
                           max(d$age, na.rm = TRUE)),
             age_missing = sum(is.na(d$age)), missing_handling = "complete case per model; stage unknown as own level",
             km_cutoff = "median of the score in the analysed oral patients")
}))))
fwrite(remark, file.path(tdir, "O_remark_profile.csv"))
print(flow); print(remark[, .(subgroup, endpoint, patients, events, fu = round(median_followup_months, 1), stage, age_missing)])

## =========================================================================
## O2. Cox models in subgroups (programmes and genes)
## =========================================================================
models <- list(univariable = character(), main = "stage", sensitivity_age = c("stage", "age"))
cox_fit <- function(d, adj) {
  d <- droplevels(d[complete.cases(d[, c("time", "event", "x", adj), with = FALSE])])
  f <- coxph(as.formula(paste("Surv(time, event) ~", paste(c("x", adj), collapse = " + "))), data = d)
  s <- summary(f)$coefficients["x", ]
  zph <- tryCatch(cox.zph(f)$table, error = function(e) NULL)
  data.table(n = f$n, events = f$nevent, HR = exp(s[["coef"]]), lo = exp(s[["coef"]] - qnorm(0.975) * s[["se(coef)"]]),
             hi = exp(s[["coef"]] + qnorm(0.975) * s[["se(coef)"]]), logHR = s[["coef"]], p = s[["Pr(>|z|)"]],
             ph_p_term = if (is.null(zph)) NA_real_ else zph["x", "p"],
             ph_p_global = if (is.null(zph)) NA_real_ else zph["GLOBAL", "p"])
}
run_group <- function(g, values) rbindlist(lapply(eps, function(ep) {
  d <- copy(dd[[g]][[ep]])
  rbindlist(lapply(names(values), function(v) {
    xv <- values[[v]][match(d$sample, names(values[[v]]))]
    stopifnot(!anyNA(xv))
    d[, x := as.numeric(scale(xv))]
    rbindlist(lapply(names(models), function(m) cbind(data.table(subgroup = g, feature = v, endpoint = ep, model = m),
                                                     cox_fit(d, models[[m]]))))
  }))
}))
prog_vals <- setNames(lapply(set_info$set, function(s) setNames(sc[[s]], sc$sample)), set_info$set)
prog <- rbindlist(lapply(names(groups), run_group, values = prog_vals))
prog <- merge(set_info, prog, by.x = "set", by.y = "feature")
prog[model == "main", FDR := p.adjust(p, "BH"), by = .(subgroup, endpoint, role)]
setorder(prog, subgroup, endpoint, role, programme, set, model)
fwrite(prog, file.path(tdir, "O_programme_cox.csv"))
print(prog[model == "main", .(subgroup, endpoint, set, n, events, HR = round(HR, 2), lo = round(lo, 2), hi = round(hi, 2),
                              p = signif(p, 2), FDR = signif(FDR, 2), ph = signif(ph_p_term, 2))])

## genes (oral only)
elig <- fread(file.path(t07a, "D_recurrence_ranked.csv"))$gene
panel <- fread(file.path(paths$step06, "tables", "D_panel_final.csv"))[set %in% c("core", "reserve", "optional (Arm B)")]
glist <- data.table(gene = c(elig, "MYBL2", "E2F4", "HIF1A", "EPAS1", panel$gene),
                    group = c(rep("07a eligible gene", length(elig)), rep("regulator", 4), paste("step 06", panel$set)))
glist <- glist[!duplicated(gene)]
tested <- intersect(glist$gene, rownames(tc$lcpm))
gcox <- run_group("oral", setNames(lapply(tested, function(g) tc$lcpm[g, ]), tested))
setnames(gcox, "feature", "gene")
gcox[model == "main", FDR := p.adjust(p, "BH"), by = endpoint]
gcox <- merge(glist, gcox, by = "gene")
setorder(gcox, endpoint, model, p)
fwrite(gcox, file.path(tdir, "O_gene_cox.csv"))
fwrite(glist[!gene %in% tested], file.path(tdir, "O_gene_not_tested.csv"))
net <- gcox[gene %in% c("NFIC", "FOXM1", "HMGB3", "MYBL2") & model %in% c("main", "sensitivity_age")]
fwrite(net[order(gene, endpoint, model)], file.path(tdir, "O_network_genes.csv"))
message(sprintf("O2. genes listed %d, tested %d", nrow(glist), length(tested)))

## =========================================================================
## O3. Site interaction, all-site: Surv ~ score * oral + stage
## =========================================================================
inter <- rbindlist(lapply(eps, function(ep) rbindlist(lapply(set_info$set, function(s) {
  d <- copy(dd[["all"]][[ep]])[!is.na(oral)]
  d[, x := as.numeric(scale(get(s)))]
  f <- coxph(Surv(time, event) ~ x * oral + stage, data = d)
  b <- coef(f); V <- vcov(f); k <- "x:oraloral"
  est <- function(w) { e <- sum(w * b[names(w)]); se <- sqrt(drop(t(w) %*% V[names(w), names(w)] %*% w))
                       c(HR = exp(e), lo = exp(e - qnorm(0.975) * se), hi = exp(e + qnorm(0.975) * se),
                         p = 2 * pnorm(-abs(e / se))) }
  i <- est(setNames(1, k)); no <- est(setNames(1, "x")); o <- est(setNames(c(1, 1), c("x", k)))
  data.table(set = s, endpoint = ep, n = f$n, events = f$nevent, n_oral = sum(d$oral == "oral"),
             interaction_HR = i[["HR"]], interaction_lo = i[["lo"]], interaction_hi = i[["hi"]], interaction_p = i[["p"]],
             HR_oral = o[["HR"]], HR_oral_lo = o[["lo"]], HR_oral_hi = o[["hi"]],
             HR_nonoral = no[["HR"]], HR_nonoral_lo = no[["lo"]], HR_nonoral_hi = no[["hi"]])
}))))
inter <- merge(set_info, inter, by = "set")
setorder(inter, endpoint, role, programme, set)
fwrite(inter, file.path(tdir, "O_site_interaction.csv"))
print(inter[, .(endpoint, set, int_HR = round(interaction_HR, 2), int_p = signif(interaction_p, 2),
                HR_oral = round(HR_oral, 2), HR_nonoral = round(HR_nonoral, 2))])

## =========================================================================
## O4. Kaplan-Meier, median split in oral (primary sets)
## =========================================================================
km <- list(); lr <- list()
for (ep in eps) for (s in set_info[role == "primary", set]) {
  d <- copy(dd[["oral"]][[ep]])[, grp := factor(fifelse(get(s) > median(get(s)), "high", "low"), c("low", "high"))]
  f <- survfit(Surv(time, event) ~ grp, data = d)
  km[[paste(ep, s)]] <- rbind(data.table(endpoint = ep, set = s, group = c("low", "high"), time = 0, surv = 1, n_risk = NA_real_),
                              data.table(endpoint = ep, set = s, group = rep(sub("grp=", "", names(f$strata)), f$strata),
                                         time = f$time, surv = f$surv, n_risk = f$n.risk))
  t <- survdiff(Surv(time, event) ~ grp, data = d)
  lr[[paste(ep, s)]] <- data.table(endpoint = ep, set = s, n_low = sum(d$grp == "low"), n_high = sum(d$grp == "high"),
                                   events_low = sum(d[grp == "low", event]), events_high = sum(d[grp == "high", event]),
                                   logrank_p = pchisq(t$chisq, 1, lower.tail = FALSE))
}
fwrite(rbindlist(km), file.path(tdir, "O_km_curves.csv"))
fwrite(rbindlist(lr), file.path(tdir, "O_km_logrank.csv"))

## =========================================================================
## O5. Oral side by side: TCGA oral vs GSE41613 (step 07b); no pooling
## =========================================================================
ep_map <- data.table(endpoint = c("PFI", "OS"), gse_endpoint = c("DSS", "OS"), level = c("primary", "secondary"))
tco <- prog[subgroup == "oral" & model %in% c("main", "sensitivity_age")]
gse <- fread(file.path(t07b, "V_programme_cox.csv"))[cohort == "GSE41613" & model == "multivariable" & endpoint %in% c("DSS", "OS")]
fmt <- function(HR, lo, hi, p) sprintf("%.2f (%.2f-%.2f), p %.2g", HR, lo, hi, p)
sbs_p <- merge(dcast(tco[, .(set, programme, role, endpoint, model, txt = fmt(HR, lo, hi, p))], set + programme + role + endpoint ~ model,
                     value.var = "txt"), ep_map, by = "endpoint")
sbs_p <- merge(sbs_p, gse[, .(set, gse_endpoint = endpoint, GSE41613 = fmt(HR, lo, hi, p), GSE41613_HR = HR)],
               by = c("set", "gse_endpoint"), all.x = TRUE)
sbs_p <- merge(sbs_p, tco[model == "main", .(set, endpoint, TCGA_oral_HR = HR)], by = c("set", "endpoint"))
setnames(sbs_p, c("main", "sensitivity_age"), c("TCGA_oral_stage", "TCGA_oral_stage_age"))
sbs_p[, same_direction := sign(log(TCGA_oral_HR)) == sign(log(GSE41613_HR))]
setorder(sbs_p, -level, role, programme, set)
fwrite(sbs_p, file.path(tdir, "O_side_by_side_programmes.csv"))

gg <- fread(file.path(t07b, "V_gene_cox.csv"))[cohort == "GSE41613" & model == "multivariable" & endpoint %in% c("DSS", "OS")]
tg <- gcox[model == "main"]
sbs_g <- merge(tg[, .(gene, group, endpoint, TCGA_oral = fmt(HR, lo, hi, p), TCGA_oral_HR = HR, TCGA_oral_FDR = FDR)], ep_map, by = "endpoint")
sbs_g <- merge(sbs_g, gg[, .(gene, gse_endpoint = endpoint, GSE41613 = fmt(HR, lo, hi, p), GSE41613_HR = HR, GSE41613_FDR = FDR)],
               by = c("gene", "gse_endpoint"), all.x = TRUE)
sbs_g[is.na(GSE41613), GSE41613 := fifelse(gene %in% c("HIF1A", "EPAS1"), "not run in step 07b", "not on GPL570")]
sbs_g[, same_direction := sign(log(TCGA_oral_HR)) == sign(log(GSE41613_HR))]
setorder(sbs_g, -level, group, gene)
fwrite(sbs_g, file.path(tdir, "O_side_by_side_genes.csv"))
print(sbs_p[level == "primary", .(set, TCGA_oral_stage, GSE41613, same_direction)])
