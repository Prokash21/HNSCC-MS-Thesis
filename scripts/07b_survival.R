## Step 07b (part B): survival models in GSE65858 and GSE41613, validation verdict
##
## Purpose: Cox models of programme scores (univariable and multivariable with
##   age + stage (+ site)), proportional hazards checks, median-split
##   Kaplan-Meier data, REMARK profile per cohort, the pre-specified validation
##   verdict for the TCGA-prognostic sets, a side-by-side table with the step 07a
##   TCGA results, and descriptive per-gene Cox models.
## Inputs:  data/processed/v2/step07b_geo.rds, results/step07a_recurrence/tables/{A_gene_sets.csv,
##          S_programme_cox.csv, D_recurrence_ranked.csv, D_regulators_ranked.csv},
##          results/step06_qPCR/tables/D_panel_final.csv
## Outputs: results/step07b_validation/tables/V_*.csv
## Date: 2026-10-28

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})
source("config/paths.R")

tdir <- file.path(paths$step07b, "tables")
t07a <- file.path(paths$step07a, "tables")
gs <- fread(file.path(t07a, "A_gene_sets.csv"))
co <- readRDS(paths$step07b_geo)
tcga_prognostic <- c("Buffa2010_hypoxia", "Gavish2023_MP2_G1S", "Gavish2023_MP1_G2M")
p_valid <- 0.05

## cohort design: endpoints (first = primary), multivariable covariates, time unit
design <- list(
  GSE65858 = list(eps = c(PFS = "primary", OS = "secondary"), adj = c("age", "stage", "site"), days_per_unit = 1),
  GSE41613 = list(eps = c(DSS = "primary", OS = "secondary", DSS_unk = "sensitivity"), adj = c("age", "stage"),
                  days_per_unit = 30.44))
models <- function(k) list(univariable = character(), multivariable = design[[k]]$adj)

## =========================================================================
## B1. Analysis data per cohort and endpoint, patient flow, REMARK profile
## =========================================================================
ep_data <- function(k, ep) {
  d <- merge(co[[k]]$clin, co[[k]]$scores, by = "gsm")
  d[, `:=`(time = get(paste0(ep, ".time")), event = get(ep))]
  d[, months := time * design[[k]]$days_per_unit / 30.44]
  d[!is.na(time) & time > 0 & !is.na(event)]
}
dd <- list()
for (k in names(design)) for (ep in names(design[[k]]$eps)) dd[[k]][[ep]] <- ep_data(k, ep)

flow <- rbindlist(lapply(names(design), function(k) {
  rbind(cbind(cohort = k, co[[k]]$flow),
        rbindlist(lapply(names(design[[k]]$eps), function(ep) {
          d <- dd[[k]][[ep]]
          data.table(cohort = k, step = c(sprintf("analysed for %s (time > 0, event known)", ep),
                                          sprintf("  complete covariates for the %s multivariable model", ep)),
                     n = c(nrow(d), sum(complete.cases(d[, design[[k]]$adj, with = FALSE]))))
        })))
}))
fwrite(flow, file.path(tdir, "V_patient_flow.csv"))

remark <- rbindlist(lapply(names(design), function(k) rbindlist(lapply(names(design[[k]]$eps), function(ep) {
  d <- dd[[k]][[ep]]
  rk <- summary(survfit(Surv(months, 1 - event) ~ 1, data = d))$table
  tab <- function(v) paste(sprintf("%s %d", names(table(v, useNA = "ifany")), table(v, useNA = "ifany")), collapse = "; ")
  data.table(cohort = k, endpoint = ep, role = design[[k]]$eps[[ep]], patients = nrow(d), events = sum(d$event),
             median_followup_months = rk[["median"]], max_followup_months = max(d$months),
             age = if (k == "GSE41613") tab(d$age_group) else sprintf("median %.1f (range %.1f-%.1f)", median(d$age),
                                                                       min(d$age), max(d$age)),
             stage = tab(d$stage), site = tab(d$site), sex = tab(d$sex), treatment = tab(d$treatment),
             missing_covariates = paste(sprintf("%s %d", design[[k]]$adj,
                                                sapply(design[[k]]$adj, function(v) sum(is.na(d[[v]])))), collapse = "; "),
             missing_handling = "complete case per model",
             km_cutoff = "median of the score in the analysed patients")
}))))
fwrite(remark, file.path(tdir, "V_remark_profile.csv"))
print(flow); print(remark[, .(cohort, endpoint, patients, events, fu = round(median_followup_months, 1), missing_covariates)])

## =========================================================================
## B2. Cox models
## =========================================================================
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
## values: named list of named vectors (names = gsm)
run_all <- function(k, values) rbindlist(lapply(names(design[[k]]$eps), function(ep) {
  d <- copy(dd[[k]][[ep]])
  rbindlist(lapply(names(values), function(v) {
    xv <- values[[v]][match(d$gsm, names(values[[v]]))]
    stopifnot(!anyNA(xv))
    d[, x := as.numeric(scale(xv))]
    mm <- models(k)
    rbindlist(lapply(names(mm), function(m) cbind(data.table(cohort = k, feature = v, endpoint = ep,
                                                             endpoint_role = design[[k]]$eps[[ep]], model = m),
                                                  cox_fit(d, mm[[m]]))))
  }))
}))

## programme scores
set_info <- unique(gs[, .(set, programme, role)])
prog <- rbindlist(lapply(names(design), function(k) {
  sc <- co[[k]]$scores
  run_all(k, setNames(lapply(set_info$set, function(s) setNames(sc[[s]], sc$gsm)), set_info$set))
}))
prog <- merge(set_info, prog, by.x = "set", by.y = "feature")
prog[model == "multivariable", FDR := p.adjust(p, "BH"), by = .(cohort, endpoint, role)]
setorder(prog, cohort, endpoint, role, programme, set, model)
fwrite(prog, file.path(tdir, "V_programme_cox.csv"))
print(prog[model == "multivariable", .(cohort, endpoint, set, n, events, HR = round(HR, 2), lo = round(lo, 2),
                                       hi = round(hi, 2), p = signif(p, 2), FDR = signif(FDR, 2), ph = signif(ph_p_term, 2))])

## =========================================================================
## B3. Pre-specified validation verdict (primary endpoint, multivariable)
## =========================================================================
pv <- prog[model == "multivariable" & endpoint_role == "primary" & role == "primary"]
val <- dcast(pv, set + programme ~ cohort, value.var = c("HR", "p"))
val[, tcga_prognostic := set %in% tcga_prognostic]
val[, `:=`(hr_above_1_both = HR_GSE65858 > 1 & HR_GSE41613 > 1,
           p_below_005_any = p_GSE65858 < p_valid | p_GSE41613 < p_valid)]
val[, verdict := fifelse(!tcga_prognostic, "no verdict (not prognostic in TCGA)",
                         fifelse(hr_above_1_both & p_below_005_any, "validated",
                                 paste("not validated:", fifelse(!hr_above_1_both & !p_below_005_any, "direction and significance",
                                                                 fifelse(!hr_above_1_both, "direction", "significance")))))]
setorder(val, -tcga_prognostic, programme, set)
fwrite(val, file.path(tdir, "V_validation_verdict.csv"))
print(val)

## side by side: TCGA (step 07a, PFI) and both cohorts (primary endpoint)
tc <- fread(file.path(t07a, "S_programme_cox.csv"))[endpoint == "PFI" & model %in% c("multivariable", "sensitivity_age")]
tc[, `:=`(cohort = "TCGA", endpoint = "PFI", model = fifelse(model == "multivariable", "multivariable (stage + site)",
                                                            "sensitivity (stage + site + age)"))]
geo <- prog[endpoint_role == "primary"]
geo[, model := fifelse(model == "multivariable", sprintf("multivariable (%s)", paste(design[[cohort[1]]]$adj, collapse = " + ")),
                       "univariable"), by = cohort]
sbs <- rbind(tc[, .(cohort, set, programme, role, endpoint, model, n, events, HR, lo, hi, p, FDR, ph_p_term)],
             geo[, .(cohort, set, programme, role, endpoint, model, n, events, HR, lo, hi, p, FDR, ph_p_term)])
fwrite(sbs, file.path(tdir, "V_side_by_side.csv"))

## =========================================================================
## B4. Genes (descriptive): 7 eligible 07a genes, MYBL2, E2F4, step 06 panel/reserve/CA9
## =========================================================================
elig <- fread(file.path(t07a, "D_recurrence_ranked.csv"))$gene
tfs <- fread(file.path(t07a, "D_regulators_ranked.csv"))[shortlist == TRUE, TF]
panel <- fread(file.path(paths$step06, "tables", "D_panel_final.csv"))[set %in% c("core", "reserve", "optional (Arm B)")]
glist <- data.table(gene = c(elig, tfs, panel$gene),
                    group = c(rep("07a eligible gene", length(elig)), rep("07a shortlisted regulator", length(tfs)),
                              paste("step 06", panel$set)))
glist <- glist[!duplicated(gene)]
gcox <- rbindlist(lapply(names(design), function(k) {
  m <- co[[k]]$expr
  tested <- intersect(glist$gene, rownames(m))
  run_all(k, setNames(lapply(tested, function(g) m[g, ]), tested))
}))
setnames(gcox, "feature", "gene")
gcox[model == "multivariable", FDR := p.adjust(p, "BH"), by = .(cohort, endpoint)]
gcox <- merge(glist, gcox, by = "gene")
setorder(gcox, cohort, endpoint, model, p)
fwrite(gcox, file.path(tdir, "V_gene_cox.csv"))
nt <- rbindlist(lapply(names(design), function(k) glist[!gene %in% rownames(co[[k]]$expr)][, cohort := k]))
fwrite(nt, file.path(tdir, "V_gene_not_tested.csv"))
message(sprintf("B4. genes listed %d; not on platform: %s", nrow(glist),
                if (nrow(nt)) paste(nt$cohort, nt$gene, collapse = ", ") else "none"))

## =========================================================================
## B5. Kaplan-Meier, median split of each primary-set score (descriptive)
## =========================================================================
km <- list(); lr <- list()
for (k in names(design)) for (ep in names(design[[k]]$eps)) for (s in set_info[role == "primary", set]) {
  d <- copy(dd[[k]][[ep]])[, grp := factor(fifelse(get(s) > median(get(s)), "high", "low"), c("low", "high"))]
  f <- survfit(Surv(months, event) ~ grp, data = d)
  km[[paste(k, ep, s)]] <- rbind(data.table(cohort = k, endpoint = ep, set = s, group = c("low", "high"), months = 0, surv = 1,
                                            n_risk = NA_real_),
                                 data.table(cohort = k, endpoint = ep, set = s, group = rep(sub("grp=", "", names(f$strata)), f$strata),
                                            months = f$time, surv = f$surv, n_risk = f$n.risk))
  t <- survdiff(Surv(months, event) ~ grp, data = d)
  lr[[paste(k, ep, s)]] <- data.table(cohort = k, endpoint = ep, set = s, n_low = sum(d$grp == "low"), n_high = sum(d$grp == "high"),
                                      events_low = sum(d[grp == "low", event]), events_high = sum(d[grp == "high", event]),
                                      logrank_p = pchisq(t$chisq, 1, lower.tail = FALSE))
}
fwrite(rbindlist(km), file.path(tdir, "V_km_curves.csv"))
fwrite(rbindlist(lr), file.path(tdir, "V_km_logrank.csv"))
