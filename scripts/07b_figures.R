## Step 07b (part C): figures
##
## Purpose: forest plot of programme HRs across TCGA and the two GEO cohorts;
##   median-split Kaplan-Meier curves per cohort for the TCGA-prognostic sets;
##   REMARK profile tables per cohort.
## Inputs:  results/step07b_validation/tables/ (V_* tables)
## Outputs: results/step07b_validation/Fig35-Fig37 (png, pdf)
## Date: 2026-10-28

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(gridExtra)
})
source("config/paths.R")

out <- paths$step07b
rd <- function(f) fread(file.path(out, "tables", f))
base <- 7; ink <- "#0b0b0b"; ink2 <- "#52514e"; rule <- "#c3c2b7"
theme_fig <- function() theme_classic(base_size = base) +
  theme(text = element_text(colour = ink), axis.text = element_text(colour = ink2),
        axis.line = element_line(linewidth = 0.3, colour = rule), axis.ticks = element_line(linewidth = 0.3, colour = rule),
        plot.title = element_text(size = base + 1, face = "bold"), plot.subtitle = element_text(size = base - 0.5, colour = ink2),
        strip.background = element_blank(), strip.text = element_text(face = "bold"), legend.key.size = unit(3, "mm"))
save_fig <- function(p, name, w, h) {
  ggsave(file.path(out, paste0(name, ".png")), p, width = w, height = h, units = "mm", dpi = 300, bg = "white",
         device = ragg::agg_png)
  ggsave(file.path(out, paste0(name, ".pdf")), p, width = w, height = h, units = "mm", bg = "white")
}
set_order <- c("Toustrup2011_hypoxia", "Buffa2010_hypoxia", "Gavish2023_MP6_hypoxia", "Puram2017_pEMT",
               "Gavish2023_MP13_EMTII", "Gavish2023_MP2_G1S", "Gavish2023_MP1_G2M", "Tirosh2016_G1S", "Tirosh2016_G2M")
set_lab <- setNames(c("Toustrup 2011 hypoxia", "Buffa 2010 hypoxia", "Gavish MP6 hypoxia", "Puram 2017 p-EMT",
                      "Gavish MP13 EMT-II", "Gavish MP2 G1/S", "Gavish MP1 G2/M", "Tirosh 2016 G1/S", "Tirosh 2016 G2/M"), set_order)
tcga_prognostic <- c("Buffa2010_hypoxia", "Gavish2023_MP2_G1S", "Gavish2023_MP1_G2M")
coh_lab <- c(TCGA = "TCGA (PFI)", GSE65858 = "GSE65858 (PFS)", GSE41613 = "GSE41613 (DSS)")
coh_col <- setNames(c("#0b0b0b", "#2a6f97", "#c0392b"), coh_lab)

## Fig 35: forest plot, multivariable HR on the primary endpoint, TCGA and both cohorts
sb <- rd("V_side_by_side.csv")[grepl("^multivariable", model)]
vd <- rd("V_validation_verdict.csv")
sb[, `:=`(cohort_lab = factor(coh_lab[cohort], rev(coh_lab)),
          lab = factor(paste0(set_lab[set], ifelse(role == "concordance", " (conc.)", "")),
                       rev(paste0(set_lab[set_order], ifelse(set_order %in% sb[role == "concordance", unique(set)], " (conc.)", "")))))]
sb[, programme := factor(programme, c("hypoxia", "p-EMT", "proliferation"))]
vtxt <- vd[, .(lab = factor(set_lab[set], levels(sb$lab)), programme = factor(programme, levels(sb$programme)),
               txt = fifelse(tcga_prognostic, sub("not validated: ", "not validated (", sub("(direction|significance)$", "\\1)", verdict)),
                             "no verdict"))]
f35 <- ggplot(sb, aes(HR, lab, colour = cohort_lab)) +
  geom_vline(xintercept = 1, linewidth = 0.3, colour = rule) +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.65), size = 0.18, linewidth = 0.3) +
  geom_text(data = vtxt, aes(x = 4.6, y = lab, label = txt), inherit.aes = FALSE, hjust = 1, size = 1.8, colour = ink2) +
  facet_grid(programme ~ ., scales = "free_y", space = "free_y") +
  scale_x_log10(limits = c(0.5, 4.6), breaks = c(0.5, 0.75, 1, 1.5, 2, 3, 4)) +
  scale_colour_manual(values = coh_col, breaks = coh_lab, labels = coh_lab, name = NULL) +
  labs(x = "hazard ratio per SD of programme score (log scale)", y = NULL,
       title = "Step 07a programmes in two independent HPV- cohorts",
       subtitle = paste0("multivariable Cox, primary endpoint. TCGA HPV- (n = 415): score + stage + site (step 07a);\n",
                         "GSE65858 HPV16 DNA- (n = 176): score + age + stage + site; GSE41613 HPV- OSCC (n = 97): score + age + stage.\n",
                         "Verdict (TCGA-prognostic sets): HR > 1 in both GEO cohorts and p < 0.05 in at least one; conc. = concordance set")) +
  theme_fig() + theme(legend.position = "bottom")
save_fig(f35, "Fig35_validation_forest", 180, 120)

## Fig 36: Kaplan-Meier, median split, TCGA-prognostic sets, primary endpoint per cohort
km <- rd("V_km_curves.csv")[set %in% tcga_prognostic & endpoint %in% c("PFS", "DSS")]
lr <- rd("V_km_logrank.csv")[set %in% tcga_prognostic & endpoint %in% c("PFS", "DSS")]
for (d in list(km, lr)) d[, `:=`(lab = factor(set_lab[set], set_lab[tcga_prognostic]),
                                 cohort_lab = factor(coh_lab[cohort], coh_lab[c("GSE65858", "GSE41613")]))]
km[, group := factor(group, c("low", "high"))]
f36 <- ggplot(km, aes(months, surv, colour = group)) + geom_step(linewidth = 0.35) +
  geom_text(data = lr, aes(x = Inf, y = 1, label = sprintf("log-rank p = %.2g\nlow %d/%d, high %d/%d events/n",
                                                           logrank_p, events_low, n_low, events_high, n_high)),
            inherit.aes = FALSE, hjust = 1.05, vjust = 1, size = 1.7, colour = ink2) +
  facet_grid(cohort_lab ~ lab) +
  scale_colour_manual(values = c(low = "#2a6f97", high = "#c0392b"), name = "score (median split)") +
  scale_y_continuous(limits = c(0, 1)) + scale_x_continuous(breaks = seq(0, 84, 24)) +
  labs(x = "months", y = "event-free probability", title = "GEO cohorts: Kaplan-Meier by TCGA-prognostic programme score",
       subtitle = "pre-specified median split; descriptive only; GSE65858 PFS, GSE41613 disease-specific survival") +
  theme_fig() + theme(legend.position = "bottom")
save_fig(f36, "Fig36_validation_km", 180, 105)

## Fig 37: REMARK profile tables per cohort
flow <- rd("V_patient_flow.csv"); rp <- rd("V_remark_profile.csv")
pcm <- rd("V_programme_cox.csv")[role == "primary" & endpoint_role == "primary"]
tt <- ttheme_minimal(base_size = 5.5, core = list(fg_params = list(hjust = 0, x = 0.02, col = ink)),
                     colhead = list(fg_params = list(hjust = 0, x = 0.02, fontface = "bold", col = ink)))
tab1 <- flow[, .(Cohort = cohort, Step = trimws(step), Patients = n)]
tab2 <- rp[, .(Cohort = cohort, Endpoint = sprintf("%s (%s)", endpoint, role), Patients = patients, Events = events,
               `Median FU (months)` = sprintf("%.1f", median_followup_months), Age = sub("range ", "", age), Stage = stage,
               Site = site, Missing = missing_covariates)]
tab3 <- dcast(pcm[model %in% c("univariable", "multivariable")], cohort + set + endpoint ~ model,
              value.var = c("HR", "lo", "hi"))[, .(
  Cohort = cohort, Set = set_lab[set], Endpoint = endpoint,
  `Univariable HR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", HR_univariable, lo_univariable, hi_univariable),
  `Multivariable HR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", HR_multivariable, lo_multivariable, hi_multivariable))]
tab3 <- merge(tab3, pcm[model == "multivariable", .(Cohort = cohort, Set = set_lab[set], Endpoint = endpoint,
                                                    `PH test p (score)` = sprintf("%.2g", ph_p_term))], by = c("Cohort", "Set", "Endpoint"))
setorder(tab3, Cohort, Set)
notes <- paste("Cut-off rule: median of each score in the analysed patients (pre-specified; no optimal cut-point search).",
               "Missing covariates: complete case per model (GSE65858: 3 patients without site are dropped from multivariable models).",
               "GSE65858 times in days (assumed; converted to months for display); GSE41613 times in months (Lohavanichbutr 2013).",
               "GSE41613 age is available only in 4 groups, entered as an ordinal score; DSS_unk = unknown-cause deaths counted as events.",
               sep = "\n")
f37 <- wrap_elements(tableGrob(tab1, rows = NULL, theme = tt)) / wrap_elements(tableGrob(tab2, rows = NULL, theme = tt)) /
  wrap_elements(tableGrob(tab3, rows = NULL, theme = tt)) +
  plot_layout(heights = c(1.15, 0.5, 1.1)) +
  plot_annotation(title = "REMARK profile: GEO validation cohorts", caption = notes,
                  theme = theme(plot.title = element_text(size = base + 1, face = "bold"),
                                plot.caption = element_text(size = base - 1, hjust = 0, colour = ink2)))
save_fig(f37, "Fig37_validation_remark", 250, 215)
