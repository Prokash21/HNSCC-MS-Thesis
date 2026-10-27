## Step 07a (part F): figures
##
## Purpose: Bill relapse vs primary programme scores; TCGA Cox forest plots
##   (PFI and OS); median-split Kaplan-Meier curves; REMARK profile table;
##   eligible-gene ranking.
## Inputs:  results/step07a_recurrence/tables/ (C_*, S_*, D_* tables)
## Outputs: results/step07a_recurrence/Fig30-Fig34 (png, pdf)
## Date: 2026-10-28

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(gridExtra)
})
source("config/paths.R")

out <- paths$step07a
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
role_col <- c(primary = "#0b0b0b", concordance = "#8a8a84")
grp_col <- c(primary = "#2a6f97", relapse = "#c0392b")

## Fig 30: Bill relapse vs primary, patient-mean single-cell scores
sc <- rd("C_bill_programme_scores.csv")
te <- rd("C_bill_programme_tests.csv")[analysis == "all relapse"]
sc[, lab := factor(paste0(set_lab[set], ifelse(role == "concordance", " (conc.)", "")),
                   paste0(set_lab[set_order], ifelse(set_order %in% sc[role == "concordance", unique(set)], " (conc.)", "")))]
te[, lab := factor(paste0(set_lab[set], ifelse(role == "concordance", " (conc.)", "")), levels(sc$lab))]
n_p <- te$n_primary[1]; n_r <- te$n_relapse[1]
f30 <- ggplot(sc, aes(group, score_single_cell, colour = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, linewidth = 0.3, colour = ink2) +
  geom_point(aes(shape = sample_type), position = position_jitter(width = 0.12, seed = 1), size = 1.1) +
  geom_text(data = te, aes(x = 1.5, y = Inf, label = sprintf("Wilcoxon p = %.2g\nFDR = %.2g", wilcox_p, FDR)),
            inherit.aes = FALSE, vjust = 1.1, size = 1.8, colour = ink2) +
  facet_wrap(~ lab, scales = "free_y", nrow = 3) +
  scale_colour_manual(values = grp_col, guide = "none") +
  scale_shape_manual(values = c(primary = 16, recurrence = 17, metastasis = 15), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.35))) +
  labs(x = NULL, y = "patient mean single-cell score (Tirosh 2016)",
       title = "Bill 2023 HPV-: programme scores, relapse vs primary",
       subtitle = sprintf("malignant cells; %d primary vs %d relapse patients (8 locoregional recurrences, 4 metastases); conc. = concordance set",
                          n_p, n_r)) +
  theme_fig() + theme(legend.position = "bottom")
save_fig(f30, "Fig30_bill_programme_scores", 170, 150)

## Fig 31: TCGA Cox forest plots (PFI and OS)
pc <- rd("S_programme_cox.csv")[model %in% c("univariable", "multivariable")]
pc[, `:=`(lab = factor(set_lab[set], rev(set_lab[set_order])), endpoint = factor(endpoint, c("PFI", "OS")),
          model = factor(model, c("multivariable", "univariable")))]
fdr_lab <- pc[model == "multivariable", .(lab, endpoint, txt = sprintf("HR %.2f (%.2f-%.2f), FDR %.2g", HR, lo, hi, FDR))]
f31 <- ggplot(pc, aes(HR, lab)) +
  geom_vline(xintercept = 1, linewidth = 0.3, colour = rule) +
  geom_pointrange(aes(xmin = lo, xmax = hi, colour = role, shape = model), position = position_dodge(width = 0.55),
                  size = 0.2, linewidth = 0.3) +
  geom_text(data = fdr_lab, aes(x = 1.62, y = lab, label = txt), hjust = 0, size = 1.8, colour = ink2) +
  facet_wrap(~ endpoint) +
  scale_x_log10(limits = c(0.75, 2.6), breaks = c(0.8, 1, 1.25, 1.5)) +
  scale_colour_manual(values = role_col, name = NULL) +
  scale_shape_manual(values = c(multivariable = 16, univariable = 1), name = NULL) +
  labs(x = "hazard ratio per SD of programme score", y = NULL,
       title = "TCGA HPV-: programme scores and outcome",
       subtitle = "Cox; multivariable = score + stage + site; labels: multivariable HR (95% CI), BH FDR within primary or concordance sets") +
  theme_fig() + theme(legend.position = "bottom")
save_fig(f31, "Fig31_tcga_cox_forest", 180, 85)

## Fig 32: Kaplan-Meier, median split of primary-set scores
km <- rd("S_km_curves.csv"); lr <- rd("S_km_logrank.csv")
km[, `:=`(lab = factor(set_lab[set], set_lab[set_order]), months = time / 30.44, group = factor(group, c("low", "high")))]
lr[, lab := factor(set_lab[set], set_lab[set_order])]
f32 <- ggplot(km, aes(months, surv, colour = group)) + geom_step(linewidth = 0.35) +
  geom_text(data = lr, aes(x = Inf, y = 1, label = sprintf("log-rank p = %.2g", logrank_p)), inherit.aes = FALSE,
            hjust = 1.05, vjust = 1, size = 1.8, colour = ink2) +
  facet_grid(endpoint ~ lab) +
  scale_colour_manual(values = c(low = "#2a6f97", high = "#c0392b"), name = "score (median split)") +
  scale_y_continuous(limits = c(0, 1)) + scale_x_continuous(breaks = c(0, 60, 120, 180)) +
  labs(x = "months", y = "event-free probability", title = "TCGA HPV-: Kaplan-Meier by programme score",
       subtitle = "pre-specified median split; descriptive only") +
  theme_fig() + theme(legend.position = "bottom")
save_fig(f32, "Fig32_tcga_km", 180, 85)

## Fig 33: REMARK profile table
flow <- rd("S_patient_flow.csv"); rp <- rd("S_remark_profile.csv")
pcm <- rd("S_programme_cox.csv")[role == "primary"]
tab1 <- flow[, .(Step = trimws(step), Patients = n)]
tab2 <- rp[, .(Endpoint = endpoint, Patients = patients, Events = events,
               `Median follow-up (months)` = sprintf("%.1f", median_followup_months),
               `Stage I-II / III / IV / unknown` = sprintf("%d / %d / %d / %d", stage_I_II, stage_III, stage_IV, stage_unknown),
               `Site oral / larynx / other` = sprintf("%d / %d / %d", site_oral, site_larynx, site_other))]
tab3 <- dcast(pcm[model %in% c("univariable", "multivariable")],
              set + endpoint ~ model, value.var = c("HR", "lo", "hi"))[, .(
  Set = set_lab[set], Endpoint = endpoint,
  `Univariable HR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", HR_univariable, lo_univariable, hi_univariable),
  `Multivariable HR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", HR_multivariable, lo_multivariable, hi_multivariable))]
tab3 <- merge(tab3, pcm[model == "multivariable", .(Set = set_lab[set], Endpoint = endpoint,
                                                    `PH test p (score)` = sprintf("%.2g", ph_p_term))], by = c("Set", "Endpoint"))
setorder(tab3, Endpoint, Set)
tt <- ttheme_minimal(base_size = 6, core = list(fg_params = list(hjust = 0, x = 0.02, col = ink)),
                     colhead = list(fg_params = list(hjust = 0, x = 0.02, fontface = "bold", col = ink)))
notes <- paste("Cut-off rule: median of each score in the analysed patients (pre-specified; no optimal cut-point search).",
               "Missing covariates: stage missing or [Discrepancy] coded as 'unknown'; site and age complete.",
               "Endpoints: TCGA-CDR PFI (primary) and OS (secondary); scores standardised (HR per SD).", sep = "\n")
f33 <- wrap_elements(tableGrob(tab1, rows = NULL, theme = tt)) / wrap_elements(tableGrob(tab2, rows = NULL, theme = tt)) /
  wrap_elements(tableGrob(tab3, rows = NULL, theme = tt)) +
  plot_layout(heights = c(1, 0.45, 1.4)) +
  plot_annotation(title = "REMARK profile: TCGA HPV- programme scores", caption = notes,
                  theme = theme(plot.title = element_text(size = base + 1, face = "bold"),
                                plot.caption = element_text(size = base - 1, hjust = 0, colour = ink2)))
save_fig(f33, "Fig33_remark_profile", 180, 190)

## Fig 34: eligible-gene ranking and the eligibility screen
el <- rd("D_recurrence_eligible_all.csv")[c1_programme_prognostic == TRUE]
rk <- rd("D_recurrence_ranked.csv")
el[, cls := fcase(eligible, "eligible", c2_gene_PFI, "gene PFI passes, other criterion fails", default = "gene PFI fails")]
pa <- ggplot(el[!is.na(PFI_p)], aes(log2(PFI_HR), -log10(PFI_p))) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.3, colour = rule) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = rule) +
  geom_point(aes(colour = cls), size = 0.9) +
  ggrepel::geom_text_repel(data = el[eligible == TRUE], aes(label = gene), size = 1.9, segment.size = 0.2) +
  scale_colour_manual(values = c(eligible = "#c0392b", `gene PFI passes, other criterion fails` = "#e69f00",
                                 `gene PFI fails` = "grey70"), name = NULL) +
  labs(x = "log2 HR per SD (PFI, multivariable)", y = expression(-log[10]~p),
       title = "a  Genes of prognostic primary sets",
       subtitle = sprintf("%d genes in Buffa hypoxia, Gavish G1/S, Gavish G2/M; %d tested in TCGA; %d eligible",
                          nrow(el), sum(!is.na(el$PFI_p)), sum(el$eligible))) +
  theme_fig() + theme(legend.position = "bottom") + guides(colour = guide_legend(nrow = 2))
lng <- melt(rk[, .(gene, rank, shortlist, `log HR PFI` = pr_logHR, `Bill log2FC` = pr_bill, `tumour share` = pr_share,
                   `log TCGA TPM` = pr_TPM, `mean (rank score)` = rank_score)],
            id.vars = c("gene", "rank", "shortlist"), variable.name = "component", value.name = "pct")
lng[, gene_lab := factor(sprintf("%d. %s%s", rank, gene, ifelse(shortlist, " *", "")),
                         rev(rk[, sprintf("%d. %s%s", rank, gene, ifelse(shortlist, " *", ""))]))]
pb <- ggplot(lng, aes(component, gene_lab, fill = pct)) + geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", pct)), size = 1.8) +
  scale_fill_gradient(low = "#f4f1ea", high = "#c0392b", limits = c(0, 1), breaks = c(0, 0.5, 1), name = "percentile") +
  labs(x = NULL, y = NULL, title = "b  Ranking of eligible genes", subtitle = "* = suggested shortlist (max 2 per programme, size 4)") +
  theme_fig() + theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom")
save_fig(pa + pb + plot_layout(widths = c(1.1, 1)), "Fig34_recurrence_ranking", 180, 95)
message("figures done")
