## Step 07c (part B): figures (exploratory oral subgroup)
##
## Purpose: forest plot of programme HRs in TCGA oral cavity vs larynx vs all
##   sites with the site-interaction p; median-split Kaplan-Meier curves in oral
##   for the TCGA-prognostic sets; oral side by side, TCGA oral vs GSE41613.
## Inputs:  results/step07c_oral/tables/ (O_* tables), results/step07a_recurrence/tables/S_programme_cox.csv
## Outputs: results/step07c_oral/Fig38-Fig40 (png, pdf)
## Date: 2026-10-28

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})
source("config/paths.R")

out <- paths$step07c
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
set_lab <- setNames(c("Toustrup 2011 hypoxia", "Buffa 2010 hypoxia", "Gavish MP6 hypoxia (conc.)", "Puram 2017 p-EMT",
                      "Gavish MP13 EMT-II (conc.)", "Gavish MP2 G1/S", "Gavish MP1 G2/M", "Tirosh 2016 G1/S (conc.)",
                      "Tirosh 2016 G2/M (conc.)"), set_order)
tcga_prognostic <- c("Buffa2010_hypoxia", "Gavish2023_MP2_G1S", "Gavish2023_MP1_G2M")
prog_levels <- c("hypoxia", "p-EMT", "proliferation")
exploratory <- "EXPLORATORY, post hoc subgroup analysis"

## Fig 38: forest, all sites (step 07a) vs oral vs larynx, PFI and OS, stage-adjusted
pc <- rd("O_programme_cox.csv")[model == "main"]
all <- fread(file.path(paths$step07a, "tables", "S_programme_cox.csv"))[model == "multivariable"]
grp_lab <- c(all = "all sites, n = 415 (stage + site; step 07a)", oral = "oral cavity, n = 277 (stage)",
             larynx = "larynx, n = 106 (stage)")
fd <- rbind(all[, .(group = "all", set, programme, endpoint, HR, lo, hi)],
            pc[, .(group = subgroup, set, programme, endpoint, HR, lo, hi)])
fd[, `:=`(group = factor(grp_lab[group], rev(grp_lab)), lab = factor(set_lab[set], rev(set_lab[set_order])),
          programme = factor(programme, prog_levels), endpoint = factor(endpoint, c("PFI", "OS")))]
it <- rd("O_site_interaction.csv")
it[, `:=`(lab = factor(set_lab[set], levels(fd$lab)), programme = factor(programme, prog_levels),
          endpoint = factor(endpoint, levels(fd$endpoint)), txt = sprintf("int. p %.2g", interaction_p))]
f38 <- ggplot(fd, aes(HR, lab, colour = group)) +
  geom_vline(xintercept = 1, linewidth = 0.3, colour = rule) +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.65), size = 0.16, linewidth = 0.3) +
  geom_text(data = it, aes(x = 2.6, y = lab, label = txt), inherit.aes = FALSE, hjust = 1, size = 1.7, colour = ink2) +
  facet_grid(programme ~ endpoint, scales = "free_y", space = "free_y") +
  scale_x_log10(limits = c(0.55, 2.6), breaks = c(0.6, 0.8, 1, 1.25, 1.6, 2)) +
  scale_colour_manual(values = setNames(c("#0b0b0b", "#c0392b", "#2a6f97"), grp_lab), breaks = unname(grp_lab), name = NULL) +
  labs(x = "hazard ratio per SD of programme score (log scale)", y = NULL,
       title = "TCGA HPV-: programme scores by site",
       subtitle = paste0(exploratory, ". Cox per subgroup: score + stage; all sites: score + stage + site (step 07a).\n",
                         "int. p = score x oral-vs-other interaction in all-site data (score * oral + stage); interaction tests have low power")) +
  theme_fig() + theme(legend.position = "bottom")
save_fig(f38, "Fig38_oral_site_forest", 190, 130)

## Fig 39: Kaplan-Meier, oral, TCGA-prognostic sets, median split
km <- rd("O_km_curves.csv")[set %in% tcga_prognostic]; lr <- rd("O_km_logrank.csv")[set %in% tcga_prognostic]
for (d in list(km, lr)) d[, `:=`(lab = factor(set_lab[set], set_lab[tcga_prognostic]), endpoint = factor(endpoint, c("PFI", "OS")))]
km[, `:=`(months = time / 30.44, group = factor(group, c("low", "high")))]
f39 <- ggplot(km, aes(months, surv, colour = group)) + geom_step(linewidth = 0.35) +
  geom_text(data = lr, aes(x = Inf, y = 1, label = sprintf("log-rank p = %.2g\nlow %d/%d, high %d/%d events/n",
                                                           logrank_p, events_low, n_low, events_high, n_high)),
            inherit.aes = FALSE, hjust = 1.05, vjust = 1, size = 1.7, colour = ink2) +
  facet_grid(endpoint ~ lab) +
  scale_colour_manual(values = c(low = "#2a6f97", high = "#c0392b"), name = "score (median split within oral)") +
  scale_y_continuous(limits = c(0, 1)) + scale_x_continuous(breaks = c(0, 60, 120, 180)) +
  labs(x = "months", y = "event-free probability", title = "TCGA HPV- oral cavity: Kaplan-Meier by programme score",
       subtitle = paste0(exploratory, "; n = 277; pre-specified median split; descriptive only")) +
  theme_fig() + theme(legend.position = "bottom")
save_fig(f39, "Fig39_oral_km", 180, 100)

## Fig 40: oral side by side, TCGA oral (PFI, stage) vs GSE41613 (DSS, age + stage)
coh_lab <- c(tcga = "TCGA HPV- oral, PFI (n = 277; score + stage)", gse = "GSE41613 HPV- OSCC, DSS (n = 97; score + age + stage; step 07b)")
coh_col <- setNames(c("#c0392b", "#7a4b94"), coh_lab)
sp <- rd("O_side_by_side_programmes.csv")[level == "primary"]
pv <- rbind(pc[subgroup == "oral" & endpoint == "PFI", .(cohort = "tcga", set, programme, HR, lo, hi)],
            fread(file.path(paths$step07b, "tables", "V_programme_cox.csv"))[cohort == "GSE41613" & model == "multivariable" &
                                                                              endpoint == "DSS", .(cohort = "gse", set, programme, HR, lo, hi)])
pv[, `:=`(cohort = factor(coh_lab[cohort], rev(coh_lab)), lab = factor(set_lab[set], rev(set_lab[set_order])),
          programme = factor(programme, prog_levels))]
pa <- ggplot(pv, aes(HR, lab, colour = cohort)) +
  geom_vline(xintercept = 1, linewidth = 0.3, colour = rule) +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.6), size = 0.16, linewidth = 0.3) +
  facet_grid(programme ~ ., scales = "free_y", space = "free_y") +
  scale_x_log10(limits = c(0.6, 4.4), breaks = c(0.75, 1, 1.5, 2, 3, 4)) +
  scale_colour_manual(values = coh_col, breaks = unname(coh_lab), name = NULL) +
  labs(x = "HR per SD (log scale)", y = NULL, title = "Programmes") +
  theme_fig() + theme(legend.position = "none")
gt <- rd("O_gene_cox.csv")[model == "main" & endpoint == "PFI", .(cohort = "tcga", gene, group, HR, lo, hi)]
gg <- fread(file.path(paths$step07b, "tables", "V_gene_cox.csv"))[cohort == "GSE41613" & model == "multivariable" & endpoint == "DSS",
                                                                  .(cohort = "gse", gene, HR, lo, hi)]
gg[, group := gt$group[match(gene, gt$gene)]]
gv <- rbind(gt, gg[, .(cohort, gene, group, HR, lo, hi)])
gv[, group := factor(fcase(group == "07a eligible gene", "07a eligible", group == "regulator", "regulators",
                           group == "step 06 core", "step 06 core", default = "reserve / CA9"),
                     c("07a eligible", "regulators", "step 06 core", "reserve / CA9"))]
gv[, `:=`(cohort = factor(coh_lab[cohort], rev(coh_lab)), gene = factor(gene, rev(sort(unique(gene)))))]
pb <- ggplot(gv, aes(HR, gene, colour = cohort)) +
  geom_vline(xintercept = 1, linewidth = 0.3, colour = rule) +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.6), size = 0.12, linewidth = 0.25) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  scale_x_log10(breaks = c(0.5, 0.75, 1, 1.5, 2, 3)) +
  scale_colour_manual(values = coh_col, breaks = unname(coh_lab), name = NULL) +
  labs(x = "HR per SD (log scale)", y = NULL, title = "Genes",
       caption = "HIF1A and EPAS1 were not run in step 07b (TCGA only); CASC9 is not on GPL570. No pooling.") +
  guides(colour = guide_legend(nrow = 2)) +
  theme_fig() + theme(legend.position = "bottom", strip.text.y = element_text(angle = 0, hjust = 0),
                      plot.caption = element_text(size = base - 1.5, colour = ink2, hjust = 0))
f40 <- (pa | pb) + plot_layout(widths = c(1, 1.1)) +
  plot_annotation(title = "Oral cavity HPV- cancer: TCGA oral subgroup vs GSE41613",
                  subtitle = paste0(exploratory, "; primary endpoint of each cohort; the two endpoints differ (PFI vs disease-specific survival)"),
                  theme = theme(plot.title = element_text(size = base + 1, face = "bold"),
                                plot.subtitle = element_text(size = base - 0.5, colour = ink2)))
save_fig(f40, "Fig40_oral_tcga_vs_gse41613", 190, 160)
