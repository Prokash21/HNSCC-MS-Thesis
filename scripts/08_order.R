## Step 08 (parts B-C): primer report, new-gene expression check, ordering file, primer map
##
## Purpose: expression check of the 5 newly designed genes (step 04b C2 log2FC
##   and FDR, TCGA paired log2FC, malignant pseudobulk CPM, TCGA paired tumour
##   TPM and expression class as in step 06, tumour-cell share); pseudogene and
##   gDNA flags (step 06 rules); the final primer ordering file (CSV and XLSX,
##   one row per oligo, optional genes on a separate sheet); Fig51 primer map
##   of the new genes (style of Fig29).
## Inputs:  results/step08_primers/tables/{A_final_gene_list.csv, B_primers_final.csv, B_primer_map.csv},
##          results/step05c_gene_atlas/tables/B_gene_atlas_matrix.csv, results/step04b_DE/tables/{D_C2_HPVneg_allpatients,
##          E_TCGA_HPVneg_paired_tumour_vs_normal}.csv, data/processed/v2/{tcga_hnsc_paired_counts, tcga_hnsc_paired_tpm}.rds,
##          data/raw/step01_metadata/cbioportal_hnsc_tcga_pan_can_atlas_2018_patient.csv, HGNC complete set
## Outputs: results/step08_primers/{FINAL_primer_order.csv, FINAL_primer_order.xlsx, Fig51_primer_map_new (png, pdf)},
##          results/step08_primers/tables/{C_new_gene_expression.csv, C_primer_flags.csv}
## Date: 2026-10-28

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(openxlsx)
})
source("config/paths.R")
source("scripts/lib_de.R")

out <- paths$step08
tdir <- file.path(out, "tables")
rd <- function(f) fread(file.path(tdir, f))
fdr_cut <- 0.05; tpm_min <- 1; tpm_ok <- 10; min_cpm <- 10
gl <- rd("A_final_gene_list.csv")
pr <- rd("B_primers_final.csv")
hg <- read_hgnc(paths$hgnc)
down_genes <- c("KRT13", "NFIC")

## =========================================================================
## B1. Expression check (all genes; reported for the 5 new genes)
## =========================================================================
at <- fread(file.path(paths$step05c, "tables", "B_gene_atlas_matrix.csv"))
tc <- readRDS(paths$tcga_counts)
cb <- fread(paths$cbio_pancan, select = c("patientId", "SUBTYPE"))
ts <- copy(tc$samples)[, hpv := cb$SUBTYPE[match(case, cb$patientId)]][hpv == "HNSC_HPV-"]
tpm <- readRDS(paths$tcga_tpm)$tpm
tum <- ts[sample_type == "Primary Tumor", sample]; nor <- ts[sample_type != "Primary Tumor", sample]
ex <- data.table(gene = gl$gene)
ex <- merge(ex, at[, .(gene, C2_log2FC = lfc_C2, C2_FDR = fdr_C2, TCGA_log2FC = lfc_TCGA, TCGA_FDR = fdr_TCGA,
                       malignant_pb_CPM = sc_malignant, normal_pb_CPM = sc_normal, tumour_share, specificity)], by = "gene", all.x = TRUE)
ex[, `:=`(TCGA_TPM_tumour = apply(tpm[gene, tum, drop = FALSE], 1, median), TCGA_TPM_normal = apply(tpm[gene, nor, drop = FALSE], 1, median))]
ex[, TPM_high := fifelse(gene %in% down_genes, TCGA_TPM_normal, TCGA_TPM_tumour)]
ex[, expression_class := fcase(TPM_high < tpm_min, "very low (<1 TPM)", TPM_high < tpm_ok, "low (1-10 TPM)", default = "moderate/high (>=10 TPM)")]
ex[, detect_CPM := fifelse(gene %in% down_genes, normal_pb_CPM, malignant_pb_CPM)]
ex[, low_expresser := TPM_high < tpm_ok | detect_CPM < min_cpm]
ex <- merge(gl[, .(gene, group, primer_source)], ex, by = "gene", sort = FALSE)
fwrite(ex, file.path(tdir, "C_expression_check.csv"))
newg <- ex[primer_source %like% "new"]
fwrite(newg, file.path(tdir, "C_new_gene_expression.csv"))
print(newg[, .(gene, C2 = round(C2_log2FC, 2), C2_FDR = signif(C2_FDR, 2), TCGA = round(TCGA_log2FC, 2), mal_CPM = round(malignant_pb_CPM, 1),
               TPM_tum = round(TCGA_TPM_tumour, 1), expression_class, share = round(tumour_share, 2), low_expresser)])

## =========================================================================
## B2. Flags (step 06 rules) and ordering file
## =========================================================================
hg_pseudo <- hg[locus_group == "pseudogene", symbol]
fl <- pr[, .(gene, junction_status, isoform_specific, main_isoforms_detected, pseudogene_products, readthrough_products)]
fl[, hgnc_pseudogenes := sapply(gene, function(g) paste(grep(paste0("^", g, "P[0-9]+[A-Z]?$"), hg_pseudo, value = TRUE), collapse = " "))]
fl[, flags := mapply(function(js, iso, mid, pp, rt, hp) {
  f <- c(if (grepl("one exon|single-exon", js)) "no junction / single exon (gDNA can amplify)",
         if (grepl("different exons", js)) "intron-flanking, no junction primer",
         if (pp != "none") paste("PREDICTED pseudogene product:", pp),
         if (hp != "") { v <- strsplit(hp, " ")[[1]]
                         paste0("HGNC pseudogenes: ", if (length(v) > 4) sprintf("%d (%s, ...)", length(v), paste(v[1:3], collapse = ", ")) else paste(v, collapse = ", ")) },
         if (rt != "none") "readthrough transcript also amplified",
         if (iso) paste("detects", mid, "main isoforms"))
  if (length(f)) paste(f, collapse = "; ") else "none"
}, junction_status, isoform_specific, main_isoforms_detected, pseudogene_products, readthrough_products, hgnc_pseudogenes)]
fl[, DNase_noRT := fifelse(grepl("one exon|single-exon|different exons", junction_status) | hgnc_pseudogenes != "" | pseudogene_products != "none",
                           "essential (DNase + no-RT)", "recommended (standard DNase + no-RT per sample)")]
fwrite(fl, file.path(tdir, "C_primer_flags.csv"))
ord <- merge(gl[, .(gene, group, order_set, idx = .I)], pr, by = c("gene", "group"))
ord <- merge(ord, fl[, .(gene, flags, DNase_noRT)], by = "gene")
setorder(ord, idx)
olig <- rbind(ord[, .(idx, gene, group, order_set, oligo = paste0(gene, "_F"), sequence = forward, Tm = f_tm, GC = f_gc, product_bp,
                      canonical_transcript = transcript, MANE = mane, junction_status, flags, DNase_noRT, primer_source)],
              ord[, .(idx, gene, group, order_set, oligo = paste0(gene, "_R"), sequence = reverse, Tm = r_tm, GC = r_gc, product_bp,
                      canonical_transcript = transcript, MANE = mane, junction_status, flags, DNase_noRT, primer_source)])
setorder(olig, idx, oligo)
olig[, `:=`(length = nchar(sequence), scale = "25 nmol", purification = "desalted", `Primer-BLAST status` = "")]
cols <- c("gene", "group", "oligo", "sequence", "length", "Tm", "GC", "product_bp", "canonical_transcript", "MANE", "junction_status",
          "flags", "DNase_noRT", "primer_source", "scale", "purification", "Primer-BLAST status")
setnames(olig, c("sequence", "Tm", "GC"), c("sequence_5to3", "Tm_C", "GC_pct"))
cols <- sub("^sequence$", "sequence_5to3", sub("^Tm$", "Tm_C", sub("^GC$", "GC_pct", cols)))
fwrite(olig[, c("order_set", cols), with = FALSE], file.path(out, "FINAL_primer_order.csv"))
wb <- createWorkbook()
hs <- createStyle(textDecoration = "bold", fgFill = "#e8e4d8", border = "bottom")
for (sh in list(c("core_22_assays", "core (22 assays)"), c("optional", "optional (order if budget allows)"))) {
  addWorksheet(wb, sh[1])
  d <- olig[order_set == sh[2], ..cols]
  writeData(wb, sh[1], d, headerStyle = hs)
  setColWidths(wb, sh[1], cols = seq_along(cols), widths = "auto")
  freezePane(wb, sh[1], firstRow = TRUE)
}
addWorksheet(wb, "notes")
writeData(wb, "notes", data.table(note = c(
  "One row per oligo, 5'->3'. Scale 25 nmol, desalted.",
  "Core sheet: 22 assays (44 oligos). Optional sheet: NFIC and ARPC1A, order only if budget allows.",
  "Primers designed with Primer3 (primer3-py 2.4.0) on the Ensembl 116 canonical (MANE Select where available) transcript;",
  "  70-150 bp product, Tm 58-62 C, junction-spanning first; checked against the Ensembl 116 transcriptome (cDNA + ncRNA).",
  "Fill 'Primer-BLAST status' after checking each pair in NCBI Primer-BLAST (RefSeq mRNA and genome databases).",
  "DNase treatment and a no-RT control are recommended for every sample; 'essential' marks assays without a junction-spanning",
  "  primer or with HGNC pseudogenes (gDNA or pseudogene signal possible).",
  sprintf("Generated %s from scripts/08_order.R.", Sys.Date()))))
saveWorkbook(wb, file.path(out, "FINAL_primer_order.xlsx"), overwrite = TRUE)
print(olig[, .N, by = order_set])
print(fl[flags != "none", .(gene, flags, DNase_noRT)])

## =========================================================================
## Fig 51: primer map of the 5 new genes (style of Fig29)
## =========================================================================
base <- 7; ink <- "#0b0b0b"; ink2 <- "#52514e"; rule <- "#c3c2b7"
theme_fig <- function() theme_classic(base_size = base) +
  theme(text = element_text(colour = ink), axis.text = element_text(colour = ink2),
        axis.line = element_line(linewidth = 0.3, colour = rule), axis.ticks = element_line(linewidth = 0.3, colour = rule),
        plot.title = element_text(size = base + 1, face = "bold"), plot.subtitle = element_text(size = base - 0.5, colour = ink2),
        strip.background = element_blank(), strip.text = element_text(face = "bold"), legend.key.size = unit(3, "mm"))
mp <- rd("B_primer_map.csv")
pg <- newg$gene
flank <- 120
rows <- lapply(pg, function(g) {
  m <- mp[gene == g]; f <- m[kind == "forward"]; r <- m[kind == "reverse"]; tx <- m[kind == "transcript"]
  lo <- max(0, f$start - flank); hi <- min(tx$end, r$end + flank)
  e <- m[kind == "exon" & end > lo & start < hi][, `:=`(start = pmax(start, lo), end = pmin(end, hi))]
  list(ex = e[, .(gene, label, x0 = start - f$start, x1 = end - f$start)],
       pr = data.table(gene = g, kind = c("forward", "reverse"), x0 = c(0, r$start - f$start), x1 = c(f$end - f$start, r$end - f$start)),
       amp = data.table(gene = g, x0 = 0, x1 = r$end - f$start))
})
exn <- rbindlist(lapply(rows, `[[`, "ex")); prr <- rbindlist(lapply(rows, `[[`, "pr")); amp <- rbindlist(lapply(rows, `[[`, "amp"))
strip <- setNames(sprintf("%s (%d bp, %s; %s)", pg, pr$product_bp[match(pg, pr$gene)], pr$main_isoforms_detected[match(pg, pr$gene)],
                          pr$transcript[match(pg, pr$gene)]), pg)
for (d in list(exn, prr, amp)) d[, gene := factor(strip[gene], levels = strip)]
exn[, shade := factor(as.integer(label) %% 2)]
f51 <- ggplot() +
  geom_rect(data = exn, aes(xmin = x0, xmax = x1, ymin = -0.25, ymax = 0.25, fill = shade), colour = "white", linewidth = 0.4) +
  geom_text(data = exn, aes(x = (x0 + x1) / 2, y = 0, label = paste0("exon ", label)), size = (base - 1.5) / .pt, colour = ink2) +
  geom_segment(data = amp, aes(x = x0, xend = x1, y = -0.42, yend = -0.42), colour = "#e34948", linewidth = 1.2) +
  geom_segment(data = prr[kind == "forward"], aes(x = x0, xend = x1, y = 0.4, yend = 0.4), colour = ink, linewidth = 0.5,
               arrow = arrow(length = unit(1.2, "mm"), type = "closed")) +
  geom_segment(data = prr[kind == "reverse"], aes(x = x1, xend = x0, y = 0.4, yend = 0.4), colour = ink, linewidth = 0.5,
               arrow = arrow(length = unit(1.2, "mm"), type = "closed")) +
  facet_wrap(~ gene, ncol = 3, scales = "free_x") +
  scale_fill_manual(values = c("0" = "#d6d5cf", "1" = "#c3c2b7"), guide = "none") +
  scale_y_continuous(limits = c(-0.55, 0.55)) +
  labs(title = "qPCR primer positions on the canonical transcript (new step 08 genes)",
       subtitle = "arrows = primers; red = amplicon; boxes = exons; strip: product size, main isoforms detected, canonical transcript",
       x = "position relative to the forward primer (bp)", y = NULL) +
  theme_fig() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
                      strip.text = element_text(hjust = 0))
ggsave(file.path(out, "Fig51_primer_map_new.png"), f51, width = 180, height = 25 + 22 * ceiling(length(pg) / 3), units = "mm", dpi = 300,
       bg = "white", device = ragg::agg_png)
ggsave(file.path(out, "Fig51_primer_map_new.pdf"), f51, width = 180, height = 25 + 22 * ceiling(length(pg) / 3), units = "mm", bg = "white")
