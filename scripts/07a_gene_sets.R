## Step 07a (part A): pre-specified programme gene sets
##
## Purpose: read the hypoxia, p-EMT and proliferation gene sets from the
##   papers' own tables, map symbols to current HGNC symbols, and report the
##   overlap of each set with the step 06 panel, reserve and CA9.
## Inputs:  data/raw/step01_metadata/papers/ (Buffa 2010 Table S5, Gavish 2023
##          Supplementary Table "Cancer MPs", Tirosh 2016 Table S5; downloaded)
##          config/gene_sets/ (Toustrup 2011 Table S2 and Puram 2017 Table S7,
##          transcribed from the PDF tables), HGNC complete set,
##          GPL570 annotation (resolves an ambiguous Buffa symbol by its probe ID),
##          results/step06_qPCR/tables/D_panel_final.csv
## Outputs: results/step07a_recurrence/tables/A_gene_sets.csv, A_set_summary.csv,
##          A_symbol_changes.csv, A_panel_overlap.csv
## Date: 2026-10-28

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})
source("config/paths.R")
source("scripts/lib_de.R")

tdir <- file.path(paths$step07a, "tables")
dir.create(tdir, recursive = TRUE, showWarnings = FALSE)
pp <- function(...) file.path(paths$papers_dir, ...)
hg <- read_hgnc(paths$hgnc)

read_txt <- function(f) { x <- trimws(readLines(file.path(paths$gene_sets_dir, f))); x[x != "" & !startsWith(x, "#")] }

## Buffa 2010, Table S5 "Common Hypoxia Signature Genes" (header on row 5)
buffa_t <- as.data.table(read_excel(pp("Buffa 2010, Br J Cancer", "41416_2010_BF6605450_MOESM10_ESM.xls"), sheet = 1, skip = 4))
buffa_t <- buffa_t[!is.na(`HUGO Symbol`), .(probe = `Affymetrix Probeset ID`, symbol = `HUGO Symbol`)]
buffa <- buffa_t$symbol

## Gavish 2023, Supplementary Table, sheet "Cancer MPs" (one column per MP, 50 genes)
gv <- read_excel(pp("Gavish2023", "41586_2023_6130_MOESM6_ESM.xlsx"), sheet = "Cancer MPs")
gv_col <- function(id) { j <- grep(sprintf("^%s\\s", id), names(gv)); stopifnot(length(j) == 1); list(name = names(gv)[j], genes = na.omit(gv[[j]])) }
mp1 <- gv_col("MP1"); mp2 <- gv_col("MP2"); mp6 <- gv_col("MP6"); mp13 <- gv_col("MP13")

## Tirosh 2016, Table S5 "Cell cycle gene-sets": phase-specific G1/S (column 1) and G2/M (column 2)
ti <- as.data.frame(read_excel(pp("Tirosh 2016, Science", "aad0501_table_s5.xlsx"), sheet = 1, col_names = FALSE,
                               .name_repair = "minimal"))
h <- which(ti[[1]] == "G1/S" & ti[[2]] == "G2/M")
stopifnot(length(h) == 1)
tir_g1s <- na.omit(ti[[1]][-(1:h)]); tir_g2m <- na.omit(ti[[2]][-(1:h)])

sets <- list(
  list(set = "Toustrup2011_hypoxia", programme = "hypoxia", role = "primary", genes = read_txt("toustrup2011_tableS2.txt"),
       source = "Toustrup 2011 Cancer Res, Table S2, genes marked ** (final 15-gene classifier)"),
  list(set = "Buffa2010_hypoxia", programme = "hypoxia", role = "primary", genes = buffa,
       source = "Buffa 2010 Br J Cancer, Table S5 Common Hypoxia Signature Genes (MOESM10_ESM.xls)"),
  list(set = "Gavish2023_MP6_hypoxia", programme = "hypoxia", role = "concordance", genes = mp6$genes,
       source = sprintf("Gavish 2023 Nature, Supplementary Table sheet 'Cancer MPs', column '%s'", mp6$name)),
  list(set = "Puram2017_pEMT", programme = "p-EMT", role = "primary", genes = read_txt("puram2017_tableS7_pEMT.txt"),
       source = "Puram 2017 Cell, Table S7 (mmc1.pdf p. 6), column p-EMT, genes 1-100"),
  list(set = "Gavish2023_MP13_EMTII", programme = "p-EMT", role = "concordance", genes = mp13$genes,
       source = sprintf("Gavish 2023 Nature, Supplementary Table sheet 'Cancer MPs', column '%s'", mp13$name)),
  list(set = "Gavish2023_MP2_G1S", programme = "proliferation", role = "primary", genes = mp2$genes,
       source = sprintf("Gavish 2023 Nature, Supplementary Table sheet 'Cancer MPs', column '%s'", mp2$name)),
  list(set = "Gavish2023_MP1_G2M", programme = "proliferation", role = "primary", genes = mp1$genes,
       source = sprintf("Gavish 2023 Nature, Supplementary Table sheet 'Cancer MPs', column '%s'", mp1$name)),
  list(set = "Tirosh2016_G1S", programme = "proliferation", role = "concordance", genes = tir_g1s,
       source = "Tirosh 2016 Science, Table S5, phase-specific genes, G1/S"),
  list(set = "Tirosh2016_G2M", programme = "proliferation", role = "concordance", genes = tir_g2m,
       source = "Tirosh 2016 Science, Table S5, phase-specific genes, G2/M"))

gs <- rbindlist(lapply(sets, function(s) data.table(set = s$set, programme = s$programme, role = s$role, source = s$source,
                                                   order = seq_along(s$genes), original = trimws(s$genes))))
## Gavish prints some open-reading-frame symbols in upper case (e.g. C4ORF3); HGNC uses "orf"
gs[, query := sub("^(C[0-9XY]+)ORF([0-9]+)$", "\\1orf\\2", original)]
mp <- hgnc_map(unique(gs$query), hg)
gs[, `:=`(gene = mp$mapped[match(query, mp$original)], mapping = mp$how[match(query, mp$original)])]
gs[query != original & !is.na(gene), mapping := paste0("case fixed; ", mapping)]
## ambiguous symbols: resolve with the source table's own Affymetrix probe (GPL570 annotation)
amb <- gs[is.na(gene) & set == "Buffa2010_hypoxia", unique(original)]
if (length(amb)) {
  stopifnot(file.exists(paths$gpl570))
  an <- fread(cmd = paste("gzip -dc", shQuote(paths$gpl570)), skip = "ID\tGene title", select = c("ID", "Gene symbol"))
  for (a in amb) {
    pr <- buffa_t[symbol == a, probe]
    cand <- unlist(strsplit(an[ID %in% pr, `Gene symbol`], "///"))
    cand <- intersect(cand, hg$symbol)
    if (length(cand) == 1) gs[original == a & set == "Buffa2010_hypoxia",
                              `:=`(gene = cand, mapping = sprintf("ambiguous; resolved by probe %s (GPL570)", pr))]
  }
}
changes <- unique(gs[mapping != "approved", .(original, gene, mapping, sets = set)])[, .(sets = paste(sets, collapse = ";")),
                                                                                    by = .(original, gene, mapping)]
gs <- gs[!is.na(gene)][, query := NULL]
gs <- gs[!duplicated(gs[, .(set, gene)])]

panel <- fread(file.path(paths$step06, "tables", "D_panel_final.csv"))
panel <- panel[set %in% c("core", "reserve", "optional (Arm B)")]
gs[, panel_set := panel$set[match(gene, panel$gene)]]

set_sum <- gs[, .(genes_in_table = NA_integer_, genes_mapped = .N, in_panel = paste(gene[!is.na(panel_set)], collapse = ";")),
              by = .(set, programme, role, source)]
set_sum[, genes_in_table := vapply(sets, function(s) length(s$genes), 0L)[match(set, vapply(sets, `[[`, "", "set"))]]
overlap <- dcast(gs[!is.na(panel_set), .(gene, panel_set, set)], gene + panel_set ~ set, fun.aggregate = length, value.var = "set")
## Jaccard overlap between sets of the same programme (primary vs concordance)
jac <- rbindlist(lapply(split(gs, gs$programme), function(d) {
  ss <- unique(d$set); cmb <- CJ(a = ss, b = ss)[a < b]
  cmb[, `:=`(programme = d$programme[1],
             shared = mapply(function(x, y) length(intersect(d[set == x, gene], d[set == y, gene])), a, b),
             jaccard = mapply(function(x, y) length(intersect(d[set == x, gene], d[set == y, gene])) /
                                length(union(d[set == x, gene], d[set == y, gene])), a, b))]
}))

fwrite(gs, file.path(tdir, "A_gene_sets.csv"))
fwrite(set_sum, file.path(tdir, "A_set_summary.csv"))
fwrite(changes, file.path(tdir, "A_symbol_changes.csv"))
fwrite(overlap, file.path(tdir, "A_panel_overlap.csv"))
fwrite(jac, file.path(tdir, "A_set_overlap.csv"))
print(set_sum[, .(set, role, genes_in_table, genes_mapped, in_panel)])
print(changes)
print(jac)
