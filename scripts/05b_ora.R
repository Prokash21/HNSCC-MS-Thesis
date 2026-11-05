## Step 05b (part C): over-representation analysis and GO redundancy reduction
##
## Purpose: ORA (clusterProfiler enricher) of HPV- up candidates (tiers 1, 1b,
##   2) and down candidates (tier 1b down + down_replicated) against the C2
##   tested-gene universe, GO BP and Hallmark; cluster significant GO BP terms
##   (ORA and replicated GSEA terms) with simplifyEnrichment and draw the GO
##   cluster heatmaps.
## Inputs:  results/step05_TF/tables/A_candidates_final.csv, results/step04b_DE/tables/D_C2_HPVneg_allpatients.csv,
##          data/processed/v2/step05b_gene_sets.rds, data/reference/msigdb/c5.go.bp.v2024.1.Hs.json,
##          results/step05b_functional/tables/B_replicated_sets.csv
## Outputs: results/step05b_functional/tables/C_*.csv, results/step05b_functional/Fig43_go_clusters (png, pdf)
## Date: 2026-10-28

suppressPackageStartupMessages({
  library(data.table)
  library(clusterProfiler)
  library(simplifyEnrichment)
})
source("config/paths.R")

out <- paths$step05b
tdir <- file.path(out, "tables")
fdr_cut <- 0.05
min_terms <- 5
sets <- readRDS(paths$step05b_sets)

## =========================================================================
## C1. Gene lists and universe
## =========================================================================
cf <- fread(file.path(paths$step05, "tables", "A_candidates_final.csv"))
up <- cf[tier %in% c("1", "1b", "2") & direction == "up", unique(gene)]
down <- cf[(tier == "1b" & direction == "down") | tier == "down_replicated", unique(gene)]
universe <- fread(file.path(paths$step04b, "tables", "D_C2_HPVneg_allpatients.csv"))$gene
stopifnot(all(c(up, down) %in% universe))
fwrite(data.table(list = c("up", "down"), genes = c(length(up), length(down)), universe = length(universe),
                  definition = c("step 05 tiers 1, 1b, 2 with direction up", "step 05 tier 1b down + down_replicated")),
       file.path(tdir, "C_ora_lists.csv"))
message(sprintf("C1. up %d, down %d, universe %d", length(up), length(down), length(universe)))

## =========================================================================
## C2. ORA: GO BP and Hallmark
## =========================================================================
t2g <- function(col) rbindlist(lapply(names(sets[[col]]), function(s) data.table(term = s, gene = sets[[col]][[s]])))
tg <- list(GOBP = t2g("GOBP"), H = t2g("H"))
ora <- rbindlist(lapply(c(up = "up", down = "down"), function(l) rbindlist(lapply(names(tg), function(col) {
  e <- enricher(get(l), universe = universe, TERM2GENE = tg[[col]], minGSSize = 15, maxGSSize = 500,
                pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1)
  if (is.null(e)) return(NULL)
  r <- as.data.table(e@result)
  r[, `:=`(list = l, collection = col)]
}))))
ora <- ora[, .(list, collection, term = ID, GeneRatio, BgRatio, pvalue, FDR = p.adjust, qvalue, genes = gsub("/", ";", geneID), Count)]
setorder(ora, list, collection, pvalue)
fwrite(ora, file.path(tdir, "C_ora_all.csv"))
fwrite(ora[FDR < fdr_cut], file.path(tdir, "C_ora_significant.csv"))
print(ora[FDR < fdr_cut, .N, by = .(list, collection)])
print(ora[FDR < fdr_cut & collection == "H", .(list, term, GeneRatio, FDR = signif(FDR, 2), genes)])

## =========================================================================
## C3. GO redundancy: simplifyEnrichment (GO semantic similarity, binary cut)
## =========================================================================
js <- jsonlite::fromJSON(paths$msigdb_gobp_json, simplifyVector = FALSE)
go_id <- sapply(js, function(x) if (is.null(x$exactSource)) NA_character_ else x$exactSource)
rep_sets <- fread(file.path(tdir, "B_replicated_sets.csv"))[collection == "GOBP"]
groups <- list(ORA_up = ora[list == "up" & collection == "GOBP" & FDR < fdr_cut, term],
               ORA_down = ora[list == "down" & collection == "GOBP" & FDR < fdr_cut, term],
               GSEA_replicated_up = rep_sets[direction == "up in tumour", set],
               GSEA_replicated_down = rep_sets[direction == "down in tumour", set])
clu <- list(); info <- list()
for (g in names(groups)) {
  terms <- groups[[g]]
  ids <- go_id[terms]
  info[[g]] <- data.table(group = g, terms = length(terms), with_go_id = sum(!is.na(ids)), clustered = FALSE)
  ids <- ids[!is.na(ids)]
  if (length(ids) < min_terms) { message(sprintf("C3. %s: %d terms, clustering skipped (< %d)", g, length(ids), min_terms)); next }
  mat <- GO_similarity(unname(ids), ont = "BP")
  cl <- simplifyGO(mat, plot = FALSE)
  base <- file.path(out, paste0("Fig43_go_clusters_", g))
  title <- sprintf("%s: %d GO BP terms (simplifyEnrichment, binary cut)", gsub("_", " ", g), nrow(mat))
  ragg::agg_png(paste0(base, ".png"), width = 2400, height = 1500, res = 220)
  ht_clusters(mat, cl$cluster, use_raster = FALSE, column_title = title)
  dev.off()
  pdf(paste0(base, ".pdf"), width = 11, height = 7)
  ht_clusters(mat, cl$cluster, use_raster = FALSE, column_title = title)
  dev.off()
  cl <- as.data.table(cl)
  cl[, `:=`(group = g, term = names(ids)[match(id, ids)])]
  clu[[g]] <- cl
  info[[g]][, `:=`(clustered = TRUE, clusters = uniqueN(cl$cluster))]
  message(sprintf("C3. %s: %d terms in %d clusters", g, nrow(cl), uniqueN(cl$cluster)))
}
fwrite(rbindlist(info, fill = TRUE), file.path(tdir, "C_go_cluster_summary.csv"))
if (length(clu)) {
  cl_all <- rbindlist(clu, fill = TRUE)
  fwrite(cl_all[order(group, cluster)], file.path(tdir, "C_go_clusters.csv"))
  ## cluster label: the terms of each cluster, largest clusters first
  lab <- cl_all[, .(n_terms = .N, terms = paste(head(sub("^GOBP_", "", term), 6), collapse = "; ")), by = .(group, cluster)]
  fwrite(lab[order(group, -n_terms)], file.path(tdir, "C_go_cluster_labels.csv"))
  print(lab[order(group, -n_terms)][, head(.SD, 6), by = group])
}
