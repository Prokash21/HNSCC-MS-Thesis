## Transcription-factor (TF) activity in tumour cells: HPV+ vs HPV-, and
## tumour vs normal squamous. decoupleR univariate linear model (ULM) with the
## CollecTRI TF-target network.
##
## Two complementary views, as the strategy note recommends:
##   1. contrast level: ULM on the edgeR statistics of each DGE comparison
##   2. patient level: ULM on each patient's tumour pseudobulk, then
##      HPV+ vs HPV- tested per TF with the patient as the unit (study-adjusted)
## A TF is called robust when both views agree (FDR < 0.05, same direction).
## Positive controls for HPV+: E2F family, TFDP1, FOXM1, MYBL2 (E7 -> RB -> E2F).
##
## Inputs:  data/processed/HNSCatlas/dge_results.rds (scripts/06)
##          CollecTRI from OmniPath (downloaded once to data/reference/collectri/)
## Outputs: HNSCC-HPV/results/06_TF_activity/ (Fig10 + tables/)
## Run from the repo root:  Rscript scripts/12_tf_activity.R

suppressPackageStartupMessages({
  library(edgeR)
  library(decoupleR)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

out_dir <- "HNSCC-HPV/results/06_TF_activity"
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
res <- readRDS("data/processed/HNSCatlas/dge_results.rds")
controls <- c("E2F1", "E2F2", "E2F3", "TFDP1", "FOXM1", "MYBL2")

## ---- CollecTRI network ----------------------------------------------------
## same processing as decoupleR::get_collectri(split_complexes = FALSE):
## NF-kB and AP-1 complexes become "NFKB" and "AP1"; inhibition -> mor -1
net_file <- "data/reference/collectri/collectri_omnipath.tsv"
if (!file.exists(net_file)) {
  dir.create(dirname(net_file), recursive = TRUE, showWarnings = FALSE)
  download.file(paste0("https://omnipathdb.org/interactions?datasets=collectri&genesymbols=yes",
                       "&organisms=9606&loops=yes&license=academic&format=tsv"), net_file)
}
raw <- read.delim(net_file)
src <- raw$source_genesymbol
cplx <- grepl("_", src)
src[cplx & grepl("NFKB|REL", src)] <- "NFKB"
src[cplx & grepl("JUN|FOS|ATF", src)] <- "AP1"
inhibits <- as.logical(raw$is_inhibition) & !as.logical(raw$is_stimulation)
net <- data.frame(source = src, target = raw$target_genesymbol, mor = ifelse(inhibits, -1, 1))
net <- net[!grepl("_", net$source) & !duplicated(net[, c("source", "target")]), ]
message(sprintf("CollecTRI: %d TFs, %d interactions", length(unique(net$source)), nrow(net)))

ulm <- function(mat) {
  run_ulm(mat = mat, net = net, .source = "source", .target = "target", .mor = "mor", minsize = 5)
}

## ---- 1. contrast level ----------------------------------------------------
## signed edgeR statistic: sign(log2FC) * sqrt(F) (QL F-test, 1 df ~ t);
## F can be a tiny negative number from rounding, so floor it at 0
stat_of <- function(tab) setNames(sign(tab$logFC) * sqrt(pmax(tab$F, 0)), tab$gene)
contrasts <- list(B  = stat_of(res$B),  A1 = stat_of(res$A1), A2 = stat_of(res$A2))
contrast_tf <- do.call(rbind, lapply(names(contrasts), function(nm) {
  s <- contrasts[[nm]]
  r <- as.data.frame(ulm(matrix(s, ncol = 1, dimnames = list(names(s), nm))))
  r$FDR <- p.adjust(r$p_value, "BH")
  r
}))
contrast_tf <- contrast_tf[, c("condition", "source", "score", "p_value", "FDR")]
names(contrast_tf) <- c("contrast", "TF", "score", "p_value", "FDR")

## ---- 2. patient level (tumour pseudobulk, HPV+ vs HPV-) --------------------
lc <- cpm(res$yB, log = TRUE, prior.count = 2)
z  <- t(scale(t(lc)))                                    # gene z-scores across patients
act <- as.data.frame(ulm(z))
act_mat <- tapply(act$score, list(act$source, act$condition), identity)
act_mat <- act_mat[, res$sB$sample]
hpv   <- factor(res$sB$HPV, levels = c("HPV-", "HPV+"))
study <- factor(res$sB$Study)
pt <- do.call(rbind, lapply(rownames(act_mat), function(tf) {
  f <- summary(lm(act_mat[tf, ] ~ study + hpv))$coefficients["hpvHPV+", ]
  data.frame(TF = tf, diff_HPVpos = f[["Estimate"]], t = f[["t value"]], p_value = f[["Pr(>|t|)"]])
}))
pt$FDR <- p.adjust(pt$p_value, "BH")

## ---- combine and save -------------------------------------------------------
hb <- contrast_tf[contrast_tf$contrast == "B", ]
tf_hpv <- merge(hb[, c("TF", "score", "FDR")], pt[, c("TF", "diff_HPVpos", "t", "FDR")],
                by = "TF", suffixes = c("_contrast", "_patient"))
tf_hpv$robust <- tf_hpv$FDR_contrast < 0.05 & tf_hpv$FDR_patient < 0.05 &
                 sign(tf_hpv$score) == sign(tf_hpv$t)
tf_hpv$higher_in <- ifelse(tf_hpv$score > 0, "HPV+", "HPV-")
tf_hpv <- tf_hpv[order(-tf_hpv$robust, tf_hpv$FDR_contrast), ]
write.csv(tf_hpv, file.path(out_dir, "tables", "TF_activity_HPVpos_vs_HPVneg.csv"), row.names = FALSE)
write.csv(contrast_tf, file.path(out_dir, "tables", "TF_activity_all_contrasts.csv"), row.names = FALSE)
write.csv(data.frame(TF = rownames(act_mat), act_mat, check.names = FALSE),
          file.path(out_dir, "tables", "TF_activity_per_patient.csv"), row.names = FALSE)

message(sprintf("HPV+ vs HPV-: %d TFs robust (both views FDR < 0.05): %d higher in HPV+, %d in HPV-",
                sum(tf_hpv$robust), sum(tf_hpv$robust & tf_hpv$higher_in == "HPV+"),
                sum(tf_hpv$robust & tf_hpv$higher_in == "HPV-")))
message("Top robust HPV+: ", paste(head(tf_hpv$TF[tf_hpv$robust & tf_hpv$higher_in == "HPV+"], 15), collapse = ", "))
message("Top robust HPV-: ", paste(head(tf_hpv$TF[tf_hpv$robust & tf_hpv$higher_in == "HPV-"], 15), collapse = ", "))
message("\nPositive controls (expected higher in HPV+):")
print(tf_hpv[tf_hpv$TF %in% controls, c("TF", "score", "FDR_contrast", "t", "FDR_patient", "robust")],
      row.names = FALSE, digits = 2)

## ---- Figure 10 -------------------------------------------------------------
base_size <- 7
ink <- "#0b0b0b"; ink2 <- "#52514e"; rule <- "#c3c2b7"; grid_col <- "#e1e0d9"
hpv_cols <- c("HPV+" = "#2a78d6", "HPV-" = "#eb6834")
col_tumor <- "#e34948"; col_norm <- "#4a3aa7"
theme_fig <- function() {
  theme_classic(base_size = base_size) +
    theme(text = element_text(colour = ink), axis.text = element_text(colour = ink2),
          axis.line = element_line(linewidth = 0.3, colour = rule),
          axis.ticks = element_line(linewidth = 0.3, colour = rule),
          plot.title = element_text(size = base_size + 1, face = "bold"),
          plot.subtitle = element_text(size = base_size - 0.5, colour = ink2),
          plot.tag = element_text(size = base_size + 3, face = "bold"),
          strip.background = element_blank(),
          strip.text = element_text(size = base_size - 0.5, face = "bold"),
          legend.text = element_text(size = base_size - 0.5), legend.key.size = unit(3, "mm"),
          plot.margin = margin(4, 6, 4, 4))
}

## a: top HPV-associated TFs (contrast score), robust ones solid
top_n <- 12
top <- rbind(head(tf_hpv[tf_hpv$robust & tf_hpv$higher_in == "HPV+", ], top_n),
             head(tf_hpv[tf_hpv$robust & tf_hpv$higher_in == "HPV-", ], top_n))
top$label <- ifelse(top$TF %in% controls, paste0(top$TF, " *"), top$TF)
top$label <- factor(top$label, levels = top$label[order(top$score)])
pa <- ggplot(top, aes(score, label, fill = higher_in)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = rule) +
  scale_fill_manual(values = hpv_cols, name = "Higher in") +
  labs(title = "TF activity: HPV+ vs HPV- tumour cells",
       subtitle = sprintf("top robust TFs (both views FDR < 0.05; %d in total); * = positive control",
                          sum(tf_hpv$robust)),
       x = "Activity score (ULM on DGE statistics)", y = NULL) +
  theme_fig() + theme(axis.text.y = element_text(colour = ink))

## b: positive controls per patient
ctl <- intersect(controls, rownames(act_mat))
dfb <- data.frame(TF = rep(ctl, times = ncol(act_mat)),
                  HPV = factor(rep(res$sB$HPV, each = length(ctl)), levels = names(hpv_cols)),
                  score = as.vector(act_mat[ctl, ]))
dfb$TF <- factor(dfb$TF, levels = ctl)
pb <- ggplot(dfb, aes(HPV, score)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, colour = ink2, linewidth = 0.3, fill = NA) +
  geom_point(aes(colour = HPV), position = position_jitter(width = 0.15, seed = 1),
             size = 0.8, stroke = 0, shape = 16) +
  facet_wrap(~ TF, nrow = 1) +
  scale_colour_manual(values = hpv_cols, guide = "none") +
  labs(title = "Positive controls per patient (E7 -> RB -> E2F axis)",
       subtitle = "expected higher in HPV+; one point = one patient's tumour cells",
       x = NULL, y = "TF activity") +
  theme_fig()

## c: tumour vs normal squamous (A1): top TFs either way
a1 <- contrast_tf[contrast_tf$contrast == "A1" & contrast_tf$FDR < 0.05, ]
a1 <- rbind(head(a1[a1$score > 0, ][order(-a1$score[a1$score > 0]), ], top_n),
            head(a1[a1$score < 0, ][order(a1$score[a1$score < 0]), ], top_n))
a1$dir <- ifelse(a1$score > 0, "Tumour", "Normal squamous")
a1$TF <- factor(a1$TF, levels = a1$TF[order(a1$score)])
pc <- ggplot(a1, aes(score, TF, fill = dir)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = rule) +
  scale_fill_manual(values = c(Tumour = col_tumor, "Normal squamous" = col_norm), name = "Higher in") +
  labs(title = "TF activity: tumour vs normal squamous",
       subtitle = "all paired patients (A1), FDR < 0.05", x = "Activity score", y = NULL) +
  theme_fig() + theme(axis.text.y = element_text(colour = ink))

fig10 <- (pa | pc) / pb + plot_layout(heights = c(1.6, 1)) + plot_annotation(tag_levels = "a")
for (ext in c("png", "pdf")) {
  ggsave(file.path(out_dir, paste0("Fig10_TF_activity.", ext)), fig10, width = 180, height = 170,
         units = "mm", dpi = 600, bg = "white", device = if (ext == "png") ragg::agg_png else "pdf")
}
message("Wrote ", file.path(out_dir, "Fig10_TF_activity.png/.pdf"))
