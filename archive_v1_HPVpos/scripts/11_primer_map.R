## Primer map: where each chosen qPCR pair sits on its transcript (zoomed to the
## amplicon), with exon boundaries, so the exon-junction design can be checked
## at a glance.
##
## Input:  HNSCC-HPV/results/05_primers/tables/ (scripts/10_primer_design.py)
## Output: HNSCC-HPV/results/05_primers/Fig9_primer_map.png/.pdf
## Run from the repo root:  Rscript scripts/11_primer_map.R

suppressPackageStartupMessages({
  library(ggplot2)
})

in_dir <- "HNSCC-HPV/results/05_primers"
mp <- read.csv(file.path(in_dir, "tables", "primer_map.csv"))
os <- read.csv(file.path(in_dir, "tables", "primer_order_sheet.csv"))
flank <- 120

base_size <- 7
ink <- "#0b0b0b"; ink2 <- "#52514e"; rule <- "#c3c2b7"

rows <- lapply(os$gene, function(g) {
  m  <- mp[mp$gene == g, ]
  f  <- m[m$kind == "forward", ]
  r  <- m[m$kind == "reverse", ]
  tx <- m[m$kind == "transcript", ]
  lo <- max(0, f$start - flank)
  hi <- min(tx$end, r$end + flank)
  ex <- m[m$kind == "exon" & m$end > lo & m$start < hi, ]
  ex$start <- pmax(ex$start, lo)
  ex$end   <- pmin(ex$end, hi)
  list(ex = transform(ex, x0 = start - f$start, x1 = end - f$start),
       pr = data.frame(gene = g, kind = c("forward", "reverse"),
                       x0 = c(0, r$start - f$start), x1 = c(f$end - f$start, r$end - f$start)),
       amp = data.frame(gene = g, x0 = 0, x1 = r$end - f$start))
})
ex  <- do.call(rbind, lapply(rows, `[[`, "ex"))
pr  <- do.call(rbind, lapply(rows, `[[`, "pr"))
amp <- do.call(rbind, lapply(rows, `[[`, "amp"))
strip <- setNames(sprintf("%s  (%s, %d bp)", os$gene, os$role, os$product_bp), os$gene)
for (d in c("ex", "pr", "amp")) {
  x <- get(d)
  x$gene <- factor(strip[x$gene], levels = strip)
  assign(d, x)
}
ex$shade <- factor(as.integer(ex$label) %% 2)

p <- ggplot() +
  geom_rect(data = ex, aes(xmin = x0, xmax = x1, ymin = -0.25, ymax = 0.25, fill = shade),
            colour = "white", linewidth = 0.4) +
  geom_text(data = ex, aes(x = (x0 + x1) / 2, y = 0, label = paste0("exon ", label)),
            size = (base_size - 1.5) / .pt, colour = ink2) +
  geom_segment(data = amp, aes(x = x0, xend = x1, y = -0.42, yend = -0.42),
               colour = "#e34948", linewidth = 1.2) +
  geom_segment(data = pr[pr$kind == "forward", ], aes(x = x0, xend = x1, y = 0.4, yend = 0.4),
               colour = ink, linewidth = 0.5, arrow = arrow(length = unit(1.2, "mm"), type = "closed")) +
  geom_segment(data = pr[pr$kind == "reverse", ], aes(x = x1, xend = x0, y = 0.4, yend = 0.4),
               colour = ink, linewidth = 0.5, arrow = arrow(length = unit(1.2, "mm"), type = "closed")) +
  facet_wrap(~ gene, ncol = 2, scales = "free_x") +
  scale_fill_manual(values = c("0" = "#d6d5cf", "1" = "#c3c2b7"), guide = "none") +
  scale_y_continuous(limits = c(-0.55, 0.55)) +
  labs(title = "qPCR primer positions on the canonical transcript",
       subtitle = "arrows = forward / reverse primer; red = amplicon; boxes = exons (a primer crossing a box edge spans an exon junction)",
       x = "Position relative to the forward primer (bp)", y = NULL) +
  theme_classic(base_size = base_size) +
  theme(text = element_text(colour = ink), axis.text.x = element_text(colour = ink2),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        axis.line.x = element_line(linewidth = 0.3, colour = rule),
        strip.background = element_blank(), strip.text = element_text(face = "bold", hjust = 0),
        plot.title = element_text(size = base_size + 1, face = "bold"),
        plot.subtitle = element_text(size = base_size - 0.5, colour = ink2))

for (ext in c("png", "pdf")) {
  ggsave(file.path(in_dir, paste0("Fig9_primer_map.", ext)), p, width = 180, height = 150,
         units = "mm", dpi = 600, bg = "white", device = if (ext == "png") ragg::agg_png else "pdf")
}
message("Wrote ", file.path(in_dir, "Fig9_primer_map.png/.pdf"))
