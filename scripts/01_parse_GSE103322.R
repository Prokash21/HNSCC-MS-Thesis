## Parse GSE103322 (Puram et al. 2017, Cell) raw supplementary file into
## analysis-ready objects: a Seurat object plus flat CSV/matrix exports.
##
## Source file layout (GSE103322_HNSCC_all_data.txt.gz):
##   col 1      = row label (blank for header row, else gene symbol)
##   cols 2..N  = cells
##   rows 1     = cell IDs (header)
##   row 2      = "processed by Maxima enzyme" (all 1, technical flag)
##   row 3      = "Lymph node" (0/1)
##   row 4      = "classified as cancer cell" (0/1)
##   row 5      = "classified as non-cancer cells" (0/1)
##   row 6      = "non-cancer cell type" (string, or 0 if cancer/unresolved)
##   rows 7..end = genes, values = log2(TPM/10 + 1)

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
})

raw_path <- "data/raw/GSE103322/GSE103322_HNSCC_all_data.txt.gz"
out_dir  <- "data/processed/GSE103322"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading raw file (this can take a minute)...")
dt <- fread(raw_path, header = FALSE, sep = "\t", quote = "", data.table = TRUE)

cell_ids <- as.character(unlist(dt[1, -1]))
lymph_node        <- as.integer(unlist(dt[3, -1]))
classified_cancer <- as.integer(unlist(dt[4, -1]))
classified_normal <- as.integer(unlist(dt[5, -1]))
noncancer_type    <- as.character(unlist(dt[6, -1]))

## derive one clean cell-type label per cell
cell_type <- ifelse(
  classified_cancer == 1, "malignant",
  ifelse(classified_normal == 1,
         paste0("non-malignant:", gsub("^-", "", noncancer_type)),
         "unresolved")
)

patient_id <- sub("_.*$", "", cell_ids)

meta <- data.table(
  cell_id           = cell_ids,
  patient_id        = patient_id,
  lymph_node        = lymph_node,
  classified_cancer = classified_cancer,
  classified_normal = classified_normal,
  noncancer_type    = noncancer_type,
  cell_type         = cell_type
)

message("Cell type breakdown:")
print(table(meta$cell_type))

## expression block: rows 7..end, col 1 = gene symbol, cols 2..N = values
expr_dt <- dt[7:nrow(dt)]
genes <- gsub("^'|'$", "", as.character(expr_dt[[1]]))
expr_mat <- as.matrix(expr_dt[, -1, with = FALSE])
mode(expr_mat) <- "numeric"
rownames(expr_mat) <- genes
colnames(expr_mat) <- cell_ids

## drop duplicate gene symbols (keep first occurrence) and empty gene names
dup <- duplicated(rownames(expr_mat)) | rownames(expr_mat) == ""
if (any(dup)) {
  message(sprintf("Dropping %d duplicate/empty gene rows", sum(dup)))
  expr_mat <- expr_mat[!dup, ]
}

message(sprintf("Expression matrix: %d genes x %d cells", nrow(expr_mat), ncol(expr_mat)))

## --- save flat exports (portable across R/Python) ---
fwrite(meta, file.path(out_dir, "cell_metadata.csv"))

expr_sparse <- Matrix(expr_mat, sparse = TRUE)
saveRDS(expr_sparse, file.path(out_dir, "expr_log2tpm_matrix.rds"))

## gene x cell log2(TPM/10+1) matrix as gzipped CSV for Python/pySCENIC use
fwrite(
  data.table(gene = rownames(expr_mat), expr_mat),
  file.path(out_dir, "expr_log2tpm_matrix.csv.gz")
)

## --- build a Seurat object for downstream R analysis ---
## Seurat expects non-negative "counts"-like data in `counts` slot; since this
## is already log2(TPM/10+1) (not raw counts), store it in both slots and skip
## NormalizeData — flag this clearly for downstream steps.
meta_df <- as.data.frame(meta)
rownames(meta_df) <- meta_df$cell_id

seu <- CreateSeuratObject(
  counts   = expr_sparse,
  meta.data = meta_df,
  project  = "GSE103322_HNSCC"
)
seu <- SetAssayData(seu, layer = "data", new.data = expr_sparse)

saveRDS(seu, file.path(out_dir, "GSE103322_seurat.rds"))

message("Done. Outputs written to ", out_dir)
message(" - cell_metadata.csv")
message(" - expr_log2tpm_matrix.rds  (sparse dgCMatrix, genes x cells)")
message(" - expr_log2tpm_matrix.csv.gz (portable, genes x cells)")
message(" - GSE103322_seurat.rds  (Seurat object; counts/data slots hold log2(TPM/10+1), NOT raw counts)")
