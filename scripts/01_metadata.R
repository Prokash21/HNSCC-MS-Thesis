## Step 01: clinical metadata and data availability
##
## Purpose: per-patient clinical metadata (site, subsite, HPV and how it was
##   tested, stage, smoking, treatment) for the 78 HNSCatlas patients; normal
##   and leukoplakia samples of GSE181919; sample summaries of GSE172577 and
##   GSE225331; TCGA-HNSC site x HPV counts, normals and outcome data; and the
##   feasibility of site-specific analysis.
## Inputs:  archive_v1_HPVpos/results/01_atlas/tables/patient_metadata.csv
##          GEO series matrices and small supplementary metadata files
##          (downloaded to data/raw/step01_metadata/, each < 50 MB)
##          cBioPortal (hnsc_tcga_pan_can_atlas_2018, hnsc_tcga) and GDC APIs
##          Supplementary tables of Kurten 2021 and Puram 2023, transcribed below
## Outputs: results/step01_metadata/
##            patient_clinical.csv          one row per atlas patient
##            choi_GSE181919_samples.csv    all 37 GSE181919 samples
##            external_scrna_samples.csv    GSE172577 and GSE225331 samples
##            dataset_summary.csv           single-cell datasets overview
##            tcga_hnsc_site_hpv.csv        TCGA-HNSC by site/subsite and HPV
##            site_feasibility.csv          single-cell vs TCGA counts per site
## Date: 2026-10-11
## Run from the repo root:  Rscript scripts/01_metadata.R

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})
source("config/paths.R")

raw_dir <- paths$step01_raw
out_dir <- paths$step01
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
max_bytes <- 50e6

## ---- download helper: skip existing files, refuse anything > 50 MB ---------
geo_url <- function(gse, sub, file) {
  sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%snnn/%s/%s/%s",
          substr(gse, 1, nchar(gse) - 3), gse, sub, file)
}
fetch <- function(url, dest = file.path(raw_dir, basename(url))) {
  if (file.exists(dest)) return(dest)
  h <- curlGetHeaders(url)
  len <- as.numeric(sub(".*:\\s*", "", grep("^content-length", h, ignore.case = TRUE, value = TRUE)))
  if (length(len) && max(len) > max_bytes) stop("Refusing to download ", url, " (", max(len), " bytes)")
  download.file(url, dest, mode = "wb", quiet = TRUE)
  dest
}
## size of a remote file without downloading it
remote_size <- function(url) {
  h <- curlGetHeaders(url)
  len <- as.numeric(sub(".*:\\s*", "", grep("^content-length", h, ignore.case = TRUE, value = TRUE)))
  if (length(len)) max(len) else NA_real_
}

## series matrix -> one row per GSM with every characteristics field
read_series <- function(gse) {
  f <- fetch(geo_url(gse, "matrix", paste0(gse, "_series_matrix.txt.gz")))
  l <- readLines(gzfile(f))
  s <- grep("^!Sample_(title|geo_accession|source_name_ch1|characteristics_ch1)\t", l, value = TRUE)
  sp <- strsplit(s, "\t")
  m <- do.call(cbind, lapply(sp, function(v) gsub("\"", "", v[-1])))
  key <- sub("^!Sample_", "", vapply(sp, `[`, "", 1))
  ch <- key == "characteristics_ch1"
  key[ch] <- tolower(gsub("[^A-Za-z]+", "_", sub(":.*", "", m[1, ch])))
  m[, ch] <- sub("^[^:]*:\\s*", "", m[, ch])
  colnames(m) <- make.unique(key)
  as.data.table(m)
}

## =========================================================================
## GOAL A: clinical metadata for the 78 atlas patients
## =========================================================================
pm <- fread(paths$patient_metadata)

## ---- Bill 2023, GSE234933 (atlas BHN<x> = Bill HN<x>) ---------------------
bill_f <- fetch(geo_url("GSE234933", "suppl", "GSE234933_MGH_HNSCC_sample_annotation.txt.gz"))
bill <- fread(bill_f)
setnames(bill, c("sample", "sex", "age", "hpv", "smoking", "site", "specimen"))
bill[, patient_id := paste0("B", sample)]
bill_src <- "GEO GSE234933 suppl. GSE234933_MGH_HNSCC_sample_annotation.txt.gz (Sex, Age, HPV-Status, Smoking history, Original anatomic site, Anatomic location of scRNA-seq specimen)"

## ---- Choi 2023, GSE181919 (atlas P<n> = Choi patient P<n>, sample C<nn>) ---
## header row lacks a name for the barcode column
choi_cells <- fread(fetch(geo_url("GSE181919", "suppl", "GSE181919_Barcode_metadata.txt.gz")),
                    skip = 1, header = FALSE,
                    col.names = c("barcode", "patient.id", "sample.id", "Gender", "Age",
                                  "tissue.type", "subsite", "hpv", "cell.type"))
choi_gsm <- read_series("GSE181919")
choi_gsm[, sample_id := sub("^.* ", "", title)]
choi_src <- "GEO GSE181919 sample characteristics (subsite, disease stage, hpv status) + suppl. GSE181919_Barcode_metadata.txt.gz (patient.id, Gender, Age)"

## ---- Kurten 2021, GSE164690: no clinical fields in GEO ---------------------
## Transcribed from Kurten et al. 2021, Nat Commun 12:7338, Supplementary
## Information, Supplementary Table 1 (HPV by p16 IHC; TNM AJCC 8th, clinical
## except patients 5 and 8 = pathological). GEO: all 18 treatment-naive.
kurten <- fread(text = "
patient_id,sex,age,smoking,alcohol,site,T,N,M,hpv
HN01,M,70-79,Yes,No,Oral cavity,T4a,N2b,M0,Neg
HN02,F,60-69,No,No,Oral cavity,T3,N2a,M0,Neg
HN03,M,80-89,No,No,Oral cavity,T4a,N0,M0,Neg
HN04,M,50-59,Yes,Yes,Oral cavity,T3,N1,M0,Neg
HN05,F,50-59,Yes,Yes,Oral cavity,pT3,pN3b,pM0,Neg
HN06,M,30-39,Yes,Yes,Oral cavity,T3,N0,M0,Neg
HN07,M,60-69,Yes,Yes,Larynx,T3,N0,M0,Neg
HN08,F,70-79,Yes,Yes,Oral cavity,pT1,pN0,pM0,Neg
HN09,F,70-79,Yes,Yes,Oral cavity,T3,N2b,M0,Neg
HN10,M,50-59,No,Yes,Oral cavity,T3,N0,M0,Neg
HN11,M,80-89,No,No,Oral cavity,T2,N0,M0,Neg
HN12,M,50-59,Yes,Yes,Oropharynx,T2,N1,M0,Pos
HN13,M,70-79,No,No,Oropharynx,T2,N0,M0,Pos
HN14,M,50-59,Yes,Yes,Oropharynx,T1,N1,M0,Pos
HN15,F,60-69,Yes,NA,Oral cavity,T2,N0,M0,Neg
HN16,M,40-49,Yes,Yes,Oropharynx,T2,N1,M0,Pos
HN17,M,50-59,Yes,Yes,Oropharynx,T1,N1,M0,Pos
HN18,M,50-59,Yes,Yes,Oropharynx,T2,N2,M0,Pos")
kurten_src <- "Kurten 2021 Nat Commun (doi:10.1038/s41467-021-27619-4) Supplementary Table 1; treatment-naive per GEO GSE164690 overall design"

## ---- Puram 2023, GSE182227: site from GEO; rest from the paper -------------
## Transcribed from Puram et al. 2023, Nat Genet 55:640, Supplementary Tables
## S1 (clinical) and S2 (pathology, p16 IHC). Patient IDs OPSCC<n> = OP<n>.
puram <- fread(text = "
patient_id,age,sex,subsite,clin_stage,hpv_related,prior_treatment,smoking,alcohol,p16,path_stage,histology
OP4,69,Male,Base of tongue,T2N1M0,Yes,No,Former,Yes,P16+,T2N1,NK SCC
OP5,58,Male,Base of tongue,T3N1M0,Yes,No,Current,No,P16+,T2N2,NK SCC
OP6,60,Female,Tonsil,T3N0M0,Yes,No,Former,Yes,P16+,T4N0,NK SCC
OP8,69,Male,Base of tongue,T3N0M0,No,Yes,Former,No,P16+ (HPV-),Biopsy,NK SCC
OP9,50,Male,Tonsil,T2N1M0,Yes,No,Never,Yes,P16+,T3N2,NK SCC
OP10,63,Male,Base of tongue,T2N2bM0,No,No,Current,No,P16-,T2N0,K SCC
OP12,81,Male,Tonsil,T2N1M0,No,No,Current,No,P16-,T2N2b,K SCC
OP13,78,Male,Base of tongue,T2N0M0,Yes,Yes,Former,No,P16+,Biopsy,NK SCC
OP14,60,Male,Base of tongue,T2N1M0,Yes,No,Never,Yes,P16+,T2N1,NK SCC
OP16,64,Male,Soft palate,T1N1M0,No,No,Former,Yes,P16-,T1N1,K SCC
OP17,49,Male,Tonsil,T2N3M0,Yes,No,Never,Yes,P16+,T2N2,NK SCC
OP19,67,Female,Tonsil,T2N1M0,No,No,Never,No,P16-,T2N2,K SCC
OP20,55,Male,Base of tongue,T1N1M0,Yes,No,Current,No,P16+,T1N1,NK SCC
OP33,56,Male,Tonsil,T1N1M0,Yes,No,Former,Yes,P16+,T2N2,NK SCC
OP34,74,Male,Tonsil,T2N1M0,Yes,No,Never,No,P16+,T2N1,NK SCC
OP35,57,Male,Base of tongue,T1N1M0,Yes,No,Former,Yes,P16+,T2N1,NK SCC")
puram_gsm <- read_series("GSE182227")
puram_src <- "GEO GSE182227 sample characteristics (source name Oropharynx, tissue, hpv status) + Puram 2023 Nat Genet (doi:10.1038/s41588-023-01357-3) Supplementary Tables S1-S2"
## which GEO samples each atlas patient's cells come from (barcode prefixes)
puram_bc <- fread(fetch(geo_url("GSE182227", "suppl", "GSE182227_barcodes.txt.gz")), header = FALSE)
puram_bc[, sample := sub("_[ACGT]+-1$", "", V1)]
puram_n <- puram_bc[, .N, by = sample]

## ---- assemble one row per atlas patient ------------------------------------
site_group <- function(x) {
  fifelse(grepl("^Oral", x), "Oral cavity",
  fifelse(grepl("Oropharynx|Base of tongue|Tonsil|Soft palate", x), "Oropharynx",
  fifelse(x == "Hypopharynx", "Hypopharynx",
  fifelse(x == "Larynx", "Larynx",
  fifelse(x == "Larynx/Hypopharynx", "Larynx/Hypopharynx",
  fifelse(x == "Nasopharynx", "Nasopharynx", NA_character_))))))
}
tnm <- function(s, part) {
  r <- regmatches(s, regexpr(paste0(part, "[0-9X][a-c]?"), s))
  out <- rep(NA_character_, length(s)); out[grepl(paste0(part, "[0-9X]"), s)] <- r; out
}

rows <- list()
for (i in seq_len(nrow(pm))) {
  p <- pm$Patient[i]
  r <- list(patient_id = p, study = pm$Study[i], hpv_atlas = pm$HPV[i], sex_atlas = pm$Sex[i],
            tumour_cells = pm$Tumor_cells[i], nonmalignant_epi_cells = pm$NonMalignant_epithelial_cells[i],
            source_sample_ids = NA, sex = NA, age = NA, site = NA, subsite = NA,
            hpv_source = NA, hpv_method = NA, T = NA, N = NA, M = NA, stage_overall = NA,
            grade_histology = NA, smoking = NA, alcohol = NA, prior_treatment = NA,
            specimen = NA, recurrence_or_survival = NA, source_of_info = NA, notes = NA)
  if (startsWith(p, "BHN")) {
    b <- bill[patient_id == p]
    r[c("source_sample_ids", "sex", "age", "site", "smoking", "specimen", "source_of_info")] <-
      list(b$sample, b$sex, as.character(b$age), b$site, b$smoking, b$specimen, bill_src)
    r$hpv_source <- ifelse(b$hpv == "Positive", "HPV+", "HPV-")
    r$hpv_method <- "UNCLEAR (not stated in GEO annotation)"
    if (b$site == "Nasopharynx") r$notes <- "Nasopharynx primary"
    if (b$site == "Unknown") r$notes <- "Unknown primary of the head and neck"
  } else if (grepl("^P[0-9]", p)) {
    c1 <- unique(choi_cells[patient.id == p & tissue.type == "CA", .(sample.id, Gender, Age, subsite, hpv)])
    g <- choi_gsm[sample_id == c1$sample.id]
    r[c("source_sample_ids", "sex", "age", "hpv_source", "specimen", "source_of_info")] <-
      list(paste0(c1$sample.id, " (", g$geo_accession, ")"), c1$Gender, as.character(c1$Age), c1$hpv,
           "Primary cancer (CA)", choi_src)
    r$site <- c(OC = "Oral cavity", OP = "Oropharynx", HP = "Hypopharynx")[[c1$subsite]]
    r$hpv_method <- "UNCLEAR (not stated in GEO or paper main text)"
    other <- unique(choi_cells[patient.id == p & tissue.type != "CA", sample.id])
    if (length(other)) r$notes <- paste("Also in GEO, not in atlas:", paste(other, collapse = ", "))
  } else if (startsWith(p, "HN")) {
    k <- kurten[patient_id == p]
    r[c("sex", "age", "site", "smoking", "alcohol", "T", "N", "M", "source_of_info")] <-
      list(k$sex, k$age, k$site, k$smoking, k$alcohol, k$T, k$N, k$M, kurten_src)
    r$source_sample_ids <- paste(p, c("CD45+", "CD45-"), collapse = "; ")
    r$hpv_source <- ifelse(k$hpv == "Pos", "HPV+", "HPV-")
    r$hpv_method <- "p16 IHC"
    r$prior_treatment <- "No (treatment-naive)"
    r$specimen <- "Primary tumour, FACS CD45+/CD45- fractions"
  } else if (startsWith(p, "OP")) {
    o <- puram[patient_id == p]
    smp <- puram_n[startsWith(sample, p) & !grepl(paste0("^", p, "[0-9]"), sample)]
    r[c("sex", "age", "subsite", "smoking", "alcohol", "prior_treatment", "source_of_info")] <-
      list(o$sex, as.character(o$age), o$subsite, o$smoking, o$alcohol, o$prior_treatment, puram_src)
    r$site <- "Oropharynx"
    r$source_sample_ids <- paste(sprintf("%s (%d cells)", smp$sample, smp$N), collapse = "; ")
    r$hpv_source <- ifelse(o$hpv_related == "Yes", "HPV+", "HPV-")
    r$hpv_method <- paste0("p16 IHC: ", o$p16, "; HPV-related per paper (viral transcripts)")
    r[c("T", "N", "M")] <- list(tnm(o$clin_stage, "T"), tnm(o$clin_stage, "N"), tnm(o$clin_stage, "M"))
    r$grade_histology <- o$histology
    r$specimen <- "Primary tumour"
    if (any(grepl("norm", smp$sample))) r$notes <- "Atlas patient includes adjacent-normal sample cells"
    if (any(grepl("LN", smp$sample))) r$notes <- paste(r$notes, "and metastatic lymph node cells")
    if (o$prior_treatment == "Yes") r$notes <- paste(na.omit(c(r$notes, "Prior HNC treatment per Table S1 (GEO says treatment-naive)")), collapse = "; ")
  }
  rows[[i]] <- r
}
clin <- rbindlist(lapply(rows, function(r) lapply(r, function(v) if (is.null(v) || length(v) == 0) NA else as.character(v))))
clin[, site_group := site_group(fifelse(is.na(subsite), site, subsite))]
clin[site_group == "Oropharynx" & is.na(subsite) & !is.na(site), subsite := NA]
clin[, hpv_match := hpv_atlas == hpv_source]
clin[, sex_match := substr(sex_atlas, 1, 1) == substr(sex, 1, 1)]
clin[, recurrence_or_survival := "NA (not reported per patient in GEO or the tables used)"]
setcolorder(clin, c("patient_id", "study", "source_sample_ids", "sex", "age", "site_group", "site", "subsite",
                    "hpv_atlas", "hpv_source", "hpv_match", "hpv_method", "T", "N", "M", "stage_overall",
                    "grade_histology", "smoking", "alcohol", "prior_treatment", "specimen",
                    "recurrence_or_survival", "tumour_cells", "nonmalignant_epi_cells", "sex_atlas",
                    "sex_match", "notes", "source_of_info"))
setorder(clin, study, patient_id)
fwrite(clin, file.path(out_dir, "patient_clinical.csv"))

message("== Goal A ==")
message("HPV mismatches atlas vs source: ", sum(!clin$hpv_match, na.rm = TRUE),
        " | sex mismatches: ", sum(!clin$sex_match, na.rm = TRUE))
print(clin[, .N, by = .(hpv_atlas, site_group)][order(hpv_atlas, site_group)])
clin[, tumour_cells := as.integer(tumour_cells)]
print(clin[tumour_cells >= 50, .N, by = .(hpv_atlas, site_group, subsite)][order(hpv_atlas, site_group)])

## Bill samples not in the atlas (recurrences, metastases)
bill[, in_atlas := patient_id %in% pm$Patient]
message("\nBill 2023 samples by specimen type and atlas membership:")
print(bill[, .N, by = .(specimen = sub(" \\(.*", "", specimen), hpv, in_atlas)][order(-in_atlas, specimen)])

## =========================================================================
## GOAL B: GSE181919 (Choi) samples, incl. normal (NL) and leukoplakia (LP)
## =========================================================================
ct <- dcast(choi_cells, sample.id ~ cell.type, fun.aggregate = length)
ct[, total_cells := rowSums(.SD), .SDcols = -1]
choi_s <- unique(choi_cells[, .(sample_id = sample.id, patient_id = patient.id, tissue_type = tissue.type,
                                subsite, hpv, sex = Gender, age = Age)])
choi_s <- merge(choi_s, ct[, .(sample_id = sample.id, total_cells, epithelial_cells = Epithelial.cells,
                               malignant_cells = Malignant.cells)], by = "sample_id")
choi_s <- merge(choi_s, choi_gsm[, .(sample_id, gsm = geo_accession)], by = "sample_id", all.x = TRUE)
choi_s[, tissue_type := c(CA = "Primary cancer", LN = "Lymph node metastasis", LP = "Leukoplakia",
                          NL = "Normal mucosa")[tissue_type]]
choi_s[, site := c(OC = "Oral cavity", OP = "Oropharynx", HP = "Hypopharynx")[subsite]]
choi_s[, patient_in_atlas := patient_id %in% pm$Patient]
choi_s[, sample_in_atlas := tissue_type == "Primary cancer" & patient_in_atlas]
## atlas totals per Choi patient equal their CA sample totals?
chk <- merge(choi_s[tissue_type == "Primary cancer", .(patient_id, total_cells)],
             pm[, .(patient_id = Patient, Total_cells)], by = "patient_id")
choi_s[, subsite := NULL]
setcolorder(choi_s, c("sample_id", "gsm", "patient_id", "tissue_type", "site", "hpv", "sex", "age"))
setorder(choi_s, tissue_type, patient_id)
fwrite(choi_s, file.path(out_dir, "choi_GSE181919_samples.csv"))
message("\n== Goal B ==")
message("Atlas Choi patients whose atlas cell total equals their CA sample: ",
        sum(chk$total_cells == chk$Total_cells), " of ", nrow(chk))
print(choi_s[tissue_type %in% c("Normal mucosa", "Leukoplakia")])
choi_files <- data.table(file = c("GSE181919_UMI_counts.txt.gz", "GSE181919_Barcode_metadata.txt.gz"))
choi_files[, bytes := vapply(file, function(f) remote_size(geo_url("GSE181919", "suppl", f)), 0)]
print(choi_files)

## =========================================================================
## GOAL C: GSE172577 (Peng) and GSE225331 (Quah) from GEO metadata only
## =========================================================================
peng <- read_series("GSE172577")
peng_fl <- fread(fetch(geo_url("GSE172577", "suppl", "filelist.txt"), file.path(raw_dir, "GSE172577_filelist.txt")))
peng_fl[, gsm := sub("_.*", "", Name)]
peng_sz <- peng_fl[Type != "TAR", .(bytes = sum(Size), files = paste(sub("^GSM[0-9]+_[^_]+_", "", Name), collapse = "+")), by = gsm]
peng_s <- merge(peng[, .(gsm = geo_accession, sample = title, tissue, treatment)], peng_sz, by = "gsm")
peng_s[, `:=`(dataset = "GSE172577 (Peng 2021)", patient_id = sample, site = "Oral cavity (OSCC, subsite NA)",
              hpv = NA_character_, sample_type = "Tumour", format = "10x mtx (raw UMI counts)")]

quah <- read_series("GSE225331")
quah[, line := sub(" \\[.*", "", title)]
quah_s <- quah[, .(n_samples = .N), by = .(line, patient_id, primary_site, stage, grade, site_of_origin, cell_line)]
quah_s[, `:=`(dataset = "GSE225331 (Quah 2023)", hpv = NA_character_,
              sample_type = paste(cell_line, "-", site_of_origin),
              format = "RData (good_counts: counts; good_tpm: TPM)")]
quah_files <- data.table(file = c("GSE225331_good_counts.RData.gz", "GSE225331_good_tpm.RData.gz",
                                  "GSE225331_good_cellinfo.RData.gz"))
quah_files[, bytes := vapply(file, function(f) remote_size(geo_url("GSE225331", "suppl", f)), 0)]

ext <- rbind(
  peng_s[, .(dataset, sample = gsm, patient_id, site, hpv, sample_type, n_cells_or_wells = NA_integer_,
             format, bytes, notes = paste(tissue, treatment, sep = "; "))],
  quah_s[, .(dataset, sample = line, patient_id, site = paste0("Oral cavity (", primary_site, ")"), hpv,
             sample_type, n_cells_or_wells = n_samples, format, bytes = NA_real_,
             notes = paste("Stage", stage, "| grade", grade))])
fwrite(ext, file.path(out_dir, "external_scrna_samples.csv"))
message("\n== Goal C ==")
print(ext[, .(dataset, sample, patient_id, site, sample_type, n_cells_or_wells, bytes)])
print(quah_files)
message("GSE172577_RAW.tar bytes: ", remote_size(geo_url("GSE172577", "suppl", "GSE172577_RAW.tar")))

zen <- fromJSON("https://zenodo.org/api/records/17880437")
message("Zenodo 17880437: ", zen$metadata$title)
print(data.table(file = zen$files$key, GB = round(zen$files$size / 1e9, 2)))

summ <- data.table(
  dataset = c("HNSCatlas (4 studies)", "GSE181919 Choi NL/LP (not in atlas)", "GSE172577 Peng",
              "GSE225331 Quah", "Kroehling/Monti HPV- atlas (Zenodo 17880437)"),
  hpv_neg_patients = c(sum(pm$HPV == "HPV-"), choi_s[tissue_type %in% c("Normal mucosa", "Leukoplakia") & hpv == "HPV-", uniqueN(patient_id)],
                       NA, NA, NA),
  hpv_pos_patients = c(sum(pm$HPV == "HPV+"), choi_s[tissue_type %in% c("Normal mucosa", "Leukoplakia") & hpv == "HPV+", uniqueN(patient_id)],
                       NA, NA, 0),
  tumour_cells = c(sum(pm$Tumor_cells), 0, NA, sum(quah_s$n_samples), NA),
  sites = c("OC, OP, larynx/hypopharynx, HP, NP, unknown", "OC, OP, HP", "Oral cavity", "Oral cavity (tongue, RMT)", "UNCLEAR"),
  normal_samples = c("No (within-tumour epithelium; Puram OP33-35 normals merged in)", "Yes (9 NL, 4 LP)", "No", "No", "UNCLEAR"),
  raw_counts = c("Yes", "Yes (UMI text matrix)", "Yes (10x mtx)", "Yes (RData counts)", "UNCLEAR"),
  in_atlas = c("Yes", "No (patients yes for 9 NL + LP15)", "No", "No", "No"))
fwrite(summ, file.path(out_dir, "dataset_summary.csv"))

## =========================================================================
## GOAL D: TCGA-HNSC site x HPV, normals, outcome data; site feasibility
## =========================================================================
cbio <- function(study) {
  f <- file.path(raw_dir, paste0("cbioportal_", study, "_patient.csv"))
  if (!file.exists(f)) {
    d <- as.data.table(fromJSON(sprintf(
      "https://www.cbioportal.org/api/studies/%s/clinical-data?clinicalDataType=PATIENT&projection=SUMMARY&pageSize=100000", study)))
    fwrite(dcast(d, patientId ~ clinicalAttributeId, value.var = "value"), f)
  }
  fread(f)
}
pan <- cbio("hnsc_tcga_pan_can_atlas_2018")   # SUBTYPE (HPV), PFS = PFI, DFS = DFI
leg <- cbio("hnsc_tcga")                      # PRIMARY_SITE_PATIENT (subsite), p16, ISH

nf <- file.path(raw_dir, "gdc_tcga_hnsc_normal_rnaseq_cases.txt")
if (!file.exists(nf)) {
  flt <- list(op = "and", content = list(
    list(op = "in", content = list(field = "cases.project.project_id", value = list("TCGA-HNSC"))),
    list(op = "in", content = list(field = "data_type", value = list("Gene Expression Quantification"))),
    list(op = "in", content = list(field = "analysis.workflow_type", value = list("STAR - Counts"))),
    list(op = "in", content = list(field = "cases.samples.sample_type", value = list("Solid Tissue Normal")))))
  u <- paste0("https://api.gdc.cancer.gov/files?size=2000&format=json&fields=cases.submitter_id&filters=",
              URLencode(as.character(toJSON(flt, auto_unbox = TRUE)), reserved = TRUE))
  hits <- fromJSON(u, simplifyVector = FALSE)$data$hits
  writeLines(unique(vapply(hits, function(h) h$cases[[1]]$submitter_id, "")), nf)
}
normals <- readLines(nf)

tcga <- merge(pan[, .(patientId, SUBTYPE, PFS_STATUS, DFS_STATUS, ICD_O_3_SITE)],
              leg[, .(patientId, subsite = PRIMARY_SITE_PATIENT, HPV_STATUS_P16, HPV_STATUS_ISH)],
              by = "patientId", all.x = TRUE)
tcga[, hpv := fifelse(SUBTYPE == "HNSC_HPV+", "HPV+", fifelse(SUBTYPE == "HNSC_HPV-", "HPV-", NA_character_))]
tcga[, site_group := fcase(
  subsite %in% c("Oral Tongue", "Floor of mouth", "Buccal Mucosa", "Alveolar Ridge", "Hard Palate", "Lip", "Oral Cavity"), "Oral cavity",
  subsite %in% c("Base of tongue", "Tonsil", "Oropharynx"), "Oropharynx",
  subsite == "Hypopharynx", "Hypopharynx",
  subsite == "Larynx", "Larynx",
  default = NA_character_)]
tcga[, has_normal := patientId %in% normals]
tab_t <- tcga[, .(tumours = .N,
                  hpv_neg = sum(hpv == "HPV-", na.rm = TRUE), hpv_pos = sum(hpv == "HPV+", na.rm = TRUE),
                  hpv_na = sum(is.na(hpv)),
                  normals_rnaseq = sum(has_normal),
                  pfi_available = sum(!is.na(PFS_STATUS) & PFS_STATUS != ""),
                  pfi_events = sum(grepl("^1", PFS_STATUS)),
                  dfi_available = sum(!is.na(DFS_STATUS) & DFS_STATUS != ""),
                  dfi_events = sum(grepl("^1", DFS_STATUS))),
              by = .(site_group, subsite)][order(site_group, subsite)]
fwrite(tab_t, file.path(out_dir, "tcga_hnsc_site_hpv.csv"))
message("\n== Goal D: TCGA ==")
print(tab_t)
message("TCGA patients: ", nrow(tcga), " | normal RNA-seq cases: ", length(normals),
        " (", sum(tcga$has_normal), " matched to clinical rows)")

## single-cell vs TCGA per site group
clin[, tumour_cells := as.integer(tumour_cells)]
sc <- clin[, .(sc_hpv_neg_all = sum(hpv_atlas == "HPV-"),
               sc_hpv_neg_50 = sum(hpv_atlas == "HPV-" & tumour_cells >= 50),
               sc_hpv_pos_all = sum(hpv_atlas == "HPV+"),
               sc_hpv_pos_50 = sum(hpv_atlas == "HPV+" & tumour_cells >= 50)),
           by = .(site_group = fcoalesce(site_group, "Unknown"))]
tg <- tcga[, .(tcga_hpv_neg = sum(hpv == "HPV-", na.rm = TRUE), tcga_hpv_pos = sum(hpv == "HPV+", na.rm = TRUE),
               tcga_normals = sum(has_normal), tcga_pfi_events = sum(grepl("^1", PFS_STATUS))),
           by = .(site_group = fcoalesce(site_group, "Unknown"))]
feas <- merge(sc, tg, by = "site_group", all = TRUE)
feas[, sc_feasible_hpv_neg := fcoalesce(sc_hpv_neg_50, 0L) >= 5]
fwrite(feas, file.path(out_dir, "site_feasibility.csv"))
print(feas)
message("\nDone. Tables in ", out_dir)
