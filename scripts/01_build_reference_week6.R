# scripts/01_build_reference_week6.R
# Purpose: Build a Week-6 (HE6W) reference Seurat object from GSE106118 merged counts
# Inputs:
#   data/reference/GSE106118_barcode_information.txt.gz
#   data/reference/GSE106118_UMI_count_merge.txt.gz
# Output:
#   data/subsets/ref6_week6.rds
# Notes: No raw FASTQs here; this uses GEO-provided processed matrices.

suppressPackageStartupMessages({
  library(data.table)
  library(Seurat)
  library(Matrix)
})

# ---- Config (use relative paths so the repo is portable) ---------------------
meta_path   <- "data/reference/GSE106118_barcode_information.txt.gz"
counts_path <- "data/reference/GSE106118_UMI_count_merge.txt.gz"

if (!file.exists(meta_path) || !file.exists(counts_path)) {
  stop("Reference files not found. Place them under data/reference/ (see data/README_data.md).")
}

# ---- Load --------------------------------------------------------------------
metadata <- fread(meta_path, showProgress = TRUE)
counts   <- fread(counts_path, showProgress = TRUE)

# Sanity checks
req_cols <- c("cell", "processed_file")
missing  <- setdiff(req_cols, names(metadata))
if (length(missing)) stop("Metadata missing required column(s): ", paste(missing, collapse=", "))

# ---- Select Week-6 cells -----------------------------------------------------
wk6_meta <- metadata[grepl("HE6W", processed_file) | grepl("HE6W", cell)]
if (nrow(wk6_meta) == 0) stop("No HE6W rows found. Check metadata$processed_file / $cell patterns.")

# Standardize gene column name in counts
gene_col <- intersect(names(counts), c("gene","Gene","GENE","GENEID","GeneID","symbol","SYMBOL","V1"))
if (!length(gene_col)) gene_col <- names(counts)[1]
setnames(counts, gene_col[1], "gene")

# Match Week-6 cells to counts columns
candidate_cols <- intersect(names(counts), wk6_meta$cell)
if (!length(candidate_cols)) {
  stop("No candidate Week-6 cell columns matched counts. Inspect naming of metadata$cell vs counts colnames.")
}

# Subset to gene + Week-6 cells
week6_counts <- counts[, c("gene", candidate_cols), with = FALSE]

# ---- QC (single pass) --------------------------------------------------------
libsizes <- colSums(week6_counts[, -1])
meanZero <- colMeans(week6_counts[, -1] == 0)

# Basic thresholds (tune as needed)
keep <- libsizes > max(2000, quantile(libsizes, 0.05, na.rm = TRUE)) & meanZero < 0.98
if (!any(keep)) stop("All Week-6 cells filtered out by QC thresholds; relax criteria or inspect data distributions.")

week6_counts_qc <- week6_counts[, c("gene", names(keep)[keep]), with = FALSE]

# Minimal diagnostics
message("Week-6 cells (pre-QC): ", ncol(week6_counts)-1,
        " | post-QC: ", sum(keep),
        " | genes: ", nrow(week6_counts_qc))

# ---- Seurat object from counts table ----------------------------------------
dashify <- function(x) {
  x <- gsub("_","-", x, fixed=TRUE)
  x <- gsub("\\s+","", x)
  make.unique(x)
}

CreateSeuratFromCounts <- function(counts_df, project = "Week6Ref",
                                   min.cells = 3, min.features = 200) {
  stopifnot(ncol(counts_df) >= 2)
  genes_orig <- counts_df[[1]]
  mat <- as.matrix(data.frame(row.names = genes_orig, counts_df[,-1, drop=FALSE]))
  rownames(mat) <- dashify(rownames(mat))
  mat <- Matrix(mat, sparse = TRUE)
  obj <- CreateSeuratObject(mat, project = project, min.cells = min.cells, min.features = min.features)
  obj
}

ref6 <- CreateSeuratFromCounts(week6_counts_qc, project = "Week6Ref")
ref6 <- NormalizeData(ref6, verbose = FALSE)

# Optional: add simple region label from cell IDs (LA/LV/RA/RV) if present
cell_ids <- colnames(ref6)
region <- sub("^HE6W.*_(LA|LV|RA|RV).*", "\\1", cell_ids)
region[!region %in% c("LA","LV","RA","RV")] <- NA
ref6$region <- region

# ---- Save --------------------------------------------------------------------
dir.create("data/subsets", showWarnings = FALSE, recursive = TRUE)
saveRDS(ref6, file = "data/subsets/ref6_week6.rds")

# Reproducibility record
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
writeLines(c(capture.output(sessionInfo())),
           "results/tables/sessionInfo_ref6.txt")

message("Saved: data/subsets/ref6_week6.rds")
