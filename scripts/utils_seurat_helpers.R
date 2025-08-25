# scripts/utils_seurat_helpers.R
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(data.table) })

dashify <- function(x) {
  x <- gsub("_","-", x, fixed = TRUE)
  x <- gsub("\\s+","", x)
  make.unique(x)
}

# Build a Seurat object from a 10x-style folder (matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz)
CreateSeuratFrom10x <- function(dir, project = basename(dir),
                                min.cells = 3, min.features = 200) {
  mtx   <- file.path(dir, "matrix.mtx.gz")
  feats <- file.path(dir, "features.tsv.gz")
  bc    <- file.path(dir, "barcodes.tsv.gz")

  stopifnot(file.exists(mtx), file.exists(feats), file.exists(bc))

  feat_df <- fread(feats, header = FALSE)
  if (ncol(feat_df) < 2) stop("features.tsv.gz should have at least 2 columns (feature_id, symbol).")
  colnames(feat_df)[1:2] <- c("feature_id","symbol")

  raw <- ReadMtx(mtx = mtx, features = feats, cells = bc)  # dgCMatrix
  genes_orig <- rownames(raw)
  rownames(raw) <- dashify(genes_orig)

  obj <- CreateSeuratObject(raw, project = project, min.cells = min.cells, min.features = min.features)
  obj$condition <- project
  obj@misc$gene_symbol_map <- data.frame(original = genes_orig, seurat = rownames(raw), stringsAsFactors = FALSE)
  obj@misc$feature_map     <- feat_df
  obj@misc$source_dir      <- dir
  obj
}
