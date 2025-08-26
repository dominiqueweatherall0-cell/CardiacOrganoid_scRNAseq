# scripts/03_gate_endothelial.R
# Purpose: Gate endothelial cells (ECs) from all experimental 10x objects + Week6 reference
# Outputs:
#   data/subsets/{ctrl_ec,emm1_ec,emm2_ec,mm_ec,ref6_ec}.rds
#   results/tables/gating_summary.csv
# Notes:
#   - Portable (relative) paths
#   - Robust mapping of Ensembl-like IDs -> gene symbols if feature_map is available
#   - Uses counts > 0 across a curated marker set; threshold is tunable

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
})

# ---- Config (edit if your layout differs) ------------------------------------
paths <- list(
  ref6_rds = "data/subsets/ref6_week6.rds",
  tenx     = list(
    Control = "data/10x/Control",
    EMM1    = "data/10x/EMM1",
    EMM2    = "data/10x/EMM2",
    MM      = "data/10x/MM"
  ),
  out_dir_subsets = "data/subsets",
  out_dir_tables  = "results/tables"
)

dir.create(paths$out_dir_subsets, recursive = TRUE, showWarnings = FALSE)
dir.create(paths$out_dir_tables,  recursive = TRUE, showWarnings = FALSE)

# ---- Helpers -----------------------------------------------------------------

dashify <- function(x) {
  x <- gsub("_", "-", x, fixed = TRUE)
  x <- gsub("\\s+", "", x)
  make.unique(x)
}

# Minimal 10x loader that also stashes a feature_map if features.tsv.gz has (ID, symbol)
CreateSeuratFrom10x <- function(dir, project = basename(dir),
                                min.cells = 3, min.features = 200) {
  mtx   <- file.path(dir, "matrix.mtx.gz")
  feats <- file.path(dir, "features.tsv.gz")
  bc    <- file.path(dir, "barcodes.tsv.gz")
  if (!file.exists(mtx) || !file.exists(feats) || !file.exists(bc)) {
    stop("Missing 10x files in: ", dir)
  }
  raw <- ReadMtx(mtx = mtx, features = feats, cells = bc)  # dgCMatrix
  genes_orig <- rownames(raw)
  genes_new  <- dashify(genes_orig)
  rownames(raw) <- genes_new

  obj <- CreateSeuratObject(raw, project = project,
                            min.cells = min.cells, min.features = min.features)
  obj$condition <- project

  # Try to read features.tsv.gz to keep an ID↔symbol map (if present)
  feat_df <- tryCatch({
    fread(feats, header = FALSE, nThread = 1)
  }, error = function(e) NULL)

  if (!is.null(feat_df) && ncol(feat_df) >= 2) {
    colnames(feat_df)[1:2] <- c("feature_id", "symbol")
    obj@misc$feature_map <- feat_df
  }
  obj@misc$gene_symbol_map <- data.frame(original = genes_orig, seurat = genes_new, stringsAsFactors = FALSE)
  obj@misc$source_dir <- dir
  obj
}

# Main gating function
gate_ec <- function(o, min_markers = 3, normalize_first = TRUE) {
  stopifnot("RNA" %in% names(o@assays))
  DefaultAssay(o) <- "RNA"

  # Curated EC marker set (with common synonyms; compared case- & dash-insensitive)
  marker_syms <- c("PECAM1","CD31","CDH5","VE-CADHERIN","KDR","VEGFR2","FLT1",
                   "ESAM","VWF","PLVAP","ENG","CD34")
  norm_key <- function(x) gsub("-", "", toupper(x))
  marker_norm <- unique(norm_key(marker_syms))

  # If rownames look like IDs and a feature_map is present, map to symbols
  rn <- rownames(o)
  looks_id <- grepl("^(ENSG\\d+|AC\\d+|AL\\d+|AP\\d+|RP\\d+|LINC\\d+\\.?\\d*)$", rn)
  if (any(looks_id) && !is.null(o@misc$feature_map)) {
    fmap <- o@misc$feature_map
    id_to_sym <- setNames(fmap$symbol, fmap$feature_id)

    # Original 10x IDs in creation order (if kept)
    orig_ids <- if (!is.null(o@misc$gene_symbol_map)) o@misc$gene_symbol_map$original else rn
    new_syms <- id_to_sym[orig_ids]

    rn2 <- rn
    ok  <- !is.na(new_syms) & nzchar(new_syms)
    rn2[ok] <- dashify(new_syms[ok])

    if (any(duplicated(rn2))) {
      # collapse duplicates by summing counts and rebuild
      mat <- GetAssayData(o, slot = "counts")
      dup_groups <- split(seq_along(rn2), rn2)
      summed <- do.call(rbind, lapply(dup_groups, function(ix) Matrix::colSums(mat[ix, , drop=FALSE])))
      summed <- as(summed, "dgCMatrix")
      o <- CreateSeuratObject(summed, meta.data = o@meta.data, project = o@project.name)
    } else {
      mat <- GetAssayData(o, slot = "counts")
      rownames(mat) <- rn2
      o[["RNA"]]@counts <- mat
      o[["RNA"]]@data   <- new("dgCMatrix")  # recompute below if asked
    }
  }

  # Determine EC genes present
  rn_norm <- norm_key(rownames(o))
  present <- rn_norm %in% marker_norm
  if (!any(present)) stop("No endothelial markers found after mapping/norming.")
  ec_genes <- rownames(o)[present]

  if (normalize_first) {
    o <- NormalizeData(o, verbose = FALSE)
  }

  # Score: count how many EC markers are detected per cell (counts > 0)
  mat <- GetAssayData(o, slot = "counts")
  o$ECscore <- Matrix::colSums(mat[ec_genes, , drop=FALSE] > 0)

  subset(o, subset = ECscore >= min_markers)
}

# ---- Load or build objects ---------------------------------------------------

message("Loading Week6 reference...")
if (!file.exists(paths$ref6_rds)) stop("Missing ref6 RDS: ", paths$ref6_rds)
ref6 <- readRDS(paths$ref6_rds)

build_or_skip <- function(dir, label) {
  rds_out <- file.path(paths$out_dir_subsets, paste0(tolower(label), ".rds"))
  if (file.exists(rds_out)) {
    readRDS(rds_out)
  } else {
    obj <- CreateSeuratFrom10x(dir, project = label)
    saveRDS(obj, rds_out)
    obj
  }
}

message("Loading/creating 10x objects...")
ctrl <- build_or_skip(paths$tenx$Control, "Control")
emm1 <- build_or_skip(paths$tenx$EMM1,    "EMM1")
emm2 <- build_or_skip(paths$tenx$EMM2,    "EMM2")
mm   <- build_or_skip(paths$tenx$MM,      "MM")

# ---- Gate endothelial cells --------------------------------------------------

gate_and_save <- function(obj, name, min_markers = 3) {
  message(sprintf("Gating ECs: %s (min_markers = %d)", name, min_markers))
  ec <- gate_ec(obj, min_markers = min_markers, normalize_first = TRUE)
  out <- file.path(paths$out_dir_subsets, paste0(tolower(name), "_ec.rds"))
  saveRDS(ec, out)
  data.frame(
    object     = name,
    n_cells_in = ncol(obj),
    n_cells_ec = ncol(ec),
    frac_ec    = round(ncol(ec) / ncol(obj), 4),
    stringsAsFactors = FALSE
  )
}

summary_df <- rbind(
  gate_and_save(ctrl, "ctrl",  min_markers = 3),
  gate_and_save(emm1, "emm1",  min_markers = 3),
  gate_and_save(emm2, "emm2",  min_markers = 3),
  gate_and_save(mm,   "mm",    min_markers = 3),
  gate_and_save(ref6, "ref6",  min_markers = 3)
)

fwrite(summary_df, file.path(paths$out_dir_tables, "gating_summary.csv"))
message("Done. Wrote gated objects to data/subsets/ and summary to results/tables/gating_summary.csv")essage("Saved EC objects and session info.")
