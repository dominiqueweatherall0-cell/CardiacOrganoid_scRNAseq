# scripts/03_gate_endothelial.R
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(data.table) })
# load inputs produced by script 02
ctrl <- readRDS("data/subsets/ctrl_raw.rds")
emm1 <- readRDS("data/subsets/emm1_raw.rds")
emm2 <- readRDS("data/subsets/emm2_raw.rds")
mm   <- readRDS("data/subsets/mm_raw.rds")

dashify <- function(x){ x <- gsub("_","-", x, fixed=TRUE); x <- gsub("\\s+","", x); make.unique(x) }

gate_ec <- function(o, collapse_duplicates = FALSE, ec_min_detected = 3){
  DefaultAssay(o) <- "RNA"

  # Endothelial panel
  marker_syms <- c("PECAM1","CD31","CDH5","VE-CADHERIN","KDR","VEGFR2","FLT1","ESAM","VWF","PLVAP","ENG","CD34")
  normalize_syms <- function(x) gsub("-", "", toupper(x))
  marker_syms <- unique(normalize_syms(marker_syms))

  # Map Ensembl-like IDs to symbols if feature_map is available
  rn <- rownames(o)
  looks_id <- grepl("^(ENSG\\d+|AC\\d+|AL\\d+|AP\\d+|RP\\d+|LINC\\d+\\.?\\d*)$", rn)
  if (any(looks_id) && !is.null(o@misc$feature_map)) {
    fmap <- o@misc$feature_map
    id_to_sym <- setNames(fmap$symbol, fmap$feature_id)
    rn2 <- rn
    new_syms <- id_to_sym[o@misc$gene_symbol_map$original]
    ok <- !is.na(new_syms) & nzchar(new_syms)
    rn2[ok] <- dashify(new_syms[ok])

    if (collapse_duplicates && any(duplicated(rn2))) {
      mat <- GetAssayData(o, assay="RNA", slot="counts")
      dup_groups <- split(seq_along(rn2), rn2)
      summed <- do.call(rbind, lapply(dup_groups, function(ix) Matrix::colSums(mat[ix,,drop=FALSE])))
      o <- CreateSeuratObject(as(summed,"dgCMatrix"), meta.data=o@meta.data, project=o@project.name)
    } else if (!all(rn == rn2)) {
      mat <- GetAssayData(o, assay="RNA", slot="counts")
      rownames(mat) <- rn2
      o[["RNA"]]@counts <- mat
      o[["RNA"]]@data   <- new("dgCMatrix")  # recomputed after NormalizeData
    }
  }

  o <- NormalizeData(o, verbose = FALSE)

  present <- normalize_syms(rownames(o)) %in% marker_syms
  if (!any(present)) stop("No endothelial markers found after mapping. Inspect rownames(o)[1:20].")
  ec_genes <- rownames(o)[present]

  mat <- GetAssayData(o, assay="RNA", slot="counts")
  o$ECscore <- Matrix::colSums(mat[ec_genes, , drop = FALSE] > 0)

  subset(o, subset = ECscore >= ec_min_detected)
}

message("Gating EC populations…")
ctrl_ec <- gate_ec(ctrl);  message("  Control ✓")
emm1_ec <- gate_ec(emm1);  message("  EMM1 ✓")
emm2_ec <- gate_ec(emm2);  message("  EMM2 ✓")
mm_ec   <- gate_ec(mm);    message("  MM ✓")

dir.create("data/subsets", recursive = TRUE, showWarnings = FALSE)
saveRDS(ctrl_ec, "data/subsets/ctrl_ec.rds")
saveRDS(emm1_ec,"data/subsets/emm1_ec.rds")
saveRDS(emm2_ec,"data/subsets/emm2_ec.rds")
saveRDS(mm_ec,  "data/subsets/mm_ec.rds")

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
writeLines(c(capture.output(sessionInfo())), "results/tables/sessionInfo_scRNA_gate.txt")

message("Saved EC objects and session info.")
