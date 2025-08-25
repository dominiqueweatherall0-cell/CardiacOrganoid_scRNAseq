# scripts/02_create_objects_10x.R
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(data.table) })
source("scripts/utils_seurat_helpers.R")

# ---- set relative folders for the four 10x datasets ----
CTRL_DIR <- "data/raw/Control"
EMM1_DIR <- "data/raw/EMM1"
EMM2_DIR <- "data/raw/EMM2"
MM_DIR   <- "data/raw/MM"

stopifnot(dir.exists(CTRL_DIR), dir.exists(EMM1_DIR), dir.exists(EMM2_DIR), dir.exists(MM_DIR))

message("Creating Seurat objects from 10x matrices…")
ctrl <- CreateSeuratFrom10x(CTRL_DIR, "Control"); message("  Control ✓")
emm1 <- CreateSeuratFrom10x(EMM1_DIR, "EMM1");     message("  EMM1 ✓")
emm2 <- CreateSeuratFrom10x(EMM2_DIR, "EMM2");     message("  EMM2 ✓")
mm   <- CreateSeuratFrom10x(MM_DIR,   "MM");       message("  MM ✓")

dir.create("data/subsets", recursive = TRUE, showWarnings = FALSE)
saveRDS(ctrl, "data/subsets/ctrl_raw.rds")
saveRDS(emm1,"data/subsets/emm1_raw.rds")
saveRDS(emm2,"data/subsets/emm2_raw.rds")
saveRDS(mm,  "data/subsets/mm_raw.rds")

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
writeLines(c(capture.output(sessionInfo())), "results/tables/sessionInfo_create_objects.txt")

message("Saved raw objects to data/subsets/")
