# CardiacOrganoid_scRNAseq
Independent single-cell RNA-seq analysis of human heart organoids (Aguirre lab GSE218582 + reference GSE106118): QC → integration → clustering → markers → pathways → vascularization hypothesis.
Quality control was applied to both the reference dataset (Week-6 fetal heart cells, GSE106118) and the experimental cardiac organoid datasets (Control, EMM1, EMM2, MM; GSE218582).

For the Week-6 reference, cells were filtered using thresholds on library size and dropout rates (libsizes > 2000 and meanZero < 0.98), ensuring removal of low-quality or sparsely expressed cells.

For the experimental organoid datasets, Seurat’s default QC thresholds were applied (min.features = 200, min.cells = 3) to exclude low-quality cells and genes with insufficient coverage. In addition, endothelial identity was enforced using a biological gating step, requiring expression of at least three canonical endothelial markers (PECAM1, CDH5, VWF, KDR, etc.). This strategy ensured that downstream analyses focused specifically on bona fide endothelial cells while minimizing technical noise.
