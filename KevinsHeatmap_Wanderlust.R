suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(slingshot)
  library(ComplexHeatmap)
  library(circlize)
  library(Matrix)
  library(grid)
  library(magicBatch)   # devtools::install_github("kbrulois/magicBatch")
  library(reticulate)   # optional but helpful for finding the python path
})

# ============================================================
# GOAL
# - Top: "Waves" + Cluster blocks (split by populations, titles)
# - Body: smooth pseudotime heatmap (BINNED)
# - Genes ordered to look like a diagonal (early peaks at top, late at bottom)
# - Entropy + CellDivision shown as HEATMAP STRIPS (same z-score/color scale as body)
#
# THIS VERSION ORDERS CELLS BY: Wanderlust pseudotime (imported from CSV)
# Optional: still run Slingshot ONLY to get weights -> entropy (can be skipped)
# ============================================================

# -----------------------------
# USER SETTINGS
# -----------------------------
SEURAT_RDS <- "annotated_mln_wt_res_0.4_noLEC.rds"

keep_idents <- c("CapEC", "CapEC2", "PCV", "HEV", "Artery")

# Expression for heatmap will come from this assay/layer (we overwrite to MAGIC later)
ASSAY_USE <- "originalexp"
LAYER_USE <- "data"

# Wanderlust pseudotime file exported from MATLAB
WANDERLUST_CSV <- "wanderlust_pseudotime.csv"   # columns: cell,pseudotime

# MAGIC python env (Windows)
CONDA_ENV <- "C:/Users/ghumm/miniconda3/envs/seuratextend"
python_path <- file.path(CONDA_ENV, "python.exe")

# Optional: compute entropy from Slingshot weights (TRUE/FALSE)
COMPUTE_ENTROPY_FROM_SLINGSHOT <- TRUE

# Heatmap binning
NBINS_PER_CLUSTER <- 40

# Smoothing controls
SPAR_GENE_SMOOTH <- 0.8
SPAR_WAVE <- 0.99

# -----------------------------
# Helper: z-score + clamp
# -----------------------------
z_clamp <- function(x, clamp = 2) {
  z <- as.numeric(scale(x))
  z[is.na(z)] <- 0
  z[z >  clamp] <-  clamp
  z[z < -clamp] <- -clamp
  z
}

# Shannon entropy from lineage weights (per cell)
entropy_from_weights <- function(w_mat, eps = 1e-12) {
  # w_mat: cells x lineages
  w <- w_mat
  w[is.na(w)] <- 0
  rs <- rowSums(w)
  rs[rs == 0] <- 1
  p <- w / rs
  # H = -sum(p log p)
  H <- -rowSums(p * log(p + eps))
  as.numeric(H)
}

# -----------------------------
# 0) Load & keep only labeled EC clusters
# -----------------------------
seu <- readRDS(SEURAT_RDS)

seu_ec <- subset(seu, idents = keep_idents)
stopifnot(all(keep_idents %in% levels(Idents(seu_ec)) | table(Idents(seu_ec)) == 0))
print(table(Idents(seu_ec)))

# -----------------------------
# 1) Assay/layer
# -----------------------------
DefaultAssay(seu_ec) <- ASSAY_USE

# -----------------------------
# 2) Genes
# -----------------------------
hev_markers    <- c("St6gal1","Glycam1","Chst4","Ccl21a","Madcam1","Fut7")
art_markers    <- c("Bmx","Depp1","Gja4","Gja5","Gkn3","Sox17")
CRP_markers    <- c("Angpt2","Apln","Esm1","Mcam","Nid2","Pdgfb","Pgf","Vim")
PCV_markers    <- c("Ackr1","Icam1","Nr2f2","Sele","Selp","Vcam1","Vwf")
CapEC_marker   <- c("Cdh13","Emcn","Gja1","Gpihbp1","Ly6a","Ly6c1","Podxl","Ramp3")
CapEC2_marker  <- c("Atf3","Cxcl1","Egr1","Fos","Fosb","Jun","Junb","Jund","Nfkbia","Nr4a1","Sgk1")
CapEC1_marker  <- c("Col4a1","Col4a2","Hlx","Id1","Id3","Igfbp3")
art_vn_markers <- c("Edn1","Eln","Fbln5","Foxn3","Klf4","Ltbp4","Ptprb")

genes_use <- unique(c(
  hev_markers, art_markers, CRP_markers, PCV_markers,
  CapEC_marker, CapEC1_marker, CapEC2_marker, art_vn_markers
))
genes_use <- intersect(genes_use, rownames(seu_ec))
stopifnot(length(genes_use) >= 10)
message("Genes used: ", length(genes_use))

# -----------------------------
# 2.25) MAGIC imputation (magicBatch) -> create a "magic" assay
# -----------------------------
stopifnot(dir.exists(CONDA_ENV))
stopifnot(file.exists(python_path))

# Make conda env DLLs discoverable on Windows
new_path <- paste(
  file.path(CONDA_ENV, "Library", "bin"),
  file.path(CONDA_ENV, "Library", "usr", "bin"),
  file.path(CONDA_ENV, "Scripts"),
  CONDA_ENV,
  Sys.getenv("PATH"),
  sep = ";"
)
Sys.setenv(PATH = new_path)

dir.create("C:/temp_magic", showWarnings = FALSE)
Sys.setenv(TMPDIR="C:/temp_magic", TEMP="C:/temp_magic", TMP="C:/temp_magic")

# Sanity check python modules
system2(python_path, c("-c", shQuote("import numpy,pandas,scipy,sklearn,magic; print('OK from R')")),
        stdout = TRUE, stderr = TRUE)

genes_magic <- genes_use

expr_gc <- GetAssayData(seu_ec, assay = ASSAY_USE, layer = LAYER_USE)[genes_magic, , drop = FALSE]
expr_gc <- as.matrix(expr_gc)

# magicBatch expects CELLS x GENES
expr_cells_genes <- t(expr_gc)

# MAR matrix input: use MNN_corrected embedding if present
rd_use <- if ("MNN_corrected" %in% names(seu_ec@reductions)) "MNN_corrected" else "pca"
mar_mat_input <- Embeddings(seu_ec, reduction = rd_use)
mar_mat_input <- mar_mat_input[, seq_len(min(30, ncol(mar_mat_input))), drop = FALSE]
mar_mat_input <- mar_mat_input[rownames(expr_cells_genes), , drop = FALSE]

MAGIC_out <- magicBatch::magicBatch(
  data = expr_cells_genes,
  mar_mat_input = mar_mat_input,
  t_param = 6,
  python_command = python_path
)

imputed_cells_genes <- MAGIC_out$imputed_data[[1]]
stopifnot(all(dim(imputed_cells_genes) == dim(expr_cells_genes)))

# Convert back to GENES x CELLS
imputed_gc <- t(imputed_cells_genes)
rownames(imputed_gc) <- colnames(expr_cells_genes)
colnames(imputed_gc) <- rownames(expr_cells_genes)

seu_ec[["magic"]] <- CreateAssayObject(data = imputed_gc)

# Switch downstream code to use MAGIC for the heatmap expression
ASSAY_USE <- "magic"
LAYER_USE <- "data"
DefaultAssay(seu_ec) <- ASSAY_USE

message("MAGIC assay created: using ASSAY_USE='magic' for the heatmap expression.")

# -----------------------------
# 2.5) Cell division / cell-cycle scoring (robust)
# -----------------------------
assay_before <- DefaultAssay(seu_ec)
DefaultAssay(seu_ec) <- "originalexp"   # score on non-MAGIC (change if your raw assay differs)

data("cc.genes.updated.2019", package = "Seurat")
cc <- cc.genes.updated.2019

toMouseStyle <- function(x) paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
s0  <- toMouseStyle(cc$s.genes)
g0  <- toMouseStyle(cc$g2m.genes)

m.s.genes   <- Seurat::CaseMatch(search = s0, match = rownames(seu_ec))
m.g2m.genes <- Seurat::CaseMatch(search = g0, match = rownames(seu_ec))

cat("Matched mouse S genes:", length(m.s.genes), "\n")
cat("Matched mouse G2M genes:", length(m.g2m.genes), "\n")

seu_ec <- Seurat::CellCycleScoring(
  seu_ec,
  s.features = m.s.genes,
  g2m.features = m.g2m.genes,
  set.ident = FALSE
)
seu_ec$CellDivision_raw <- seu_ec$S.Score + seu_ec$G2M.Score

DefaultAssay(seu_ec) <- assay_before

# -----------------------------
# 3) PSEUDOTIME: load Wanderlust pseudotime (REPLACES Slingshot ordering)
# -----------------------------
stopifnot(file.exists(WANDERLUST_CSV))
wl <- read.csv(WANDERLUST_CSV, stringsAsFactors = FALSE)
stopifnot(all(c("cell","pseudotime") %in% colnames(wl)))

wanderlust_pt <- setNames(wl$pseudotime, wl$cell)

# Align to this Seurat object's cell names
pt_all <- wanderlust_pt[colnames(seu_ec)]
LINEAGE_COL <- "Wanderlust"

# sanity
message("Wanderlust pseudotime matched cells: ", sum(!is.na(pt_all)), " / ", length(pt_all))
stopifnot(sum(!is.na(pt_all)) > 50)

# -----------------------------
# 3b) OPTIONAL: run Slingshot ONLY to get weights -> entropy strip
#      (heatmap ordering still uses Wanderlust!)
# -----------------------------
entropy_cell <- rep(NA_real_, ncol(seu_ec))
names(entropy_cell) <- colnames(seu_ec)
w_mat <- NULL

if (COMPUTE_ENTROPY_FROM_SLINGSHOT) {
  
  rd_base <- if ("MNN_corrected" %in% names(seu_ec@reductions)) "MNN_corrected" else
    if ("umap.mnn" %in% names(seu_ec@reductions)) "umap.mnn" else
      if ("pca" %in% names(seu_ec@reductions)) "pca" else NULL
  stopifnot(!is.null(rd_base))
  
  X <- Embeddings(seu_ec, reduction = rd_base)
  X <- X[, seq_len(min(10, ncol(X))), drop = FALSE]
  
  sce <- as.SingleCellExperiment(seu_ec)
  sce$cluster <- as.character(Idents(seu_ec))
  reducedDim(sce, "SLING") <- X
  
  start_clus <- "CapEC"
  if (!(start_clus %in% sce$cluster)) start_clus <- unique(sce$cluster)[1]
  
  sce <- slingshot(sce, clusterLabels = "cluster", reducedDim = "SLING", start.clus = start_clus)
  
  w_mat <- as.matrix(slingCurveWeights(sce))  # cells x lineages
  rownames(w_mat) <- colnames(seu_ec)
  
  entropy_cell <- entropy_from_weights(w_mat)
  names(entropy_cell) <- colnames(seu_ec)
  
  message("Computed entropy from Slingshot curve weights (ordering still Wanderlust).")
} else {
  message("Skipping Slingshot entropy; Entropy strip will be flat/zero unless you provide entropy_cell.")
  entropy_cell[is.na(entropy_cell)] <- 0
}

# Store Entropy / CellDivision in meta.data
seu_ec$Entropy <- entropy_cell[colnames(seu_ec)]
seu_ec$Entropy_z <- z_clamp(seu_ec$Entropy)
seu_ec$CellDivision <- seu_ec$CellDivision_raw
seu_ec$CellDivision_z <- z_clamp(seu_ec$CellDivision)

# -----------------------------
# 4) Keep only pseudotime cells and order them by Wanderlust
# -----------------------------
cells_keep <- names(pt_all)[!is.na(pt_all)]
pt_keep    <- pt_all[cells_keep]
seg_keep   <- factor(as.character(Idents(seu_ec)[cells_keep]), levels = keep_idents)

ord_cells <- order(pt_keep)
cells_ord <- cells_keep[ord_cells]
pt_ord    <- pt_keep[ord_cells]
seg_ord   <- seg_keep[ord_cells]

entropy_ord <- seu_ec$Entropy[cells_ord]
div_ord     <- seu_ec$CellDivision[cells_ord]

# -----------------------------
# 5) Bin per population (EQUAL WIDTH PANELS)
# -----------------------------
mat <- GetAssayData(seu_ec, assay = ASSAY_USE, layer = LAYER_USE)[genes_use, cells_ord, drop = FALSE]
mat <- as.matrix(mat)

mat_bin_list <- list()
bin_cluster <- character(0)
bin_entropy <- numeric(0)
bin_div     <- numeric(0)

for (cl in keep_idents) {
  idx <- which(seg_ord == cl)
  if (length(idx) < 10) next
  
  cuts_cl <- cut(seq_along(idx), breaks = NBINS_PER_CLUSTER, labels = FALSE)
  
  cl_bin <- sapply(seq_len(NBINS_PER_CLUSTER), function(b) {
    rowMeans(mat[, idx[cuts_cl == b], drop = FALSE])
  })
  colnames(cl_bin) <- paste0(cl, "_B", seq_len(ncol(cl_bin)))
  mat_bin_list[[cl]] <- cl_bin
  
  bin_cluster <- c(bin_cluster, rep(cl, NBINS_PER_CLUSTER))
  
  bin_entropy <- c(bin_entropy, sapply(seq_len(NBINS_PER_CLUSTER), function(b) {
    median(entropy_ord[idx[cuts_cl == b]], na.rm = TRUE)
  }))
  bin_div <- c(bin_div, sapply(seq_len(NBINS_PER_CLUSTER), function(b) {
    median(div_ord[idx[cuts_cl == b]], na.rm = TRUE)
  }))
}

mat_bin <- do.call(cbind, mat_bin_list)
bin_cluster <- factor(bin_cluster, levels = keep_idents)

stopifnot(ncol(mat_bin) == length(bin_cluster))
stopifnot(length(bin_entropy) == ncol(mat_bin))
stopifnot(length(bin_div) == ncol(mat_bin))

# -----------------------------
# 6) Z-score the binned gene matrix + smooth across bins
# -----------------------------
mat_z <- t(scale(t(mat_bin)))
mat_z[is.na(mat_z)] <- 0
mat_z[mat_z > 2] <- 2
mat_z[mat_z < -2] <- -2

mat_z <- t(apply(mat_z, 1, function(x) smooth.spline(x = seq_along(x), y = x, spar = SPAR_GENE_SMOOTH)$y))
mat_z[is.na(mat_z)] <- 0
mat_z[mat_z > 2] <- 2
mat_z[mat_z < -2] <- -2

# -----------------------------
# 7) Gene ordering for a diagonal across POPULATIONS
# -----------------------------
get_peak_within <- function(g) {
  peaks <- numeric(length(keep_idents))
  for (i in seq_along(keep_idents)) {
    cl <- keep_idents[i]
    cols <- which(bin_cluster == cl)
    if (length(cols) == 0) { peaks[i] <- NA; next }
    v <- mat_z[g, cols]
    peaks[i] <- which.max(v) / length(cols)  # scaled 0-1 within that panel
  }
  peaks
}

peak_mat <- t(sapply(rownames(mat_z), get_peak_within))
colnames(peak_mat) <- keep_idents

panel_mean <- sapply(keep_idents, function(cl) {
  cols <- which(bin_cluster == cl)
  rowMeans(mat_z[, cols, drop = FALSE])
})
colnames(panel_mean) <- keep_idents

best_panel <- apply(panel_mean, 1, function(x) keep_idents[which.max(x)])
panel_index <- match(best_panel, keep_idents)

within_peak <- peak_mat[cbind(rownames(peak_mat), best_panel)]
within_peak[is.na(within_peak)] <- 0.5

global_rank <- panel_index + within_peak
best_strength <- apply(panel_mean, 1, max)

gene_order <- order(global_rank, -best_strength)

genes_ord <- rownames(mat_z)[gene_order]
mat_z_ord <- mat_z[genes_ord, , drop = FALSE]

# -----------------------------
# 8) Entropy + CellDivision strips use SAME scoring as body (z across bins)
# -----------------------------
entropy_z <- z_clamp(bin_entropy)
div_z     <- z_clamp(bin_div)

# -----------------------------
# 9) Colors + top annotation (Waves + Cluster + Entropy strip + CellDivision strip)
# -----------------------------
col_fun <- colorRamp2(
  c(-2, -0.5, 0, 1, 2),
  c("black", "#3b007a", "#7a00cc", "#ffae00", "#ffffcc")
)

seg_cols <- setNames(
  c("#4FD1C5", "#63B3ED", "#B794F4", "#F6AD55", "#FBB6CE"),
  keep_idents
)

# Waves: smooth membership of each population across bins
seg_levels <- levels(bin_cluster)
S <- sapply(seg_levels, function(s) as.numeric(as.character(bin_cluster) == s))
x <- seq_len(ncol(mat_z_ord))

S_smooth <- sapply(seq_along(seg_levels), function(j) {
  y <- S[, j]
  fit <- smooth.spline(x = x, y = y, spar = SPAR_WAVE)
  yhat <- predict(fit, x)$y
  pmin(1, pmax(0, yhat))
})
colnames(S_smooth) <- seg_levels

mx <- apply(S_smooth, 2, max)
mx[mx == 0] <- 1
S_smooth <- sweep(S_smooth, 2, mx, "/")

wave_gp <- grid::gpar(col = unname(seg_cols[seg_levels]), lwd = 2)

top_ha <- HeatmapAnnotation(
  Waves = anno_lines(S_smooth, gp = wave_gp, ylim = c(0, 1), axis = FALSE),
  Cluster = bin_cluster,
  Entropy = anno_simple(entropy_z, col = col_fun),
  CellDivision = anno_simple(div_z, col = col_fun),
  col = list(Cluster = seg_cols),
  annotation_name_side = "left",
  annotation_height = unit.c(
    unit(14, "mm"),  # waves
    unit(4, "mm"),   # cluster bar
    unit(4, "mm"),   # entropy strip
    unit(4, "mm")    # cell division strip
  )
)

# -----------------------------
# 10) Heatmap
# -----------------------------
ht <- ComplexHeatmap::Heatmap(
  mat_z_ord,
  name = "z",
  top_annotation = top_ha,
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_split = bin_cluster,
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_names_gp = grid::gpar(fontsize = 9),
  border = FALSE
)

ComplexHeatmap::draw(ht)

grid::grid.text(
  paste0("MLN WT EC trajectory heatmap (", LINEAGE_COL, ")"),
  x = unit(0.5, "npc"),
  y = unit(2.0, "npc"),
  gp = grid::gpar(fontsize = 12, fontface = "bold")
)

# -----------------------------
# 11) Debug messages
# -----------------------------
message("Pseudotime range (Wanderlust, kept cells): ",
        paste0(signif(range(pt_keep, na.rm = TRUE), 4), collapse = " to "))

message("Entropy range (cell-level, non-NA): ",
        paste0(signif(range(entropy_cell[!is.na(entropy_cell)]), 4), collapse = " to "))

message("Entropy range (bin-level): ",
        paste0(signif(range(bin_entropy, na.rm = TRUE), 4), collapse = " to "))

message("CellDivision range (raw): ",
        paste0(signif(range(seu_ec$CellDivision_raw, na.rm = TRUE), 4), collapse = " to "))

if (!is.null(w_mat)) {
  print(table(rowSums(w_mat > 0.05, na.rm = TRUE)))
}

system2(python_path, "--version", stdout = TRUE, stderr = TRUE)