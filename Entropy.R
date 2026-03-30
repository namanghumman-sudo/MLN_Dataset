###############################################################################
# FULL END-TO-END SCRIPT
# - Loads annotated Seurat object
# - Keeps EC clusters
# - (Optional) MAGIC imputation -> makes "magic" assay used for heatmap
# - Runs Slingshot on a stable embedding (MNN/PCA) and ALSO stores UMAP for plotting
# - Computes Entropy from slingCurveWeights (NOT left undefined)
# - Computes CellDivision from CellCycleScoring (or fallback MKI67/PCNA)
# - Makes Kevin-style equal-width binned pseudotime heatmap with
#   Waves + Cluster bar + Entropy strip + CellDivision strip
###############################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(slingshot)
  library(ComplexHeatmap)
  library(circlize)
  library(Matrix)
  library(grid)
  library(magicBatch)   # devtools::install_github("kbrulois/magicBatch")
})

set.seed(1234)

# =========================
# USER SETTINGS
# =========================
RDS_IN <- "annotated_mln_wt_res_0.4_noLEC.rds"  # your annotated object
keep_idents <- c("CapEC", "CapEC2", "PCV", "HEV", "Artery")

# Expression source assay/layer BEFORE MAGIC
ASSAY_INPUT <- "originalexp"
LAYER_INPUT <- "data"

# Use MAGIC?
DO_MAGIC <- TRUE

# Your working python for magicBatch
PYTHON_PATH <- "C:/Users/ghumm/miniconda3/envs/seuratextend/python.exe"

# MAGIC parameters
MAGIC_T <- 6

# Slingshot settings
START_CLUS <- "CapEC"
END_CLUS   <- c("Artery", "HEV")  # encourages branching

# Heatmap bins per cluster (equal-width panels)
NBINS_PER_CLUSTER <- 40

# =========================
# 0) Load + subset
# =========================
seu <- readRDS(RDS_IN)
seu_ec <- subset(seu, idents = keep_idents)
DefaultAssay(seu_ec) <- ASSAY_INPUT

message("Cells per cluster:")
print(table(Idents(seu_ec)))

# =========================
# 1) Genes for heatmap
# =========================
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
message("Heatmap genes used: ", length(genes_use))

# =========================
# 2) (Optional) MAGIC -> create "magic" assay used for heatmap
# =========================
ASSAY_USE <- ASSAY_INPUT
LAYER_USE <- LAYER_INPUT

if (DO_MAGIC) {
  stopifnot(file.exists(PYTHON_PATH))
  
  # expression genes x cells from input assay/layer
  expr_gc <- GetAssayData(seu_ec, assay = ASSAY_INPUT, layer = LAYER_INPUT)[genes_use, , drop = FALSE]
  expr_gc <- as.matrix(expr_gc)
  
  # MAGIC expects cells x genes
  expr_cells_genes <- t(expr_gc)
  
  # mar_mat_input: use MNN_corrected if available else pca
  rd_use <- if ("MNN_corrected" %in% names(seu_ec@reductions)) "MNN_corrected" else "pca"
  mar_mat_input <- Embeddings(seu_ec, reduction = rd_use)
  mar_mat_input <- mar_mat_input[, seq_len(min(30, ncol(mar_mat_input))), drop = FALSE]
  mar_mat_input <- mar_mat_input[rownames(expr_cells_genes), , drop = FALSE]
  
  message("Running MAGIC (t=", MAGIC_T, ") using embedding: ", rd_use)
  MAGIC_out <- magicBatch::magicBatch(
    data = expr_cells_genes,
    mar_mat_input = mar_mat_input,
    t_param = MAGIC_T,
    python_command = PYTHON_PATH
  )
  
  imputed_cells_genes <- MAGIC_out$imputed_data[[1]]
  stopifnot(all(dim(imputed_cells_genes) == dim(expr_cells_genes)))
  
  # back to genes x cells
  imputed_gc <- t(imputed_cells_genes)
  rownames(imputed_gc) <- colnames(expr_cells_genes)
  colnames(imputed_gc) <- rownames(expr_cells_genes)
  
  seu_ec[["magic"]] <- CreateAssayObject(data = imputed_gc)
  DefaultAssay(seu_ec) <- "magic"
  ASSAY_USE <- "magic"
  LAYER_USE <- "data"
  
  message("MAGIC assay created. Heatmap will use ASSAY_USE='magic'.")
} else {
  message("MAGIC disabled. Heatmap will use ASSAY_USE='", ASSAY_USE, "'.")
}

# =========================
# 3) CellDivision (CellCycleScoring or fallback)
# =========================
s.genes <- g2m.genes <- character(0)

if (exists("cc.genes.updated.2019", where = asNamespace("Seurat"), inherits = FALSE)) {
  cc <- get("cc.genes.updated.2019", envir = asNamespace("Seurat"))
  s.genes  <- intersect(cc$s.genes, rownames(seu_ec))
  g2m.genes <- intersect(cc$g2m.genes, rownames(seu_ec))
} else if (exists("cc.genes", where = asNamespace("Seurat"), inherits = FALSE)) {
  cc <- get("cc.genes", envir = asNamespace("Seurat"))
  s.genes  <- intersect(cc$s.genes, rownames(seu_ec))
  g2m.genes <- intersect(cc$g2m.genes, rownames(seu_ec))
}

if (length(s.genes) >= 10 && length(g2m.genes) >= 10) {
  seu_ec <- CellCycleScoring(seu_ec, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
  seu_ec$CellDivision_raw <- seu_ec$S.Score + seu_ec$G2M.Score
  message("CellDivision_raw computed from CellCycleScoring.")
} else {
  genes_prolif <- intersect(c("Mki67","MKI67","Pcna","PCNA"), rownames(seu_ec))
  if (length(genes_prolif) >= 1) {
    expr_pick <- GetAssayData(seu_ec, assay = ASSAY_USE, layer = LAYER_USE)
    pv <- colMeans(as.matrix(expr_pick[genes_prolif, , drop = FALSE]), na.rm = TRUE)
    seu_ec$CellDivision_raw <- as.numeric(scale(pv))
    message("CellDivision_raw fallback from: ", paste(genes_prolif, collapse = ", "))
  } else {
    seu_ec$CellDivision_raw <- 0
    warning("No cell-cycle genes and no MKI67/PCNA. CellDivision_raw set to 0.")
  }
}

# =========================
# 4) Slingshot (fit on MNN/PCA) + store UMAP for plotting
# =========================
sce <- as.SingleCellExperiment(seu_ec)
sce$cluster <- factor(as.character(Idents(seu_ec)), levels = keep_idents)

# Put fitting embedding into reducedDim(sce, "SLING")
if ("MNN_corrected" %in% names(seu_ec@reductions)) {
  reducedDim(sce, "SLING") <- Embeddings(seu_ec, "MNN_corrected")
  sling_embed_name <- "MNN_corrected"
} else if ("pca" %in% names(seu_ec@reductions)) {
  reducedDim(sce, "SLING") <- Embeddings(seu_ec, "pca")
  sling_embed_name <- "pca"
} else {
  stop("Need 'MNN_corrected' or 'pca' in seu_ec@reductions to run slingshot.")
}

# limit dims for stability
reducedDim(sce, "SLING") <- reducedDim(sce, "SLING")[, 1:min(30, ncol(reducedDim(sce, "SLING"))), drop = FALSE]

# Add UMAP coords for plotting (if present)
if ("umap.mnn" %in% names(seu_ec@reductions)) {
  reducedDim(sce, "UMAP") <- Embeddings(seu_ec, "umap.mnn")[colnames(seu_ec), , drop = FALSE]
} else if ("UMAP" %in% names(seu_ec@reductions)) {
  reducedDim(sce, "UMAP") <- Embeddings(seu_ec, "UMAP")[colnames(seu_ec), , drop = FALSE]
}

message("Running slingshot on: ", sling_embed_name, " (stored as reducedDim 'SLING')")
sce <- slingshot(
  sce,
  clusterLabels = "cluster",
  reducedDim    = "SLING",
  start.clus    = START_CLUS,
  end.clus      = END_CLUS
)

pt_mat <- as.matrix(slingPseudotime(sce))     # cells x lineages
w_mat  <- as.matrix(slingCurveWeights(sce))   # cells x lineages

message("Lineages: ", ncol(pt_mat))
print(slingLineages(sce))
print(colSums(!is.na(pt_mat)))
print(table(rowSums(w_mat > 0.05, na.rm = TRUE)))

# Choose lineage with most assigned cells for ordering
n_assigned <- colSums(!is.na(pt_mat))
best_idx <- which.max(n_assigned)

colnames(pt_mat) <- paste0("slingPseudotime_", seq_len(ncol(pt_mat)))
LINEAGE_COL <- paste0("slingPseudotime_", best_idx)

# Store pseudotime into Seurat meta.data
pt_mat <- pt_mat[colnames(seu_ec), , drop = FALSE]
seu_ec@meta.data[, colnames(pt_mat)] <- pt_mat

pt_all <- pt_mat[, LINEAGE_COL]
stopifnot(sum(!is.na(pt_all)) > 50)

# =========================
# 5) Entropy from curve weights (THIS is what was missing/misaligned before)
# =========================
w <- w_mat[colnames(seu_ec), , drop = FALSE]
w[is.na(w)] <- 0

rs <- rowSums(w)
rs[rs == 0] <- NA
p <- w / rs

entropy_cell <- -rowSums(ifelse(p > 0, p * log(p), 0), na.rm = TRUE)
seu_ec$Entropy <- entropy_cell

# (optional) normalized entropy 0..1
seu_ec$Entropy_norm <- entropy_cell / log(ncol(p))

message("Entropy summary:")
print(summary(seu_ec$Entropy))
message("Entropy in cells with >=2 lineages weight>0.05:")
both <- rowSums(p > 0.05, na.rm=TRUE) >= 2
print(summary(seu_ec$Entropy[both]))

# Save weights too (optional)
seu_ec@misc$slingshot <- list(SlingPseudotime = pt_mat, SlingCurveWeights = w_mat)

# Z-score/clamp helper
z_clamp_cell <- function(x) {
  z <- as.numeric(scale(x))
  z[is.na(z)] <- 0
  z[z > 2] <- 2
  z[z < -2] <- -2
  z
}

seu_ec$Entropy_z      <- z_clamp_cell(seu_ec$Entropy)
seu_ec$CellDivision   <- seu_ec$CellDivision_raw
seu_ec$CellDivision_z <- z_clamp_cell(seu_ec$CellDivision)

message("CellDivision_raw range: ", paste(range(seu_ec$CellDivision_raw, na.rm=TRUE), collapse=" to "))

# =========================
# 6) Plot Slingshot (optional)
# =========================
if ("UMAP" %in% reducedDimNames(sce)) {
  um <- reducedDim(sce, "UMAP")
  cl <- sce$cluster
  
  plot(um[,1], um[,2],
       col = as.numeric(factor(cl)), pch = 16, cex = 0.5,
       xlab = "UMAP_1", ylab = "UMAP_2",
       main = "UMAP with Slingshot curves")
  lines(SlingshotDataSet(sce), lwd = 2, col = "black")
  
  centroids <- aggregate(um, by = list(cluster = cl), FUN = median)
  text(centroids[,2], centroids[,3], labels = centroids$cluster, cex = 0.9, font = 2)
} else {
  X <- reducedDim(sce, "SLING")
  cl <- sce$cluster
  plot(X[,1], X[,2],
       col = as.numeric(factor(cl)), pch = 16, cex = 0.5,
       xlab = "SLING_1", ylab = "SLING_2",
       main = "Slingshot curves (fit space)")
  lines(SlingshotDataSet(sce), lwd = 2, col = "black")
}

# =========================
# 7) Order cells by pseudotime + bin per cluster (equal-width)
# =========================
cells_keep <- names(pt_all)[!is.na(pt_all)]
pt_keep    <- pt_all[cells_keep]
seg_keep   <- factor(as.character(Idents(seu_ec)[cells_keep]), levels = keep_idents)

ord_cells <- order(pt_keep)
cells_ord <- cells_keep[ord_cells]
seg_ord   <- seg_keep[ord_cells]

entropy_ord <- seu_ec$Entropy[cells_ord]
div_ord     <- seu_ec$CellDivision_raw[cells_ord]

# expression genes x ordered cells from chosen assay/layer
mat <- GetAssayData(seu_ec, assay = ASSAY_USE, layer = LAYER_USE)[genes_use, cells_ord, drop = FALSE]
mat <- as.matrix(mat)

mat_bin_list <- list()
bin_cluster  <- character(0)
bin_entropy  <- numeric(0)
bin_div      <- numeric(0)

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

# =========================
# 8) Z-score + smooth matrix, then order genes diagonally
# =========================
mat_z <- t(scale(t(mat_bin)))
mat_z[is.na(mat_z)] <- 0
mat_z[mat_z > 2] <- 2
mat_z[mat_z < -2] <- -2

mat_z <- t(apply(mat_z, 1, function(x) smooth.spline(x = seq_along(x), y = x, spar = 0.8)$y))
mat_z[is.na(mat_z)] <- 0
mat_z[mat_z > 2] <- 2
mat_z[mat_z < -2] <- -2

get_peak_within <- function(g) {
  peaks <- numeric(length(keep_idents))
  for (i in seq_along(keep_idents)) {
    cl <- keep_idents[i]
    cols <- which(bin_cluster == cl)
    if (length(cols) == 0) { peaks[i] <- NA; next }
    v <- mat_z[g, cols]
    peaks[i] <- which.max(v) / length(cols)
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
mat_z_ord <- mat_z[rownames(mat_z)[gene_order], , drop = FALSE]

# =========================
# 9) Entropy + CellDivision tracks (z-score across bins)
# =========================
z_clamp_bins <- function(x) {
  z <- as.numeric(scale(x))
  z[is.na(z)] <- 0
  z[z > 2] <- 2
  z[z < -2] <- -2
  z
}
entropy_z <- z_clamp_bins(bin_entropy)
div_z     <- z_clamp_bins(bin_div)

message("Entropy bin range: ", paste(range(bin_entropy, na.rm=TRUE), collapse=" to "))
message("Div bin range: ", paste(range(bin_div, na.rm=TRUE), collapse=" to "))

# =========================
# 10) Heatmap annotations + draw
# =========================
col_fun <- colorRamp2(
  c(-2, -0.5, 0, 1, 2),
  c("black", "#3b007a", "#7a00cc", "#ffae00", "#ffffcc")
)

seg_cols <- setNames(
  c("#4FD1C5", "#63B3ED", "#B794F4", "#F6AD55", "#FBB6CE"),
  keep_idents
)

# Waves
seg_levels <- levels(bin_cluster)
S <- sapply(seg_levels, function(s) as.numeric(as.character(bin_cluster) == s))
x <- seq_len(ncol(mat_z_ord))
SPAR_WAVE <- 0.99
S_smooth <- sapply(seq_along(seg_levels), function(j) {
  y <- S[, j]
  fit <- smooth.spline(x = x, y = y, spar = SPAR_WAVE)
  yhat <- predict(fit, x)$y
  pmin(1, pmax(0, yhat))
})
colnames(S_smooth) <- seg_levels
mx <- apply(S_smooth, 2, max); mx[mx == 0] <- 1
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
    unit(14, "mm"),
    unit(4, "mm"),
    unit(4, "mm"),
    unit(4, "mm")
  )
)

ht <- Heatmap(
  mat_z_ord,
  name = "z",
  top_annotation = top_ha,
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_split = bin_cluster,
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 9),
  border = FALSE
)

draw(ht)

grid::grid.text(
  paste0("MLN WT EC trajectory heatmap (", LINEAGE_COL, ", fit=", sling_embed_name,
         ", MAGIC=", DO_MAGIC, ")"),
  x = unit(0.5, "npc"),
  y = unit(2.0, "npc"),
  gp = gpar(fontsize = 12, fontface = "bold")
)

# =========================
# 11) Save object for reuse (optional)
# =========================
saveRDS(seu_ec, file = "seu_ec_with_magic_slingshot_entropy_division.rds")
message("Saved: seu_ec_with_magic_slingshot_entropy_division.rds")