################################################################################
## MATLAB tSpace integration end-to-end (Windows)
################################################################################

################################################################################
## 0) Libraries
################################################################################
library(Seurat)
library(ggplot2)
library(dplyr)
library(SingleCellExperiment)
library(ComplexHeatmap)
library(Matrix)
library(batchelor)
library(DelayedMatrixStats)
library(circlize)
library(viridis)

set.seed(1234)

################################################################################
## 1) Load object
################################################################################
clean_obj <- readRDS("C:/Users/ghumm/Downloads/mlnwt_fastmnn.rds")

seu <- as.Seurat(
  x = clean_obj,
  counts = "counts",
  data   = "logcounts"
)
DefaultAssay(seu) <- "originalexp"

################################################################################
## 2) Remove contaminating populations (example: LEC filter)
################################################################################
lec         <- c("Prox1","Lyve1","Pdpn")
lymphocytes <- c("Ptprc","Cd52")
fib_ret     <- c("Pdpn","Ccl19","Pdgfra")

modules <- list(
  LEC        = lec,
  Lymphocyte = lymphocytes,
  FRC        = fib_ret
)

seu <- AddModuleScore(seu, features = modules, name = names(modules))
colnames(seu@meta.data)[grep("^LEC", colnames(seu@meta.data))] <- "LEC_Score"

LEC_cutoff <- quantile(seu$LEC_Score, 0.95, na.rm = TRUE)
seu <- subset(seu, subset = LEC_Score < LEC_cutoff)

################################################################################
## 3) Highly variable genes + remove CC genes
################################################################################
DefaultAssay(seu) <- "originalexp"   # make sure you're scoring the right assay

data("cc.genes.updated.2019", package="Seurat")

s.genes   <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

# Keep only genes present in your object
s.use   <- intersect(s.genes, rownames(seu))
g2m.use <- intersect(g2m.genes, rownames(seu))

length(s.use)
length(g2m.use)

# If gene case is different (e.g. all lowercase), try case-matching:
if (length(s.use) < 10 || length(g2m.use) < 10) {
  rn <- rownames(seu)
  
  s.use   <- rn[match(toupper(s.genes),   toupper(rn), nomatch = 0)]
  g2m.use <- rn[match(toupper(g2m.genes), toupper(rn), nomatch = 0)]
  
  s.use   <- unique(s.use[s.use != ""])
  g2m.use <- unique(g2m.use[g2m.use != ""])
}

if (length(s.use) < 10 || length(g2m.use) < 10) {
  stop("Still not enough CC genes found. Your gene IDs may be Ensembl IDs or a different naming scheme.")
}

seu <- CellCycleScoring(
  seu,
  s.features = s.use,
  g2m.features = g2m.use
)

seu$Dividing <- (seu$S.Score + seu$G2M.Score) > 0
################################################################################
## 4) Check MNN reduction exists
################################################################################
if (!("MNN_corrected" %in% names(seu@reductions))) {
  stop("MNN_corrected reduction not found in seu@reductions.")
}

################################################################################
## 5) Remove cell-cycle effects via fastMNN alignment
################################################################################

sce <- as.SingleCellExperiment(seu)
reducedDim(sce, "MNN") <- Embeddings(seu, "MNN_corrected")

sce_aligned <- fastMNN(
  sce,
  batch = seu$Dividing,
  d = 30
)

seu[["MNN_cc"]] <- CreateDimReducObject(
  embeddings = reducedDim(sce_aligned, "corrected"),
  key = "MNNCC_",
  assay = DefaultAssay(seu)
)

################################################################################
## 6) Clustering (for annotation)
################################################################################
seu <- FindNeighbors(seu, reduction="MNN_cc", dims=1:30)
seu <- FindClusters(seu, resolution=0.4)

################################################################################
## 7) Rename EC clusters (EDIT as needed)
################################################################################
seu <- RenameIdents(seu,
                    "0"="CapEC",
                    "1"="CapEC2",
                    "2"="Artery",
                    "3"="PCV",
                    "4"="HEV")

keep_idents <- c("CapEC","CapEC2","PCV","HEV","Artery")
seu_ec <- subset(seu, idents=keep_idents)

################################################################################
## 8) MATLAB tSpace: export -> run MATLAB -> import
################################################################################
work_dir <- "C:/Users/ghumm/Downloads/cyt3-master/cyt3-master/src/wanderlust"
matlab_exe <- "C:/Users/ghumm/OneDrive/Desktop/Research/Single_Cell/bin/matlab.exe"   # change if needed

out_all <- file.path(work_dir, "tspace_allData.csv")
out_te  <- file.path(work_dir, "tspace_tExplain.csv")
out_pe  <- file.path(work_dir, "tspace_pExplain.csv")
in_csv  <- file.path(work_dir, "seurat_export.csv")

cmd <- paste0(
  "cd('", gsub("\\\\", "/", work_dir), "'); ",
  "setenv('path2tSpaceInput','", gsub("\\\\", "/", in_csv), "'); ",
  "setenv('path2tSpaceOutput','", gsub("\\\\", "/", out_all), "'); ",
  "setenv('path2tSpaceOutput2','", gsub("\\\\", "/", out_te), "'); ",
  "setenv('path2tSpaceOutput3','", gsub("\\\\", "/", out_pe), "'); ",
  "[allData,tspacem,megaMat,tExplain,pExplain]=tspace_ml(30,30,1,10,0,5,30); ",
  "exit"
)

system2(
  matlab_exe,
  args = c("-batch", shQuote(cmd)),
  stdout = TRUE,
  stderr = TRUE
)

file.exists(out_all)
file.exists(out_te)
file.exists(out_pe)

## 8a) Export embedding matrix for MATLAB (cells x dims)
tspace_input <- Embeddings(seu_ec, "MNN_cc")
list.files(work_dir, pattern = "\\.csv$", full.names = TRUE)
# Must match cell order
stopifnot(identical(rownames(tspace_input), colnames(seu_ec)))

## 8b) Define file paths
work_dir <- "C:/Users/ghumm/Downloads/cyt3-master/cyt3-master/src/wanderlust"
dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)

in_csv  <- file.path(work_dir, "seurat_export.csv")
out_all <- file.path(work_dir, "tspace_allData.csv")
out_te  <- file.path(work_dir, "tspace_tExplain.csv")
out_pe  <- file.path(work_dir, "tspace_pExplain.csv")

write.csv(tspace_input, in_csv, row.names = FALSE, quote = FALSE)

## 8c) Call MATLAB (EDIT matlab_exe to match your install)
matlab_exe <- "C:/Users/ghumm/OneDrive/Desktop/Research/Single_Cell/bin/matlab.exe"
work_dir   <- "C:/Users/ghumm/Downloads/cyt3-master/cyt3-master/src/wanderlust"
out_all    <- file.path(work_dir, "tspace_allData.csv")

# Save exact cell IDs used for MATLAB input
write.table(
  colnames(seu_ec),
  file = file.path(work_dir, "cell_ids.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

args <- c(
  "-batch",
  sprintf("cd('%s'); run_wanderlust_export", work_dir)
)

out <- tryCatch(
  system2(matlab_exe, args = args, stdout = TRUE, stderr = TRUE),
  error = function(e) e
)

cat(paste(out, collapse = "\n"))

if (!file.exists(out_all)) {
  stop("MATLAB finished but tspace_allData.csv not found.")
}

allData <- as.matrix(read.csv(out_all, header = FALSE))
matlab_cells <- scan(file.path(work_dir, "cell_ids.txt"), what = character(), quiet = TRUE)

# Derive component counts
p_dim <- ncol(tspace_input)
nPC_p <- if (p_dim > 40) 20 else max(2, round(p_dim / 2))

numPop <- 10
nPC_t  <- if (numPop > 40) 20 else max(2, round(numPop / 2))

col_tScore <- (3 + nPC_p):(2 + nPC_p + nPC_t)
tPC <- allData[, col_tScore, drop = FALSE]

cat("dim(tPC):", dim(tPC), "\n")
cat("length(matlab_cells):", length(matlab_cells), "\n")
cat("ncol(seu_ec):", ncol(seu_ec), "\n")

if (nrow(tPC) != length(matlab_cells)) {
  stop("Mismatch: MATLAB output rows do not match exported cell IDs.")
}

# Reorder/subset Seurat object to exact MATLAB cells
seu_ec_tspace <- subset(seu_ec, cells = matlab_cells)
seu_ec_tspace <- seu_ec_tspace[, matlab_cells]

if (ncol(seu_ec_tspace) != nrow(tPC)) {
  stop("Mismatch after subsetting: Seurat cells do not match MATLAB output.")
}

rownames(tPC) <- matlab_cells
colnames(tPC) <- paste0("tPC_", seq_len(ncol(tPC)))
seu_ec_tspace[["tPC"]] <- CreateDimReducObject(
  embeddings = tPC,
  key = "tPC_",
  assay = DefaultAssay(seu_ec_tspace)
)
################################################################################
## 10) Plot tSpace embedding
################################################################################
names(seu_ec_tspace@reductions)

seu_ec <- seu_ec_tspace

p1 <- DimPlot(
  seu_ec,
  reduction = "tPC",
  group.by = "ident",
  label = TRUE,
  pt.size = 1.5
)

p1

pt <- Embeddings(seu_ec_tspace, "tPC")[, 1]
pt <- (pt - min(pt)) / (max(pt) - min(pt))
seu_ec_tspace$MLN_tSpace_pseudotime <- pt

p2 <- FeaturePlot(
  seu_ec_tspace,
  reduction = "tPC",
  features = "MLN_tSpace_pseudotime",
  pt.size = 1.5
) +
  scale_color_viridis_c()

p1
p2

df <- data.frame(
  Embeddings(seu_ec_tspace, "tPC")[, 1:2],
  pt = seu_ec_tspace$MLN_tSpace_pseudotime
)

df <- df[order(df$pt), ]

plot(
  df$tPC_1, df$tPC_2,
  col = scales::alpha("grey30", 0.4),
  pch = 16,
  cex = 0.7,
  xlab = "tPC_1",
  ylab = "tPC_2"
)

ss <- smooth.spline(x = df$tPC_1, y = df$tPC_2, spar = 0.7)
lines(ss$x, ss$y, col = "red", lwd = 3)

emb <- as.data.frame(Embeddings(seu_ec_tspace, "tPC")[, 1:2])
emb$ident <- Idents(seu_ec_tspace)

centers <- aggregate(cbind(tPC_1, tPC_2) ~ ident, data = emb, FUN = median)

centers$ident <- factor(
  centers$ident,
  levels = c("Artery", "CapEC", "CapEC2", "PCV", "HEV")
)

centers <- centers[order(centers$ident), ]

plot(
  emb$tPC_1, emb$tPC_2,
  col = scales::alpha("grey30", 0.35),
  pch = 16,
  cex = 0.7,
  xlab = "tPC_1",
  ylab = "tPC_2"
)

points(centers$tPC_1, centers$tPC_2, col = "red", pch = 16, cex = 2)
lines(centers$tPC_1, centers$tPC_2, col = "red", lwd = 3)

text(centers$tPC_1, centers$tPC_2, labels = centers$ident, pos = 3)