################################################################################
## Unified MLN endothelial workflow for either:
##  - mlnwt_fastmnn.rds
##  - magicbatch_sub_seu(.rds)
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(Matrix)
  library(SingleCellExperiment)
  library(slingshot)
  library(batchelor)
})

set.seed(1234)

################################################################################
## 1) Load object robustly
################################################################################
infile <- "C:/Users/ghumm/Downloads/magicbatch_sub_seu"   # or mlnwt_fastmnn.rds
base_dir <- "MLN_outputs"
dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)

obj <- readRDS(infile)

# Convert only if needed
if (inherits(obj, "Seurat")) {
  seu <- obj
} else if (inherits(obj, "SingleCellExperiment")) {
  assay_names <- SummarizedExperiment::assayNames(obj)
  
  use_counts <- if ("counts" %in% assay_names) "counts" else assay_names[1]
  use_data   <- if ("logcounts" %in% assay_names) "logcounts" else use_counts
  
  seu <- as.Seurat(obj, counts = use_counts, data = use_data)
} else {
  stop("Object is neither Seurat nor SingleCellExperiment.")
}

message("Loaded object class: ", class(obj)[1])
message("Assays: ", paste(names(seu@assays), collapse = ", "))
message("Reductions: ", paste(names(seu@reductions), collapse = ", "))

################################################################################
## 2) Pick assay robustly
################################################################################
assay_names <- names(seu@assays)

preferred_assay <- if ("originalexp" %in% assay_names) {
  "originalexp"
} else if ("RNA" %in% assay_names) {
  "RNA"
} else {
  DefaultAssay(seu)
}

DefaultAssay(seu) <- preferred_assay
message("Using assay: ", preferred_assay)

################################################################################
## 3) Remove contaminants
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

# robust renaming
meta_cols <- colnames(seu@meta.data)
lec_col   <- grep("^LEC[0-9]+$", meta_cols, value = TRUE)[1]
lym_col   <- grep("^Lymphocyte[0-9]+$", meta_cols, value = TRUE)[1]
frc_col   <- grep("^FRC[0-9]+$", meta_cols, value = TRUE)[1]

if (!is.na(lec_col)) colnames(seu@meta.data)[match(lec_col, colnames(seu@meta.data))] <- "LEC_Score"
if (!is.na(lym_col)) colnames(seu@meta.data)[match(lym_col, colnames(seu@meta.data))] <- "Lymphocyte_Score"
if (!is.na(frc_col)) colnames(seu@meta.data)[match(frc_col, colnames(seu@meta.data))] <- "FRC_Score"

if (!"LEC_Score" %in% colnames(seu@meta.data)) {
  stop("LEC_Score not created.")
}

LEC_cutoff <- as.numeric(quantile(seu$LEC_Score, 0.95, na.rm = TRUE))
seu <- subset(seu, subset = LEC_Score < LEC_cutoff)

message("Cells after LEC filter: ", ncol(seu))

################################################################################
## 4) Optional cell-cycle correction only if starting from pre-CC object
################################################################################
# Use this for mlnwt_fastmnn.rds if you want to recreate MNN_cc.
# Skip automatically if MNN_cc already exists.

if (!"MNN_cc" %in% names(seu@reductions) && "MNN_corrected" %in% names(seu@reductions)) {
  message("Computing cell-cycle corrected MNN_cc from MNN_corrected...")
  
  data("cc.genes.updated.2019", package = "Seurat")
  s.genes   <- cc.genes.updated.2019$s.genes
  g2m.genes <- cc.genes.updated.2019$g2m.genes
  
  rn <- rownames(seu)
  s.use   <- intersect(s.genes, rn)
  g2m.use <- intersect(g2m.genes, rn)
  
  if (length(s.use) < 10 || length(g2m.use) < 10) {
    s.use   <- rn[match(toupper(s.genes),   toupper(rn), nomatch = 0)]
    g2m.use <- rn[match(toupper(g2m.genes), toupper(rn), nomatch = 0)]
    s.use   <- unique(s.use[s.use != ""])
    g2m.use <- unique(g2m.use[g2m.use != ""])
  }
  
  if (length(s.use) >= 10 && length(g2m.use) >= 10) {
    seu <- CellCycleScoring(seu, s.features = s.use, g2m.features = g2m.use)
    seu$Dividing <- (seu$S.Score + seu$G2M.Score) > 0
    
    sce <- as.SingleCellExperiment(seu)
    reducedDim(sce, "MNN") <- Embeddings(seu, "MNN_corrected")
    
    sce_aligned <- fastMNN(
      sce,
      batch = seu$Dividing,
      d = min(30, ncol(Embeddings(seu, "MNN_corrected")))
    )
    
    seu[["MNN_cc"]] <- CreateDimReducObject(
      embeddings = reducedDim(sce_aligned, "corrected"),
      key = "MNNCC_",
      assay = DefaultAssay(seu)
    )
  } else {
    message("Not enough CC genes detected; skipping CC correction.")
  }
}

################################################################################
## 5) Choose best non-UMAP reduction
################################################################################
candidate_reductions <- c(
  "MNN_cc",
  "MNN_corrected",
  "MNN_CORRECTED",
  "harmony",
  "integrated_pca",
  "rpca",
  "cca",
  "pca"
)

avail_red <- names(seu@reductions)
use_red <- candidate_reductions[candidate_reductions %in% avail_red][1]

if (is.na(use_red) || is.null(use_red)) {
  message("No suitable reduction found. Computing PCA.")
  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(seu, verbose = FALSE)
  seu <- ScaleData(seu, verbose = FALSE)
  seu <- RunPCA(seu, npcs = 30, verbose = FALSE)
  use_red <- "pca"
}

message("Using reduction: ", use_red)
dims_use <- 1:min(30, ncol(Embeddings(seu, use_red)))

################################################################################
## 6) Clustering and UMAP
################################################################################
umap_name <- paste0("umap.", tolower(use_red))

seu <- FindNeighbors(seu, reduction = use_red, dims = dims_use, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.4, verbose = FALSE)
seu <- RunUMAP(seu, reduction = use_red, dims = dims_use,
               reduction.name = umap_name, verbose = FALSE)

print(DimPlot(seu, reduction = umap_name, label = TRUE) + ggtitle(paste("UMAP from", use_red)))

saveRDS(seu, file = file.path(base_dir, "mln_clustered_unified.rds"))
################################################################################
## 7) Rename clusters after inspection
################################################################################
Idents(seu) <- "seurat_clusters"
print(levels(Idents(seu)))
print(table(Idents(seu)))

# EDIT this mapping to match your current cluster IDs
seu <- RenameIdents(seu,
                    "0" = "Apq7+ CapEC",
                    "1" = "CapEC2",
                    "2" = "Artery",
                    "3" = "PCV",
                    "4" = "CapEC3",
                    "5" = "HEV",
                    "6" = "Activated EC")

seu$celltype <- as.character(Idents(seu))

DimPlot(seu, reduction = umap_name, group.by = "celltype", label = TRUE) +
  ggtitle("Annotated EC populations")
################################################################################
## 8) Subset EC states for MATLAB tSpace / Wanderlust
################################################################################
Idents(seu) <- "celltype"

keep_idents <- c("Apq7+ CapEC", "CapEC2", "CapEC3", "Artery", "PCV", "HEV")
seu_ec <- subset(seu, idents = keep_idents)

Idents(seu_ec) <- "celltype"
print(table(Idents(seu_ec)))

################################################################################
## 9) Choose embedding for MATLAB export
################################################################################
avail_red <- names(seu_ec@reductions)
message("Available reductions: ", paste(avail_red, collapse = ", "))

candidate_reductions <- c(
  "MNN_cc",
  "MNN_corrected",
  "MNN_CORRECTED",
  "mnn",
  "fastmnn",
  "harmony",
  "integrated_pca",
  "rpca",
  "cca",
  "pca"
)

red_use <- candidate_reductions[candidate_reductions %in% avail_red][1]

if (is.na(red_use) || is.null(red_use) || length(red_use) == 0) {
  message("No suitable reduction found. Computing PCA...")
  DefaultAssay(seu_ec) <- if ("originalexp" %in% Assays(seu_ec)) "originalexp" else DefaultAssay(seu_ec)
  seu_ec <- NormalizeData(seu_ec, verbose = FALSE)
  seu_ec <- FindVariableFeatures(seu_ec, verbose = FALSE)
  seu_ec <- ScaleData(seu_ec, verbose = FALSE)
  seu_ec <- RunPCA(seu_ec, npcs = 50, verbose = FALSE)
  red_use <- "pca"
}

message("Using reduction for MATLAB export: ", red_use)

################################################################################
## 10) Shared paths
################################################################################
work_dir <- normalizePath(
  "C:/Users/ghumm/Downloads/cyt3-master/cyt3-master/src/wanderlust",
  winslash = "/",
  mustWork = TRUE
)

matlab_exe <- normalizePath(
  "C:/Users/ghumm/OneDrive/Desktop/Research/Single_Cell/bin/matlab.exe",
  winslash = "/",
  mustWork = TRUE
)

################################################################################
## 11) Wanderlust helper
################################################################################
run_wanderlust_from_seurat <- function(
    seu_ec,
    red_use,
    umap_name,
    base_dir,
    work_dir,
    matlab_exe,
    wait_seconds = 10
) {
  in_csv    <- file.path(work_dir, "seurat_export.csv")
  cell_txt  <- file.path(work_dir, "cell_ids.txt")
  out_csv   <- file.path(work_dir, "wanderlust_pseudotime.csv")
  log_file  <- file.path(work_dir, "matlab_wanderlust_log.txt")
  wrapper_m <- file.path(work_dir, "run_wanderlust_batch.m")
  
  for (f in c(in_csv, cell_txt, out_csv, log_file, wrapper_m)) {
    if (file.exists(f)) file.remove(f)
  }
  
  tspace_input <- Embeddings(seu_ec, red_use)
  
  if (is.null(rownames(tspace_input))) {
    stop("Embeddings matrix has no rownames.")
  }
  
  write.csv(tspace_input, in_csv, row.names = FALSE, quote = FALSE)
  write.table(
    rownames(tspace_input),
    file = cell_txt,
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
  stopifnot(
    nrow(tspace_input) == length(scan(cell_txt, what = character(), quiet = TRUE))
  )
  
  message("Wrote Wanderlust input matrix to: ", in_csv)
  message("Wrote Wanderlust cell IDs to: ", cell_txt)
  
  wrapper_code <- c(
    "try",
    "    cd(fileparts(mfilename('fullpath')));",
    "    addpath(genpath(pwd));",
    "    run('run_wanderlust_export.m');",
    "catch ME",
    "    disp(getReport(ME, 'extended'));",
    "    exit(1);",
    "end",
    "exit(0);"
  )
  writeLines(wrapper_code, wrapper_m)
  
  status <- system2(
    matlab_exe,
    args = c("-batch", "run_wanderlust_batch"),
    stdout = log_file,
    stderr = log_file
  )
  
  cat("Wanderlust MATLAB exit status:", status, "\n")
  
  if (file.exists(log_file)) {
    cat(readLines(log_file), sep = "\n")
  }
  
  if (!identical(status, 0L)) {
    stop("MATLAB Wanderlust failed. Check log: ", log_file)
  }
  
  found <- FALSE
  for (i in seq_len(wait_seconds * 2)) {
    if (file.exists(out_csv)) {
      found <- TRUE
      break
    }
    Sys.sleep(0.5)
  }
  
  if (!found) {
    stop("MATLAB finished but wanderlust_pseudotime.csv was not found.\nCheck log: ", log_file)
  }
  
  wl <- read.csv(out_csv, header = TRUE, stringsAsFactors = FALSE)
  
  if (!all(c("cell", "pseudotime") %in% colnames(wl))) {
    stop("wanderlust_pseudotime.csv must contain columns 'cell' and 'pseudotime'")
  }
  
  seu_ec$Wanderlust_pseudotime <- NA_real_
  
  idx <- match(wl$cell, colnames(seu_ec))
  if (any(is.na(idx))) {
    missing_cells <- wl$cell[is.na(idx)]
    stop(
      "Some Wanderlust output cells were not found in seu_ec. Example: ",
      paste(head(missing_cells, 10), collapse = ", ")
    )
  }
  
  seu_ec$Wanderlust_pseudotime[idx] <- wl$pseudotime
  
  p_wl <- FeaturePlot(
    seu_ec,
    reduction = umap_name,
    features = "Wanderlust_pseudotime",
    pt.size = 1.2
  ) + ggtitle("Wanderlust pseudotime")
  
  print(p_wl)
  
  saveRDS(seu_ec, file = file.path(base_dir, "mln_wanderlust_unified.rds"))
  message("Saved Wanderlust-updated object to: ", file.path(base_dir, "mln_wanderlust_unified.rds"))
  
  seu_ec
}

################################################################################
## 12) tSpace helper
################################################################################
run_tspace_from_seurat <- function(
    seu_ec,
    red_use,
    work_dir,
    matlab_exe,
    numPop = 10
) {
  in_csv   <- file.path(work_dir, "seurat_export.csv")
  out_all  <- file.path(work_dir, "tspace_allData.csv")
  out_te   <- file.path(work_dir, "tspace_tExplain.csv")
  out_pe   <- file.path(work_dir, "tspace_pExplain.csv")
  log_file <- file.path(work_dir, "matlab_tspace_log.txt")
  cell_txt <- file.path(work_dir, "cell_ids.txt")
  
  for (f in c(in_csv, out_all, out_te, out_pe, log_file, cell_txt)) {
    if (file.exists(f)) file.remove(f)
  }
  
  tspace_input <- Embeddings(seu_ec, red_use)
  
  if (is.null(rownames(tspace_input))) {
    stop("Embeddings matrix has no rownames.")
  }
  
  write.csv(tspace_input, in_csv, row.names = FALSE, quote = FALSE)
  write.table(
    rownames(tspace_input),
    file = cell_txt,
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
  stopifnot(
    nrow(tspace_input) == length(scan(cell_txt, what = character(), quiet = TRUE))
  )
  
  message("Wrote tSpace input matrix to: ", in_csv)
  message("Wrote tSpace cell IDs to: ", cell_txt)
  
  cmd <- paste0(
    "cd('", gsub("\\\\", "/", work_dir), "'); ",
    "setenv('path2tSpaceInput','", gsub("\\\\", "/", in_csv), "'); ",
    "setenv('path2tSpaceOutput','", gsub("\\\\", "/", out_all), "'); ",
    "setenv('path2tSpaceOutput2','", gsub("\\\\", "/", out_te), "'); ",
    "setenv('path2tSpaceOutput3','", gsub("\\\\", "/", out_pe), "'); ",
    "[allData,tspacem,megaMat,tExplain,pExplain]=tspace_ml(30,30,1,10,0,5,30); ",
    "exit"
  )
  
  out <- tryCatch(
    system2(
      matlab_exe,
      args = c("-batch", shQuote(cmd)),
      stdout = TRUE,
      stderr = TRUE
    ),
    error = function(e) e
  )
  
  writeLines(as.character(out), log_file)
  
  cat("---- tSpace MATLAB console output ----\n")
  cat(paste(out, collapse = "\n"), "\n")
  
  status <- attr(out, "status")
  cat("---- tSpace MATLAB exit status ----\n")
  print(status)
  
  Sys.sleep(1)
  
  if (!file.exists(out_all)) {
    stop("MATLAB finished but tspace_allData.csv was not found.\nCheck log: ", log_file)
  }
  
  allData <- as.matrix(read.csv(out_all, header = FALSE))
  matlab_cells <- scan(cell_txt, what = character(), quiet = TRUE)
  
  p_dim <- ncol(tspace_input)
  nPC_p <- if (p_dim > 40) 20 else max(2, round(p_dim / 2))
  nPC_t <- if (numPop > 40) 20 else max(2, round(numPop / 2))
  
  col_tScore <- (3 + nPC_p):(2 + nPC_p + nPC_t)
  tPC <- allData[, col_tScore, drop = FALSE]
  
  if (nrow(tPC) != length(matlab_cells)) {
    stop("Mismatch: tSpace output rows do not match exported cell IDs.")
  }
  
  seu_ec_tspace <- subset(seu_ec, cells = matlab_cells)
  seu_ec_tspace <- seu_ec_tspace[, matlab_cells]
  
  if (ncol(seu_ec_tspace) != nrow(tPC)) {
    stop("Mismatch after subsetting: Seurat cells do not match tSpace output.")
  }
  
  rownames(tPC) <- matlab_cells
  colnames(tPC) <- paste0("tPC_", seq_len(ncol(tPC)))
  
  seu_ec_tspace[["tPC"]] <- CreateDimReducObject(
    embeddings = tPC,
    key = "tPC_",
    assay = DefaultAssay(seu_ec_tspace)
  )
  
  pt <- Embeddings(seu_ec_tspace, "tPC")[, 1]
  pt <- (pt - min(pt)) / (max(pt) - min(pt))
  seu_ec_tspace$MLN_tSpace_pseudotime <- pt
  
  seu_ec_tspace
}

################################################################################
## 13) Run Wanderlust
################################################################################
seu_ec <- run_wanderlust_from_seurat(
  seu_ec = seu_ec,
  red_use = red_use,
  umap_name = umap_name,
  base_dir = base_dir,
  work_dir = work_dir,
  matlab_exe = matlab_exe
)

################################################################################
## 14) Run tSpace
################################################################################
seu_ec_tspace <- run_tspace_from_seurat(
  seu_ec = seu_ec,
  red_use = red_use,
  work_dir = work_dir,
  matlab_exe = matlab_exe
)

################################################################################
## 15) Plot both Wanderlust and tSpace
################################################################################
#pairs(Embeddings(seu_ec_tspace, "tPC")[,1:5])

p_cells <- DimPlot(
  seu_ec,
  reduction = umap_name,
  group.by = "celltype",
  label = TRUE,
  pt.size = 1.2
) + ggtitle("EC populations")

p_wl <- FeaturePlot(
  seu_ec,
  reduction = umap_name,
  features = "Wanderlust_pseudotime",
  pt.size = 1.2
) + ggtitle("Wanderlust pseudotime")

p_tspace <- DimPlot(
  seu_ec_tspace,
  reduction = "tPC",
  dims = c(1, 3),
  group.by = "celltype",
  label = TRUE,
  pt.size = 1.5
) + ggtitle("tSpace embedding (tPC_1 vs tPC_3)")

p_tspace_pt <- FeaturePlot(
  seu_ec_tspace,
  reduction = "tPC",
  dims = c(1, 3),
  features = "MLN_tSpace_pseudotime",
  pt.size = 1.5
) + 
  scale_color_viridis_c() +
  ggtitle("tSpace pseudotime (tPC_1 vs tPC_3)")

print(p_cells)
print(p_wl)
print(p_tspace)
print(p_tspace_pt)

################################################################################
## 16) Optional tSpace trajectory-like overlays
################################################################################
df <- data.frame(
  tPC_1 = Embeddings(seu_ec_tspace, "tPC")[, 1],
  tPC_3 = Embeddings(seu_ec_tspace, "tPC")[, 3],
  pt = seu_ec_tspace$MLN_tSpace_pseudotime
)
df <- df[order(df$pt), ]

plot(
  df$tPC_1, df$tPC_3,
  col = scales::alpha("grey30", 0.4),
  pch = 16,
  cex = 0.7,
  xlab = "tPC_1",
  ylab = "tPC_3",
  main = "tSpace with smoothed pseudotime curve"
)

ss <- smooth.spline(x = df$tPC_1, y = df$tPC_3, spar = 0.7)
lines(ss$x, ss$y, col = "red", lwd = 3)

emb <- data.frame(
  tPC_1 = Embeddings(seu_ec_tspace, "tPC")[, 1],
  tPC_3 = Embeddings(seu_ec_tspace, "tPC")[, 3],
  ident = Idents(seu_ec_tspace)
)

centers <- aggregate(cbind(tPC_1, tPC_3) ~ ident, data = emb, FUN = median)

desired_levels <- c("Artery", "Apq7+ CapEC", "CapEC2", "CapEC3", "PCV", "HEV")
centers$ident <- factor(centers$ident, levels = desired_levels)
centers <- centers[order(centers$ident), ]

plot(
  emb$tPC_1, emb$tPC_3,
  col = scales::alpha("grey30", 0.35),
  pch = 16,
  cex = 0.7,
  xlab = "tPC_1",
  ylab = "tPC_3",
  main = "tSpace cluster-center path"
)

points(centers$tPC_1, centers$tPC_3, col = "red", pch = 16, cex = 2)
lines(centers$tPC_1, centers$tPC_3, col = "red", lwd = 3)
text(centers$tPC_1, centers$tPC_3, labels = centers$ident, pos = 3)

################################################################################
## 17) Save final objects
################################################################################
saveRDS(seu_ec, file = file.path(base_dir, "mln_wanderlust_unified.rds"))
saveRDS(seu_ec_tspace, file = file.path(base_dir, "mln_tspace_unified.rds"))

message("Saved Wanderlust object: ", file.path(base_dir, "mln_wanderlust_unified.rds"))
message("Saved tSpace object: ", file.path(base_dir, "mln_tspace_unified.rds"))