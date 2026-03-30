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