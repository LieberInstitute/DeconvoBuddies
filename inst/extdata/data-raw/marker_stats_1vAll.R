## code to prepare `marker_stats_1vAll` dataset goes here
## Use spatialLIBD to fetch the snRNA-seq dataset
sce_path_zip <- fetch_deconvo_data("sce")

## unzip and load the data
sce_path <- unzip(sce_path_zip, exdir = tempdir())
sce <- HDF5Array::loadHDF5SummarizedExperiment(
  file.path(tempdir(), "sce_DLPFC_annotated")
)

## exclude Ambiguous cell type
sce <- sce[, sce$cellType_broad_hc != "Ambiguous"]
sce$cellType_broad_hc <- droplevels(sce$cellType_broad_hc)

## We're going to subset to the first 5k genes to save memory
sce <- sce[seq_len(5000), ]

## check the final dimensions of the dataset
dim(sce)

marker_stats_1vAll <- findMarkers_1vAll(
  sce = sce, # sce is the SingleCellExperiment with our data
  assay_name = "counts",
  cellType_col = "cellType_broad_hc", # column in colData with cell type info
  mod = "~BrNum" # Control for donor stored in "BrNum" with mod
)

# lobstr::obj_size(marker_stats_1vAll)
# 3.47 MB

usethis::use_data(marker_stats_1vAll, overwrite = TRUE)
