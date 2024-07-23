# A script describing the steps involved in making the data object(s). This 
# includes where the original data were downloaded from, pre-processing, and
# how the final R object was made. Include a description of any steps performed
# outside of R with third party software. Output of the script should be files
# on disk ready to be pushed to S3. If data are to be hosted on a personal web 
# site instead of S3, this file should explain any manipulation of the data 
# prior to hosting on the web site. For data hosted on a public web site with 
# no prior manipulation this file is not needed.

#### sce_test ####
library(SingleCellExperiment)
library(here)

sce_path_zip <- fetch_deconvo_data("sce")
sce_path <- unzip(sce_path_zip, exdir = tempdir())
sce <- HDF5Array::loadHDF5SummarizedExperiment(
  file.path(tempdir(), "sce_DLPFC_annotated"))

# lobstr::obj_size(sce)
# 172.28 MB

## exclude Ambiguous cell type
sce <- sce[,sce$cellType_broad_hc != "Ambiguous"]
sce$cellType_broad_hc <- droplevels(sce$cellType_broad_hc)

dim(sce)
## Check the broad cell type distribution 
table(sce$cellType_broad_hc)

## subset genes
## select 557 genes Mean Ratio > 2 
## https://github.com/LieberInstitute/Human_DLPFC_Deconvolution/blob/main/processed-data/08_bulk_deconvolution/markers_MeanRatio_over2.txt
markers <- scan(here("data-raw", "markers_MeanRatio_over2.txt" ), what="", sep="\n")
sce_DLPFC_example <- sce[rowData(sce)$gene_id %in% markers,]

## select 10k random nuc
sce_DLPFC_example <- sce_DLPFC_example[,sample(colnames(sce_DLPFC_example), 10000)]
table(sce_DLPFC_example$cellType_broad_hc)
# Astro EndoMural     Micro     Oligo       OPC     Excit     Inhib 
# 738       425       297      1877       344      4411      1908 

lobstr::obj_size(sce_DLPFC_example)
# 21.36 MB

class(sce_DLPFC_example)
# [1] "SingleCellExperiment"
# attr(,"package")
# [1] "SingleCellExperiment"
typeof(sce_DLPFC_example)
typeof(counts(sce_DLPFC_example))
# [1] "S4"

## Drop Reduced Dims
reducedDim(sce_DLPFC_example, "GLMPCA_approx") <- NULL
reducedDim(sce_DLPFC_example, "TSNE") <- NULL
reducedDim(sce_DLPFC_example, "UMAP") <- NULL
reducedDim(sce_DLPFC_example, "HARMONY") <- NULL
lobstr::obj_size(sce_DLPFC_example)
# 4.78 MB

counts(sce_DLPFC_example) <- NULL
logcounts(sce_DLPFC_example) <- as.matrix(logcounts(sce_DLPFC_example))

# sce_DLPFC_example <- realize(sce_DLPFC_example)

lobstr::obj_size(sce_DLPFC_example)
# 58.54 MB



# lobstr::obj_size(sce_DLPFC_example)

## output zipped is 13.9 MB
# HDF5Array::saveHDF5SummarizedExperiment(sce_DLPFC_example,
#                                         dir = ("~/Desktop/sce_DLPFC_example"),
#                                         replace = TRUE)
# test <- loadHDF5SummarizedExperiment(("~/Desktop/sce_DLPFC_example"))
# counts(test)

## Rdata file is 4.3 MB
save(sce_DLPFC_example, file = "~/Desktop/sce_DLPFC_example.Rdata")
# 
# load("~/Downloads/sce_DLPFC_example.Rdata", verbose = TRUE)
# SingleCellExperiment::logcounts(sce_DLPFC_example)[1:5,1:5]
# class(assays(sce_DLPFC_example)$logcounts)

