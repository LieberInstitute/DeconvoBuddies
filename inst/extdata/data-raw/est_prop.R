## code to prepare `est_prop` dataset goes here

# library("DeconvoBuddies")
library("SingleCellExperiment")
library("Biobase")
library("BisqueRNA")
library("tidyverse")

#### load data ####
## use fetch deconvo data to load rse_gene
rse_gene <- fetch_deconvo_data("rse_gene")
rownames(rse_gene) <- rowData(rse_gene)$Symbol

## Use spatialLIBD to fetch the snRNA-seq dataset used in this project
sce_path_zip <- fetch_deconvo_data("sce")

## unzip and load the data
sce_path <- unzip(sce_path_zip, exdir = tempdir())
sce <- HDF5Array::loadHDF5SummarizedExperiment(
  file.path(tempdir(), "sce_DLPFC_annotated")
)


## exclude Ambiguous cell type
sce <- sce[, sce$cellType_broad_hc != "Ambiguous"]
sce$cellType_broad_hc <- droplevels(sce$cellType_broad_hc)

# table(sce$cellType_broad_hc)

#### get marker genes ####
# calculate the Mean Ratio of genes for each cell type
marker_stats <- get_mean_ratio(sce,
                               cellType_col = "cellType_broad_hc",
                               gene_ensembl = "gene_id",
                               gene_name = "gene_name"
)


marker_genes <- marker_stats |>
  filter(MeanRatio.rank <= 25 & gene %in% rownames(rse_gene))

# check how many genes for each cell type (some genes are not in both datasets)
marker_genes |> count(cellType.target)

# create a vector of marker genes to subset data before deconvolution
marker_genes <- marker_genes |> pull(gene)

## convert bulk data to Expression set, sub-setting to marker genes
## include sample ID
exp_set_bulk <- Biobase::ExpressionSet(
  assayData = assays(rse_gene[marker_genes, ])$counts,
  phenoData = AnnotatedDataFrame(
    as.data.frame(colData(rse_gene))[c("SAMPLE_ID")]
  )
)

#### Run Bisque ####
## convert snRNA-seq data to Expression set, sub-setting to marker genes
## include cell type and donor information
exp_set_sce <- Biobase::ExpressionSet(
  assayData = as.matrix(assays(sce[marker_genes, ])$counts),
  phenoData = AnnotatedDataFrame(
    as.data.frame(colData(sce))[, c("cellType_broad_hc", "BrNum")]
  )
)

## check for nuclei with 0 marker expression
zero_cell_filter <- colSums(exprs(exp_set_sce)) != 0
message("Exclude ", sum(!zero_cell_filter), " cells")

exp_set_sce <- exp_set_sce[, zero_cell_filter]

## Run Bisque with bulk and single cell ExpressionSet inputs
est_prop <- ReferenceBasedDecomposition(
  bulk.eset = exp_set_bulk,
  sc.eset = exp_set_sce,
  cell.types = "cellType_broad_hc",
  subject.names = "BrNum",
  use.overlap = FALSE
)

est_prop <- t(est_prop$bulk.props)

#### Save data ####
setequal(rownames(est_prop), colnames(rse_gene))
usethis::use_data(est_prop, overwrite = TRUE)
