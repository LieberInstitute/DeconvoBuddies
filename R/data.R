#' Test Estimated Cell Type Proportions
#'
#' A test dataset of estiamted protions for 5 cell types over 100 samples.
#'
#' @format A data.frame object.
#' @source <https://github.com/LieberInstitute/brainstorm/blob/master/data-raw/est_prop.R>
"est_prop"

#' Markers stats from sce.test
#'
#' A tibble containing the marker stats from `get_mean_ratio2` for `sce.test`.
#'
#' @format tibble.
#' @source <https://github.com/LieberInstitute/brainstorm/blob/master/data-raw/marker_test.R>
"marker_test"

#' Test bulk rse dataset
#'
#' TODO Description.
#'
#' @format A SummarizedExperiment object.
#' @source <https://github.com/LieberInstitute/brainstorm/blob/master/data-raw/genotyped.R>
"rse_bulk_test"

#' TODO genotyped
#'
#' TODO Description.
#'
#' @format A SingleCellExperiment object.
#' @source <https://github.com/LieberInstitute/brainstorm/blob/master/data-raw/sce.test.R>
"sce.test"

#' Toy SCE object for testing
#'
#' TODO Description.
#'
#' @format A SingleCellExperiment object.
#' @source <https://github.com/LieberInstitute/brainstorm/blob/master/data-raw/sce_ab.R>
"sce_ab"

#' Cell Type Proportions estimated from RNAScope
#' 
#' Cell type proportion estimates from high quality images from Deconvolution 
#' Benchmark (TODO cite paper)
#' 
#' `SAMPLE_ID` : DLPFC Tissue block + RNAScope combination. 
#' `Sample` : DLFPC Tissue block (Donor BrNum + DLPFC position).   
#' `Combo` : RNAScope probe combination, either "Circle" marking cell types Astro
#' Endo, Inhib, or "Star" marking Excit, Micro, and OligoOPC. 
#' `cell_type` : The cell type measured. 
#' `Confidence` : Image confidence, this dataset has been filtered to the high & Ok confidence images. 
#' `n_cell` : the number of cells counted for the Sample and cell type. 
#' `prop` : the calculated cell type proportion from n_cell. 
#' `n_cell_sn` : number of nuclei in the corresponding snRNA-seq data. 
#' `prop_sn` : cell type proportion from the snRNA-seq data.
#'
#' @format A data.frame object.
#' @source <https://github.com/LieberInstitute/brainstorm/blob/master/data-raw/RNAScope.R>
"RNAScope_prop"
