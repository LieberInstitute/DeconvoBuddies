#' Test Estimated Cell Type Proportions
#'
#' A test dataset of estiamted protions for 5 cell types over 100 samples.
#' 
#' 11.62 kB
#'
#' @format A data.frame object.
#' @source <https://github.com/LieberInstitute/brainstorm/blob/master/data-raw/est_prop.R>
"est_prop"

#' Markers stats from sce.test
#'
#' A tibble containing the marker stats from `get_mean_ratio2` for `sce.test`.
#' 
#' 402.60 kB
#' 
#' @format tibble.
#' @source <https://github.com/LieberInstitute/brainstorm/blob/master/data-raw/marker_test.R>
"marker_test"

#' Test bulk rse dataset
#'
#' A test rse_gene object, 1000 genes x 100 samples
#' 
#' 976.77 kB
#' @format A SummarizedExperiment object.
#' @source <https://github.com/LieberInstitute/brainstorm/blob/master/data-raw/genotyped.R>
"rse_bulk_test"

#' Test subset of a real snRNA-seq data
#'
#' A subset of 1000 genes and 500 nuclei from Tran, Maynard et al., Neuron 2021
#' including 8 cell types. 
#' 
#' 3.44 MB
#'
#' @format A SingleCellExperiment object.
#' @source <https://github.com/LieberInstitute/brainstorm/blob/master/data-raw/sce.test.R>
"sce.test"

#' Toy SCE object for testing
#'
#' An example sce object for testing, with two cell types A and B. 
#' 
#' Generated with DeconvoBuddies::make_test_sce()
#' 38.26 kB
#'
#' @format A SingleCellExperiment object.
#' @source <https://github.com/LieberInstitute/brainstorm/blob/master/data-raw/sce_ab.R>
"sce_ab"

#' Cell Type Proportions estimated from RNAScope
#' 
#' Cell type proportion estimates from high quality images from 
#' Huuki-Myers et al., bioRxiv, 2024, 
#' doi: <https://doi.org/10.1101/2024.02.09.579665>.
#'  
#' 11.49 kB
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
