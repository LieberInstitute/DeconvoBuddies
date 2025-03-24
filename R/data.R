#' Test Estimated Cell Type Proportions
#'
#' A test dataset of estimated proportions for 5 cell types over 100 samples.
#'
#' 11.62 kB
#'
#' @details These are the columns of the `data.frame` object:
#' * cell_A: estimated proportions for cell type A
#' * cell_B: estimated proportions for cell type B
#' * cell_C: estimated proportions for cell type C
#' * cell_D: estimated proportions for cell type D
#' * cell_E: estimated proportions for cell type E
#'
#' @examples
#' ## R Note that the `rowSums(est_prop)` is equal to 1,
#' ## with a small error tolerance.
#' data("est_prop")
#' summary(rowSums(est_prop) - 1)
#'
#' ## You can check this yourself with:
#' all(round(rowSums(est_prop), 3) == 1)
#'
#' # To view source
#' system.file("extdata", "data-raw", "est_prop.R", package = "DeconvoBuddies")
#' @format A `data.frame` object.
#' @source <https://github.com/LieberInstitute/brainstorm/blob/master/inst/extdata/data-raw/est_prop.R>
#' @usage data("est_prop")
"est_prop"

#' Markers stats from sce_DLPFC_example
#'
#' A tibble containing the marker stats from `get_mean_ratio()` for
#' `sce_DLPFC_example`.
#'
#' 402.60 kB
#'
#' @examples
#' # To view source
#' system.file("extdata", "data-raw", "marker_test.R", package = "DeconvoBuddies")
#' @format A `tibble::tibble()`. See `get_mean_ratio()` for more details on the column
#' names.
#' @source <https://github.com/LieberInstitute/brainstorm/blob/master/inst/extdata/data-raw/marker_test.R>
#' @usage data("marker_test")
"marker_test"

#' Test bulk rse dataset
#'
#' A test `rse_gene` object with data for 1000 genes across 100 samples.
#'
#' 976.77 kB
#'
#' @examples
#' # To view source
#' system.file("extdata", "data-raw", "rse_bulk_test.R", package = "DeconvoBuddies")
#' @format A [SummarizedExperiment][SummarizedExperiment::SummarizedExperiment-class] object.
#' @source <https://github.com/LieberInstitute/brainstorm/blob/master/inst/extdata/data-raw/rse_bulk_test.R>
#' @usage data("rse_bulk_test")
"rse_bulk_test"

#' Toy SCE object for testing
#'
#' An example `sce` object for testing, with two cell types A and B.
#'
#' Generated with `DeconvoBuddies::make_test_sce()`
#' 38.26 kB
#'
#' @examples
#' # To view source
#' system.file("extdata", "data-raw", "sce_ab.R", package = "DeconvoBuddies")
#' @format A [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment-class] object.
#' @source <https://github.com/LieberInstitute/brainstorm/blob/master/inst/extdata/data-raw/sce_ab.R>
#' @usage data("sce_ab")
"sce_ab"

#' Cell Type Proportions estimated from RNAScope
#'
#' Cell type proportion estimates from high quality images from
#' Huuki-Myers et al., bioRxiv, 2024,
#' doi: <https://doi.org/10.1101/2024.02.09.579665>.
#'
#' 11.49 kB
#'
#' These are the columns of the `data.frame` object:
#' * `SAMPLE_ID` : DLPFC Tissue block + RNAScope combination.
#' * `Sample` : DLFPC Tissue block (Donor BrNum + DLPFC position).
#' * `Combo` : RNAScope probe combination, either "Circle" marking cell types Astro
#' Endo, Inhib, or "Star" marking Excit, Micro, and OligoOPC.
#' * `cell_type` : The cell type measured.
#' * `Confidence` : Image confidence, this dataset has been filtered to the high & Ok confidence images.
#' * `n_cell` : the number of cells counted for the Sample and cell type.
#' * `prop` : the calculated cell type proportion from n_cell.
#' * `n_cell_sn` : number of nuclei in the corresponding snRNA-seq data.
#' * `prop_sn` : cell type proportion from the snRNA-seq data.
#' 
#' NOTE: For the RNAScope assay utilized here only 3 cell types could be measured 
#' at once. Two consecutive tissue slices were used from each brain block to 
#' measure one combination of three major cell types, then the other three 
#' (differentiated as circle and star). DAPI was used to mark the nuclei, in 
#' each tissue section, so the overall number of cells is recorded but only a 
#' fraction has a cell type label, unlabeled nuclei are classified "other". 
#' 
#' The two sections combined should get close to identifying the type of all the
#' cells, but often the combined "non-other" fractions are around 0.7 and not 1.
#' This could be due to a few reasons: the sections while as similar as 
#' possible are not the same tissue, error in the assay and/or image processing
#' not labeling all cells possible, or presence of rare cell types not marked 
#' by the selected probes. 
#'   
#' With all of this considered, we still think the RNAScope estimates are good 
#' approximations of the cell type proportions in these samples. 
#'
#' For more info check out <https://doi.org/10.1101/2024.02.09.579665> Figure 2,
#' and Methods: RNAScope/Immunofluorescence Data Generation and HALO Analysis.
#'
#' @examples
#' # To view source
#' system.file("extdata", "data-raw", "RNAScope_prop.R", package = "DeconvoBuddies")
#' @format A `data.frame` object.
#' @source <https://github.com/LieberInstitute/brainstorm/blob/master/inst/extdata/data-raw/RNAScope_prop.R>
#' @usage data("RNAScope_prop")
"RNAScope_prop"

#' 1vAll Marker Statistics example data
#'
#' A tibble contatinting the marker statistics calculated for 5k genes from 
#' DLPFC snRNA-seq dataset by `findMarkers_1vAll`.
#'
#' 3.47 MB
#'
#' @examples
#' # To view source
#' system.file("extdata", "data-raw", "marker_stats_1vAll.R", package = "DeconvoBuddies")
#' @format A `data.frame` object.
#' @source <https://github.com/LieberInstitute/brainstorm/blob/master/inst/extdata/data-raw/RNAScope_prop.R>
#' @usage data("marker_stats_1vAll")
"marker_stats_1vAll"

