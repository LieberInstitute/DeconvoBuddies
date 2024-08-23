#' Download Human DLPFC Deconvolution Data
#'
#'  This function downloads the processed data for the experiment documented
#'  at <https://github.com/LieberInstitute/Human_DLPFC_Deconvolution>.
#'  Internally, this function downloads the data from `ExperimentHub`.
#'
#'  Ee currently waiting for <https://doi.org/10.1101/2024.02.09.579665> to
#'  pass peer review at a journal, which could lead to changes requested by the
#'  peer reviewers on the processed data for this study. Thus, this function
#'  temporarily downloads the files from Dropbox using
#'  `BiocFileCache::bfcrpath()` unless the files are present already at
#'  `destdir`.
#'
#'  Note that `ExperimentHub` and `BiocFileCache` will cache the data and
#'  automatically detect if you have previously downloaded it, thus making it
#'  the preferred way to interact with the data.
#'
#'  This function is based on `spatialLIBD::fetch_data()`.
#'
#' @param type A `character(1)` specifying which file you want to download.
#' * `rse_gene`: A [RangedSummarizedExperiment-class][SummarizedExperiment::RangedSummarizedExperiment-class]
#' with 110 bulk RNA-seq samples x 21k genes. (41 MB)
#' * `sce`: A [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment-class]
#' object with Human DLPFC snRNA-seq data. 77k nuclei x 36k genes (172 MB)
#' * `sce_DLPFC_example`: An example subset of `sec`
#' [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment-class]
#' with 10k nuclei x 557 genes (49 MB)
#'
#' @param destdir The destination directory to where files will be downloaded
#' to in case the `ExperimentHub` resource is not available. If you already
#' downloaded the files, you can set this to the current path where the files
#' were previously downloaded to avoid re-downloading them.
#' @param eh An `ExperimentHub` object
#' [ExperimentHub-class][ExperimentHub::ExperimentHub-class].
#' @param bfc A `BiocFileCache` object
#' [BiocFileCache-class][BiocFileCache::BiocFileCache-class]. Used when
#' `eh` is not available.
#'
#' @return The requested object: `rse_gene` that you assign to an object
#' @export
#' @import ExperimentHub
#' @import BiocFileCache
#' @importFrom AnnotationHub query
#' @importFrom methods is
#' @importFrom spatialLIBD fetch_data
#'
#' @examples
#' ## Download the bulk RNA gene expression data
#' ## A RangedSummarizedExperiment (41.16 MB)
#'
#' if (!exists("rse-gene")) rse_gene <- fetch_deconvo_data("rse_gene")
#'
#' ## explore bulk data
#' rse_gene
#'
#' ## load example snRNA-seq data
#' ## A SingleCellExperiment (4.79 MB)
#' if (!exists("sce_DLPFC_example")) sce_DLPFC_example <- fetch_deconvo_data("sce_DLPFC_example")
#'
#' ## explore example sce data
#' sce_DLPFC_example
#'
#' ## check the logcounts
#' SingleCellExperiment::logcounts(sce_DLPFC_example)[1:5, 1:5]
#'
#' \dontrun{
#' ## download the full sce experiment object
#' sce_path_zip <- fetch_deconvo_data("sce")
#' sce_path <- unzip(sce_path_zip, exdir = tempdir())
#' sce <- HDF5Array::loadHDF5SummarizedExperiment(
#'     file.path(tempdir(), "sce_DLPFC_annotated")
#' )
#' }
fetch_deconvo_data <- function(type = c("rse_gene", "sce", "sce_DLPFC_example"),
    destdir = tempdir(),
    eh = ExperimentHub::ExperimentHub(),
    bfc = BiocFileCache::BiocFileCache()) {
    rse_gene <- sce_DLPFC_example <- NULL

    ## Choose a type among the valid options
    type <- match.arg(type)

    ## Check inputs
    stopifnot(methods::is(eh, "ExperimentHub"))

    if (type == "rse_gene") {
        tag <- "Human_DLPFC_deconvolution_bulkRNAseq_DeconvoBuddies"
        hub_title <-
            "Human_DLPFC_deconvolution_bulkRNAseq_DeconvoBuddies"

        ## While EH is not set-up
        file_name <-
            "Human_DLPFC_deconvolution_bulkRNAseq_DeconvoBuddies"
        url <-
            "https://www.dropbox.com/scl/fi/9eyg9e1r98t73wyzsuxhr/rse_gene.Rdata?rlkey=sw2djr71y954yw4o3xrmjv59b&dl=1"
    } else if (type == "sce_DLPFC_example") {
        tag <- "Human_DLPFC_deconvolution_snRNAseq_DeconvoBuddies"
        hub_title <-
            "Human_DLPFC_deconvolution_snRNAseq_DeconvoBuddies"

        ## While EH is not set-up
        file_name <-
            "Human_DLPFC_deconvolution_example_snRNAseq_DeconvoBuddies"
        url <-
            "https://www.dropbox.com/scl/fi/w9q71rh7rd36c2l50nyi2/sce_DLPFC_example.Rdata?rlkey=v3z4u8ru0d2y12zgdl1az07q9&st=1dcfqc1i&dl=1"
    } else if (type == "sce") {
        sce_path <- spatialLIBD::fetch_data("spatialDLPFC_snRNAseq")
        return(sce_path)
    } else {
        stop("Datatype Unkown")
    }

    file_path <- file.path(destdir, file_name)

    ## Use local data if present
    if (!file.exists(file_path)) {
        q <-
            AnnotationHub::query(eh,
                pattern = c(tag, hub_title)
            )

        if (length(q) == 1) {
            ## ExperimentHub has the data - will eventually upload
            # return(res)
        } else {
            ## ExperimentHub backup: download from Dropbox
            file_path <- BiocFileCache::bfcrpath(bfc, url)
        }
    }

    ## Now load the data if possible
    message(Sys.time(), " loading file ", file_path)

    if (type == "rse_gene") {
        load(file_path)
        return(rse_gene)
    } else if (type == "sce_DLPFC_example") {
        load(file_path)
        return(sce_DLPFC_example)
    } else {
        file_path
    }
}
