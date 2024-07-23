

#' Download Human DLPFC Deconvolution Data
#' 
#'  This function downloads from ExperimentHub, if ExperimentHub is not 
#'  available, this function will download the files from Dropbox using 
#'  BiocFileCache::bfcrpath() unless the files are present already at destdir. 
#'  
#'  Note that ExperimentHub and BiocFileCache will cache the data and 
#'  automatically detect if you have previously downloaded it, thus making it 
#'  the preferred way to interact with the data. 
#'  
#'  Based on spatialLIBD::fetch_data()
#'
#' @param type A `character(1)` specifying which file you want to download. 
#' `rse_gene` {RangedSummarizedExperiment} with 110 bulk RNA-seq samples
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
#' ## explore data
#' rse_gene
#' # class: RangedSummarizedExperiment 
#' # dim: 21745 110 
#' # metadata(1): SPEAQeasy_settings
#' # assays(2): counts logcounts
#' # rownames(21745): ENSG00000227232.5 ENSG00000278267.1 ... ENSG00000210195.2 ENSG00000210196.2
#' # rowData names(11): Length gencodeID ... gencodeTx passExprsCut
#' # colnames(110): 2107UNHS-0291_Br2720_Mid_Bulk 2107UNHS-0291_Br2720_Mid_Cyto ... AN00000906_Br8667_Mid_Cyto
#' # AN00000906_Br8667_Mid_Nuc
#' # colData names(78): SAMPLE_ID Sample ... diagnosis qc_class
#' 
#' ## load example snRNA-seq data
#' \dontrun{ ## TODO fix
#' ## A SingleCellExperiment (4.79 MB)
#' if (!exists("sce_DLPFC_example")) sce_DLPFC_example <- fetch_deconvo_data("sce_DLPFC_example")
#' 
#' 
#' sce_DLPFC_example
#' #class: SingleCellExperiment 
#' #dim: 557 10000 
#' #metadata(3): Samples cell_type_colors cell_type_colors_broad
#' #assays(2): counts logcounts
#' #rownames(557): GABRD PRDM16 ... AFF2 MAMLD1
#' #rowData names(7): source type ... gene_type binomial_deviance
#' #colnames(10000): 16_CGAAGTTTCGCACGAC-1 11_TTGGGATCAACCGCTG-1 ... 8_CGCATAAGTTAAACCC-1 16_AGCTACATCCCGAGAC-1
#' #colData names(32): Sample Barcode ... cellType_layer layer_annotation
#' #reducedDimNames(0):
#' #   mainExpName: NULL
#' #altExpNames(0):
#' }
#' 
#' \dontrun{
#' sce_path_zip <- fetch_deconvo_data("sce")
#' sce_path <- unzip(sce_path_zip, exdir = tempdir())
#' sce <- HDF5Array::loadHDF5SummarizedExperiment(
#'     file.path(tempdir(), "sce_DLPFC_annotated")
#' )
#' }
fetch_deconvo_data <- function(type = c("rse_gene", "sce", "sce_DLPFC_example"),
                             destdir = tempdir(),
                             eh = ExperimentHub::ExperimentHub(),
                             bfc = BiocFileCache::BiocFileCache()){
  
  rse_gene <- sce_DLPFC_example <- NULL
  
  ## Choose a type among the valid options
  type <- match.arg(type)
  
  ## Check inputs
  stopifnot(methods::is(eh, "ExperimentHub"))
  
  if(type == "rse_gene") {
    tag <- "Human_DLPFC_deconvolution_bulkRNAseq_DeconvoBuddies"
    hub_title <-
      "Human_DLPFC_deconvolution_bulkRNAseq_DeconvoBuddies"

    ## While EH is not set-up
    file_name <-
      "Human_DLPFC_deconvolution_bulkRNAseq_DeconvoBuddies"
    url <-
      "https://www.dropbox.com/scl/fi/9eyg9e1r98t73wyzsuxhr/rse_gene.Rdata?rlkey=sw2djr71y954yw4o3xrmjv59b&dl=1"
  } else if(type == "sce_DLPFC_example") {
    tag <- "Human_DLPFC_deconvolution_snRNAseq_DeconvoBuddies"
    hub_title <-
      "Human_DLPFC_deconvolution_snRNAseq_DeconvoBuddies"

    ## While EH is not set-up
    file_name <-
      "Human_DLPFC_deconvolution_example_snRNAseq_DeconvoBuddies"
    url <-
      "https://www.dropbox.com/scl/fi/w9q71rh7rd36c2l50nyi2/sce_DLPFC_example.Rdata?rlkey=v3z4u8ru0d2y12zgdl1az07q9&st=1dcfqc1i&dl=1"
  } else if(type == "sce"){
    
    sce_path <- spatialLIBD::fetch_data("spatialDLPFC_snRNAseq") 
    return(sce_path)
  }
  else {
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
  
    }else if (type == "sce_DLPFC_example") {
    
    load(file_path)
    return(sce_DLPFC_example)
  
    } else {
    file_path
    } 
  
}
