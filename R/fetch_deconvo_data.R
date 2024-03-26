

#' Download Human DLPFC Deconvolution Data
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
#' @return The requested object
#' @export
#' @import ExperimentHub
#' @import BiocFileCache
#'
#' @examples
#' rse_gene <- fetch_deconvo_data("rse_gene")
#' 
fetch_deconvo_data <- function(type = c("rse_gene"),
                             destdir = tempdir(),
                             eh = ExperimentHub::ExperimentHub(),
                             bfc = BiocFileCache::BiocFileCache()){
  
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
  
    } else {
    file_path
  }
  
}
