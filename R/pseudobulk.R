
#' Pseudobulks sce data over donor x cell type
#'
#' @param sce Single cell experiment object
#' @param cellType_col Column name on colData of the sce that denotes the celltype
#'
#' @return
#' @export
#'
#' @examples
pseudobulk <- function(sce, cellType_col = "cellType"){
  sce$pb <- paste0(sce$donor,"_",sce[[cellType_col]])

  clusIndex = rafalib::splitit(sce$pb)
  return(clusIndex)
  pbcounts <- purrr::map(clusIndex, ~rowSums(assays(sce)$counts[ ,.x]))

  # pbcounts <- sapply(clusIndex, function(ii){
  #   rowSums(assays(sce)$counts[ ,ii])
  # }
  # )
  # Compute LSFs at this level
  sizeFactors.PB.all <- scuttle::librarySizeFactors(pbcounts)
  # Normalize with these LSFs
  geneExprs.temp <- t(apply(pbcounts, 1, function(x) {log2(x/sizeFactors.PB.all + 1)}))
  rownames(geneExprs.temp) <- rowData(sce)$Symbol
  return(geneExprs.temp)
}
