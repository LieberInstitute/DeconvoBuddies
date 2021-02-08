
#' Pseudobulks sce data over donor x cell type
#'
#' @param sce Single cell experiment object
#' @param cellType_col Column name on colData of the sce that denotes the celltype
#' @param add_symbol Use Symbol from rowData as rownames
#' @return
#' @export
#'
#' @examples
#'pb_ab <- pseudobulk(sce_ab)
#'pb_test <- pseudobulk(sce.test, cellType_col = "cellType.Broad")

pseudobulk <- function(sce, cellType_col = "cellType", add_symbol = FALSE){
  sce$pb <- paste0(sce$donor,"_",sce[[cellType_col]])

  clusIndex = suppressWarnings(rafalib::splitit(sce$pb))
  # pbcounts <- purrr::map(clusIndex, ~rowSums(assays(sce)$counts[ ,.x]))

  pbcounts <- sapply(clusIndex, function(ii){
    rowSums(
      as.matrix(assays(sce)$counts[ ,ii, drop = FALSE])
      )
  }
  )

  # Compute LSFs at this level
  # sizeFactors.PB.all <- scuttle::librarySizeFactors(pbcounts)
  sizeFactors.PB.all <- colSums(pbcounts)
  # Normalize with these LSFs
  geneExprs.temp <- t(apply(pbcounts, 1, function(x) {log2(x/sizeFactors.PB.all + 1)}))

  if(add_symbol){
    rownames(geneExprs.temp) <- rowData(sce)$Symbol
  }

  return(geneExprs.temp)
}
