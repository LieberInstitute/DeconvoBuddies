
#' Pseudobulks sce data over donor x cell type
#'
#' @param sce Single cell experiment object
#' @param cell_group_cols list of column names that create pseudobulk groups
#' @param add_symbol Use Symbol from rowData as rownames
#' @return Matrix of pseudobukled logcounts
#' @export
#'
#' @examples
#' pb_ab <- pseudobulk(sce_ab)
#' head(pb_ab)
#'
#' pb_ab_region <- pseudobulk(sce_ab, cell_group_cols = c("cellType", "donor", "region"))
#' head(pb_ab_region)
#'
#' pb_test <- pseudobulk(sce.test, cell_group_cols = c("donor", "cellType.Broad"))
#' head(pb_test)
#' @importFrom rafalib splitit
#' @importFrom SummarizedExperiment colData assays rowData
pseudobulk <- function(sce, cell_group_cols = c("donor", "cellType"), add_symbol = FALSE) {

    ## check all columns exist
    stopifnot(all(cell_group_cols %in% colnames(SummarizedExperiment::colData(sce))))

    ## create pd label col
    pb <- sce[[cell_group_cols[[1]]]]
    for (c in cell_group_cols[-1]) {
        pb <- paste0(pb, "_", sce[[c]])
    }
    # message("Unique Groups: ", length(unique(pb)))
    sce$pb <- pb

    clusIndex <- suppressWarnings(rafalib::splitit(sce$pb))
    # pbcounts <- purrr::map(clusIndex, ~rowSums(assays(sce)$counts[ ,.x]))

    pbcounts <- vapply(clusIndex, function(ii) {
        rowSums(
            as.matrix(SummarizedExperiment::assays(sce)$counts[, ii, drop = FALSE])
        )
    }, double(nrow(sce)))

    # Compute LSFs at this level
    # sizeFactors.PB.all <- scuttle::librarySizeFactors(pbcounts)
    sizeFactors.PB.all <- colSums(pbcounts)
    sizeFactors.PB.all <- sizeFactors.PB.all / mean(sizeFactors.PB.all)
    # Normalize with these LSFs
    geneExprs.temp <- t(apply(pbcounts, 1, function(x) {
        log2(x / sizeFactors.PB.all + 1)
    }))

    if (add_symbol) {
        rownames(geneExprs.temp) <- rowData(sce)$Symbol
    }

    return(geneExprs.temp)
}
