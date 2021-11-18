
#' Pseudobulks sce data over donor x cell type
#'
#' @param sce Single cell experiment object
#' @param cell_group_cols list of column names that create pseudobulk groups
#' @param add_symbol Use Symbol from rowData as rownames
#' @param just_counts Whether to return just the counts matrix or a SummarizedExperiment object (Default)
#' @return Matrix of pseudobukled logcounts or SummarizedExperiment
#' @export
#'
#' @examples
#' pb_ab <- pseudobulk(sce_ab, just_counts = FALSE)
#' head(pb_ab)
#'
#' pb_ab_region <- pseudobulk(sce_ab, cell_group_cols = c("cellType", "donor", "region"))
#' head(pb_ab_region)
#'
#' pb_test <- pseudobulk(sce.test, cell_group_cols = c("donor", "cellType.Broad"))
#' head(assays(pb_test)$counts)
#' @importFrom rafalib splitit
#' @importFrom SummarizedExperiment SummarizedExperiment colData assays rowData
pseudobulk <- function(sce, cell_group_cols = c("donor", "cellType"), add_symbol = FALSE, just_counts = FALSE) {

    ## check all columns exist
    stopifnot(all(cell_group_cols %in% colnames(SummarizedExperiment::colData(sce))))

    ## create pd label col
    pb <- sce[[cell_group_cols[[1]]]]
    for (c in cell_group_cols[-1]) {
        pb <- paste0(pb, "_", sce[[c]])
    }
    message("Unique Groups: ", length(unique(pb)))
    sce$pb <- pb
    
    ## create phenotype data
    pd <- colData(sce)[, cell_group_cols]
    pd <- unique(pd)
    rownames(pd) <- unique(pb)

    clusIndex <- suppressWarnings(rafalib::splitit(sce$pb))

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
    
    if(just_counts){
        return(geneExprs.temp)
    } else{
        pd_se <- SummarizedExperiment::SummarizedExperiment(
            assays = list(counts = geneExprs.temp),
            rowData = SummarizedExperiment::rowData(sce),
            colData = pd
        )
        return(pd_se)
    }

    
}
