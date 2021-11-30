
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
#' ## Default, pseudobulk by cellType and donor
#' pb_label_ab <- pseudobulk(sce_ab, just_counts = FALSE)
#' head(pb_label_ab)
#' 
#' ## Psuedobulk by user defined columns
#' pb_label_ab_region <- pseudobulk(sce_ab, cell_group_cols = c("cellType", "donor", "region"))
#' head(pb_label_ab_region)
#'
#' ## Use the Gene Symbol as the rownames
#' pb_label_test <- pseudobulk(sce.test, cell_group_cols = c("donor", "cellType.Broad"), add_symbol = TRUE)
#' head(pb_label_test)
#' 
#' ## Only output the count matrix 
#' pb_label_test_counts <- pseudobulk(sce.test, cell_group_cols = c("donor", "cellType.Broad"), just_counts = TRUE)
#' head(pb_label_test_counts)
#' 
#' @importFrom rafalib splitit
#' @importFrom SummarizedExperiment SummarizedExperiment colData assays rowData
pseudobulk <- function(sce, cell_group_cols = c("donor", "cellType"), add_symbol = FALSE, just_counts = FALSE) {

    ## check all columns exist
    stopifnot(all(cell_group_cols %in% colnames(SummarizedExperiment::colData(sce))))

    ## create pd label col
    pb_label <- sce[[cell_group_cols[[1]]]]
    for (c in cell_group_cols[-1]) {
        pb_label <- paste0(pb_label, "_", sce[[c]])
    }
    message("Unique Groups: ", length(unique(pb_label)))
    sce$pb_label <- pb_label
    
    ## create phenotype data
    pd <- colData(sce)[, cell_group_cols]
    pd <- unique(pd)
    rownames(pd) <- unique(pb_label)

    clusIndex <- suppressWarnings(rafalib::splitit(sce$pb_label))

    pb_labelcounts <- vapply(clusIndex, function(ii) {
        rowSums(
            as.matrix(SummarizedExperiment::assays(sce)$counts[, ii, drop = FALSE])
        )
    }, double(nrow(sce)))

    # Compute Library Size Factors at this level
    sizeFactors.pb_label.all <- colSums(pb_labelcounts)
    sizeFactors.pb_label.all <- sizeFactors.pb_label.all / mean(sizeFactors.pb_label.all)
    
    # Normalize with these LSFs
    geneExprs.temp <- t(apply(pb_labelcounts, 1, function(x) {
        log2(x / sizeFactors.pb_label.all + 1)
    }))
    
    if (add_symbol) {
        stopifnot("Symbol" %in% colnames(SummarizedExperiment::rowData(sce)))
        message("Using Symbol as rownames")
        rownames(sce) <- SummarizedExperiment::rowData(sce)$Symbol
    }
    
    ## match col and row names
    colnames(geneExprs.temp) <- rownames(pd)
    rownames(geneExprs.temp) <- rownames(sce)
    
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
