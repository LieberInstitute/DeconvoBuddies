#' Calculate pw standard fold change for each gene x cell type
#'
#' @param sce single cell experiment object
#' @param cellType_col Column name on colData of the sce that denotes the celltype
#' @param assay_name Name of the assay to use for calculation
#' @param add_symbol Add the gene symbol column to the marker stats table
#'
#' @return Table of 1 vs. ALL std log fold change + p-values for each gene x cell type
#' @export
#'
#' @examples
#' markers_pw <- findMarkers_pw(sce.test)
#' head(markers_pw)
#' @importFrom purrr map2
#' @importFrom dplyr mutate
#' @importFrom scran findMarkers
findMarkers_pw <- function(sce, assay_name = "counts", cellType_col = "cellType", add_symbol = FALSE) {
    cell_types <- unique(sce[[cellType_col]])
    names(cell_types) <- cell_types

    ## Traditional t-test with design as in PB'd/limma approach
    pd <- as.data.frame(SummarizedExperiment::colData(sce))
    # message("donor" %in% colnames(pd))

    mod <- with(pd, stats::model.matrix(~donor))
    mod <- mod[, -1, drop = F] # intercept otherwise automatically dropped by `findMarkers()`

    fm.std <- scran::findMarkers(sce,
        groups = sce[[cellType_col]],
        assay.type = assay_name, design = mod, test.type = "t",
        std.lfc = TRUE,
        direction = "up", pval.type = "all", full.stats = F
    )


    ct_stats <- mapply(.combine_stats, fm.std, names(fm.std))
    all_stats <- do.call(rbind, ct_stats)
    return(all_stats)
}

.combine_stats <- function(stats, cellType.target) {
    stats$cellType.target <- cellType.target
    stats$Gene <- rownames(stats)
    rownames(stats) <- NULL
    stats <- stats[, c("Gene", "cellType.target", "summary.logFC", "p.value", "FDR")]
    return(stats)
}
