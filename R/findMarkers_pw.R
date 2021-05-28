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
#' pw <- findMarkers_pw(sce.test)
#' @importFrom purrr map
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

    # fm <- scran::findMarkers(sce, groups = sce[[cellType_col]],
    #                                  assay.type = assay_name, design=mod, test.type ="t",
    #                                  direction="up", pval.type="all", full.stats=T)
    # summary_fm <- fm[,1:3]

    fm.std <- scran::findMarkers(sce,
        groups = sce[[cellType_col]],
        assay.type = assay_name, design = mod, test.type = "t",
        std.lfc = TRUE,
        direction = "up", pval.type = "all", full.stats = T
    )

    ct_stats <- map(fm.std, names(fm.std), ~ combine_stats(.y, cell_types, .x))
    all_stats <- do.call(rbind, ct_stats)

    return(all_stats)
}

.combine_stats <- function(cellType.target, cell_types, fm) {
    cell_types.nonTarget <- cell_types[cell_types != cellType.target]
    cell_types_cols <- paste0("stats.", cell_types.nonTarget)

    stats <- do.call(rbind, as.list(fm[, cell_types_cols]))
    stats$cellType.target <- cellType.target
    stats$cellType <- rep(cell_types.nonTarget, each = nrow(fm))
    return(stats)
}
