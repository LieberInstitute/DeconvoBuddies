#' Calculate 1 vs. All standard fold change for each gene x cell type, wrapper
#' function for scran::findMarkers
#'
#' @param sce single cell experiment object
#' @param assay_name Name of the assay to use for calculation
#' @param cellType_col Column name on colData of the sce that denotes the celltype
#' @param add_symbol Add the gene symbol column to the marker stats table
#' @param mod String specifying the model used as design in findMarkers. 
#' Can be `NULL` (default) if there are no blocking terms with uninteresting 
#' factors as documented at [pairwiseTTests][scran::pairwiseTTests].
#' @param verbose Boolean choosing to print progress messages or not
#'
#' @return Tibble of 1 vs. ALL std log fold change + p-values for each gene x cell type
#' -   `gene` is the name of the gene (from rownames(`sce`)).
#' -   `logFC` the log fold change from the DE test
#' -   `log.p.value` the log of the p-value of the DE test
#' -   `log.FDR` the log of the False Discovery Rate adjusted p.value
#' -   `std.logFC` the standard logFC
#' -   `cellType.target` the cell type we're finding marker genes for
#' -   `std.logFC.rank` the rank of `std.logFC` for each cell type
#' -   `std.logFC.anno` is an annotation of the `std.logFC` value
#'    helpful for plotting.
#' @export
#'
#' @examples
#' ## load example SingleCellExperiment
#' if (!exists("sce_DLPFC_example")) sce_DLPFC_example <- fetch_deconvo_data("sce_DLPFC_example")
#' ## Explore properties of the sce object
#' sce_DLPFC_example
#'
#' ## this data contains logcounts of gene expression
#' SummarizedExperiment::assays(sce_DLPFC_example)$logcounts[1:5, 1:5]
#'
#' ## nuclei are classified in to cell types
#' table(sce_DLPFC_example$cellType_broad_hc)
#'
#' ## Get the 1vALL stats for each gene for each cell type defined in `cellType_broad_hc`
#' marker_stats_1vAll <- findMarkers_1vAll(
#'  sce = sce_DLPFC_example,
#'  assay_name = "logcounts",
#'  cellType_col = "cellType_broad_hc",
#'  mod = "~BrNum"
#' )
#' 
#' head(marker_stats_1vAll)
#'
#' @importFrom purrr map
#' @importFrom dplyr mutate
#' @importFrom scran findMarkers
#' @importFrom tibble rownames_to_column as_tibble add_column
findMarkers_1vAll <- function(sce, 
                              assay_name = "counts", 
                              cellType_col = "cellType", 
                              add_symbol = FALSE,
                              mod = NULL, 
                              verbose = TRUE) {
    # RCMD Fix
    gene <- rank_marker <- cellType.target <- std.logFC <- rowData <- Symbol <- NULL

    cell_types <- unique(sce[[cellType_col]])
    names(cell_types) <- cell_types

    ## Traditional t-test with design as in PB'd/limma approach
    pd <- as.data.frame(SummarizedExperiment::colData(sce))
    # message("donor" %in% colnames(pd))

    if (!is.null(mod)) {
        mod <- with(pd, stats::model.matrix(as.formula(mod)))
        mod <- mod[, -1, drop = FALSE] # intercept otherwise automatically dropped by `findMarkers()`
    }

    markers.t.1vAll <- map(cell_types, function(x) {
        if (verbose) message(Sys.time(), " - Find markers for: ", x)

        sce$contrast <- ifelse(sce[[cellType_col]] == x, 1, 0)

        fm <- scran::findMarkers(sce,
            groups = sce$contrast,
            assay.type = assay_name, design = mod, test.type = "t",
            direction = "up", pval.type = "all", full.stats = TRUE
        )
        fm <- fm[[2]]$stats.0

        fm.std <- scran::findMarkers(sce,
            groups = sce$contrast,
            assay.type = assay_name, design = mod, test.type = "t",
            std.lfc = TRUE,
            direction = "up", pval.type = "all", full.stats = TRUE
        )
        fm.std <- fm.std[[2]]$stats.0
        colnames(fm.std)[[1]] <- "std.logFC"

        return(cbind(fm, fm.std[, 1, drop = FALSE]))
    })

    if (verbose) message("Building Table - ", Sys.time())
    markers.t.1vAll.table <- do.call("rbind", markers.t.1vAll) |>
        as.data.frame() |>
        tibble::rownames_to_column("gene") |>
        mutate(gene = gsub("\\.\\d+", "", gene)) |>
        tibble::as_tibble() |>
        tibble::add_column(cellType.target = rep(cell_types, each = nrow(sce))) |>
        dplyr::group_by(cellType.target) |>
        mutate(
            std.logFC.rank = dplyr::row_number(),
            std.logFC.anno = paste0(" std logFC = ", round(std.logFC, 3))
        )

    if (add_symbol) {
        markers.t.1vAll.table$Symbol <- rowData(sce)[markers.t.1vAll.table$gene, ]$Symbol
        markers.t.1vAll.table <- markers.t.1vAll.table |>
            mutate(feature_marker = paste0(stringr::str_pad(rank_marker, 4, "left"), ": ", Symbol))
    }

    if (verbose) message("** Done! **\n")
    return(markers.t.1vAll.table)
}
