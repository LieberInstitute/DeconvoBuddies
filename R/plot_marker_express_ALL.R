#' Plot the top marker genes for ALL cell types
#'
#' This function plots the top n marker genes for a all cell types based off of
#' the `stats` table from `get_mean_ratio()` in a multi-page PDF file.
#' The gene expression is plotted as violin plot with `plot_gene_express()` and
#' adds annotations to each plot.
#'
#' @param pdf_fn A `character()` of the pdf filename to plot to, if `NULL` returns all plots
#' @inheritParams plot_marker_express
#'
#' @return A PDF file with violin plots for the expression of top marker genes
#' for all cell types.
#' @export
#'
#' @examples
#' #' ## Fetch sce example data
#' if (!exists("sce_DLPFC_example")) sce_DLPFC_example <- fetch_deconvo_data("sce_DLPFC_example")
#'
#' ## load example marker stats
#' data("marker_test")
#'
#' # Plot marker gene expression to PDF, one page per cell type in stats
#' pdf_file <- tempfile("test_marker_expression_ALL", fileext = ".pdf")
#'
#' plot_marker_express_ALL(
#'     sce_DLPFC_example,
#'     cellType_col = "cellType_broad_hc",
#'     stat = marker_test,
#'     pdf_fn = pdf_file
#' )
#'
#' if (interactive()) browseURL(pdf_file)
#' @family expression plotting functions
#' @importFrom ggplot2 ggplot geom_violin geom_text facet_wrap stat_summary
#' @importFrom SummarizedExperiment colData
#' @importFrom purrr map
plot_marker_express_ALL <- function(sce,
    stats,
    pdf_fn = NULL,
    n_genes = 10,
    rank_col = "MeanRatio.rank",
    anno_col = "MeanRatio.anno",
    gene_col = "gene",
    cellType_col = "cellType",
    color_pal = NULL,
    plot_points = FALSE) {
    stopifnot(cellType_col %in% colnames(colData(sce)))

    if (is.factor(sce[[cellType_col]])) {
        cell_types <- levels(sce[[cellType_col]])
    } else {
        cell_types <- unique(sce[[cellType_col]])
    }

    if (!all(cell_types %in% stats$cellType.target)) {
        # missing <- cell_types[!cell_types %in% stats$cellType.target]
        stop("Stats is missing cell types, check you're using the correct marker stats data and cellType_col")
    }

    marker_plots <- purrr::map(
        cell_types,
        ~ plot_marker_express(
            sce = sce,
            stats = stats,
            cell_type = .x,
            n_genes = n_genes,
            rank_col = rank_col,
            anno_col = anno_col,
            gene_col = gene_col,
            cellType_col = cellType_col,
            color_pal = color_pal,
            plot_points = plot_points
        )
    )

    if (is.null(pdf_fn)) {
        return(marker_plots)
    } else {
        grDevices::pdf(pdf_fn)
        print(marker_plots)
        grDevices::dev.off()
    }
}
