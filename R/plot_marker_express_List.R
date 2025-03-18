#' Plot a nested list of genes as a multi-page pdf
#'
#' This function plots a nested list of genes as a multi-page PDF, one for each
#' sub list. A use case is plotting known marker genes for multiple cell types
#' over cell type clusters with unknown identities.
#'
#' @param gene_list A named `list()` of `character()` vectors containing the
#' names of genes to plot.
#' @param gene_name_col The `character(1)` name of `rowData()` matching the
#' gene name from `gene_list`.
#' @inheritParams plot_marker_express_ALL
#'
#' @return A PDF file with violin plots for the expression of top marker genes
#' for all cell types.
#'
#' @export
#'
#' @examples
#' ## Fetch sce example data
#' if (!exists("sce_DLPFC_example")) sce_DLPFC_example <- fetch_deconvo_data("sce_DLPFC_example")
#'
#' ## Create list-of-lists of genes to plot, names of sub-list become title of page
#' my_gene_list <- list(Inhib = c("GAD2", "SAMD5"), Astro = c("RGS20", "PRDM16"))
#'
#' # Return a list of plots
#' plots <- plot_marker_express_List(
#'     sce_DLPFC_example,
#'     gene_list = my_gene_list,
#'     cellType_col = "cellType_broad_hc"
#' )
#'
#' print(plots[[1]])
#'
#' # Plot marker gene expression to PDF, one page per cell type in stats
#' pdf_file <- tempfile("test_marker_expression_List", fileext = ".pdf")
#'
#' plot_marker_express_List(
#'     sce_DLPFC_example,
#'     gene_list = my_gene_list,
#'     pdf_fn = pdf_file,
#'     cellType_col = "cellType_broad_hc"
#' )
#'
#' if (interactive()) browseURL(pdf_file)
#'
#' @family expression plotting functions
#' @importFrom ggplot2 ggplot geom_violin geom_text facet_wrap stat_summary
#' @importFrom SummarizedExperiment colData
plot_marker_express_List <- function(
        sce,
        gene_list,
        pdf_fn = NULL,
        cellType_col = "cellType",
        gene_name_col = "gene_name",
        color_pal = NULL,
        plot_points = FALSE) {
    stopifnot(cellType_col %in% colnames(colData(sce)))

    if (!identical(rownames(sce), SummarizedExperiment::rowData(sce)[[gene_name_col]])) {
        message("Using ", gene_name_col, " as gene names")
        rownames(sce) <- SummarizedExperiment::rowData(sce)[[gene_name_col]]
    }

    marker_plots <- purrr::map2(
        gene_list, names(gene_list),
        ~ plot_gene_express(
            sce = sce,
            genes = .x,
            cat = cellType_col,
            color_pal = color_pal,
            plot_points = plot_points,
            title = .y
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
