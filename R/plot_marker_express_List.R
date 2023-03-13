#' Plot a nested list of genes as a multi-page pdf
#'
#' This function plots a nested list of genes as a multi-page PDF, one for each sub list.
#' A use case is plotting known marker genes for multiple cell types over cell type
#' clusters with unknown identities.
#'
#' @param sce [SummarizedExperiment-class][SummarizedExperiment::SummarizedExperiment-class] object
#' @param gene_list A named list containing the names of genes to plot
#' @param pdf_fn A `character()` of the pdf file name to plot to
#' @param cellType_col The `character()` name of `colData` column containing cell type for `sce` data
#' @param gene_name_col The `character()` name of `rowData` matching the gene name from `gene_list`
#' @param color_pal  A named `character(1)` vector that contains a color pallet matching the values in `cellType_col`.
#' @param plot_points A logical indicating whether to plot points over the violin,
#' defaults to `FALSE` as these often become over plotted and quite large (especially when saved as PDF)
#'
#' @return A pdf with violin plots for the expression of top marker genes for all cell types
#'
#' @export
#'
#' @examples
#' marker_list <- list(Oligo = c("RNF220", "CLIC4"), Astro = c("PRDM16", "F3"))
#' # plot_marker_express_List(sce.test, gene_list = marker_list, pdf_fn = "./plots/test_marker_expression_List.pdf")
#' # plot_marker_express_List(sce.test, gene_list = marker_list, cellType_col = "cellType.Broad", pdf_fn = "./plots/test_marker_expression_List.pdf")
#' # plot_marker_express_List(sce.test, gene_list = marker_list, pdf_fn = "./plots/test_marker_expression_List.pdf", plot_points = TRUE)
#'
#' @family expression plotting functions
#' @importFrom ggplot2 ggplot geom_violin geom_text facet_wrap stat_summary
plot_marker_express_List <- function(
        sce,
        gene_list,
        pdf_fn = "marker_expression.pdf",
        cellType_col = "cellType",
        gene_name_col = "Symbol",
        color_pal = NULL,
        plot_points = FALSE) {
    stopifnot(cellType_col %in% colnames(colData(sce)))

    if (!identical(rownames(sce), rowData(sce)[[gene_name_col]])) {
        message("Using ", gene_name_col, " as gene names")
        rownames(sce) <- rowData(sce)[[gene_name_col]]
    }

    grDevices::pdf(pdf_fn)
    message("Plotting Genes for...")

    for (i in seq(length(gene_list))) {
        n <- names(gene_list)[[i]]
        message("\t", n)
        print(plot_gene_express(
            sce = sce,
            genes = gene_list[[i]],
            cat = cellType_col,
            color_pal = color_pal,
            plot_points = plot_points,
            title = n
        ))
    }

    message("Done!")
    grDevices::dev.off()
}
