#' Plot the gene expression of a list of genes in a SCE object
#'
#' This function plots the expression of one or more genes as a violin plot,
#' over a user defined category, typically a cell type annotation.
#'
#' @param sce [SummarizedExperiment-class][SummarizedExperiment::SummarizedExperiment-class] object
#' @param genes  A `list()` of `character(1)` specifying the genes to plot, this should match the format of `rownames(sce)`
#' @param assay_name A `character(1)` specifying the name of the
#' [assay()][SummarizedExperiment::SummarizedExperiment-class] in the
#' `sce` object to use to rank expression values. Defaults to `logcounts` since
#' it typically contains the normalized expression values.
#'
#' @param cat  A `character(1)` specifying the name of the categorical variable
#' to group the cells or nuclei by. Defaults to `cellType`.
#' @param color_pal  A named `character(1)` vector that contains a color pallet matching the `cat` values.
#' @param title A `character(1)` to title the plot
#' @param plot_points A logical indicating whether to plot points over the violin,
#' defaults to `FALSE` as these often become overplotted and quite large (especially when saved as PDF)
#' @param ncol = Number of columns for the facet in the final plot. Defaults to 2.
#'
#' @return A `ggplot()` violin plot for selected genes
#' @export
#'
#' @examples
#' ## Using Symbol as rownames makes this more human readable
#' #rownames(sce.test) <- SummarizedExperiment::rowData(sce.test)$Symbol
#' plot_gene_express(sce = sce_ab, genes = c("G-D1_A"))
#' #plot_gene_express(sce = sce.test, genes = c("RNF220", "CSF3R"))
#' #plot_gene_express(sce = sce.test, genes = c("RNF220", "CSF3R"), plot_points = TRUE)
#' #plot_gene_express(sce = sce.test, assay_name = "counts", genes = c("RNF220", "CSF3R"))
#' #plot_gene_express(sce = sce.test, assay_name = "counts", genes = c("RNF220", "CSF3R"), title = "Inhib Markers")
#'
#' @family expression plotting functions
#'
plot_gene_express <- function(sce,
    genes,
    assay_name = "logcounts",
    cat = "cellType",
    color_pal = NULL,
    title = NULL,
    plot_points = FALSE,
    ncol = 2) {
    stopifnot(any(genes %in% rownames(sce)))

    if (!cat %in% colnames(colData(sce))) {
        message("ERROR '", cat, "' is not a column name in colData(sce), check that `cat` matches this sce")
        stop()
    }

    stopifnot(assay_name %in% SummarizedExperiment::assayNames(sce))

    value <- median <- NULL

    cat_df <- as.data.frame(colData(sce))[, cat, drop = FALSE]
    expression_long <- reshape2::melt(as.matrix(SummarizedExperiment::assays(sce)[[assay_name]][genes, , drop = FALSE]))

    cat <- cat_df[expression_long$Var2, ]
    expression_long <- cbind(expression_long, cat)

    expression_violin <- ggplot(data = expression_long, aes(x = cat, y = value)) +
        # ggplot2::geom_violin(aes(fill = cat), scale = "width") +
        ggplot2::facet_wrap(~Var1, ncol = ncol) +
        ggplot2::labs(
            y = paste0("Expression (", assay_name, ")"),
            title = title
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            legend.position = "None", axis.title.x = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
            strip.text.x = ggplot2::element_text(face = "italic")
        )

    if (plot_points) {
        expression_violin <- expression_violin +
            ggplot2::geom_violin(aes(color = cat), scale = "width") +
            ggplot2::geom_jitter(aes(color = cat),
                position = ggplot2::position_jitter(seed = 1, width = 0.2), size = .5
            ) +
            ggplot2::stat_summary(
                fun = median,
                geom = "crossbar",
                width = 0.3
            )

        if (!is.null(color_pal)) expression_violin <- expression_violin + ggplot2::scale_color_manual(values = color_pal)

        return(expression_violin)
    }

    expression_violin <- expression_violin +
        ggplot2::geom_violin(aes(fill = cat), scale = "width") +
        ggplot2::stat_summary(
            fun = median,
            geom = "crossbar",
            width = 0.3
        )

    if (!is.null(color_pal)) expression_violin <- expression_violin + ggplot2::scale_fill_manual(values = color_pal)

    # expression_violin
    return(expression_violin)
}
