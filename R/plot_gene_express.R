#' Plot the gene expression of a list of genes in a SCE object
#'
#' This function plots the expression of one or more genes as a violin plot,
#' over a user defined category, typically a cell type annotation. The plots are
#' made using `ggplot2`.
#'
#' @param sce A
#' [SummarizedExperiment-class][SummarizedExperiment::SummarizedExperiment-class]
#' object or one inheriting it.
#' @param genes  A `character()` vector specifying the genes to plot,
#' this should match the format of `rownames(sce)`.
#' @param assay_name A `character(1)` specifying the name of the
#' [assay()][SummarizedExperiment::SummarizedExperiment-class] in the
#' `sce` object to use to rank expression values. Defaults to `logcounts` since
#' it typically contains the normalized expression values.
#' @param category  A `character(1)` specifying the name of the categorical
#' variable to group the cells or nuclei by. Defaults to `cellType`.
#' @param color_pal  A named `character(1)` vector that contains a color palette
#' matching the `category` values.
#' @param title A `character(1)` to title the plot.
#' @param plot_points A `logical(1)` indicating whether to plot points over the
#' violin, defaults to `FALSE` as these often become over plotted and quite
#' large (especially when saved as PDF).
#' @param ncol An `integer(1)` specifying the number of columns for the facet in
#' the final plot. Defaults to 2.
#' @param plot_type A `character(1)` specifying whether to plot a 'violin' 
#' (default) or 'boxplot'.
#' @param free_y `logical(1)` indicating whether to use "free" y-axis between 
#' genes (relevant to `facet_wrap`).
#'
#' @return A `ggplot()` violin plot for selected genes.
#' @export
#'
#' @examples
#' ## Using Symbol as rownames makes this more human readable
#' data("sce_ab")
#' plot_gene_express(sce = sce_ab, genes = c("G-D1_A"))
#' plot_gene_express(sce = sce_ab, genes = c("G-D1_A"), plot_type = "boxplot")
#'
#' # Access example data
#' if (!exists("sce_DLPFC_example")) sce_DLPFC_example <- fetch_deconvo_data("sce_DLPFC_example")
#'
#' ## plot expression of two genes
#' plot_gene_express(
#'     sce = sce_DLPFC_example,
#'     category = "cellType_broad_hc",
#'     genes = c("GAD2", "CD22")
#' )
#' 
#' ## plot as boxplot
#' plot_gene_express(
#'     sce = sce_DLPFC_example,
#'     category = "cellType_broad_hc",
#'     genes = c("GAD2", "CD22"),
#'     plot_type = "boxplot"
#' )
#'
#' ## plot points - note this creates large images and is easy to over plot
#' plot_gene_express(
#'     sce = sce_DLPFC_example,
#'     category = "cellType_broad_hc",
#'     genes = c("GAD2", "CD22"), 
#'     plot_points = TRUE
#' )
#' 
#' ## with boxplot
#' plot_gene_express(
#'     sce = sce_DLPFC_example,
#'     category = "cellType_broad_hc",
#'     genes = c("GAD2", "CD22"), 
#'     plot_points = TRUE,
#'     plot_type = "boxplot"
#' )
#' 
#' ## Use free y-axis between genes
#' plot_gene_express(
#'     sce = sce_DLPFC_example,
#'     category = "cellType_broad_hc",
#'     genes = c("GAD2", "CD22"), 
#'     plot_points = TRUE,
#'     plot_type = "boxplot",
#'     free_y = FALSE
#' )
#'
#' ## Add title
#' plot_gene_express(
#'     sce = sce_DLPFC_example,
#'     category = "cellType_broad_hc",
#'     genes = c("GAD2", "CD22"),
#'     title = "My Genes"
#' )
#'
#'## Add color pallet
#'my_cell_colors <- create_cell_colors(cell_types = levels(sce_DLPFC_example$cellType_broad_hc))
#'
#'plot_gene_express(
#'     sce = sce_DLPFC_example,
#'     category = "cellType_broad_hc",
#'     genes = c("GAD2", "CD22"),
#'     color_pal =  my_cell_colors,
#'     plot_type = "boxplot",
#'     plot_points = TRUE
#' )
#'
#' @family expression plotting functions
#'
plot_gene_express <- function(
        sce,
        genes,
        assay_name = "logcounts",
        category = "cellType",
        color_pal = NULL,
        title = NULL,
        plot_points = FALSE,
        ncol = 2,
        plot_type = c("violin", "boxplot"),
        free_y = FALSE) {
  ##check inputs
  stopifnot(any(genes %in% rownames(sce)))
  plot_type <- match.arg(plot_type)
  facet_scales <- ifelse(free_y, "free_y", "fixed")

    if (!category %in% colnames(colData(sce))) {
        stop(
            category,
            "' is not a column name in colData(sce), check that `category` matches this sce",
            call. = FALSE
        )
    }

    stopifnot(assay_name %in% SummarizedExperiment::assayNames(sce))

    value <- median <- NULL
    
    ## reformat data
    category_df <- as.data.frame(colData(sce))[, category, drop = FALSE]
    expression_long <- reshape2::melt(as.matrix(SummarizedExperiment::assays(sce)[[assay_name]][genes, , drop = FALSE]))

    category <- category_df[expression_long$Var2, ]
    expression_long <- cbind(expression_long, category)

    ## create plot
    expression_plot <- ggplot(data = expression_long, aes(x = category, y = value)) +
        ggplot2::facet_wrap(~Var1, ncol = ncol, scales = facet_scales) +
        ggplot2::labs(
            y = paste0("Expression (", assay_name, ")"),
            title = title
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            legend.position = "None",
            axis.title.x = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
            strip.text.x = ggplot2::element_text(face = "italic")
        )

    if (plot_points) {
      if(plot_type == "violin"){
        
        ##violin plot w/ points
        expression_plot <- expression_plot +
          ggplot2::geom_violin(aes(color = category), scale = "width") +
          ggplot2::geom_jitter(aes(color = category),
                               position = ggplot2::position_jitter(seed = 1, width = 0.2), 
                               size = 0.5,
                               alpha = 0.5
          ) +
          ggplot2::stat_summary(
            fun = mean,
            geom = "crossbar",
            width = 0.3
          )
      } else if(plot_type == "boxplot"){
        expression_plot <- expression_plot +
          ggplot2::geom_boxplot(aes(color = category,), outlier.shape = NA) +
          ggplot2::geom_jitter(aes(color = category),
                               position = ggplot2::position_jitter(seed = 1, width = 0.2), 
                               size = .5,
                               alpha = 0.5)
      }
        

        if (!is.null(color_pal)) expression_plot <- expression_plot + ggplot2::scale_color_manual(values = color_pal)

        return(expression_plot)
    }
    
    if(plot_type == "violin"){
    ##violin plot w/o points
    expression_plot <- expression_plot +
        ggplot2::geom_violin(aes(fill = category), scale = "width") +
        ggplot2::stat_summary(
            fun = mean,
            geom = "crossbar",
            width = 0.3
        )
    } else if(plot_type == "boxplot"){
      expression_plot <- expression_plot +
        ggplot2::geom_boxplot(aes(fill = category)) 
        
    }
    if (!is.null(color_pal)) expression_plot <- expression_plot + ggplot2::scale_fill_manual(values = color_pal)

    # expression_plot
    return(expression_plot)
}
