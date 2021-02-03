
#' create plotExpression plots for top n genes for each cell type
#'
#' @param sce sce object that contains logcounts assay
#' @param stats table generated by get_mean_ratio and/or findMarkers_1vAll
#' @param n_genes number of markers you'd like to plot
#' @param rank_col name of column to rank genes by
#' @param anno_col name of column containing annotation
#' @param cellType_col came of colData column containing cell type for sce data
#'
#' @return plotExpression style violin plot for selected marker genes
#' @export
#'
#' @examples
plot_marker_express <- function(sce, stats, cell_type ,n_genes, rank_col, anno_col, cellType_col =  "cellType"){
  title = paste(cell_type, "Top", n_genes ,"markers")
  # message(title)

  max_digits <- nchar(n_genes)

  cell_stats <- stats %>%
    dplyr::rename(rank = rank_col, anno = anno_col) %>%
    dplyr::filter(cellType.target == cell_type,
                  rank <= n_genes) %>%
    mutate(Feature = paste0(stringr::str_pad(rank, max_digits, "left"), ": ", Symbol))

  marker_sce <- sce[cell_stats$gene,]
  rownames(marker_sce) <- cell_stats$Feature

  pe <- scater::plotExpression(marker_sce, exprs_values = "logcounts", features = cell_stats$Feature,
                       x= cellType_col, colour_by= cellType_col, point_alpha=0.5, point_size=.2, ncol=5,
                       add_legend=F) +
    # scale_color_manual(values = cell_colors)+
    ggplot2::stat_summary(fun = mean, fun.min = mean, fun.max = mean,
                 geom = "crossbar", width = 0.3) +
    ggplot2::geom_text(data = cell_stats, ggplot2::aes(x = -Inf, y = Inf, label = anno),
              vjust = "inward", hjust = "inward",size = 2.5)+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::ggtitle(label=title)

  return(pe)
}