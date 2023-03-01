

#' Plot the gene expression of a list of genes in a SCE object 
#'
#' @param sce [SummarizedExperiment-class][SummarizedExperiment::SummarizedExperiment-class] object
#' @param genes  A `list()` of `character(1)` specifying the names of genes to plot, 
#' @param assay A `character(1)` specifying the name of the
#' [assay()][SummarizedExperiment::SummarizedExperiment-class] in the
#' `sce` object to use to rank expression values. Defaults to `logcounts` since
#' it typically contains the normalized expression values.
#' 
#' @param cat  A `character(1)` specifying the name of the categorical variable 
#' to group the cells or nuclei by. Defaults to `cellType`.
#' @param color_pal  A named `character(1)` vector that contains a color pallet matching the `cat` values.
#' @param title A `character(1)` to title the plot 
#' @param points A logical indicating whether to plot points over the violin,
#' defaults to `FALSE` as these often become overplotted and quite large (especially when saved as PDF)
#'
#' @return A `ggplot()` violin plot for selected genes
#' @export
#'
#' @examples
#' ## Using Symbol as rownames makes this more human readable
#' rownames(sce.test) <- rowData(sce.test)$Symbol
#' plot_gene_express(sce = sce.test, genes = c("F3"))
#' plot_gene_express(sce = sce.test, genes = c("RNF220","CSF3R"))
#' plot_gene_express(sce = sce.test, genes = c("RNF220","CSF3R"), points = TRUE)
#' plot_gene_express(sce = sce.test, assay = "counts", genes = c("RNF220","CSF3R"))
#' plot_gene_express(sce = sce.test, assay = "counts", genes = c("RNF220","CSF3R"), title = "Inhib Markers")
#' 
plot_gene_express <- function(sce, genes, assay = "logcounts", cat = "cellType", color_pal = NULL, title = NULL, points = FALSE){
  
  stopifnot(any(genes %in% rownames(sce)))
  stopifnot(cat %in% colnames(colData(sce)))
  stopifnot(assay %in% SummarizedExperiment::assayNames(sce))
  
  cat_df <- as.data.frame(colData(sce))[,cat, drop = FALSE]
  expression_long <- reshape2::melt(as.matrix(assays(sce)[[assay]][genes,,drop=FALSE])) 
  
  cat <- cat_df[expression_long$Var2,]
  expression_long <- cbind(expression_long, cat)
  
  expression_violin <- ggplot(data = expression_long, aes(x = cat, y = value)) +
    # ggplot2::geom_violin(aes(fill = cat), scale = "width") +
    ggplot2::facet_wrap(~Var1, ncol = 2)+
    ggplot2::labs(y = paste0("Expression (", assay,")"),
                  title = title) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "None",axis.title.x=ggplot2::element_blank(),
                   axis.text.x=ggplot2::element_text(angle=90,hjust=1),
                   strip.text.x = ggplot2::element_text(face = "italic")) 
  
  if(points){
    message("plotting points")
    
    expression_violin <- expression_violin +
      ggplot2::geom_violin(aes(color = cat), scale = "width") +
      ggplot2::geom_jitter(aes(color = cat),
                           position = ggplot2::position_jitter(seed = 1, width = 0.2), size = .5) +
      ggplot2::stat_summary(fun = median, 
                            geom = "crossbar", 
                            width = 0.3)
    
    if(!is.null(color_pal)) expression_violin <- expression_violin + scale_color_manual(values = color_pal)
    
    return(expression_violin)
  }
  
  expression_violin <- expression_violin + 
    ggplot2::geom_violin(aes(fill = cat), scale = "width") +
    ggplot2::stat_summary(fun = median, 
                          geom = "crossbar", 
                          width = 0.3)
  
  if(!is.null(color_pal)) expression_violin <- expression_violin + scale_fill_manual(values = color_pal)
  
  # expression_violin
  return(expression_violin)
  
}
