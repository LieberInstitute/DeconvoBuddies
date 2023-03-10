#' Plot the nested-listed genes as a multi-page pdf
#'
#' This function plots a nested list of genes as a multipage PDF, one for each sub list. 
#' A use case is plotting known marker genes for multiple cell types over cell type
#' clusters with unknown identities. 
#'
#' @param sce [SummarizedExperiment-class][SummarizedExperiment::SummarizedExperiment-class] object
#' @param gene_list A named list contating the names of genes to plot
#' @param pdf_fn A `character()` of the pdf filename to plot to
#' @param cellType_col The `character()` name of colData column containing cell type for sce data
#' @param gene_name_col The `character()` name of rowData matching the gene name from `gene_list`
#' @return A pdf with violin plots for the expression of top marker genes for all cell types
#' @export
#'
#' @examples
#' marker_list <- list(Oligo = c("RNF220","CLIC4"), Astro = c("PRDM16","F3"))
#' 
#' plot_marker_express_List(sce.test, gene_list = marker_list, pdf_fn = "./plots/test_marker_expression_List.pdf")
#' @family expression plotting functions
#' @importFrom ggplot2 ggplot geom_violin geom_text facet_wrap stat_summary
plot_marker_express_List <- function(sce,
                                    stats,
                                    pdf_fn = "marker_expression.pdf",
                                    gene_list,
                                    cellType_col = "cellType",
                                    gene_name_col = "Symbol",
                                    color_pal = NULL,
                                    plot_points = FALSE) {
  
  stopifnot(cellType_col %in% colnames(colData(sce)))

  if(!identical(rownames(sce), rowData(sce)[[gene_name_col]])){
    message("Using ", gene_name_col," as gene names")
    rownames(sce) <- rowData(sce)[[gene_name_col]] 
  }
  
  
  grDevices::pdf(pdf_fn)
  message("Plotting Genes for...")
  for (ml in gene_list) {
    message("\t", ml)
    print(plot_gene_express(sce = sce, genes = ml))
  }
  message("Done!")
  grDevices::dev.off()
}
