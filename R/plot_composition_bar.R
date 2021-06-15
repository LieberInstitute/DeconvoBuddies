
#' Create barplot of average cell type composition
#'
#' @param prop_long data.frame of cell type portions in long form
#' @param x_col category to divide samples by
#' @param prop_col name of column containing proportion values
#' @param ct_col name of column containing cell type names
#' @param add_text Add rounded proption value to bars
#'
#' @return
#' @export
#'
#' @examples
#' 
#' pd <- SummarizedExperiment::colData(rse_bulk_test) %>%
#' as.data.frame()
#' 
#' est_prop_long <- est_prop %>% 
#' tibble::rownames_to_column("RNum") %>%
#' tidyr::pivot_longer(!RNum, names_to = "cell_type", values_to ="prop") %>%
#' left_join(pd %>% select(RNum, Dx))
#' 
#' plot_composition_bar(est_prop_long)
#' plot_composition_bar(est_prop_long, x_col = "Dx")
#' plot_composition_bar(est_prop_long, x_col = "RNum", add_text = FALSE)
#' @importFrom dplyr rename group_by summarise mutate arrange
#' @importFrom ggplot2 ggplot geom_bar geom_text aes theme element_text
plot_composition_bar <- function(prop_long, 
                                 x_col = "ALL", 
                                 prop_col = "prop", 
                                 ct_col = "cell_type",
                                 add_text = TRUE){
  
  # ct_col <- dplyr::enquo(ct_col)
  cell_type <- prop <- mean_prop <- x_cat <- anno_y <- NULL
  
  mean_prop <- prop_long %>%
    dplyr::mutate(ALL = "ALL") %>%
    dplyr::rename(cell_type = ct_col, prop = prop_col, x_cat = x_col) %>%
    dplyr::group_by(cell_type, x_cat) %>%
    dplyr::summarise(mean_prop = mean(prop))  %>%
    dplyr::arrange(cell_type) %>%
    dplyr::group_by(x_cat) %>%
    dplyr::mutate(anno_y = (1 - cumsum(mean_prop)) + (mean_prop *.5))
  
  comp_barplot <- ggplot2::ggplot(data = mean_prop,
                                  ggplot2::aes(x = x_cat, y = mean_prop, fill = cell_type)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::labs(x = x_cat, y = "Mean Proportion", fill ='Cell Type') +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1))
  
  if(add_text){
    comp_barplot <- comp_barplot + ggplot2::geom_text(ggplot2::aes(y = anno_y, label = format(round(mean_prop,3),3)))
  }
  
  return(comp_barplot)
}
