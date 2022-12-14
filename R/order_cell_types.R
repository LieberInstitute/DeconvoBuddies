
#' Order Cell Types Glia, then Neurons, then Ambig/drop
#'
#' @param cell_types 
#'
#' @return
#' @export
#'
#' @examples
order_cell_types <- function(cell_types) {
  neun_mask <- grepl("Excit|Inhib", cell_types)
  drop_mask <- grepl("Ambiguous|drop", cell_types)
  
  neun_ct <- cell_types[neun_mask]
  glia_ct <- cell_types[!(neun_mask|drop_mask)]
  drop_ct <- cell_types[drop_mask]
  
  ordered_cell_types <- c(sort(glia_ct), sort(neun_ct), sort(drop_ct))
  return(ordered_cell_types)
}
