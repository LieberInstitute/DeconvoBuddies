#' Create a Cell Color Pallet for Deconvolution
#'
#' @param cell_types list of cell types
#' @param pallet Choice of base pallet "classic", "gg", or "tableau"
#' @param split delineating character in cell names
#' @param preview plot preview of colors
#'
#' @return
#' @export
#'
#' @examples
#' create_cell_colors(pallet = "classic")
#' create_cell_colors(pallet = "classic", preview = TRUE)
#' create_cell_colors(pallet = "tableau", preview = TRUE)
#' create_cell_colors(cell_types = c("A.1","A.2", "B.1", "C"), split = "\\.", pallet = "gg",preview = TRUE)
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom rafalib splitit
#' @importFrom purrr map2
create_cell_colors <- function(cell_types = c("Inhib", "Excit","Astro","Micro","Oligo","OPC"),
                               pallet = c("classic", "gg", "tableau"),
                               split = NA,
                               preview = FALSE) {

  pallet <- match.arg(pallet)

  base_cell_types <- unique(ss(cell_types, pattern = split))

  cell_colors <- list()

  if(pallet == "classic"){
    cell_colors <- RColorBrewer::brewer.pal(n = length(base_cell_types), name = "Set1")
  }else if(pallet == "gg"){
    cell_colors <- gg_color_hue(length(base_cell_types))
  }else if(pallet == "tableau"){
    cell_colors <- tableau20[1:length(base_cell_types)]
  }

  names(cell_colors) <- base_cell_types

  if(!identical(base_cell_types, cell_types)){
    split_cell_types <- cell_types[!cell_types %in% base_cell_types]
    base_split <- rafalib::splitit(ss(split_cell_types, split))

    split_scale_colors <- purrr::map2(names(base_split), base_split,
                                     ~.scale_cell_colors(cell_colors[[.x]],
                                                         split_cell_types[.y]))

    split_scale_colors <- unlist(split_scale_colors)
    cell_colors <- c(cell_colors, split_scale_colors)
  }

  if(preview){
    par(las=2) # make label text perpendicular to axis
    par(mar=c(5,8,4,2)) # increase y-axis margin.
    barplot(rep(1, length(cell_colors)),
            col=cell_colors,
            horiz=TRUE,
            axes=FALSE,
            names.arg = names(cell_colors))
  }

  return(cell_colors)

}

.scale_cell_colors <- function(color, cell_types){

  n_ct <- length(cell_types)
  scale_colors <- grDevices::colorRampPalette(c(color, "white"))(n_ct + 1)
  scale_colors <- head(scale_colors, n_ct)
  names(scale_colors)  <- cell_types

  return(scale_colors)
}
