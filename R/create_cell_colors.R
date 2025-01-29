#' Create a cell color pallet for plots
#'
#' This function returns a `character()` vector with valid R colors for a given
#' input `character()` of unique cell types. These were colors that have been
#' useful in our experience.
#'
#' @param cell_types A `character()` vector listing unique cell types.
#' @param pallet_name A `character(1)` indicating choice of included pallets:
#' 
#' * `"classic"`: classic set of 8 cell type colors from LIBD, checked for 
#' visability and color blind accessibility. Default pallet. 
#' * `"gg"` : mimic colors automatically picked by ggplot. 
#' * `"tableau"` : 20 distinct colors from tableau color pallet, good for 
#' large number of cell type. 
#' 
#' @param pallet A `character()` vector listing user provided color pallet. If 
#' provided, overrides pallet selection with pallet_name.
#' @param split delineating `character(1)` after which suffixes will be ignored.
#' This is useful for cases when say `A.1` and `A.2` are both to be considered
#' subtypes of cell type `A` (here `split = "\\."`).
#' @param preview A `logical(1)` indicating whether to make a plot to preview
#' the colors.
#'
#' @return A named `character()` vector of R and hex color values compatible
#' with `ggplot2:scale_color_manual()`.
#' @export
#'
#' @examples
#' ## create cell colors with included pallets
#' create_cell_colors(pallet_name = "classic")
#' create_cell_colors(pallet_name = "classic", preview = TRUE)
#' create_cell_colors(pallet_name = "tableau", preview = TRUE)
#' 
#' ## use custom colors
#' my_colors <- c("darkorchid4", "deeppink4", "aquamarine3", "darkolivegreen1")
#' create_cell_colors(cell_type = c("A", "B", "C", "D"), 
#'                    pallet = my_colors, 
#'                    preview = TRUE)
#'                    
#' ## use Rcolor brewer
#' create_cell_colors(cell_type = c("A", "B", "C"), 
#'                    pallet = RColorBrewer::brewer.pal(n = 3, name = "Set1"),
#'                    previe = TRUE)
#' 
#' ## Options for subtype handling
#' ## Provide unique colors for cell subtypes (DEFAULT)
#' create_cell_colors(
#'     cell_types = c("A.1", "A.2", "B.1", "C", "D"),
#'     pallet_name = "classic",
#'     preview = TRUE
#' )
#'
#' ## Provide gradient colors for A.1 and A.2 by using the "split" argument
#' ## returns a base cell type color & subtype colors
#' create_cell_colors(
#'     cell_types = c("A.1", "A.2", "B.1", "C", "D"),
#'     split = "\\.",
#'     pallet_name = "classic",
#'     preview = TRUE
#' )
#' 
#' ## try with custom colors
#' create_cell_colors(
#'     cell_types = c("A.1", "A.2", "B.1", "C", "D"),
#'     split = "\\.",
#'     pallet = my_colors,
#'     preview = TRUE
#' )
#' 
#' @importFrom grDevices colorRampPalette
#' @importFrom rafalib splitit
#' @importFrom purrr map2
#' @importFrom graphics barplot par
#' @importFrom grDevices hcl
#' @importFrom utils head
create_cell_colors <- function(
    cell_types = c("Astro", 
                   "Micro", 
                   "Endo", 
                   "Oligo", 
                   "OPC", 
                   "Excit", 
                   "Inhib", 
                   "Other"),
    pallet_name = c("classic",
                    "gg", 
                    "tableau"),
    pallet = NULL,
    split = NA,
    preview = FALSE) {
  
  ## check number of cell types
  stopifnot(length(cell_types) > 0)
  stopifnot(is.character(cell_types))
  
  base_cell_types <- unique(ss(cell_types, pattern = split))
  nct <- length(base_cell_types)
  # if (nct < 3) stop("Need 3 or more base cell types")
  cell_colors <- list()
  
  ## check pallet selection
  if(is.null(pallet_name) & is.null(pallet)){
    stop("must select a pallet_name or provide custom pallet")
    
  } else if(!is.null(pallet)){ ## use custom pallet
    stopifnot(is.character(pallet))
    cell_colors = pallet
    message(sprintf("Creating custom pallet for %d cell types", nct))
    
  } else { ## use user provided pallet
    pallet_name <- match.arg(pallet_name)
    message(sprintf("Creating %s pallet for %d cell types", pallet_name, nct))
    
    if (pallet_name == "gg") {
      cell_colors <- gg_color_hue(nct)
    } else if (pallet_name == "tableau") {
      cell_colors <- tableau20[seq(nct)]
    } else if (pallet_name == "classic"){
      cell_colors = c("#3BB273",
                      "#FF56AF",
                      "#663894",
                      "#F57A00",
                      "#D2B037",
                      "#247FBC",
                      "#E83E38",
                      "#4E586A")
    }
  }
  
  #### match cell types and colors ####
  if(length(cell_colors) < nct){ ## error if not enough colors
    stop(sprintf("more cell types (%d) than colors in pallet (%d)", nct, length(cell_colors)))
    
  } else if(length(cell_colors) > nct) { ## subset large pallet
    # message(sprintf("more colors (%d) than cell types (%d), using first (%d) colors", length(cell_colors), nct, nct))
    cell_colors <- cell_colors[seq(nct)]
  }
  ## assign cell types to colors
  names(cell_colors) <- base_cell_types
  
  ## handle cell subtype gradient
  if (!identical(base_cell_types, cell_types)) {
    split_cell_types <- cell_types[!cell_types %in% base_cell_types]
    base_split <- rafalib::splitit(ss(split_cell_types, split))
    
    split_scale_colors <- purrr::map2(
      names(base_split), base_split,
      ~ .scale_cell_colors(
        cell_colors[[.x]],
        split_cell_types[.y]
      )
    )
    message(sprintf("Creating subtype gradients for %d base cell types", length(split_scale_colors)))
    split_scale_colors <- unlist(split_scale_colors)
    cell_colors <- c(cell_colors, split_scale_colors)
  }
  
  ## plot preview
  if (preview) {
    par(las = 2) # make label text perpendicular to axis
    par(mar = c(5, 8, 4, 2)) # increase y-axis margin.
    bp <- barplot(rep(1, length(cell_colors)),
            col = cell_colors,
            horiz = TRUE,
            axes = FALSE,
            names.arg = names(cell_colors)
    )
    text(y = bp, x = rep(.5, length(cell_colors)), cell_colors)
  }
  
  return(cell_colors)
}

.scale_cell_colors <- function(color, cell_types) {
  n_ct <- length(cell_types)
  scale_colors <- grDevices::colorRampPalette(c(color, "white"))(n_ct + 1)
  scale_colors <- utils::head(scale_colors, n_ct)
  names(scale_colors) <- cell_types
  
  return(scale_colors)
}
