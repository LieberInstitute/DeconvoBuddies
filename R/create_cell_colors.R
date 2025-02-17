#' Create a cell color palette for plots
#'
#' This function returns a `character()` vector with valid R colors for a given
#' input `character()` of unique cell types. These were colors that have been
#' useful in our experience.
#'
#' @param cell_types A `character()` vector listing unique cell types.
#' @param palette_name A `character(1)` indicating choice of included palettes:
#' 
#' * `"classic"`: classic set of 8 cell type colors from LIBD, checked for 
#' visability and color blind accessibility. Default palette. 
#' * `"gg"` : mimic colors automatically picked by ggplot. 
#' * `"tableau"` : 20 distinct colors from tableau color palette, good for 
#' large number of cell type. 
#' 
#' @param palette A `character()` vector listing user provided color palette. If 
#' provided, overrides palette selection with palette_name.
#' @param split delineating `character(1)` after which suffixes will be ignored.
#' This is useful for cases when say `A.1` and `A.2` are both to be considered
#' fine subtypes of broad cell type `A` (here `split = "\\."`). When used the 
#' function returns a nested list of borad and fine cell types. 
#' @param preview A `logical(1)` indicating whether to make a plot to preview
#' the colors.
#'
#' @return A named `character()` vector of R and hex color values compatible
#' with `ggplot2:scale_color_manual()`.
#' @export
#'
#' @examples
#' ## create cell colors with included palettes
#' create_cell_colors(palette_name = "classic")
#' create_cell_colors(palette_name = "classic", preview = TRUE)
#' create_cell_colors(palette_name = "tableau", preview = TRUE)
#' 
#' ## use custom colors
#' my_colors <- c("darkorchid4", "deeppink4", "aquamarine3", "darkolivegreen1")
#' create_cell_colors(cell_type = c("A", "B", "C", "D"), 
#'                    palette = my_colors, 
#'                    preview = TRUE)
#'                    
#' ## use Rcolor brewer
#' create_cell_colors(cell_type = c("A", "B", "C"), 
#'                    palette = RColorBrewer::brewer.pal(n = 3, name = "Set1"),
#'                    previe = TRUE)
#' 
#' ## Options for subtype handling
#' ## Provide unique colors for cell subtypes (DEFAULT) - returns one level list
#' create_cell_colors(
#'     cell_types = c("A.1", "A.2", "B.1", "C", "D"),
#'     palette_name = "classic",
#'     preview = FALSE
#' )
#'
#' ## Provide gradient colors for A.1 and A.2 by using the "split" argument
#' ## returns a nested list with broad & fine cell type colors, fine cell types
#' ## are gradient with the top level matching the broad cell type
#' create_cell_colors(
#'     cell_types = c("A.1", "A.2", "B.1", "C", "D"),
#'     split = "\\.",
#'     palette_name = "classic",
#'     preview = TRUE
#' )
#' 
#' ## try with custom colors
#' create_cell_colors(
#'     cell_types = c("A.1", "A.2", "B.1", "C", "D"),
#'     split = "\\.",
#'     palette = my_colors,
#'     preview = TRUE
#' )
#' 
#' @importFrom grDevices colorRampPalette
#' @importFrom rafalib splitit
#' @importFrom purrr map2
#' @importFrom graphics barplot par text
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
    palette_name = c("classic",
                    "gg", 
                    "tableau"),
    palette = NULL,
    split = NA,
    preview = FALSE) {
  
  ## check number of cell types
  stopifnot(length(cell_types) > 0)
  stopifnot(is.character(cell_types))
  
  broad_cell_types <- unique(ss(cell_types, pattern = split))
  nct <- length(broad_cell_types)
  cell_colors <- list()
  
  ## check palette selection
  if(is.null(palette_name) & is.null(palette)){
    stop("must select a palette_name or provide custom palette")
    
  } else if(!is.null(palette)){ ## use custom palette
    stopifnot(is.character(palette))
    cell_colors = palette
    message(sprintf("Creating custom palette for %d broad cell types", nct))
    
  } else { ## use user provided palette
    palette_name <- match.arg(palette_name)
    message(sprintf("Creating %s palette for %d broad cell types", palette_name, nct))
    
    if (palette_name == "gg") {
      cell_colors <- gg_color_hue(nct)
    } else if (palette_name == "tableau") {
      cell_colors <- tableau20[seq(nct)]
    } else if (palette_name == "classic"){
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
    stop(sprintf("more cell types (%d) than colors in palette (%d)", nct, length(cell_colors)))
    
  } else if(length(cell_colors) > nct) { ## subset large palette
    # message(sprintf("more colors (%d) than cell types (%d), using first (%d) colors", length(cell_colors), nct, nct))
    cell_colors <- cell_colors[seq(nct)]
  }
  ## assign cell types to colors
  names(cell_colors) <- broad_cell_types
  
  ## handle cell subtype gradient
  if (!identical(broad_cell_types, cell_types)) {
    split_cell_types <- cell_types[!cell_types %in% broad_cell_types]
    broad_split <- rafalib::splitit(ss(split_cell_types, split))
    
    split_scale_colors <- purrr::map2(
      names(broad_split), broad_split,
      ~ .scale_cell_colors(
        cell_colors[[.x]],
        split_cell_types[.y]
      )
    )
    message(sprintf("Creating fine cell type gradients for %d cell types", length(split_scale_colors)))
    split_scale_colors <- unlist(split_scale_colors)
    cell_colors <- c(cell_colors, split_scale_colors)
    cell_colors <- list(broad = cell_colors[broad_cell_types],
                        fine = cell_colors[cell_types])
    plot_cell_colors <- cell_colors$fine
  } else {
    plot_cell_colors <- cell_colors
  }
  
  ## plot preview
  if (preview) {
    plot_cell_colors <- rev(plot_cell_colors) # flip order
    par(las = 2) # make label text perpendicular to axis
    par(mar = c(5, 8, 4, 2)) # increase y-axis margin.
    bp <- barplot(rep(1, length(plot_cell_colors)),
            col = plot_cell_colors,
            horiz = TRUE,
            axes = FALSE,
            names.arg = names(plot_cell_colors)
    )
    text(y = bp, x = rep(.5, length(cell_colors)), plot_cell_colors)
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
