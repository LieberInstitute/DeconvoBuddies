#' Create a cell color pallet for plots
#'
#' This function returns a `character()` vector with valid R colors for a given
#' input `character()` of unique cell types. These were colors that have been
#' useful in our experience.
#'
#' @param cell_types A `character()` vector listing unique cell types.
#' @param pallet Choice of base pallet `"classic"`, `"gg"`, or `"tableau"`.
#' @param split delineating `character(1)` after which suffixes will be ignored.
#' This is useful for cases when say `A.1` and `A.2` are both to be considered
#' as cell type `A` (here `split = "\\."`).
#' @param preview A `logical(1)` indicating whether to make a plot to preview
#' the colors.
#'
#' @return A named `character()` vector of R and hex color values compatible
#' with `ggplot2:scale_color_manual()`.
#' @export
#'
#' @examples
#' create_cell_colors(pallet = "classic")
#' create_cell_colors(pallet = "classic", preview = TRUE)
#' create_cell_colors(pallet = "tableau", preview = TRUE)
#'
#' ## Consider A.1 and A.2 as two different cell types (default)
#' create_cell_colors(
#'     cell_types = c("A.1", "A.2", "B.1", "C", "D"),
#'     pallet = "gg",
#'     preview = TRUE
#' )
#'
#' ## Consider A.1 and A.2 as cell type A by using the "split" argument
#' create_cell_colors(
#'     cell_types = c("A.1", "A.2", "B.1", "C"),
#'     split = "\\.",
#'     pallet = "gg",
#'     preview = TRUE
#' )
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom rafalib splitit
#' @importFrom purrr map2
#' @importFrom graphics barplot par
#' @importFrom grDevices hcl
#' @importFrom utils head
create_cell_colors <- function(cell_types = c("Astro", "Micro", "Oligo", "OPC", "Inhib", "Excit"),
    pallet = c("classic", "gg", "tableau"),
    split = NA,
    preview = FALSE) {
    pallet <- match.arg(pallet)

    base_cell_types <- unique(ss(cell_types, pattern = split))
    nct <- length(base_cell_types)
    if (nct < 3) stop("Need 3 or more base cell types")

    cell_colors <- list()

    if (pallet == "classic") {
        cell_colors <- RColorBrewer::brewer.pal(n = nct, name = "Set1")
        cell_colors <- c(seq(3, nct), cell_colors[c(1, 2)])
    } else if (pallet == "gg") {
        cell_colors <- gg_color_hue(nct)
    } else if (pallet == "tableau") {
        cell_colors <- tableau20[seq(nct)]
    }

    names(cell_colors) <- base_cell_types

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

        split_scale_colors <- unlist(split_scale_colors)
        cell_colors <- c(cell_colors, split_scale_colors)
    }

    if (preview) {
        par(las = 2) # make label text perpendicular to axis
        par(mar = c(5, 8, 4, 2)) # increase y-axis margin.
        barplot(rep(1, length(cell_colors)),
            col = cell_colors,
            horiz = TRUE,
            axes = FALSE,
            names.arg = names(cell_colors)
        )
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
