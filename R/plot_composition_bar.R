#' Create barplot of average cell type composition
#'
#' Given a long formatted `data.frame`, this function creates a barplot for
#' the average cell type composition among a set of samples (donors) using
#' `ggplot2`.
#'
#' @param prop_long A `data.frame` of cell type portions in long form
#' @param sample_col A `character(1)` specifying the name of column in
#' `prop_long` that identifies samples.
#' @param x_col A `character(1)` specifying the name of column in
#' `prop_long` that specifies the category to divide samples by.
#' @param prop_col A `character(1)` specifying the name of column in
#' `prop_long` that contains proportion values.
#' @param ct_col A `character(1)` specifying the name of column in
#' `prop_long` containing cell type names.
#' @param add_text A `logical(1)` determining whether to add the rounded
#' proportion value to the bars.
#' @param min_prop_text A `numeric(1)` specifying the minimum proportion to
#' display text. Values greater than (>) `min_prop_text` will be displayed.
#'
#' @return A stacked barplot `ggplot2` object representing the mean proportion
#' of cell types for each group.
#' @export
#'
#' @examples
#' # Load example data
#' data("rse_bulk_test")
#' data("est_prop")
#'
#' # extract relevant colData from the example RangedSummarizedExperiment object
#' pd <- SummarizedExperiment::colData(rse_bulk_test) |>
#'     as.data.frame()
#'
#' # combine with the example estimated proportions in a long style table
#' est_prop_long <- est_prop |>
#'     tibble::rownames_to_column("RNum") |>
#'     tidyr::pivot_longer(!RNum, names_to = "cell_type", values_to = "prop") |>
#'     dplyr::inner_join(pd |> dplyr::select(RNum, Dx))
#'
#' est_prop_long
#'
#' # Create composition bar plots
#' # Mean composition of all samples
#' plot_composition_bar(est_prop_long)
#'
#' # Mean composition by Dx
#' plot_composition_bar(est_prop_long, x_col = "Dx")
#'
#' # control minimum value of text to add
#' plot_composition_bar(est_prop_long, x_col = "Dx", min_prop_text = 0.1)
#'
#' # plot all samples, then facet by Dx
#' plot_composition_bar(est_prop_long, x_col = "RNum", add_text = FALSE) +
#'     ggplot2::facet_wrap(~Dx, scales = "free_x")
#'
#' @importFrom dplyr rename group_by summarise mutate arrange
#' @importFrom ggplot2 ggplot geom_bar geom_text aes theme element_text
plot_composition_bar <- function(prop_long,
    sample_col = "RNum",
    x_col = "ALL",
    prop_col = "prop",
    ct_col = "cell_type",
    add_text = TRUE,
    min_prop_text = 0) {
    x_cat <- cell_type <- anno_y <- NULL

    # ct_col <- dplyr::enquo(ct_col)
    mean_prop <- .get_cat_prop(prop_long, sample_col, x_col, prop_col, ct_col) |>
        dplyr::ungroup()

    if (x_col == "ALL") x_col <- NULL

    comp_barplot <- ggplot2::ggplot(
        data = mean_prop,
        ggplot2::aes(x = x_cat, y = mean_prop, fill = cell_type)
    ) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::labs(x = x_col, y = "Mean Proportion", fill = "Cell Type") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))

    if (add_text) {
        comp_barplot <- comp_barplot +
            ggplot2::geom_text(
                ggplot2::aes(
                    y = mean_prop,
                    label = ifelse(mean_prop > min_prop_text,
                        format(round(mean_prop, 3), 3),
                        ""
                    )
                ),
                position = ggplot2::position_stack(vjust = 0.5)
            )
    }

    return(comp_barplot)
}


.get_cat_prop <- function(prop_long,
    sample_col = "RNum",
    x_col = "ALL",
    prop_col = "prop",
    ct_col = "cell_type") {
    cell_type <- prop <- mean_prop <- x_cat <- anno_y <- sum_prop <- n <- NULL

    prop_long <- prop_long |>
        dplyr::mutate(
            ALL = "ALL",
            sample = !!as.symbol(sample_col)
        ) |>
        dplyr::rename(cell_type = ct_col, prop = prop_col, x_cat = x_col)

    n_sample <- prop_long |>
        dplyr::group_by(x_cat) |>
        dplyr::summarise(n = length(unique(sample)))

    cat_prop <- prop_long |>
        dplyr::group_by(cell_type, x_cat) |>
        dplyr::mutate(sum_prop = sum(prop)) |>
        dplyr::slice(1) |>
        dplyr::left_join(n_sample, by = "x_cat") |>
        dplyr::mutate(mean_prop = sum_prop / n) |>
        dplyr::arrange(cell_type) |>
        dplyr::group_by(x_cat)

    return(cat_prop)
}
