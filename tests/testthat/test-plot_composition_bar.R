# Load example data
data("rse_bulk_test")
data("est_prop")

pd <- SummarizedExperiment::colData(rse_bulk_test) |>
    as.data.frame()

est_prop_long <- est_prop |>
    tibble::rownames_to_column("RNum") |>
    tidyr::pivot_longer(!RNum, names_to = "cell_type", values_to = "prop") |>
    dplyr::inner_join(pd |> dplyr::select(RNum, Dx), by = dplyr::join_by(RNum))

test_plot <- plot_composition_bar(est_prop_long)
test_plot_dx <- plot_composition_bar(est_prop_long, x_col = "Dx")

test_that("Expected plot labels", {
    ## Overall plot
    expect_identical(test_plot$labels$y, "Mean Proportion")
    expect_identical(test_plot$labels$fill, "Cell Type")
    ## Dx plot
    expect_identical(test_plot_dx$labels$y, "Mean Proportion")
    expect_identical(test_plot_dx$labels$x, "Dx")
})



# test_that("Error with Bad Input",{
#   expect_error(plot_composition_bar(est_prop_long, x_col = "NOT_THERE"))
# }
# )
