if (!exists("sce_DLPFC_example")) sce_DLPFC_example <- fetch_deconvo_data("sce_DLPFC_example")

data("marker_test")

## Example plot
test_plot <- plot_marker_express(
    sce = sce_DLPFC_example,
    stat = marker_test,
    cellType_col = "cellType_broad_hc",
    cell_type = "Astro",
    gene_col = "gene"
)

test_that("Expected labels", {
    expect_identical(test_plot$labels$title, "Astro Top 4 Markers")
    expect_identical(test_plot$labels$y, "Expression (logcounts)")
})

test_that("Error with bad input", {
    expect_error(plot_marker_express(sce = sce_DLPFC_example, cellType_col = "not_there", stat = marker_test, cell_type = "Astro"))
    expect_error(plot_marker_express(sce = sce_DLPFC_example, stat = marker_test, cell_type = "not_there"))
})
