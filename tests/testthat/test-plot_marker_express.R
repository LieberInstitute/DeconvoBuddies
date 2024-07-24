if (!exists("sce_DLPFC_example")) sce_DLPFC_example <- fetch_deconvo_data("sce_DLPFC_example")

test_that("Data checks", {
    expect_error(plot_marker_express(sce = sce_DLPFC_example, cellType_col = "not_there", stat = marker_test, cell_type = "Astro"))
    expect_error(plot_marker_express(sce = sce_DLPFC_example, stat = marker_test, cell_type = "not_there"))
})
