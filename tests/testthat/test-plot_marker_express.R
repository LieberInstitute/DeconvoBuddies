test_that("Data checks", {
  expect_error(plot_marker_express(sce = sce.test,cellType_col = "not_there", stat = marker_test, cell_type = "Astro"))
  expect_error(plot_marker_express(sce = sce.test, stat = marker_test, cell_type = "not_there"))
})
