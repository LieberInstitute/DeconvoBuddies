## load Example data
if (!exists("sce_DLPFC_example")) sce_DLPFC_example <- fetch_deconvo_data("sce_DLPFC_example")
data("marker_test")

test_cell_types <- levels(marker_test$cellType.target)

test_plot_list <- plot_marker_express_ALL(
  sce_DLPFC_example,
  cellType_col = "cellType_broad_hc",
  stat = marker_test
)

test_that("Return list of plots",{
  expect_equal(length(test_plot_list), length(test_cell_types))
  expect_equal(unlist(lapply(lapply(test_plot_list,"[[","labels"),"[[","title")),
               paste(test_cell_types, "Top 10 Markers"))
})
