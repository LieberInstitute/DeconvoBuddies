mrt_ab <- suppressMessages(get_mean_ratio(sce_ab, cellType_col = "cellType"))

test_that("Means all zero", {
    expect_equal(sum(mrt_ab$mean.2nd), 0)
})

if (!exists("sce_DLPFC_example")) sce_DLPFC_example <- fetch_deconvo_data("sce_DLPFC_example")

test_that("Warn for <10 cells",{
  expect_warning(get_mean_ratio(sce_DLPFC_example, cellType_col = "cellType_hc"))
})
