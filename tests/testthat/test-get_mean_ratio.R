mrt_ab <- suppressMessages(get_mean_ratio(sce_ab, cellType_col = "cellType"))

test_that("Means all zero", {
    expect_equal(sum(mrt_ab$mean.2nd), 0)
})
