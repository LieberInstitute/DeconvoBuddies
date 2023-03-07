mrt_ab <- suppressMessages(get_mean_ratio(sce_ab))

test_that("Means all zero", {
    expect_equal(sum(mrt_ab$mean), 0)
})
