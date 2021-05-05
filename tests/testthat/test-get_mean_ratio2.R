mrt2 <- get_mean_ratio2(sce_ab, assay_name = "logcounts", cellType_col = "cellType", add_symbol = TRUE)

test_that("Correct Dims", {
  expect_equal(nrow(mrt2), 4)
})
