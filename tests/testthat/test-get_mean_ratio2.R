if (!exists("sce_DLPFC_example")) sce_DLPFC_example <- fetch_deconvo_data("sce_DLPFC_example")

mrt <- get_mean_ratio(sce_DLPFC_example, cellType_col = "cellType_broad_hc")
mrt2 <- .get_mean_ratio2(sce_DLPFC_example, assay_name = "logcounts", cellType_col = "cellType_broad_hc", add_symbol = TRUE)

test_that("Correct Dims", {
    expect_equal(nrow(mrt2), nrow(mrt))
    expect_equal(ncol(mrt2), ncol(mrt))
})


test_that("Ratios are equal", {
    expect_equal(mrt2$ratio, mrt$MeanRatio)
})
