test_that("multiplication works", {
    mrt2 <- get_mean_ratio2(sce.test, assay_name = "logcounts", cellType_col = "cellType.Broad", add_symbol = TRUE)
})
