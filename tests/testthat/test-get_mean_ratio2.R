mrt <- suppressMessages(get_mean_ratio(sce.test, assay_name = "logcounts", cellType_col = "cellType", add_symbol = TRUE))
mrt2 <- get_mean_ratio2(sce.test, assay_name = "logcounts", cellType_col = "cellType", add_symbol = TRUE)

test_that("Correct Dims", {
    expect_equal(nrow(mrt2), nrow(mrt))
    expect_equal(ncol(mrt2) + 2, ncol(mrt))
})

r1 <- mrt %>% dplyr::filter(gene == "ENSG00000187147", cellType.target == "Oligo")
r2 <- mrt %>% dplyr::filter(gene == "ENSG00000187147", cellType.target == "Oligo")

test_that("One Ratio Matches", {
    expect_equal(r2$ratio, r1$ratio)
})
