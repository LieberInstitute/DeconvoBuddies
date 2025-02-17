test_that("make expected test sce", {
    test_sce <- make_test_sce(n_cell = 100, n_gene = 100, n_cellType = 4, n_donor = 2)
    expect_s4_class(test_sce, "SingleCellExperiment")
    expect_equal(ncol(test_sce), 100)
    expect_equal(nrow(test_sce), 100)
})
