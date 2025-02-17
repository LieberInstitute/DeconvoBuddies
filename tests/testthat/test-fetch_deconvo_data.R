test_that("ERROR for bad data values", {
    expect_error(fetch_deconvo_data(type = "BADVALUE"))
})

test_that("Download rse_gene", {
    suppressMessages(rse_gene <- fetch_deconvo_data("rse_gene"))
    ## expect to download S4 Summarized Experiment w/ 110 samples
    expect_s4_class(rse_gene, "RangedSummarizedExperiment")
    expect_equal(ncol(rse_gene), 110)
})

test_that("Download sce_DLPFC_example", {
    suppressMessages(sce_DLPFC_example <- fetch_deconvo_data("sce_DLPFC_example"))
    ## expect to download S4 SingleCellExperiment w/ 10k nuclei
    expect_s4_class(sce_DLPFC_example, "SingleCellExperiment")
    expect_equal(ncol(sce_DLPFC_example), 10000)
})
