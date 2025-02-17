## example data
if (!exists("sce_DLPFC_example")) sce_DLPFC_example <- fetch_deconvo_data("sce_DLPFC_example")
marker_stats_1vAll <- suppressMessages(
    findMarkers_1vAll(
        sce = sce_DLPFC_example,
        assay_name = "logcounts",
        cellType_col = "cellType_broad_hc",
        mod = "~BrNum"
    )
)

test_that("1vALL stats returned for each gene", {
    gene_count <- marker_stats_1vAll |> dplyr::count(cellType.target)
    testthat::expect_true(all(gene_count$n == nrow(sce_DLPFC_example)))
})
