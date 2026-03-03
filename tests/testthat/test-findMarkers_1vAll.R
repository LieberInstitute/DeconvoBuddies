## example data
if (!exists("sce_DLPFC_example")) sce_DLPFC_example <- fetch_deconvo_data("sce_DLPFC_example")
marker_stats_no_raw <- suppressMessages(
    findMarkers_1vAll(
        sce = sce_DLPFC_example,
        assay_name = "logcounts",
        cellType_col = "cellType_broad_hc",
        mod = "~BrNum",
        raw_logFC = FALSE
    )
)
marker_stats_with_raw <- suppressMessages(
    findMarkers_1vAll(
        sce = sce_DLPFC_example,
        assay_name = "logcounts",
        cellType_col = "cellType_broad_hc",
        mod = "~BrNum",
        raw_logFC = TRUE
    )
)

test_that("1vALL stats returned for each gene", {
    gene_count <- marker_stats_no_raw |> dplyr::count(cellType.target)
    testthat::expect_true(all(gene_count$n == nrow(sce_DLPFC_example)))
})

test_that(
    "1vALL stats have expected columns",
    {
        expected_cols_no_raw <- c(
            "gene", "std.logFC", "log.p.value", "log.FDR",
            "cellType.target", "std.logFC.rank", "std.logFC.anno"
        )
        expected_cols_with_raw <- c(
            "gene", "std.logFC", "log.p.value", "log.FDR", "logFC",
            "cellType.target", "std.logFC.rank", "std.logFC.anno"
        )
        expect_equal(colnames(marker_stats_no_raw), expected_cols_no_raw)
        expect_equal(colnames(marker_stats_with_raw), expected_cols_with_raw)
    }
)
