pb_ab <- pseudobulk(sce_ab, cell_group_cols = c("donor", "cellType"), just_counts = TRUE)

test_that("Correct Dim", {
    expect_equal(nrow(pb_ab), nrow(sce_ab))
    expect_equal(ncol(pb_ab), 4)
})

test_that("rowMedians are Zero", {
    expect_equal(sum(matrixStats::rowMedians(pb_ab)), 0)
})
