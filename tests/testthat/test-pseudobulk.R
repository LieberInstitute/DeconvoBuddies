pb_ab <- pseudobulk(sce_ab, cell_group_cols = c("donor", "cellType"))

test_that("Correct Dim", {
    expect_equal(nrow(pb_ab), nrow(sce_ab))
    expect_equal(ncol(pb_ab), 4)
})

test_that("rowMedians are Zero", {
    expect_equal(sum(Biobase::rowMedians(pb_ab)), 0)
})
