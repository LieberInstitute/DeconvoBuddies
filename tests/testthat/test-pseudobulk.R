test_that("rowSums of sce_ab", {
  pb_ab <- pseudobulk(sce_ab)
  ones <- rep(1, nrow(pb_ab))
  names(ones) <- row.names(pb_ab)
  expect_equal(rowSums(pb_ab), ones)
})
