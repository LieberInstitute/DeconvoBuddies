data(sce_ab)
test_plot <- plot_gene_express(sce = sce_ab, genes = c("G-D1_A"))

test_that("Error with Missing Gene", {
    expect_error(plot_gene_express(sce = sce_ab, genes = c("MISSING_GENE")))
})


test_that("Expected y-axis label", {
    expect_identical(test_plot$labels$y, "Expression (logcounts)")
})
