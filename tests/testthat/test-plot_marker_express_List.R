## load Example data
if (!exists("sce_DLPFC_example")) sce_DLPFC_example <- fetch_deconvo_data("sce_DLPFC_example")
my_gene_list <- list(Inhib = c("GAD2", "SAMD5"), Astro = c("RGS20", "PRDM16"))

test_plot_list <- plot_marker_express_List(
  sce_DLPFC_example,
  gene_list = my_gene_list,
  cellType_col = "cellType_broad_hc"
)

test_that("Return list plots matches my_gene_list",{
  expect_equal(length(test_plot_list), length(my_gene_list))
  expect_equal(unname(unlist(lapply(lapply(test_plot_list,"[[","labels"),"[[","title"))),
               names(my_gene_list))
})
