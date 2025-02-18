## code to prepare `markers_test` dataset goes here

marker_test <- get_mean_ratio(sce_DLPFC_example, cellType_col = "cellType_broad_hc", gene_name = "gene_name", gene_ensembl = "gene_id")

usethis::use_data(marker_test, overwrite = TRUE)
