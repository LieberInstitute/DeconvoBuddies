## code to prepare `markers_test` dataset goes here

marker_test <- get_mean_ratio2(sce.test)

usethis::use_data(marker_test, overwrite = TRUE)
