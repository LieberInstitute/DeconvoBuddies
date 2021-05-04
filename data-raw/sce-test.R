## code to prepare `sce.test` dataset goes here
load(here::here("data", "sce.test.Rdata"), verbose = TRUE)
usethis::use_data(sce.test, overwrite = TRUE)
