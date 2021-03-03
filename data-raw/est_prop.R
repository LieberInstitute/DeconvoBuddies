## code to prepare `est_prop` dataset goes here

n_rows <- 100
n_cols <- 5
samples <- 1:n_rows %>% stringr::str_pad(nchar(n_rows), side = "left", pad = 0) %>% paste0("R",.)
cell_types <- paste0("cell_",LETTERS[1:n_cols])

est_prop <- runif( n = n_rows * n_cols) %>% matrix(n_rows)
est_prop <- t(t(est_prop) * (2*n_cols:1))
est_prop <- est_prop/rowSums(est_prop)

rownames(est_prop) <- samples
colnames(est_prop) <- cell_types

est_prop <- as.data.frame(est_prop)

usethis::use_data(est_prop, overwrite = TRUE)
