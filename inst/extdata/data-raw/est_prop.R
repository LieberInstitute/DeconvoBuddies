## code to prepare `est_prop` dataset goes here
set.seed(567)
data("rse_bulk_test")
n_samples <- ncol(rse_bulk_test)
n_cols <- 5
cell_types <- paste0("cell_", LETTERS[seq(n_cols)])

est_prop <- matrix(runif(n = n_samples * n_cols), n_samples)
est_prop <- t(t(est_prop) * (2 * n_cols:1))
est_prop <- est_prop / rowSums(est_prop)

rownames(est_prop) <- colnames(rse_bulk_test)
colnames(est_prop) <- cell_types

est_prop <- as.data.frame(est_prop)

setequal(rownames(est_prop), colnames(rse_bulk_test))

usethis::use_data(est_prop, overwrite = TRUE)
