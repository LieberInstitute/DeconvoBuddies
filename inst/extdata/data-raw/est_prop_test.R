## code to prepare `est_prop_test` dataset goes here
set.seed(567)
data("rse_bulk_test")
n_samples <- ncol(rse_bulk_test)
n_cols <- 5
cell_types <- paste0("cell_", LETTERS[seq(n_cols)])

est_prop_test <- matrix(runif(n = n_samples * n_cols), n_samples)
est_prop_test <- t(t(est_prop_test) * (2 * n_cols:1))
est_prop_test <- est_prop_test / rowSums(est_prop_test)

rownames(est_prop_test) <- colnames(rse_bulk_test)
colnames(est_prop_test) <- cell_types

est_prop_test <- as.data.frame(est_prop_test)

setequal(rownames(est_prop_test), colnames(rse_bulk_test))
all(rowSums(est_prop_test) - 1 < 0.001)

usethis::use_data(est_prop_test, overwrite = TRUE)
