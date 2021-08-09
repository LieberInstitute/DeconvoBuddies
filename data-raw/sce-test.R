## code to prepare `sce.test` dataset goes here
load(here::here("data", "sce.test.Rdata"), verbose = TRUE)
print(object.size(sce.test), units = "auto")
dim(sce.test)
table(sce.test$cellType)
## Drop Excit 3 & 4
sce.test <- sce.test[, !sce.test$cellType %in% c("Excit.3", "Excit.4")]
sce.test$cellType <- droplevels(sce.test$cellType)

## Drop Oligo down to ~ 500 cells
oligo_filter <- !sce.test$cellType == "Oligo"
oligo_filter[!oligo_filter] <- sample(c(TRUE, rep(FALSE, 5)), sum(!oligo_filter), replace = TRUE)

sce.test <- sce.test[, oligo_filter]
## select 500 cells
sce.test <- sce.test[, sample(colnames(sce.test), 500)]

levels(sce.test$cellType) <- c("Astro", "Micro", "Oligo", "OPC", "Excit.1", "Excit.2", "Inhib.1", "Inhib.2")

print(object.size(sce.test), units = "auto")
usethis::use_data(sce.test, overwrite = TRUE)
