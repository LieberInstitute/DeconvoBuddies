## code to prepare `sce_ab` dataset goes here

pd <- S4Vectors::DataFrame(
    donor = rep(c("D1", "D2"), each = 50),
    cellType = rep(rep(c("A", "B"), each = 25), 2),
    region = rep(rep(c("Front", "Back"), each = 10), 5)
)
rownames(pd) <- paste0("S", stringr::str_pad(1:nrow(pd), width = nchar(nrow(pd)), pad = "0"))

combo <- paste0(pd$donor, "_", pd$cellType)
tc <- table(combo)

counts <- do.call(rbind, map(names(tc), ~ as.integer(.x == combo)))
rownames(counts) <- paste0("G-", names(tc))
sce_ab <- SingleCellExperiment::SingleCellExperiment(list(counts = counts),
    colData = pd
)

usethis::use_data(sce_ab, overwrite = TRUE)
