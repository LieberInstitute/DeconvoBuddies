#' Title
#'
#' @param n_cell Number of cells
#' @param n_gene Number of genes
#' @param n_cellType Number of cell types
#' @param n_donor Number of donors
#'
#' @return SingelCellExperiment object with randomly generated counts and colData
#' @export
#'
#' @examples
#' test <- make_test_sce()
#' table(test$cellType, test$donor)
#' @importFrom S4Vectors DataFrame
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment colData
#' @importFrom stats rpois
#' @importFrom stringr str_pad
make_test_sce <- function(n_cell = 100, n_gene = 100, n_cellType = 4, n_donor = 2) {
    counts <- matrix(rpois(n_cell * n_gene, lambda = 10), ncol = n_cell, nrow = n_gene)
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))

    donors <- paste0("D", 1:n_donor)
    cells <- paste0("Cell", 1:n_cellType)

    pd <- S4Vectors::DataFrame(
        donor = sample(donors, n_cell, replace = TRUE),
        cellType = sample(cells, n_cell, replace = TRUE)
    )

    width <- nchar(as.character(n_cell))
    rownames(pd) <- paste0("S", stringr::str_pad(1:n_cell, width = width, pad = "0"))
    SummarizedExperiment::colData(sce) <- pd

    rownames(sce) <- paste0("G", 1:n_gene)

    return(sce)
}
