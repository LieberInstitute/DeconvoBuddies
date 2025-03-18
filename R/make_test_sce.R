#' Simulate a
#' [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment-class] for
#' testing.
#'
#' The counts are simulated from a poisson distribution with `stats::rpois()`.
#' Use `set.seed()` if you want the results to be reproducible.
#'
#' @param n_cell An `integer(1)` specifying the number of cells.
#' @param n_gene An `integer(1)` specifying the number of genes.
#' @param n_cellType An `integer(1)` specifying the number of cell types.
#' @param n_donor An `integer(1)` specifying the number of donors.
#'
#' @return A
#' [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment-class]
#' object with randomly generated counts and `colData()`.
#' @export
#'
#' @examples
#' ## Create an example sce using default values.
#' set.seed(20240823)
#' test <- make_test_sce()
#'
#' ## Let's check the number of cells per cell type from each donor
#' addmargins(table(test$cellType, test$donor))
#' @importFrom S4Vectors DataFrame
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment colData
#' @importFrom stats rpois
#' @importFrom stringr str_pad
make_test_sce <- function(n_cell = 100,
    n_gene = 100,
    n_cellType = 4,
    n_donor = 2) {
    counts <- matrix(rpois(n_cell * n_gene, lambda = 10),
        ncol = n_cell,
        nrow = n_gene
    )
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))

    donors <- paste0("D", seq(n_donor))
    cells <- paste0("Cell", seq(n_cellType))

    pd <- S4Vectors::DataFrame(
        donor = sample(donors, n_cell, replace = TRUE),
        cellType = sample(cells, n_cell, replace = TRUE)
    )

    width <- nchar(as.character(n_cell))
    rownames(pd) <- paste0(
        "S",
        stringr::str_pad(seq_len(n_cell), width = width, pad = "0")
    )
    SummarizedExperiment::colData(sce) <- pd

    rownames(sce) <- paste0("G", seq(n_gene))

    return(sce)
}
