#' Get Mean Ratio for Each Gene x Cell Type
#'
#' @param sce Single cell experiment object
#' @param cellType_col Column name on colData of the sce that denotes the celltype
#' @param assay_name Name of the assay to use for calculation
#' @param add_symbol Add the gene symbol column to the marker stats table
#' @return Table of mean ratio for each x cell type
#' @export
#'
#' @examples
#' markers_ratio <- get_mean_ratio(sce.test)
#' head(markers_ratio)
#' @importFrom stats median
#' @importFrom dplyr select rename left_join group_by ungroup summarise slice arrange mutate
#' @importFrom tibble rownames_to_column
get_mean_ratio <- function(sce, cellType_col = "cellType", assay_name = "logcounts", add_symbol = FALSE) {
    # RCMD fix
    id <- Var1 <- Var2 <- value <- gene <- cellType <- logcounts <- median <- NULL
    cellType.target <- cellType <- Symbol <- mean.target <- ratio <- rank_ratio <- NULL

    sce_celltypes <- as.data.frame(SummarizedExperiment::colData(sce)) %>%
        dplyr::select(cellType = cellType_col) %>%
        tibble::rownames_to_column(var = "id") %>%
        mutate(id = as.character(id))

    message("nrow cellType: ", nrow(sce_celltypes))

    gene_stat <- as.matrix(SummarizedExperiment::assays(sce)[[assay_name]]) %>%
        reshape2::melt() %>%
        dplyr::rename(gene = Var1, id = Var2, logcounts = value) %>%
        mutate(id = as.character(id)) %>%
        dplyr::left_join(sce_celltypes, by = "id") %>%
        dplyr::group_by(gene, cellType) %>%
        dplyr::summarise(
            median = median(logcounts),
            mean = mean(logcounts)
        ) %>%
        dplyr::ungroup()

    message("build target stat")
    target_stat <- gene_stat %>%
        dplyr::filter(median != 0) %>%
        dplyr::rename(
            cellType.target = cellType,
            mean.target = mean,
            median.target = median
        )

    message("build mean ratio")
    mean_ratio <- target_stat %>%
        dplyr::right_join(gene_stat, by = "gene") %>%
        dplyr::filter(cellType.target != cellType) %>%
        mutate(ratio = mean.target / mean) %>%
        arrange(gene, cellType.target, ratio) %>%
        dplyr::group_by(gene, cellType.target) %>%
        dplyr::slice(1) %>%
        dplyr::group_by(cellType.target) %>%
        arrange(-ratio) %>%
        mutate(
            rank_ratio = dplyr::row_number(),
            anno_ratio = paste0(cellType.target, "/", cellType, " = ", round(ratio, 3))
        )

    if (add_symbol) {
        mean_ratio$Symbol <- SummarizedExperiment::rowData(sce)[mean_ratio$gene, ]$Symbol

        mean_ratio <- mean_ratio %>%
            dplyr::select(-median) %>%
            mutate(feature_ratio = paste0(stringr::str_pad(rank_ratio, 4, "left"), ": ", Symbol))
    }

    return(mean_ratio)
}
