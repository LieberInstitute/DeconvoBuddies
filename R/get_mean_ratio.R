#' Get Mean Ratio for Each Gene x Cell Type
#'
#' Calculate the Mean Ratio value and rank for each gene for each cell type in
#' the `sce` object, to identify effective marker genes for deconvolution.
#'
#' Note if a cell type has < 10 cells the MeanRatio results may be unstable.
#' See rational in OSCA:
#' <https://bioconductor.org/books/3.19/OSCA.multisample/multi-sample-comparisons.html#performing-the-de-analysis>.
#'
#' @param sce [SummarizedExperiment-class][SummarizedExperiment::SummarizedExperiment-class]
#' (or any derivative class) object containing single cell/nucleus gene
#' expression data.
#' @param cellType_col A `character(1)` name of the column in the
#' [colData()][SummarizedExperiment::SummarizedExperiment-class] of `sce` that
#' denotes the cell type or group of interest.
#' @param assay_name A `character(1)` specifying the name of the
#' [assay()][SummarizedExperiment::SummarizedExperiment-class] in the
#' `sce` object to use to rank expression values. Defaults to `logcounts` since
#' it typically contains the normalized expression values.
#' @param gene_ensembl A `character(1)` specifying the `rowData(sce_pseudo)`
#' column with the ENSEMBL gene IDs. This will be used by `layer_stat_cor()`.
#' @param gene_name A `character(1)` specifying the `rowData(sce_pseudo)`
#' column with the gene names (symbols).
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying how to
#' potentially parallelize key matrix operations.
#'
#' @return A `tibble::tibble()` with the `MeanRatio` values for each gene x cell
#' type.
#' * `gene` is the name of the gene (from rownames(`sce`)).
#' * `cellType.target` is the cell type we're finding marker genes for.
#' * `mean.target` is the mean expression of `gene` for `cellType.target`.
#' * `cellType.2nd` is the second highest non-target cell type.
#' * `mean.2nd` is the mean expression of `gene` for `cellType.2nd`.
#' * `MeanRatio` is the ratio of `mean.target/mean.2nd`.
#' * `MeanRatio.rank` is the rank of `MeanRatio` for the cell type.
#' * `MeanRatio.anno` is an annotation of the `MeanRatio` calculation helpful
#' for plotting.
#' * `gene_ensembl` & `gene_name` optional columns from `rowData(sce)` specified
#' by the user to add gene information.
#'
#' @export
#'
#'
#' @examples
#' ## load example SingleCellExperiment
#' if (!exists("sce_DLPFC_example")) sce_DLPFC_example <- fetch_deconvo_data("sce_DLPFC_example")
#' ## Explore properties of the sce object
#' sce_DLPFC_example
#'
#' ## this data contains logcounts of gene expression
#' SummarizedExperiment::assays(sce_DLPFC_example)$logcounts[1:5, 1:5]
#'
#' ## nuclei are classified in to cell types
#' table(sce_DLPFC_example$cellType_broad_hc)
#'
#' ## Get the mean ratio for each gene for each cell type defined in
#' ## `cellType_broad_hc`
#' get_mean_ratio(sce_DLPFC_example, cellType_col = "cellType_broad_hc")
#'
#' # Option to specify gene_name as the "Symbol" column from rowData
#' # this will be added to the marker stats output
#' SummarizedExperiment::rowData(sce_DLPFC_example)
#'
#' ## specify rowData col names for gene_name and gene_ensembl
#' get_mean_ratio(sce_DLPFC_example,
#'     cellType_col = "cellType_broad_hc",
#'     gene_name = "gene_name",
#'     gene_ensembl = "gene_id"
#' )
#'
#' @family marker gene functions
#'
#' @import dplyr
#' @importFrom purrr map
#' @importFrom purrr map2
#' @importFrom MatrixGenerics rowMedians
#' @importFrom DelayedMatrixStats rowMedians
#' @importFrom MatrixGenerics rowMeans
#' @importFrom tibble tibble
#' @importFrom BiocParallel bplapply SerialParam
get_mean_ratio <- function(sce,
    cellType_col,
    assay_name = "logcounts",
    gene_ensembl = NULL,
    gene_name = NULL,
    BPPARAM = BiocParallel::SerialParam()) {
    # RCMD Fix
    cellType.target <- NULL
    cellType <- NULL
    ratio <- NULL
    rank_ratio <- NULL
    anno_ratio <- NULL

    ## check inputs are valid
    stopifnot(cellType_col %in% colnames(colData(sce)))
    stopifnot(assay_name %in% names(SummarizedExperiment::assays(sce)))

    cell_types <- unique(sce[[cellType_col]])
    names(cell_types) <- cell_types

    ct_table <- table(sce[[cellType_col]])

    if (any(ct_table < 10)) warning("One or more cell types has < 10 cells, this may result in unstable marker genes results. Check details of get_mean_ratio() for more info")

    sce_assay <- SummarizedExperiment::assays(sce)[[assay_name]]

    ## Get mean and median expression for each gene for each cell type
    result_list = BiocParallel::bplapply(
        cell_types,
        function(cell_type, sce_assay, cellType_col) {
            result_list = list()

            assay_piece = sce_assay[, sce[[cellType_col]] == cell_type]

            result_list[['cell_means']] = tibble::tibble(
                mean = unname(MatrixGenerics::rowMeans(assay_piece)),
                cellType = cell_type,
                gene = rownames(assay_piece)
            )

            result_list[['cell_medians']] = unname(MatrixGenerics::rowMedians(assay_piece)) != 0

            return(result_list)
        },
        BPPARAM = BPPARAM,
        sce_assay = sce_assay,
        cellType_col = cellType_col
    )

    cell_means = dplyr::bind_rows(purrr::map(result_list, ~ .x[['cell_means']]))

    ## Filter and calculate ratio for each celltype. This is not parallelized as the size of each
    ## table is only dependent on the number of genes, which should make the tables small (and we
    ## avoid overhead here)
    ratio_tables <- purrr::map(
        cell_types,
        ~ .get_ratio_table(.x, cell_means, result_list[[as.character(.x)]][['cell_medians']])
    )

    ratio_tables <- dplyr::bind_rows(ratio_tables) |>
        dplyr::mutate(anno_ratio = paste0(cellType.target, "/", cellType, ": ", base::round(ratio, 3))) |>
        dplyr::rename(
            cellType.2nd = cellType,
            mean.2nd = mean,
            MeanRatio = ratio,
            MeanRatio.rank = rank_ratio,
            MeanRatio.anno = anno_ratio
        )

    ## Add gene ensemble and gene_name if specified
    if (!is.null(gene_ensembl)) {
        if (gene_ensembl %in% colnames(SummarizedExperiment::rowData(sce))) {
            ratio_tables$gene_ensembl <- SummarizedExperiment::rowData(sce)[ratio_tables$gene, ][[gene_ensembl]]
        } else {
            warning("'", gene_ensembl, "' not in col rowData, gene_ensembl not included in output")
        }
    }

    if (!is.null(gene_name)) {
        if (gene_name %in% colnames(SummarizedExperiment::rowData(sce))) {
            ratio_tables$gene_name <- SummarizedExperiment::rowData(sce)[ratio_tables$gene, ][[gene_name]]
        } else {
            warning("'", gene_name, "' not in col rowData, gene_name not included in output")
        }
    }

    return(ratio_tables)
}


.get_ratio_table <- function(x, cell_means, cell_medians) {
    # RCMD Fix
    mean.target <- NULL
    gene <- NULL
    ratio <- NULL
    cellType.target <- NULL
    cellType <- NULL

    # filter target median != 0
    # filter for target means
    target_mean <- cell_means[cell_means$cellType == x, ]
    target_mean <- target_mean[cell_medians, ]
    colnames(target_mean) <- c("mean.target", "cellType.target", "gene")

    nontarget_mean <- cell_means[cell_means$cellType != x, ]

    ratio_table <- dplyr::left_join(target_mean, nontarget_mean, by = "gene") |>
        dplyr::mutate(ratio = mean.target / mean) |>
        dplyr::group_by(gene) |>
        dplyr::arrange(ratio) |>
        dplyr::slice(1) |>
        dplyr::select(gene, cellType.target, mean.target, cellType, mean, ratio) |>
        dplyr::arrange(-ratio) |>
        dplyr::ungroup() |>
        dplyr::mutate(rank_ratio = dplyr::row_number())

    return(ratio_table)
}
