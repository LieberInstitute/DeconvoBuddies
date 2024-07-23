#' Get Mean Ratio for Each Gene x Cell Type
#'
#' Calculate the mean ratio value and rank for each gene for each cell type in the `sce`
#' object, to identify effective marker genes for deconvolution.
#'
#' Improved efficiency and ability to handle large data sets from `get_mean_ratio()`.
#'
#' @param sce [SummarizedExperiment-class][SummarizedExperiment::SummarizedExperiment-class] object
#' @param cellType_col A `character(1)` name of the column in the
#' [colData()][SummarizedExperiment::SummarizedExperiment-class] of `sce` that
#' denotes the cell type or group of interest
#' @param assay_name A `character(1)` specifying the name of the
#' [assay()][SummarizedExperiment::SummarizedExperiment-class] in the
#' `sce` object to use to rank expression values. Defaults to `logcounts` since
#' it typically contains the normalized expression values.
#' @param add_symbol a logical indicating whether the gene symbol column to the marker stats table
#'
#' @return Table of mean ratio for each x cell type
#' @export
#'
#'
#' @examples
#' #' ## load example SingleCellExperiment
#' if (!exists("sce_DLPFC_example")) sce_DLPFC_example <- fetch_deconvo_data("sce_DLPFC_example")
#' 
#' ## Get the mean ratio for each gene for each cell type defined in `cellType_broad_hc`
#' get_mean_ratio2(sce_DLPFC_example, cellType_col = "cellType_broad_hc")
#' 
#' #  gene       cellType.target mean.target cellType    mean ratio rank_ratio anno_ratio             
#' #<chr>      <fct>                 <dbl> <fct>      <dbl> <dbl>      <int> <chr>                  
#' #1 CD22       Oligo                  1.36 OPC       0.0730  18.6          1 Oligo/OPC: 18.625      
#' #2 LINC01608  Oligo                  2.39 Micro     0.142   16.8          2 Oligo/Micro: 16.829    
#' #3 FOLH1      Oligo                  1.59 OPC       0.101   15.7          3 Oligo/OPC: 15.684      
#' #4 SLC5A11    Oligo                  2.14 Micro     0.145   14.7          4 Oligo/Micro: 14.697    
#' #5 AC012494.1 Oligo                  2.42 OPC       0.169   14.3          5 Oligo/OPC: 14.282 
#' 
#' @importFrom dplyr mutate
#' @importFrom dplyr arrange
#' @importFrom purrr map
#' @importFrom purrr map2
#' @importFrom matrixStats rowMedians
get_mean_ratio2 <- function(sce, cellType_col = "cellType", assay_name = "logcounts", add_symbol = TRUE) {
    # RCMD fix
    cellType.target <- NULL
    cellType <- NULL
    ratio <- NULL

    cell_types <- unique(sce[[cellType_col]])
    names(cell_types) <- cell_types

    sce_assay <- as.matrix(SummarizedExperiment::assays(sce)[[assay_name]])

    ## Get mean expression for each gene for each cellType
    cell_means <- map(cell_types, ~ as.data.frame(base::rowMeans(sce_assay[, sce[[cellType_col]] == .x])))

    cell_means <- do.call("rbind", cell_means)
    colnames(cell_means) <- "mean"
    ## Define columns
    cell_means$cellType <- rep(cell_types, each = nrow(sce))
    cell_means$gene <- rep(rownames(sce), length(cell_types))
    # print(head(cell_means))

    ## Filter and calculate ratio for each celltype
    ratio_tables <- map(cell_types, ~ .get_ratio_table(.x, sce, sce_assay, cellType_col, cell_means))

    ratio_tables <- do.call("rbind", ratio_tables)

    # max_digits <- nchar(max(ratio_tables$ratio_tables))

    if (add_symbol) {
        ratio_tables$Symbol <- SummarizedExperiment::rowData(sce)[ratio_tables$gene, ]$Symbol
        # ratio_tables <- ratio_tables %>%
        #    mutate(ratio_anno = paste0(stringr::str_pad(rank_ratio, max_digits, "left"),": ",Symbol))
    }
    ratio_tables <- ratio_tables %>%
        mutate(anno_ratio = paste0(cellType.target, "/", cellType, ": ", base::round(ratio, 3)))

    return(ratio_tables)
}


.get_ratio_table <- function(x, sce, sce_assay, cellType_col, cell_means) {
    # RCMD Fix
    mean.target <- NULL
    gene <- NULL
    ratio <- NULL
    cellType.target <- NULL
    cellType <- NULL

    # filter target median != 0
    median_index <- matrixStats::rowMedians(sce_assay[, sce[[cellType_col]] == x]) != 0
    # message("Median == 0: ", sum(!median_index))
    # filter for target means
    target_mean <- cell_means[cell_means$cellType == x, ]
    target_mean <- target_mean[median_index, ]
    colnames(target_mean) <- c("mean.target", "cellType.target", "gene")

    nontarget_mean <- cell_means[cell_means$cellType != x, ]

    ratio_table <- dplyr::left_join(target_mean, nontarget_mean, by = "gene") %>%
        mutate(ratio = mean.target / mean) %>%
        dplyr::group_by(gene) %>%
        arrange(ratio) %>%
        dplyr::slice(1) %>%
        dplyr::select(gene, cellType.target, mean.target, cellType, mean, ratio) %>%
        arrange(-ratio) %>%
        dplyr::ungroup() %>%
        mutate(rank_ratio = dplyr::row_number())

    return(ratio_table)
}
