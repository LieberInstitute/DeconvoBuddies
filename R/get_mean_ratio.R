#' Get Mean Ratio for Each Gene x Cell Type
#'
#' Calculate the mean ratio value and rank for each gene for each cell type in the `sce`
#' object, to identify effective marker genes for deconvolution.
#'
#' Improved argument names and documentaion, but same functionalty from `get_mean_ratio2()`.
#'
#' @param sce [SummarizedExperiment-class][SummarizedExperiment::SummarizedExperiment-class] object
#' @param cellType_col A `character(1)` name of the column in the
#' [colData()][SummarizedExperiment::SummarizedExperiment-class] of `sce` that
#' denotes the cell type or group of interest
#' @param assay_name A `character(1)` specifying the name of the
#' [assay()][SummarizedExperiment::SummarizedExperiment-class] in the
#' `sce` object to use to rank expression values. Defaults to `logcounts` since
#' it typically contains the normalized expression values.
#' @param gene_ensembl A `character(1)` specifying the `rowData(sce_pseudo)`
#' column with the ENSEMBL gene IDs. This will be used by `layer_stat_cor()`.
#' @param gene_name A `character(1)` specifying the `rowData(sce_pseudo)`
#' column with the gene names (symbols).
#'
#' @return A `tibble` with the `MeanRatio` values for each gene x cell type. 
#' * `gene` is the name of the gene (from rownames(`sce`)). 
#' * `cellType.target` is the cell type we're finding marker genes for. 
#' * `mean.target` is the mean expression of `gene` for `cellType.target`. 
#' * `cellType.2nd` is the second highest non-target cell type. 
#' * `mean.2nd` is the mean expression of `gene` for `cellType.2nd`.
#' * `MeanRatio` is the ratio of `mean.target/mean.2nd`. 
#' * `MeanRatio.rank` is the rank of `MeanRatio` for the cell type. 
#' * `MeanRatio.anno` is an annotation of the `MeanRatio` calculation helpful for plotting. 
#' * `gene_ensembl` & `gene_name` optional cols spcified by the user to add gene infomation
#'  
#' @export
#'
#'
#' @examples
#' ## Get the mean ratio for each gene for each cell type defined in `cellType`
#' get_mean_ratio(sce.test, cellType_col = "cellType")
#' 
#' #   gene            cellType.target mean.target cellType.2nd mean.2nd MeanRatio MeanRatio.rank MeanRatio.anno        
#' #  <chr>           <fct>                 <dbl> <fct>           <dbl>     <dbl>          <int> <chr>                 
#' # 1 ENSG00000230615 Inhib.2               1.00  Excit.2         0.239      4.20              1 Inhib.2/Excit.2: 4.197
#' # 2 ENSG00000162512 Inhib.2               1.71  Astro           0.512      3.35              2 Inhib.2/Astro: 3.346  
#' # 3 ENSG00000137965 Inhib.2               0.950 Astro           0.413      2.30              3 Inhib.2/Astro: 2.301  
#' # 4 ENSG00000060718 Inhib.2               3.32  Astro           1.47       2.26              4 Inhib.2/Astro: 2.264  
#' # 5 ENSG00000162631 Inhib.2               3.55  Astro           1.62       2.19              5 Inhib.2/Astro: 2.19  
#'
#' # Option to specify gene_name as the "Symbol" column from rowData
#' # this will be added to the marker stats output
#' SummarizedExperiment::rowData(sce.test)
#' get_mean_ratio(sce.test, cellType_col = "cellType", gene_name = "Symbol", gene_ensembl = "gene_id")
#' # A tibble: 1,778 × 10
#' # gene            cellType.target mean.target cellType.2nd mean.2nd MeanRatio MeanRatio.rank MeanRatio.anno   gene_ensembl gene_name
#' # <chr>           <fct>                 <dbl> <fct>           <dbl>     <dbl>          <int> <chr>            <chr>        <chr>    
#' # 1 ENSG00000230615 Inhib.2               1.00  Excit.2         0.239      4.20              1 Inhib.2/Excit.2… ENSG0000023… AL139220…
#' # 2 ENSG00000162512 Inhib.2               1.71  Astro           0.512      3.35              2 Inhib.2/Astro: … ENSG0000016… SDC3     
#' # 3 ENSG00000137965 Inhib.2               0.950 Astro           0.413      2.30              3 Inhib.2/Astro: … ENSG0000013… IFI44    
#' 
#' @importFrom dplyr mutate
#' @importFrom dplyr arrange
#' @importFrom purrr map
#' @importFrom purrr map2
#' @importFrom matrixStats rowMedians
get_mean_ratio <- function(sce, 
                            cellType_col, 
                            assay_name = "logcounts", 
                            gene_ensembl = NULL,
                            gene_name = NULL
                            ) {
    # RCMD fix
    cellType.target <- NULL
    cellType <- NULL
    ratio <- NULL
    
    ## check inputs are valid
    stopifnot(cellType_col %in% colnames(colData(sce)))
    stopifnot(assay_name %in% names(SummarizedExperiment::assays(sce.test)))

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
    ratio_tables <- map(cell_types, ~ .get_ratio_table(.x, 
                                                       sce, 
                                                       sce_assay, 
                                                       cellType_col, 
                                                       cell_means))

    ratio_tables <- do.call("rbind", ratio_tables)  |>
      mutate(anno_ratio = paste0(cellType.target, "/", cellType, ": ", base::round(ratio, 3))) |>
      rename(cellType.2nd = cellType,
             mean.2nd = mean, 
             MeanRatio = ratio, 
             MeanRatio.rank = rank_ratio,
             MeanRatio.anno = anno_ratio)

    ## Add gene ensemble and gene_name if specified
    if (!is.null(gene_ensembl)) {
      if(gene_ensembl %in% colnames(SummarizedExperiment::rowData(sce))){
        ratio_tables$gene_ensembl <- SummarizedExperiment::rowData(sce)[ratio_tables$gene, ][[gene_ensembl]]
      } else {
        warning("'",gene_ensembl,"' not in col rowData, gene_ensembl not included in output" )
      }
    }
    
    if (!is.null(gene_name)) {
      if(gene_name %in% colnames(SummarizedExperiment::rowData(sce))){
        ratio_tables$gene_name <- SummarizedExperiment::rowData(sce)[ratio_tables$gene, ][[gene_name]]
      } else {
        warning("'",gene_name,"' not in col rowData, gene_name not included in output" )
      }
    }
    
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

    ratio_table <- dplyr::left_join(target_mean, nontarget_mean, by = "gene") |>
        mutate(ratio = mean.target / mean) |>
        dplyr::group_by(gene) |>
        arrange(ratio) |>
        dplyr::slice(1) |>
        dplyr::select(gene, cellType.target, mean.target, cellType, mean, ratio) |>
        arrange(-ratio) |>
        dplyr::ungroup() |>
        mutate(rank_ratio = dplyr::row_number())

    return(ratio_table)
}
