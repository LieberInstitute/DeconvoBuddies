
#' Get Mean Ratio for Each Gene x Cell Type
#'
#' @param sce Single cell experiment object
#' @param cellType_col Column name on colData of the sce that denotes the celltype
#' @param assay_name Name of assay to use for calculation
#'
#' @return Table of mean ratio for each x cell type
#' @export
#'
#' @examples
#' set.seed(127)
#' test_sce <- make_test_sce()
#' mrt <- get_mean_ratio2(test_sce)
#'
#'@importFrom dplyr mutate
#'@importFrom purrr map
#'@importFrom jaffelab ss
get_mean_ratio2 <- function(sce, cellType_col =  "cellType", assay_name = "counts", add_symbol = FALSE){

   celltypes <- unique(sce[[cellType_col]])
   names(celltypes) <- celltypes
   print(celltypes)

   sce_assay <- assays(sce)[[assay_name]]

   cell_means <- map(celltypes, ~as.data.frame(Matrix::rowMeans(sce_assay[,sce[[cellType_col]] == .x])))
   # edit cell means table
   cell_means <- do.call("rbind", cell_means)
   colnames(cell_means) <- "mean"
   cell_means$cellType = ss(rownames(cell_means),"\\.")
   cell_means$gene = ss(rownames(cell_means),"\\.",2)


   ratio_tables <- map(celltypes, function(x){
      #filter target median != 0
      median_index <- rowMedians(as.matrix(sce_assay[,sce[[cellType_col]] == x])) != 0
      # message("Median == 0: ", sum(!median_index))
      #filter for target means
      target_mean <- cell_means[cell_means$cellType == x,]
      target_mean <- target_mean[median_index,]
      colnames(target_mean) <- c("target_mean","target_cellType","gene")

      nontarget_mean <- cell_means[cell_means$cellType != x,]


      dplyr::left_join(target_mean, nontarget_mean, by = "gene") %>%
         mutate(ratio = target_mean/mean) %>%
         dplyr::group_by(gene) %>%
         arrange(ratio) %>%
         slice(1) %>%
         select(gene, target_cellType, target_mean, cellType, mean, ratio) %>%
         arrange(-ratio) %>%
         dplyr::ungroup() %>%
         mutate(ratio_rank = dplyr::row_number())

   })

   ratio_tables <- do.call("rbind", ratio_tables)

   return(ratio_tables)
}
