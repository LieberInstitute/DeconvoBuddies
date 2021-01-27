
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
#'@importFrom dplyr mutate, left_join, group_by, arrange, slice
#'@importFrom purrr map, map2, pluck
#'@importFrom jaffelab ss
get_mean_ratio2 <- function(sce, cellType_col =  "cellType", assay_name = "counts", add_symbol = FALSE){

   celltypes <- unique(sce[[cellType_col]])
   names(celltypes) <- celltypes
   sce_assay <- assays(sce)[[assay_name]]

   cell_means <- map(celltypes, ~as.data.frame(rowMeans(sce_assay[,sce[[cellType_col]] == .x])))
   # edit cell means table
   cell_means <- do.call("rbind", cell_means)
   colnames(cell_means) <- "mean"
   cell_means$cellType = ss(rownames(cell_means),"\\.")
   cell_means$gene = ss(rownames(cell_means),"\\.",2)


   ratio_tables <- map(celltypes, function(x){
      # median_index <- rowMedians(as.matrix(sce_assay[,sce[[cellType_col]] == .x])) != 0
      # temp_cell_means <- cell_means[median_index,]
      target_mean <- cell_means[cell_means$cellType == x,]
      colnames(target_mean) <- c("target_mean","target_cellType","gene")

      nontarget_mean <- cell_means[cell_means$cellType != x,]


      left_join(target_mean, nontarget_mean, by = "gene") %>%
         mutate(ratio = target_mean/mean) %>%
         group_by(gene) %>%
         arrange(ratio) %>%
         slice(1)

   })

   return(ratio_tables)
}
