
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
#'@importFrom dplyr arrange
#'@importFrom purrr map
#'@importFrom purrr map2
#'@importFrom jaffelab ss
get_mean_ratio2 <- function(sce, cellType_col =  "cellType", assay_name = "counts", add_symbol = FALSE){

   cell_types <- unique(sce[[cellType_col]])
   names(cell_types) <- cell_types

   sce_assay <- as.matrix(assays(sce)[[assay_name]])

   ## Get mean expression for each gene for each cellType
   cell_means <- map(cell_types, ~as.data.frame(base::rowMeans(sce_assay[,sce[[cellType_col]] == .x])))

   cell_means <- do.call("rbind", cell_means)
   colnames(cell_means) <- "mean"
   ## Define columns
   cell_means$cellType = rep(cell_types, each = nrow(sce))
   cell_means$gene = rep(rownames(sce),length(cell_types))
   # print(head(cell_means))

   ## Filter and calculate ratio for each celltype
   ratio_tables <- map(cell_types, function(x){
      #filter target median != 0
      median_index <- rowMedians(sce_assay[,sce[[cellType_col]] == x]) != 0
      # message("Median == 0: ", sum(!median_index))
      #filter for target means
      target_mean <- cell_means[cell_means$cellType == x,]
      target_mean <- target_mean[median_index,]
      colnames(target_mean) <- c("mean.target","cellType.target","gene")

      nontarget_mean <- cell_means[cell_means$cellType != x,]


      dplyr::left_join(target_mean, nontarget_mean, by = "gene") %>%
         mutate(ratio = mean.target/mean) %>%
         dplyr::group_by(gene) %>%
         arrange(ratio) %>%
         dplyr::slice(1) %>%
         dplyr::select(gene, cellType.target, mean.target, cellType, mean, ratio) %>%
         arrange(-ratio) %>%
         dplyr::ungroup() %>%
         mutate(rank_ratio = dplyr::row_number())

   })

   ratio_tables <- do.call("rbind", ratio_tables)
   return(ratio_tables)
}
