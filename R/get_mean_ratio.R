
#' Get Mean Ratio for Each Gene x Cell Type
#'
#' @param sce Single cell experiment object
#' @param cellType_col Column name on colData of the sce that denotes the celltype
#'
#' @return Table of mean ratio for each x cell type
#' @export
#'
#' @examples
get_mean_ratio <- function(sce, cellType_col =  "cellType"){

  sce_celltypes <- as.data.frame(colData(sce)) %>%
    select(cellType = !!sym(cellType_col)) %>%
    rownames_to_column(var = "id") %>%
    mutate(id = as.character(id))

  message("nrow cellType: ", nrow(sce_celltypes))

   gene_stat <- as.matrix(assays(sce)$logcounts) %>%
     melt() %>%
     rename(gene = Var1, id = Var2, logcounts = value) %>%
     mutate(id = as.character(id)) %>%
     left_join(sce_celltypes, by = "id")  %>%
     group_by(gene, cellType) %>%
     summarise(median_logcount = median(logcounts),
               mean_logcount = mean(logcounts)) %>%
     ungroup()

   print(object.size(gene_stat), units = "auto")

   message("build target stat")
   target_stat <- gene_stat %>%
     filter(median_logcount != 0) %>%
     rename(cellType.target = cellType,
            mean_logcount.target = mean_logcount,
            median_logcount.target = median_logcount)

   message("build mean ratio")
   mean_ratio <- target_stat %>% right_join(gene_stat, by = "gene") %>%
     filter(cellType.target != cellType) %>%
     mutate(ratio = mean_logcount.target/mean_logcount) %>%
     arrange(gene, cellType.target,ratio) %>%
     group_by(gene, cellType.target) %>%
     slice(1) %>%
     group_by(cellType.target) %>%
     arrange(-ratio) %>%
     mutate(rank_ratio = row_number(),
            anno_ratio = paste0(cellType.target,"/",cellType," = ",round(ratio, 3)))

   mean_ratio$Symbol <- rowData(sce)[mean_ratio$gene,]$Symbol

   mean_ratio <- mean_ratio %>%
      select(-median_logcount) %>%
      mutate(feature_ratio = paste0(str_pad(rank_ratio, 4, "left"),": ",Symbol))

   return(mean_ratio)
}
