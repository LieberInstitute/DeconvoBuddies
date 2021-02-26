#' Calculate pw standard fold change for each gene x cell type
#'
#' @param sce single cell experiment object
#' @param assay Name of the assay to use for calculation
#' @param cellType_col Column name on colData of the sce that denotes the celltype
#'
#' @return Table of 1 vs. ALL std log fold change + p-values for each gene x cell type
#' @export
#'
#' @examples
#' markers.amy.t.design <- findMarkers_pw(sce.test)
#' sapply(markers.amy.t.design, function(x){table(x$FDR<0.05)})
#'@importFrom purrr map
#'@importFrom dplyr mutate
findMarkers_pw <- function(sce, assay_name = "counts", cellType_col = "cellType", add_symbol = FALSE){

  cell_types <- unique(sce[[cellType_col]])
  names(cell_types) <- cell_types

  ## Traditional t-test with design as in PB'd/limma approach
  pd <- as.data.frame(SummarizedExperiment::colData(sce))
  # message("donor" %in% colnames(pd))

  mod <- with(pd, stats::model.matrix(~donor))
  mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`

  fm <- scran::findMarkers(sce, groups = sce[[cellType_col]],
                                   assay.type = assay_name, design=mod, test.type ="t",
                                   direction="up", pval.type="all", full.stats=T)

  fm.std <- scran::findMarkers(sce, groups = sce[[cellType_col]],
                               assay.type = assay_name, design=mod, test.type ="t",
                               std.lfc = TRUE,
                               direction = "up", pval.type = "all", full.stats=T)



  # markers.t.1vAll.table <- do.call("rbind",markers.t.1vAll) %>%
  #   as.data.frame() %>%
  #   tibble::rownames_to_column("gene") %>%
  #   mutate(gene = gsub("\\.\\d+","",gene)) %>%
  #   tibble::as_tibble() %>%
  #   tibble::add_column(cellType.target = rep(cell_types, each = nrow(sce))) %>%
  #   dplyr::group_by(cellType.target) %>%
  #   mutate(rank_marker = dplyr::row_number(),
  #          anno_logFC = paste0(" std logFC = ", round(std.logFC,3)))
  #
  # if(add_symbol){
  #   markers.t.1vAll.table$Symbol <- rowData(sce)[markers.t.1vAll.table$gene,]$Symbol
  #   markers.t.1vAll.table <- markers.t.1vAll.table %>%
  #     mutate(feature_marker = paste0(stringr::str_pad(rank_marker, 4, "left"),": ",Symbol))
  # }

  return(fm)
}
