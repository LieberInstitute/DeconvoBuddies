
#' Calculate 1 vs. All standard fold change for each gene x cell type
#'
#' @param sce single cell experiment object
#' @param assay Name of the assay to use for calculation
#' @param cellType_col Column name on colData of the sce that denotes the celltype
#'
#' @return Table of 1 vs. ALL std log fold change + p-values for each gene x cell type
#' @export
#'
#' @examples
findMarkers_1vAll <- function(sce, assay = "logcounts", cellType_col = "cellType"){
  ct <- levels(sce[[cellType_col]])
  names(ct) <- ct
  ## Traditional t-test with design as in PB'd/limma approach ===
  mod <- with(colData(sce), model.matrix(~ donor))
  mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`

  markers.t.1vAll <- map(ct, function(x){
    sce$contrast <- ifelse(sce[[cellType_col]]==x, 1, 0)
    fm <- findMarkers(sce, groups=sce$contrast,
                assay.type=assay, design=mod, test="t",
                direction="up", pval.type="all", full.stats=T)
    fm <- fm[[2]]$stats.0
    fm.std <- findMarkers(sce, groups=sce$contrast,
                                     assay.type="logcounts", design=mod, test="t",
                                     std.lfc=TRUE,
                                     direction="up", pval.type="all", full.stats=T)
    fm.std <- fm.std[[2]]$stats.0
    colnames(fm.std)[[1]] <- "std.logFC"
    return(cbind(fm, fm.std[,1, drop=FALSE]))
  })

  markers.t.1vAll.table <- do.call("rbind",markers.t.1vAll) %>%
    as.data.frame() %>%
    rownames_to_column("gene")%>%
    mutate(gene = gsub("\\.\\d+","",gene)) %>%
    as_tibble %>%
    add_column(cellType.target = rep(ct, each = nrow(sce))) %>%
    group_by(cellType.target) %>%
    mutate(rank_marker = row_number(),
           anno_logFC = paste0(" std logFC = ", round(std.logFC,3)))

  markers.t.1vAll.table$Symbol <- rowData(sce)[markers.t.1vAll.table$gene,]$Symbol

  markers.t.1vAll.table <- markers.t.1vAll.table %>%
  mutate(feature_marker = paste0(str_pad(rank_marker, 4, "left"),": ",Symbol))

  return(markers.t.1vAll.table)
}
