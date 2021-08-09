

#' Find the poportion of zero counts pre group for each gene
#'
#' @param sce Single cell experiment object
#' @param groups column names in the colData of sce that define groups of interest
#'
#' @return
#' @export
#'
#' @examples
#' prop_zero(sce_ab)
#' sce_ab$g2 <- paste0(sce_ab$cellType, "_", sce_ab$donor)
#' prop_zero(sce_ab, group_col = "g2")
#' prop_zero(sce.test)
#' @importFrom purrr map_dfc
prop_zero <- function(sce, group_col = "cellType") {
    all_groups <- unique(sce[[group_col]])
    names(all_groups) <- all_groups

    gene_propZero <- purrr::map_dfc(all_groups, function(g) {
        sce_temp <- sce[, sce[[group_col]] == g]
        prop_zero <- rowSums(as.matrix(assays(sce_temp)$counts == 0) / ncol(sce_temp))
        return(prop_zero)
    })

    gene_propZero <- as.data.frame(gene_propZero)
    rownames(gene_propZero) <- rownames(sce)
    return(gene_propZero)
}

#' Find the maximum proportion of zero counts between groups for each gene
#'
#' @param sce Single cell experiment object
#' @param groups column names in the colData of sce that define groups of interest
#'
#' @return
#' @export
#'
#' @examples
#' max_prop_zero(sce = sce_ab)
#' head(max_prop_zero(sce = sce.test))
max_prop_zero <- function(sce, group_col = "cellType") {
    prop_zero <- prop_zero(sce = sce, group_col = group_col)
    max_prop_zero <- apply(prop_zero, 1, max)
    return(max_prop_zero)
}
