#' Calculate 1 vs. All standard fold change for each gene x cell type, wrapper
#' function for `scran::findMarkers()`.
#'
#' For each cell type, this function computes the statistics comparing that cell
#' type (the "1") against all other cell types combined ("All").
#'
#' See <https://github.com/MarioniLab/scran/issues/57> for a more in depth
#' discussion about the standard log fold change statistics provided by
#' `scran::findMarkers()`.
#'
#' See also <https://youtu.be/IaclszgZb-g> for a LIBD
#' rstats club presentation on "Finding and interpreting marker genes in
#' sc/snRNA-seq data". The companion notes are available at
#' <https://docs.google.com/document/d/1BeMtKgE7gpmNywInndVC9o_ufopn-U2EZHB32bO7ObM/edit?usp=sharing>.
#'
#' @param sce A
#' [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment-class]
#' object.
#' @param assay_name Name of the assay to use for calculation. See
#' see `assayNames(sce)` for possible values.
#' @param cellType_col Column name on `colData(sce)` that denotes the cell type.
#' @param add_symbol A `logical(1)` indicating whether to add the gene symbol
#' column to the marker stats table.
#' @param mod A `character(1)` string specifying the model used as design in
#' `scran::findMarkers()`. It can be `NULL` (default) if there are no blocking
#' terms with uninteresting
#' factors as documented at [pairwiseTTests][scran::pairwiseTTests].
#' @param verbose A `logical(1)` choosing whether to print progress messages or
#' not.
#' @param direction A `character(1)` for the choice of direction tested for
#' gene cell type markers: `"up"` (default), `"any"`, or `"down"`. Impacts
#' p-values: if `"up"` genes with logFC < 0 will have `p.value = 1`.
#'
#' @return A `tibble::tibble()` of 1 vs. ALL standard log fold change + p-values
#' for each gene x cell type.
#' -   `gene` is the name of the gene (from `rownames(sce)`).
#' -   `logFC` the log fold change from the DE test
#' -   `log.p.value` the log of the p-value of the DE test
#' -   `log.FDR` the log of the False Discovery Rate adjusted p.value
#' -   `std.logFC` the standard logFC.
#' -   `cellType.target` the cell type we're finding marker genes for
#' -   `std.logFC.rank` the rank of `std.logFC` for each cell type
#' -   `std.logFC.anno` is an annotation of the `std.logFC` value
#'    helpful for plotting.
#' @export
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
#' ## Get the 1vALL stats for each gene for each cell type defined in
#' ## `cellType_broad_hc`
#' marker_stats_1vAll <- findMarkers_1vAll(
#'     sce = sce_DLPFC_example,
#'     assay_name = "logcounts",
#'     cellType_col = "cellType_broad_hc",
#'     mod = "~BrNum"
#' )
#'
#' ## explore output, top markers have high logFC
#' head(marker_stats_1vAll)
#'
#' @importFrom purrr map
#' @importFrom dplyr mutate
#' @importFrom scran findMarkers
#' @importFrom tibble rownames_to_column as_tibble add_column
findMarkers_1vAll <- function(sce,
    assay_name = "counts",
    cellType_col = "cellType",
    add_symbol = FALSE,
    mod = NULL,
    verbose = TRUE,
    direction = "up") {
    # RCMD Fix
    gene <- rank_marker <- cellType.target <- std.logFC <- rowData <- Symbol <- NULL

    cell_types <- unique(sce[[cellType_col]])
    names(cell_types) <- cell_types

    ## Traditional t-test with design as in PB'd/limma approach
    pd <- as.data.frame(SummarizedExperiment::colData(sce))
    # message("donor" %in% colnames(pd))

    if (verbose) message("Running 1vALL Testing for ", direction, "-regulated genes")

    if (!is.null(mod)) {
        mod <- with(pd, stats::model.matrix(as.formula(mod)))
        mod <- mod[, -1, drop = FALSE] # intercept otherwise automatically dropped by `findMarkers()`
    }

    markers.t.1vAll <- map(cell_types, function(x) {
        if (verbose) message(Sys.time(), " - Find markers for: ", x)

        sce$contrast <- ifelse(sce[[cellType_col]] == x, 1, 0)

        fm <- scran::findMarkers(sce,
            groups = sce$contrast,
            assay.type = assay_name, design = mod, test.type = "t",
            direction = direction, pval.type = "all", full.stats = TRUE
        )
        fm <- fm[[2]]$stats.0

        fm.std <- scran::findMarkers(sce,
            groups = sce$contrast,
            assay.type = assay_name, design = mod, test.type = "t",
            std.lfc = TRUE,
            direction = direction, pval.type = "all", full.stats = TRUE
        )
        fm.std <- fm.std[[2]]$stats.0
        colnames(fm.std)[[1]] <- "std.logFC"

        return(cbind(fm, fm.std[, 1, drop = FALSE]))
    })

    if (verbose) message(Sys.time(), " - Building Table")
    markers.t.1vAll.table <- do.call("rbind", markers.t.1vAll) |>
        as.data.frame() |>
        tibble::rownames_to_column("gene") |>
        mutate(gene = gsub("\\.\\d+", "", gene)) |>
        tibble::as_tibble() |>
        tibble::add_column(cellType.target = rep(cell_types, each = nrow(sce))) |>
        dplyr::group_by(cellType.target) |>
        mutate(
            std.logFC.rank = dplyr::row_number(),
            std.logFC.anno = paste0(" std logFC = ", round(std.logFC, 3))
        )

    if (add_symbol) {
        markers.t.1vAll.table$Symbol <- rowData(sce)[markers.t.1vAll.table$gene, ]$Symbol
        markers.t.1vAll.table <- markers.t.1vAll.table |>
            mutate(feature_marker = paste0(stringr::str_pad(rank_marker, 4, "left"), ": ", Symbol))
    }

    if (verbose) message("** Done! **\n")
    return(markers.t.1vAll.table)
}
