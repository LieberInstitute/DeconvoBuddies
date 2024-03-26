# A script to make the metadata.csv file located in inst/extdata of the package.
# See ?AnnotationHubData::makeAnnotationHubMetadata for a description of the 
# metadata.csv file, expected fields and data types. This 
# AnnotationHubData::makeAnnotationHubMetadata() function can be used to 
# validate the metadata.csv file before submitting the package.

library("here")

# outdir <- ""
pkgname <- "DeconvoBuddies"

meta <- data.frame(
  Title = c(
    "rse_gene"
  ),
  Description = c(
    "RangedSummarizedExperiment with bulk gene RNA expression data of Human DLPFC, generated at the Lieber Institute for Brain Development (LIBD) and available through the DeconvoBuddies Bioconductor package."
  ),BiocVersion = "3.18",
  Genome = "GRCh38",
  SourceType = "GTF",
  SourceUrl = "https://github.com/LieberInstitute/DeconvoBuddies",
  SourceVersion = "Mar 26 2023",
  Species = "Homo sapiens",
  TaxonomyId = "9606",
  Coordinate_1_based = TRUE,
  DataProvider = "Lieber Institute for Brain Development (LIBD)",
  Maintainer = "Louise Huuki-Myers <lahuuki@gmail.com>",
  RDataClass = "RangedSummarizedExperiment",
  DispatchClass = "Rda",
  RDataPath = file.path(
    pkgname,
    c(
      "rse_gene")
  ),
  Tags = "DeconvoBuddies:LIBD:human:deconvolution",
  row.names = NULL,
  stringsAsFactors = FALSE
)

write.csv(
  meta,
  file = here::here("inst", "extdata", "metadata.csv"),
  row.names = FALSE
)

## Check
if (FALSE) {
  AnnotationHubData::makeAnnotationHubMetadata(here::here(), fileName = "metadata.csv")
}
# Error in .checkValidViews(.views) : 
#   [1] Please add either ExperimentHubSoftware or AnnotationHubSoftware to biocViews list in DESCRIPTION.

