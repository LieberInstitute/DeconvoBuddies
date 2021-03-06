% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_mean_ratio2.R
\name{get_mean_ratio2}
\alias{get_mean_ratio2}
\title{Get Mean Ratio for Each Gene x Cell Type}
\usage{
get_mean_ratio2(
  sce,
  cellType_col = "cellType",
  assay_name = "logcounts",
  add_symbol = TRUE
)
}
\arguments{
\item{sce}{Single cell experiment object}

\item{cellType_col}{Column name on colData of the sce that denotes the celltype}

\item{assay_name}{Name of assay to use for calculation}

\item{add_symbol}{Add the gene symbol column to the marker stats table}
}
\value{
Table of mean ratio for each x cell type
}
\description{
Get Mean Ratio for Each Gene x Cell Type
}
\examples{
get_mean_ratio2(sce.test)
}
