
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DeconvoBuddies

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/lahuuki/DeconvoBuddies/branch/main/graph/badge.svg)](https://codecov.io/gh/lahuuki/DeconvoBuddies?branch=main)
<!-- badges: end -->

The goal of `DeconvoBuddies` is to provide helper functions for the
deconvolution process

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `DeconvoBuddies` using
from [Bioconductor](http://bioconductor.org/) the following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("DeconvoBuddies")
```

And the development version from [GitHub](https://github.com/) with:

``` r
BiocManager::install("lahuuki/DeconvoBuddies")
```

## Example

Get mean ratios for each gene x cell type

``` r
library("DeconvoBuddies")
library("dplyr")
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
## basic example code

ratios <- get_mean_ratio2(sce.test)
#> Loading required package: SingleCellExperiment
#> Loading required package: SummarizedExperiment
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: 'matrixStats'
#> The following object is masked from 'package:dplyr':
#> 
#>     count
#> 
#> Attaching package: 'MatrixGenerics'
#> The following objects are masked from 'package:matrixStats':
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#>     rowWeightedSds, rowWeightedVars
#> Loading required package: GenomicRanges
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Warning: package 'BiocGenerics' was built under R version 4.0.5
#> Loading required package: parallel
#> 
#> Attaching package: 'BiocGenerics'
#> The following objects are masked from 'package:parallel':
#> 
#>     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
#>     clusterExport, clusterMap, parApply, parCapply, parLapply,
#>     parLapplyLB, parRapply, parSapply, parSapplyLB
#> The following objects are masked from 'package:dplyr':
#> 
#>     combine, intersect, setdiff, union
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from 'package:base':
#> 
#>     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
#>     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
#>     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
#>     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
#>     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
#>     union, unique, unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: 'S4Vectors'
#> The following objects are masked from 'package:dplyr':
#> 
#>     first, rename
#> The following object is masked from 'package:base':
#> 
#>     expand.grid
#> Loading required package: IRanges
#> 
#> Attaching package: 'IRanges'
#> The following objects are masked from 'package:dplyr':
#> 
#>     collapse, desc, slice
#> Loading required package: GenomeInfoDb
#> Warning: package 'GenomeInfoDb' was built under R version 4.0.5
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: 'Biobase'
#> The following object is masked from 'package:MatrixGenerics':
#> 
#>     rowMedians
#> The following objects are masked from 'package:matrixStats':
#> 
#>     anyMissing, rowMedians
fc <- findMarkers_1vAll(sce.test)

(marker_stats <- left_join(ratios, fc, by = c("gene", "cellType.target")))
#> # A tibble: 2,504 x 15
#>    gene       cellType.target mean.target cellType  mean ratio rank_ratio Symbol
#>    <chr>      <fct>                 <dbl> <fct>    <dbl> <dbl>      <int> <chr> 
#>  1 ENSG00000… Oligo                  5.98 Astro    1.29   4.62          1 RNF220
#>  2 ENSG00000… Oligo                  2.19 OPC      0.737  2.98          2 CLIC4 
#>  3 ENSG00000… Oligo                  2.73 Micro    1.35   2.02          3 ATG4C 
#>  4 ENSG00000… Oligo                  3.36 Excit.3  1.70   1.98          4 SORT1 
#>  5 ENSG00000… Oligo                  2.83 Excit.4  1.50   1.88          5 SLC22…
#>  6 ENSG00000… Oligo                  5.16 Excit.1  2.90   1.78          6 TTLL7 
#>  7 ENSG00000… Oligo                  3.68 OPC      2.26   1.63          7 PTBP2 
#>  8 ENSG00000… Oligo                  1.52 Astro    0.944  1.61          8 PADI2 
#>  9 ENSG00000… Oligo                  1.57 Micro    1.01   1.55          9 SRRM1 
#> 10 ENSG00000… Oligo                  4.11 Excit.2  2.67   1.54         10 DNAJC6
#> # … with 2,494 more rows, and 7 more variables: ratio_anno <chr>, logFC <dbl>,
#> #   log.p.value <dbl>, log.FDR <dbl>, std.logFC <dbl>, rank_marker <int>,
#> #   anno_logFC <chr>
```

<img src="man/figures/README-plot_marker_expression-1.png" width="100%" />

## Citation

Below is the citation output from using `citation('DeconvoBuddies')` in
R. Please run this yourself to check for any updates on how to cite
**DeconvoBuddies**.

``` r
print(citation('DeconvoBuddies'), bibtex = TRUE)
#> 
#> To cite package 'DeconvoBuddies' in publications use:
#> 
#>   Lousie Huuki (2021). DeconvoBuddies: Helper Functions for LIBD
#>   Deconvolution. R package version 0.99.0.
#>   https://github.com/lahuuki/DeconvoBuddies
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {DeconvoBuddies: Helper Functions for LIBD Deconvolution},
#>     author = {Lousie Huuki},
#>     year = {2021},
#>     note = {R package version 0.99.0},
#>     url = {https://github.com/lahuuki/DeconvoBuddies},
#>   }
```

Please note that the `DeconvoBuddies` was only made possible thanks to
many other R and bioinformatics software authors, which are cited either
in the vignettes and/or the paper(s) describing this package.

## Code of Conduct

Please note that the `DeconvoBuddies` project is released with a
[Contributor Code of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

## Development tools

-   Continuous code testing is possible thanks to [GitHub
    actions](https://www.tidyverse.org/blog/2020/04/usethis-1-6-0/)
    through *[usethis](https://CRAN.R-project.org/package=usethis)*,
    *[remotes](https://CRAN.R-project.org/package=remotes)*, and
    *[rcmdcheck](https://CRAN.R-project.org/package=rcmdcheck)*
    customized to use [Bioconductor’s docker
    containers](https://www.bioconductor.org/help/docker/) and
    *[BiocCheck](https://bioconductor.org/packages/3.12/BiocCheck)*.
-   Code coverage assessment is possible thanks to
    [codecov](https://codecov.io/gh) and
    *[covr](https://CRAN.R-project.org/package=covr)*.
-   The [documentation website](http://.github.io/DeconvoBuddies) is
    automatically updated thanks to
    *[pkgdown](https://CRAN.R-project.org/package=pkgdown)*.
-   The code is styled automatically thanks to
    *[styler](https://CRAN.R-project.org/package=styler)*.
-   The documentation is formatted thanks to
    *[devtools](https://CRAN.R-project.org/package=devtools)* and
    *[roxygen2](https://CRAN.R-project.org/package=roxygen2)*.

For more details, check the `dev` directory.

This package was developed using
*[biocthis](https://bioconductor.org/packages/3.12/biocthis)*.
