# DeconvoBuddies 1.1.4

BUG FIXES

Package `BisqueRNA` is no longer a suggested dependency (no longer available on 
CRAN). 

Vignette "Deconvolution Benchmark in Human DLPFC" now loads pre-computed 
`est_prop` data instead of running Bisque deconvolution.

# DeconvoBuddies 1.1.3

NEW FEATURES

* `plot_gene_express()` now has the option to use "free_y" axis.

# DeconvoBuddies 1.1.2

NEW FEATURES

* `plot_gene_express()` now has the option to change `plot_type` to 'violin' or
'boxplot'

# DeconvoBuddies 1.1.1

BUG FIXES

* `get_mean_ratio()` was initially coercing sparse matrices into in memory
matrices. We resolved this by using `MatrixGenerics::rowMeans()` and
`MatrixGenerics::rowMedians()`. This issue was reported by @cyntsc.

# DeconvoBuddies 0.99.0

NEW FEATURES

* Initial version of `DeconvoBuddies` that introduces the _Mean Ratio_
method for identifying cluster marker genes as implemented in
`get_mean_ratio()`. This method is described in more detail at
<https://doi.org/10.1101/2024.02.09.579665>. This package also provides a
wrapper to `scran::findMarkers()` for identifying marker genes expressed in
one cluster compared to all remaining ones. See `findMarkeres_1vAll()` for more
details. Additionally, `DeconvoBuddies` provides plotting functions for
visualizing gene expression as violin plots across different clusters. These
plots are much smaller in size compared to other ones you can make with 
`scater::plotExpression()`. See `plot_gene_express()`, `plot_marker_express()`,
and related functions.
