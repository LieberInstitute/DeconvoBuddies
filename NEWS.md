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
