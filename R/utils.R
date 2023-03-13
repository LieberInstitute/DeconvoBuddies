#' @importFrom grDevices hcl

ss <- function(x, pattern, slot = 1, ...) {
    vapply(strsplit(x = x, split = pattern, ...), "[", character(1), slot)
}


gg_color_hue <- function(n) {
    hues <- seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[seq(n)]
}
