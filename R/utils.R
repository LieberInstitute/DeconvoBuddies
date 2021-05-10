

ss <- function(x, pattern, slot = 1, ...) {
  sapply(strsplit(x = x, split = pattern, ...), "[", slot)
}
