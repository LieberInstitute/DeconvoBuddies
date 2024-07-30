## build tar ball with Build > Build source Package

setwd("..")
BiocCheck::BiocCheck(
  "DeconvoBuddies_0.99.0.tar.gz",
  `quit-with-status` = FALSE,
  `no-check-R-ver` = TRUE,
  `no-check-bioc-help` = TRUE
)

## check notes in log
