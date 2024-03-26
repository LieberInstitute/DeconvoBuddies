## code to prepare `RNAScope_prop` dataset goes here

## load proportion data from Human_DLPFC_Deconvolution/processed-data/03_HALO/08_explore_proportions/HALO_cell_type_proportions.csv
RNAScope_prop <- read.csv("data-raw/HALO_cell_type_proportions.csv")
dim(RNAScope_prop)

RNAScope_prop |> dplyr::count(Confidence)

RNAScope_prop <- RNAScope_prop |> dplyr::filter(Confidence != "Low")

usethis::use_data(RNAScope_prop, overwrite = TRUE)
