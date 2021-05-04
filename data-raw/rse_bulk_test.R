## code to prepare `rse_bulk_test` dataset goes here
set.seed(1234)

nrows <- 2000
ncols <- 100
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)

gene_sym <- paste0("G", stringr::str_pad(1:nrows, nchar(nrows), side = "left", pad = 0))

rowRanges <- GRanges(rep(c("chr1", "chr2"), c(500, nrows - 500)),
    IRanges(floor(runif(nrows, 1e5, 1e6)), width = 100),
    strand = sample(c("+", "-"), nrows, TRUE),
    feature_id = sprintf("ID%03d", 1:nrows),
    Symbol = gene_sym
)

BrNums <- paste0("Br", stringr::str_pad(1:ncols, nchar(ncols), side = "left", pad = 0))
RNums_raw <- sample(1:(ncols * 10), ncols, replace = FALSE)
RNums <- paste0("R", stringr::str_pad(RNums_raw, nchar(max(RNums_raw)), side = "left", pad = 0))

colData <- DataFrame(
    RNum = RNums,
    BrNum = BrNums,
    Sex = sample(c("M", "F"), ncols, replace = TRUE),
    Dx = sample(c("Case", "Control"), ncols, replace = TRUE),
    Age = runif(ncols, min = -1, max = 100)
)

rse_bulk_test <- SummarizedExperiment(
    assays = list(counts = counts),
    rowRanges = rowRanges, colData = colData
)

colnames(rse_bulk_test) <- RNums
rownames(rse_bulk_test) <- gene_sym

usethis::use_data(rse_bulk_test, overwrite = TRUE)
