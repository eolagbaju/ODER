megadepth::install_megadepth()
if (!exists("rec_url")) {
    rec_url <- recount::download_study(
        project = "SRP012682",
        type = "samples",
        download = FALSE
    ) # file_cache is an internal function to download a bigwig file from a link
    # if the file has been downloaded recently, it will be retrieved from a cache
}
bw_path <- file_cache(rec_url[1])

ex_opt_ers <- GenomicRanges::GRanges( # this is created to not overlap
    seqnames = S4Vectors::Rle(c("chr21"), c(3)),
    ranges = IRanges::IRanges(
        start = c(5032176, 5033408, 5034717),
        end = c(5032217, 5033425, 5034756)
    )
)

bw_paths <- c(bw_path, bw_path)
col_info <- data.frame(dummy1 = rep("run1", length(bw_paths)))

example_cm <- suppressWarnings(get_count_matrix(bw_paths = bw_paths, annot_ers = ex_opt_ers, cols = col_info))
example_cm2 <- suppressWarnings(get_count_matrix(bw_paths = bw_path, annot_ers = ex_opt_ers))

test_that("get_count_matrix works", {
    expect_equal(length(bw_paths), ncol(SummarizedExperiment::assay(example_cm)))
    expect_equal(length(ex_opt_ers), nrow(SummarizedExperiment::assay(example_cm)))
    expect_equal(col_info, as.data.frame(SummarizedExperiment::colData(example_cm)))
    expect_true(methods::is(example_cm, "RangedSummarizedExperiment"))
    expect_true(methods::is(example_cm2, "RangedSummarizedExperiment"))
    expect_true(identical(SummarizedExperiment::assay(example_cm)[, 1], SummarizedExperiment::assay(example_cm2)[, 1]))
})
