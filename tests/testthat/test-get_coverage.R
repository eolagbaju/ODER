gtex_metadata <- recount::all_metadata("gtex")
gtex_metadata <- gtex_metadata %>%
    as.data.frame() %>%
    dplyr::filter(project == "SRP012682")

url <- recount::download_study(
    project = "SRP012682",
    type = "samples",
    download = FALSE
)

bw_path <- .file_cache(url[1])

bw_plus <- ODER:::.file_cache(url[58])
bw_minus <- ODER:::.file_cache(url[84])

test_strand_ers <- get_strand_ers(
    bw_pos = bw_plus, bw_neg = bw_minus, auc_raw_pos = gtex_metadata[["auc"]][58],
    auc_raw_neg = gtex_metadata[["auc"]][84], auc_target = 40e6 * 100,
    chrs = "chr21", mccs = c(5, 10), mrgs = c(10, 20), bw_chr = "chr"
)

test_that("get_chr_info works", {
    chrs <- c("chr1", "chr2", "chr3")
    test_info <- get_chr_info(chrs = chrs, genome = "hg38")

    expect_false(test_info %contain% "chr4")
    expect_error(
        get_chr_info(
            chrs = "chromononsense",
            genome = "hg38"
        ),
        "Non-UCSC chromosome format was entered"
    )
    expect_equal(nrow(test_info), length(chrs))
    expect_type(test_info[["chrom"]], "character")
    expect_type(test_info[[2]], "integer")
    expect_type(test_info[[3]], "logical")
    expect_type(test_info[[4]], "logical")
})

test_that("get_coverage works", {
    expect_error(
        get_coverage(
            bw_paths = bw_path,
            auc_raw = "Not a number",
            auc_target = 40e6 * 100,
            chrs = c("chr21", "chr22")
        ),
        "Please enter a valid number for the auc values"
    )
    expect_error(
        get_coverage(
            bw_paths = bw_path,
            auc_raw = auc_example, # auc_example is from the data folder
            auc_target = "Not a number",
            chrs = c("chr21", "chr22")
        ),
        "Please enter a valid number for the auc values"
    )
    expect_error(
        get_coverage(
            bw_paths = "",
            auc_raw = auc_example,
            auc_target = 40e6 * 100,
            chrs = c("chr21", "chr22")
        ),
        "Please check your bigwig file paths are correct"
    )
    expect_error(
        get_coverage(
            auc_raw = auc_example,
            auc_target = 40e6 * 100,
            chrs = c("chr21", "chr22")
        ),
        "bw_paths is empty"
    )
    expect_error(
        get_coverage(
            bw_paths = paste0(
                "~/.cache/BiocFileCache/25c6571687b_",
                "SRR660824_SRS389722_SRX222703_male_lung.png"
            ),
            auc_raw = auc_example,
            auc_target = 40e6 * 100,
            chrs = c("chr21", "chr22")
        ),
        "Please check your bigwig file paths are correct"
    )
    expect_type(coverage_example, "list") # coverage_example is from the data folder
    expect_type(coverage_example[["chr21"]][["meanCoverage"]], "S4")
    expect_true(methods::is(coverage_example[["chr21"]][["meanCoverage"]], "Rle"))
    expect_equal(length(c("chr21", "chr22")), length(coverage_example))
})

test_that("get_ers works", {
    test_er_df <- as.data.frame(ers_example[[1]][[1]])

    expect_error(
        get_ers(
            mccs = c(5, 10),
            mrgs = c(10, 20)
        ),
        "Coverage is missing"
    )
    expect_error(
        get_ers(
            coverage = coverage_example,
            mrgs = c(10, 20)
        ),
        "Mean Coverage Cutoff is empty"
    )
    expect_error(
        get_ers(
            coverage = coverage_example,
            mccs = c(5, 10)
        ),
        "Max Region Gap is empty"
    )

    expect_type(ers_example, "list") # ers_example is from the data folder
    expect_type(ers_example[[1]][[1]], "S4")
    expect_true(methods::is(ers_example[[1]][[1]], "GenomicRanges"))
    expect_true(methods::is(test_er_df, "data.frame"))
})

test_that("get_strand_ers works", {
    expect_true("+" %in% unlist(as.list(BiocGenerics::strand(test_strand_ers[[1]][[1]]))))
    expect_true("-" %in% unlist(as.list(BiocGenerics::strand(test_strand_ers[[1]][[1]]))))
    expect_false("*" %in% unlist(as.list(BiocGenerics::strand(test_strand_ers[[1]][[1]]))))
})
