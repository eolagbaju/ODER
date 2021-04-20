gtex_metadata <- recount::all_metadata("gtex")

gtex_metadata <- gtex_metadata %>%
    as.data.frame() %>%
    dplyr::filter(project == "SRP012682")

# obtain path to example bw on recount2
url <- recount::download_study(
    project = "SRP012682",
    type = "samples",
    download = FALSE
)

bw_path <- .file_cache(url[1])

test_chrs <- c("chr21", "chr22")
test_coverage <- get_coverage(
    bw_paths = bw_path,
    auc_raw = gtex_metadata[["auc"]][1],
    auc_target = 40e6 * 100,
    chrs = test_chrs
)

test_that("get_chr_info works", {
    chrs <- c("chr1", "chr2", "chr3")
    test_info <- get_chr_info(chrs = chrs, genome = "hg38")

    expect_false(test_info %contain% "chr4")
    expect_error(get_chr_info(chrs = "chromononsense", genome = "hg38"), "Non-UCSC chromosome format was entered")
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
            auc_raw = gtex_metadata[["auc"]][1],
            auc_target = "Not a number",
            chrs = c("chr21", "chr22")
        ),
        "Please enter a valid number for the auc values"
    )
    expect_error(
        get_coverage(
            bw_paths = "",
            auc_raw = gtex_metadata[["auc"]][1],
            auc_target = 40e6 * 100,
            chrs = c("chr21", "chr22")
        ),
        "Please check your bigwig file paths are correct"
    )
    expect_error(
        get_coverage(
            auc_raw = gtex_metadata[["auc"]][1],
            auc_target = 40e6 * 100,
            chrs = c("chr21", "chr22")
        ),
        "bw_paths is empty"
    )
    expect_error(
        get_coverage(
            bw_paths = "~/.cache/BiocFileCache/25c6571687b_SRR660824_SRS389722_SRX222703_male_lung.png",
            auc_raw = gtex_metadata[["auc"]][1],
            auc_target = 40e6 * 100,
            chrs = c("chr21", "chr22")
        ),
        "Please check your bigwig file paths are correct"
    )
    expect_type(test_coverage, "list")
    expect_type(test_coverage[["chr21"]][["meanCoverage"]], "S4") # output should be an Rle which is an S4
    expect_equal(length(test_chrs), length(test_coverage))
})

test_that("get_ers works", {
    test_ers <- get_ers(coverage = test_coverage, mccs = c(5, 10), mrgs = c(10, 20))
    expect_error(get_ers(mccs = c(5, 10), mrgs = c(10, 20)), "Coverage is missing")
    expect_error(get_ers(coverage = coverage, mrgs = c(10, 20)), "Mean Coverage Cutoff is empty")
    expect_error(get_ers(coverage = coverage, mccs = c(5, 10)), "Max Region Gap is empty")

    expect_type(test_ers, "list")
    expect_type(test_ers[[1]][[1]], "S4")
})
