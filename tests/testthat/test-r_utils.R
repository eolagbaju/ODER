

test_that("bw_check works", {
    test_bw_paths <- c("testpath1.bw", "testpath2.bq")
    test_bw_paths2 <- c("testpath3.bw", "testpath4.bw")

    expect_true(bw_check("testpath1.bw"))
    expect_false(bw_check("testpath2.bq"))
    expect_false(bw_check(test_bw_paths))
    expect_true(bw_check(test_bw_paths2))
})


test_that("chr_formatting works", {
    expect_equal(chr_formatting("chr1", "chr"), "chr1")
    expect_equal(chr_formatting("chr1", "nochr"), "1")
    expect_equal(chr_formatting("chrM", "nochr"), "MT")
})


test_that("informatting works", {
    default <- c(
        "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
        "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
        "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"
    )

    expect_equal(informatting(c("chr1", "chr2")), c("chr1", "chr2"))
    expect_equal(informatting(c("1", "chr2")), c("chr1", "chr2"))
    expect_equal(informatting(c(1, 2, "X")), c("chr1", "chr2", "chrX"))
    expect_equal(informatting(""), default)
    expect_equal(informatting("M"), "chrM")
    expect_equal(informatting("MT"), "chrM")
    expect_equal(informatting(c(1, 2, "M")), c("chr1", "chr2", "chrM"))
})
