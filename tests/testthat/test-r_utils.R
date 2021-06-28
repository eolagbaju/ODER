

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

test_that("informatting2 works", {
    default <- c(
        "1", "2", "3", "4", "5", "6", "7", "8", "9",
        "10", "11", "12", "13", "14", "15", "16", "17",
        "18", "19", "20", "21", "22", "X", "Y", "MT"
    )

    expect_equal(informatting2(c("1", "2")), c("1", "2"))
    expect_equal(informatting2(c("chr1", "chr2")), c("1", "2"))
    expect_equal(informatting2(c("chr1", "chr2", "chrX")), c(1, 2, "X"))
    expect_equal(informatting2(""), default)
    expect_equal(informatting2("chrM"), "MT")
    expect_equal(informatting2(c("chr1", "chr2", "chrM")), c(1, 2, "MT"))
})

test_that("inbetween works", {
    expect_true(inbetween(5, 4, 6))
    expect_false(inbetween(2, 4, 6))
    expect_true(inbetween(4, 4, 6))
})

test_gr1 <- GenomicRanges::GRanges("chr21", IRanges::IRanges(5026423, 5323718), strand = "+")
test_gr2 <- GenomicRanges::GRanges("chr21", IRanges::IRanges(c(13261708, 13482787), c(13482707, 13851281)), strand = "-")
test_gr3 <- GenomicRanges::GRanges("chr21", IRanges::IRanges(c(25085652, 25086122), c(25086438, 25087197)), strand = "+")

test_that("colgrs works", {
    expect_false(colgrs(test_gr1))
    expect_false(colgrs(test_gr2))
    expect_true(colgrs(test_gr3))
})

test_that("colgrs works", {
    expect_true(inv_colgrs(test_gr1))
    expect_true(inv_colgrs(test_gr2))
    expect_false(inv_colgrs(test_gr3))
})
