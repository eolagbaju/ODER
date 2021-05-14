

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
