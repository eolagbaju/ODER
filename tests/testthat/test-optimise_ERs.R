gtf_url <- "http://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz"
test_gtf_tempfile <- file.path(tempdir(), "gtf_file.gtx")
download.file(url = gtf_url, destfile = test_gtf_tempfile)

test_that("get_exons works", {
    expect_error(
        get_exons(gtf = test_gtf_tempfile, ucsc_chr = T, ignore.strand = T),
        "Please check your gtf file path"
    )
    # expect_type(coverage_example, "list")
    # expect_type(coverage_example[["chr21"]][["meanCoverage"]], "S4") # output should be an Rle which is an S4
    # methods::is(coverage_example[["chr21"]][["meanCoverage"]], "Rle")
    # expect_equal(length(c("chr21", "chr22")), length(coverage_example))
})
