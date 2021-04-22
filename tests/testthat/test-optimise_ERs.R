gtf_url <- paste0(
    "http://ftp.ensembl.org/pub/release-103/",
    "gtf/homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz"
)

gtf_path <- ODER:::.file_cache(gtf_url)

test_exons <- get_exons(
    gtf = gtf_path,
    ucsc_chr = TRUE,
    ignore.strand = TRUE
)

test_that("get_exons works", {
    expect_error(
        get_exons(
            gtf = "gtf_file.gtx",
            ucsc_chr = TRUE,
            ignore.strand = TRUE
        ),
        "Please check your gtf file path"
    )
    expect_error(
        get_exons(
            gtf = "gtf_file.gtx.gz",
            ucsc_chr = TRUE,
            ignore.strand = TRUE
        ),
        "Please check your gtf file path"
    )
    expect_true(methods::is(test_exons, "GenomicRanges"))
})
