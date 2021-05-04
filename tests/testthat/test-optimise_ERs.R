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

test_ers_delta <- get_ers_delta(
    ers = ers_example,
    opt_exons = test_exons,
    delta_fun = .delta
)

test_opt_ers <- get_opt_ers(
    ers = ers_example,
    ers_delta = test_ers_delta
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

test_that("get_ers_delta works", {
    expect_error(
        get_ers_delta(
            opt_exons = test_exons,
            delta_fun = .delta
        ),
        "No ERs were entered"
    )
    expect_error(
        get_ers_delta(
            ers = ers_example,
            delta_fun = .delta
        ),
        "No opt_exons were entered"
    )


    expect_true(methods::is(test_ers_delta, "data.frame"))
})

test_that("get_opt_ers works", {
    expect_error(
        get_opt_ers(
            ers_delta = test_ers_delta
        ),
        "No ERs were entered"
    )
    expect_error(
        get_opt_ers(
            ers = ers_example
        ),
        "No ers_delta were entered"
    )

    expect_true(methods::is(test_opt_ers, "list"))

    expect_equal(length(test_opt_ers), 3)
})
