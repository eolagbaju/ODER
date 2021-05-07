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
    ers = ers_example, # ers_example is from the data folder
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
    ers_delta_cols <- c("mcc", "mrg", "sum", "mean", "median", "n_eq_0", "propor_eq_0")

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

    expect_true(all(ifelse(colnames(test_ers_delta) == ers_delta_cols, TRUE, FALSE)))
    expect_true(methods::is(test_ers_delta, "data.frame"))

    expect_type(test_ers_delta[["mcc"]], "double")
    expect_type(test_ers_delta[["mrg"]], "double")
    expect_type(test_ers_delta[["sum"]], "integer")
    expect_type(test_ers_delta[["mean"]], "double")
    expect_type(test_ers_delta[["median"]], "double")
    expect_type(test_ers_delta[["n_eq_0"]], "integer")
    expect_type(test_ers_delta[["propor_eq_0"]], "double")
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
    expect_true(methods::is(test_opt_ers[[1]], "GenomicRanges"))
    expect_true(methods::is(test_opt_ers[[2]], "character"))
    expect_true(methods::is(test_opt_ers[[3]], "data.frame"))

    expect_equal(length(test_opt_ers[[2]]), 2)
    expect_equal(length(test_opt_ers), 3)
})
