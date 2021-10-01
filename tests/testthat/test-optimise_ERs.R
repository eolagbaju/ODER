if (!exists("gtf_path")) {
    gtf_url <- paste0(
        "http://ftp.ensembl.org/pub/release-103/",
        "gtf/homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz"
    )
    gtf_path <- file_cache(gtf_url)
}

if (!exists("gtf_grs")) {
    gtf_grs <- rtracklayer::import(gtf_path)
}

if (!exists("test_exons")) {
    test_exons <- get_exons(
        gtf = gtf_grs,
        ucsc_chr = TRUE,
        ignore.strand = TRUE
    )
}

suppressWarnings(test_3p_exons <- get_exons(gtf = gtf_grs, ucsc_chr = TRUE, ignore.strand = TRUE, biotype = "Three Prime"))
suppressWarnings(test_5p_exons <- get_exons(gtf = gtf_grs, ucsc_chr = TRUE, ignore.strand = TRUE, biotype = "Five Prime"))
suppressWarnings(test_int_exons <- get_exons(gtf = gtf_grs, ucsc_chr = TRUE, ignore.strand = TRUE, biotype = "Internal"))
suppressWarnings(test_lnc_exons <- get_exons(gtf = gtf_grs, ucsc_chr = TRUE, ignore.strand = TRUE, biotype = "lncRNA"))
suppressWarnings(test_nc_exons <- get_exons(gtf = gtf_grs, ucsc_chr = TRUE, ignore.strand = TRUE, biotype = "ncRNA"))
suppressWarnings(test_ps_exons <- get_exons(gtf = gtf_grs, ucsc_chr = TRUE, ignore.strand = TRUE, biotype = "Pseudogene"))

data(gtex_SRP012682_SRX222703_lung_ers_1, package = "ODER")

if (!exists("test_ers_delta")) {
    test_ers_delta <- get_ers_delta(
        ers = gtex_SRP012682_SRX222703_lung_ers_1, # gtex_SRP012682_SRX222703_lung_ers_1 is from the data folder
        opt_exons = test_exons,
        delta_fun = .delta
    )
}

if (!exists("test_opt_ers")) {
    test_opt_ers <- get_opt_ers(
        ers = gtex_SRP012682_SRX222703_lung_ers_1,
        ers_delta = test_ers_delta
    )
}

test_grs <- GenomicRanges::GRanges(
    seqnames = S4Vectors::Rle(c("chr1"), c(10)),
    ranges = IRanges::IRanges(
        start = c(12975, 24738, 18157, 17915, 17554, 17233, 16858, 16570, 15796, 15005),
        end = c(13052, 24891, 18354, 18061, 17728, 17368, 17055, 16723, 15947, 15038),
        names = head(letters, 10)
    ),
    strand = S4Vectors::Rle((c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
    score = 1:10,
    GC = seq(1, 0, length = 10)
)

test_deltas <- .delta(test_grs, test_exons)

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
    expect_true(methods::is(test_3p_exons, "GenomicRanges"))
    expect_true(methods::is(test_5p_exons, "GenomicRanges"))
    expect_true(methods::is(test_int_exons, "GenomicRanges"))
    expect_true(methods::is(test_lnc_exons, "GenomicRanges"))
    expect_true(methods::is(test_nc_exons, "GenomicRanges"))
    expect_true(methods::is(test_ps_exons, "GenomicRanges"))
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
            ers = gtex_SRP012682_SRX222703_lung_ers_1,
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
            ers = gtex_SRP012682_SRX222703_lung_ers_1
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

test_that(".delta works", {
    expect_equal(test_deltas[["sum"]], 66)
    expect_equal(test_deltas[["mean"]], 16.5)
    expect_equal(test_deltas[["median"]], 0)
    expect_equal(test_deltas[["n_eq_0"]], 3)
    expect_equal(test_deltas[["propor_eq_0"]], 0.75)
})
