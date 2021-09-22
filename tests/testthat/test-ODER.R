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

if (!exists("test_ers_delta")) {
    test_ers_delta <- get_ers_delta(
        ers = gtex_lung_ers_1, # gtex_lung_ers_1 is from the data folder
        opt_exons = test_exons,
        delta_fun = .delta
    )
}

if (!exists("test_opt_ers")) {
    test_opt_ers <- get_opt_ers(
        ers = gtex_lung_ers_1,
        ers_delta = test_ers_delta
    )
}

if (!exists("rec_url")) {
    rec_url <- recount::download_study(
        project = "SRP012682",
        type = "samples",
        download = FALSE
    )
}
bw_path <- ODER::file_cache(rec_url[1])


test_ODER_opt_ers <- ODER(
    bw_paths = bw_path, auc_raw = gtex_lung_auc_1,
    auc_target = 40e6 * 100, chrs = c("chr21", "chr22"),
    genome = "hg38", mccs = c(5, 10), mrgs = c(10, 20),
    gtf = gtf_path, ucsc_chr = TRUE, ignore.strand = TRUE,
    exons_no_overlap = NULL, bw_chr = "chr"
)

test_that("ODER works", {
    expect_equal(test_opt_ers, test_ODER_opt_ers)
    expect_true(methods::is(test_ODER_opt_ers, "list"))
    expect_true(methods::is(test_ODER_opt_ers[[1]], "GenomicRanges"))
    expect_true(methods::is(test_ODER_opt_ers[[2]], "character"))
    expect_true(methods::is(test_ODER_opt_ers[[3]], "data.frame"))
})
