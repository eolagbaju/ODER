megadepth::install_megadepth()
if (!exists("gtex_metadata")) {
    gtex_metadata <- recount::all_metadata("gtex")
    gtex_metadata <- gtex_metadata %>%
        as.data.frame() %>%
        dplyr::filter(project == "SRP012682")
}
if (!exists("rec_url")) {
    rec_url <- recount::download_study(
        project = "SRP012682",
        type = "samples",
        download = FALSE
    )
}
run1 <- gtex_metadata[["run"]][[84]]
run2 <- gtex_metadata[["run"]][[91]]
runs <- c(run1, run2)
col_info <- as.data.frame(runs)

bw_path <- ODER:::.file_cache(rec_url[84])
bw_path2 <- ODER:::.file_cache(rec_url[91])
bw_paths <- c(bw_path, bw_path2)
if (!exists("gtf_path")) {
    gtf_url <- "http://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz"
    gtf_path <- ODER:::.file_cache(gtf_url)
}
if (!exists("test_opt_ers84")) {
    test_opt_ers84 <- suppressWarnings(ODER(
        bw_paths = bw_path, auc_raw = gtex_metadata[["auc"]][84],
        auc_target = 40e6 * 100, chrs = c("chr21", "chr22"),
        genome = "hg38", mccs = c(2, 4, 6, 8, 10), mrgs = c(10, 20, 30),
        gtf = gtf_path, ucsc_chr = TRUE, ignore.strand = TRUE,
        exons_no_overlap = NULL, bw_chr = "chr"
    ))
}
if (!exists("test_opt_ers91")) {
    test_opt_ers91 <- suppressWarnings(ODER(
        bw_paths = bw_path, auc_raw = gtex_metadata[["auc"]][91],
        auc_target = 40e6 * 100, chrs = c("chr21", "chr22"),
        genome = "hg38", mccs = c(2, 4, 6, 8, 10), mrgs = c(10, 20, 30),
        gtf = gtf_path, ucsc_chr = TRUE, ignore.strand = TRUE,
        exons_no_overlap = NULL, bw_chr = "chr"
    ))
}

test_juncs <- SummarizedExperiment::rowRanges(dasper::junctions_example)

if (!exists("test_gstate")) {
    test_gstate <- suppressWarnings(generate_genomic_state(
        gtf = gtf_path,
        chrs_to_keep = c("21", "22"),
        ensembl = TRUE
    ))
}

aers <- suppressWarnings(annotatERs(test_opt_ers84[["opt_ers"]],
    junc_data = test_juncs,
    gtf_path = gtf_path, ensembl = TRUE, genom_state = test_gstate
))
aers2 <- suppressWarnings(annotatERs(test_opt_ers91[["opt_ers"]],
    junc_data = test_juncs,
    gtf_path = gtf_path, ensembl = TRUE, genom_state = test_gstate
))
aersl <- GenomicRanges::GRangesList(aers, aers2)

test_cm <- suppressWarnings(get_count_matrix(bw_paths = bw_paths, annot_ers = aersl, cols = col_info))
test_cm2 <- suppressWarnings(get_count_matrix(bw_paths = bw_path, annot_ers = aers))

test_that("get_count_matrix works", {
    expect_equal(length(bw_paths), ncol(SummarizedExperiment::assay(test_cm)))
    expect_equal(length(aersl[[1]]), nrow(SummarizedExperiment::assay(test_cm)))
    expect_true(methods::is(test_cm, "RangedSummarizedExperiment"))
    expect_true(identical(SummarizedExperiment::assay(test_cm)[, 1], SummarizedExperiment::assay(test_cm2)[, 1]))
})
