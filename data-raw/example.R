## code to prepare `example` dataset goes here

gtex_metadata <- recount::all_metadata("gtex")

gtex_metadata <- gtex_metadata %>%
    as.data.frame() %>%
    dplyr::filter(project == "SRP012682")

url <- recount::download_study(
    project = "SRP012682",
    type = "samples",
    download = FALSE
)

gtex_lung_auc_1 <- gtex_metadata[["auc"]][1]
bw_path <- ODER::file_cache(url[1])

gtex_lung_coverage_1 <- ODER::get_coverage(
    bw_paths = bw_path,
    auc_raw = gtex_lung_auc_1,
    auc_target = 40e6 * 100,
    chrs = c("chr21", "chr22")
)

gtex_lung_ers_1 <- ODER::get_ers(
    coverage = gtex_lung_coverage_1, mccs = c(5, 10), mrgs = c(10, 20)
)

gtf_url <- paste0(
    "http://ftp.ensembl.org/pub/release-103/gtf/",
    "homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz"
)
gtf_path <- ODER::file_cache(gtf_url)
exons_no_overlap <- get_exons(
    gtf = gtf_path, ucsc_chr = TRUE, ignore.strand = TRUE
)

gtex_lung_erdelta_1 <- get_ers_delta(
    ers = gtex_lung_ers_1, opt_exons = exons_no_overlap, delta_fun = .delta
)

usethis::use_data(gtex_lung_auc_1, overwrite = TRUE)
usethis::use_data(gtex_lung_coverage_1, overwrite = TRUE)
usethis::use_data(gtex_lung_ers_1, overwrite = TRUE)
usethis::use_data(gtex_lung_erdelta_1, overwrite = TRUE)
