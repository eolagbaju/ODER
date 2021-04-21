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

auc_example <- gtex_metadata[["auc"]][1]
bw_path <- ODER:::.file_cache(url[1])

coverage_example <- ODER::get_coverage(
    bw_paths = bw_path,
    auc_raw = auc_example,
    auc_target = 40e6 * 100,
    chrs = c("chr21", "chr22")
)

ers_example <- ODER::get_ers(coverage = coverage_example, mccs = c(5, 10), mrgs = c(10, 20))

usethis::use_data(auc_example, overwrite = TRUE)
usethis::use_data(coverage_example, overwrite = TRUE) # might have to exclude - 1.4MB
usethis::use_data(ers_example, overwrite = TRUE)
