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

gtf_url <- "http://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz"
gtf_path <- ODER:::.file_cache(gtf_url)
exons_no_overlap <- get_exons(gtf = gtf_path, ucsc_chr = TRUE, ignore.strand = TRUE)

ers_delta_example <- get_ers_delta(ers = ers_example, opt_exons = exons_no_overlap, delta_fun = .delta)

# gtex_url <- "https://storage.googleapis.com/gtex_analysis_v6p/rna_seq_data/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct.gz"
# gtex_path <- ODER:::.file_cache(gtex_url)
# library(data.table)
# gtex_data <- fread(gtex_path)


usethis::use_data(auc_example, overwrite = TRUE)
usethis::use_data(coverage_example, overwrite = TRUE) # might have to exclude - 1.4MB
usethis::use_data(ers_example, overwrite = TRUE)
usethis::use_data(ers_delta_example, overwrite = TRUE)
