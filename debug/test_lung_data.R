gtex_metadata <- recount::all_metadata("gtex")
gtex_metadata <- gtex_metadata %>%
  as.data.frame() %>%
  dplyr::filter(project == "SRP012682")


url <- recount::download_study(
  project = "SRP012682",
  type = "samples",
  download = FALSE
)
# .file_cache is an internal function to download a bigwig file from a link if the file has been downloaded recently, it will be retrieved from a cache
bw_path <- ODER:::.file_cache(url[1])
gtf_url <- "http://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz"
gtf_path <- ODER:::.file_cache(gtf_url)

opt_ers <- ODER(
  bw_paths = bw_path, auc_raw = gtex_metadata[["auc"]][1],
  auc_target = 40e6 * 100, chrs = "",
  genome = "hg38", mccs = c(2, 4, 6, 8, 10), mrgs = c(10, 20, 30),
  gtf = gtf_path, ucsc_chr = TRUE, ignore.strand = TRUE,
  exons_no_overlap = NULL, bw_chr = "chr"
)

load(file = "/data/recount/GTEx_SRP012682/gtex_split_read_table_annotated_rda/lung_split_read_table_annotated.rda")
lung_junc <- gtex_split_read_table_annotated_only_junc_coverage

aers <- annotatERs(opt_ers[["opt_ers"]], junc_data = lung_junc,
                   gtf_path = gtf_path, ensembl = TRUE)

ref_results <- refine_ERs(aers)

refined_ers <- ref_results[[1]]

changes <- ref_results[[2]]
