gtex_metadata <- recount::all_metadata("gtex")
gtex_metadata <- gtex_metadata %>%
  as.data.frame() %>%
  dplyr::filter(project == "SRP012682")

# obtain path to example bw on recount2
url <- recount::download_study(
  project = "SRP012682",
  type = "samples",
  download = FALSE
)

# 84  91   98  106  141  217  245  253  280  301  306  310  318  328  330  375  396  488  522  547  561  588  599
bw_plus <- ODER:::.file_cache(url[58])
bw_minus <- ODER:::.file_cache(url[84])

stranded_ers2 <- get_strand_ers(bw_pos = bw_plus, bw_neg = bw_minus, auc_raw_pos = gtex_metadata[["auc"]][58], 
                               auc_raw_neg = gtex_metadata[["auc"]][84], auc_tar_pos = 40e6 * 100, 
                               auc_tar_neg = 40e6 * 100, chrs = c("chr21","chr22"), mccs = c(5,10), mrgs = c(10, 20))

stranded_ers <- get_strand_ers(bw_pos = bw_plus, bw_neg = bw_minus, auc_raw_pos = gtex_metadata[["auc"]][58], 
                               auc_raw_neg = gtex_metadata[["auc"]][84], auc_tar_pos = 40e6 * 100, 
                               auc_tar_neg = 40e6 * 100, chrs = c("chr21","chr22"), mccs = c(2, 4, 6, 8, 10), mrgs = c(10, 20, 30))

coverage_plus <- get_coverage(
  bw_paths = bw_plus,
  auc_raw = gtex_metadata[["auc"]][58],
  auc_target = 40e6 * 100, # target 40 million coverage with 100 bp length reads
  chrs = c("chr21")
)

coverage_minus <- get_coverage(
  bw_paths = bw_minus,
  auc_raw = gtex_metadata[["auc"]][84],
  auc_target = 40e6 * 100, # target 40 million coverage with 100 bp length reads
  chrs = c("chr21")
)
 
ers_plus <- get_ers(coverage = coverage_plus, mccs = c(2, 4, 6, 8, 10),mrgs = c(10, 20, 30))
ers_minus <- get_ers(coverage = coverage_minus, mccs = c(2, 4, 6, 8, 10), mrgs =  c(10, 20, 30))
# ers_plus <- get_ers(coverage = coverage_plus, mccs = c(5, 10), mrgs = c(10, 20))
# ers_minus <- get_ers(coverage = coverage_minus, mccs = c(5, 10), mrgs = c(10, 20))

sublist <- vector("list", length(ers_plus[[1]]))
ers_combi <- vector("list", length(ers_plus))

for (i in 1:length(ers_combi)){
  ers_combi[[i]] <- sublist
}

for (i in 1:length(ers_plus)){
  names(ers_combi) <- names(ers_plus)
  for (j in 1:length(ers_plus[[i]])){
    # print(names(ers_plus[[i]])[j])
    # print(j)
    names(ers_combi[[i]]) <- names(ers_plus[[j]])
    strand(ers_plus[[i]][[j]]) <- "+"
    strand(ers_minus[[i]][[j]]) <- "-"
    ers_combi[[i]][[j]] <- c(ers_plus[[i]][[j]],ers_minus[[i]][[j]])
    GenomeInfoDb::sortSeqlevels(ers_combi[[i]][[j]])
    BiocGenerics::sort(ers_combi[[i]][[j]])
  }
}
gtf_url <- "http://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz"
gtf_path <- ODER:::.file_cache(gtf_url)

exons_no_overlap <- get_exons(gtf = gtf_path, ucsc_chr = TRUE, ignore.strand = FALSE)

ers_delta <- get_ers_delta(ers = stranded_ers, opt_exons = exons_no_overlap)

opt_ers <- get_opt_ers(ers = stranded_ers, ers_delta = ers_delta)

fibro_juncs <- SummarizedExperiment::rowRanges(dasper::junctions_example)

annot_ers <- annotatERs(
  opt_ers = opt_ers[["opt_ers"]], junc_data = fibro_juncs,
  gtf_path = gtf_path, chrs_to_keep = c("21", "22"), ensembl = TRUE
)

refined_ers <- refine_ERs(annot_ers)
