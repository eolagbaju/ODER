gtex_metadata <- recount::all_metadata("gtex")
gtex_metadata <- gtex_metadata %>%
  as.data.frame() %>%
  dplyr::filter(project == "SRP012682")


url <- recount::download_study(
  project = "SRP012682",
  type = "samples",
  download = FALSE
)

run1 <- gtex_metadata[["run"]][[84]]
run2 <- gtex_metadata[["run"]][[91]]
runs <- c(run1,run2)
col_info <- as.data.frame(runs)

 
bw_path <- ODER:::.file_cache(url[84])
bw_path2 <- ODER:::.file_cache(url[91])
bw_paths <- c(bw_path,bw_path2)
gtf_url <- "http://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz"
gtf_path <- ODER:::.file_cache(gtf_url)

opt_ers <- ODER(
  bw_paths = bw_path, auc_raw = gtex_metadata[["auc"]][84],
  auc_target = 40e6 * 100, chrs = c("chr21","chr22"),
  genome = "hg38", mccs = c(2, 4, 6, 8, 10), mrgs = c(10, 20, 30),
  gtf = gtf_path, ucsc_chr = TRUE, ignore.strand = TRUE,
  exons_no_overlap = NULL, bw_chr = "chr"
)

opt_ers2 <- ODER(
  bw_paths = bw_path, auc_raw = gtex_metadata[["auc"]][91],
  auc_target = 40e6 * 100, chrs = c("chr21","chr22"),
  genome = "hg38", mccs = c(2, 4, 6, 8, 10), mrgs = c(10, 20, 30),
  gtf = gtf_path, ucsc_chr = TRUE, ignore.strand = TRUE,
  exons_no_overlap = NULL, bw_chr = "chr"
)

test_juncs <- SummarizedExperiment::rowRanges(dasper::junctions_example)

aers <- annotatERs(opt_ers[["opt_ers"]], junc_data = test_juncs,
                   gtf_path = gtf_path, ensembl = TRUE)
aers2 <- annotatERs(opt_ers2[["opt_ers"]], junc_data = test_juncs,
                   gtf_path = gtf_path, ensembl = TRUE)
aersl <- GenomicRanges::GRangesList(aers,aers2)
###### generate gene count matrix ######

setest <- get_count_matrix(bw_paths = bw_paths, annot_ers = aersl, cols = NULL)
######
dir <- tempdir()
fil <- tempfile("annotation.bed")
annot_bed1 <- rtracklayer::export.bed(aers,fil)
annot_bed1 <- rtracklayer::export.bed(aers,"annotation.bed")
annot_bed2 <- rtracklayer::export.bed(aers2,"annotation.bed")

mean_coverage1 <- megadepth::get_coverage(bigwig_file = bw_path, op = "mean", annotation = annot_bed1)
mean_coverage2 <- megadepth::get_coverage(bigwig_file = bw_path2, op = "mean", annotation = annot_bed2)

mcmat1 <- as.matrix(mcols(mean_coverage1)[["score"]])
mcmat2 <- as.matrix(mcols(mean_coverage2)[["score"]])
combimat <- cbind(mcmat1,mcmat2) # better to prebuild final matrix and then assign the columns


# chr_info <- get_chr_info(chrs = c("chr21","chr22"), genome = "hg38")
# chrseqlen1 <- chr_info[["size"]]
# names(chrseqlen1) <- c("chr21","chr22")
# GenomeInfoDb::seqlengths(test_grs) <- chrseqlen1
# GenomeInfoDb::seqlengths(test_grs2) <- chrseqlen2

# rtracklayer::export(object = test_grs, format = "bigwig", con = "test1.bw") #, con = big_file1)
# rtracklayer::export(object = test_grs2, format = "bigwig", con =  "test2.bw")

ranges <- GenomicRanges::granges(mean_coverage1)
names(combimat) <- c(1:length(ranges))
names(ranges) <- c(1:length(ranges))
rownames(col_info) <- as.character(c(1:nrow(col_info)))
se <- SummarizedExperiment::SummarizedExperiment(assays = combimat,
                                                 rowRanges = ranges,
                                                 colData = col_info)
