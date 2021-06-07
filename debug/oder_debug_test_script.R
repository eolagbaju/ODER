library("ODER")
library(magrittr)
sra_metadata <- recount::all_metadata("sra")

sra_metadata <- sra_metadata %>%
  as.data.frame() %>%
  dplyr::filter(project == "DRP000425")

# obtain path to example bw on recount2
test_url <- recount::download_study(
  project = "DRP000425",
  type = "samples",
  download = FALSE
)

# test_bw_path <- ODER:::.file_cache(test_url[1])

#https://rdrr.io/bioc/rtracklayer/man/BigWigFile.html
# which <- GenomicRanges::GRanges(c("2", "2"), IRanges(c(1, 300), c(400, 1000)))
#impbw <- rtracklayer::import(con="/data/RNA_seq_diag/mito/bw/L1901-3262F.all.bw",which = which)
# Warning message:
#   In .local(con, format, text, ...) :
#   'which' contains seqnames not known to BigWig file: chr21
 impbw <- rtracklayer::import(con="/data/RNA_seq_diag/mito/bw/L1901-3262F.all.bw")
 impbw_df <- data.frame(unique(impbw))
 impbw_chr_names <- unique(impbw_df[["seqnames"]])
# 
# test_bw_path <- "/data/RNA_seq_diag/mito/bw/L1901-3262F.all.bw"
test_bw_path <- "/home/ssethi/Data_for_ODER/GCB1.bw"

test_coverage <- get_coverage(
  bw_paths = test_bw_path,
  auc_raw = sra_metadata[["auc"]][1],
  auc_target = 40e6 * 100, # target 40 million coverage with 100 bp length reads
  chrs = c("chr21", "chr22"),
  bw_chr = "nochr"
) # defaults to chr1-22, chrX, chrY, chrM

test_ers <- get_ers(coverage = test_coverage, mccs = c(5,10), mrgs = c(10, 20))

# gtf_url <- "http://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz"
# gtf_path <- ODER:::.file_cache(gtf_url)
# exons_no_overlap <- get_exons(gtf = gtf_path, ucsc_chr = TRUE, ignore.strand = TRUE)
test_gtf <- "/home/ssethi/Data_for_ODER/Homo_sapiens.GRCh38.94.gtf"
exons_no_overlap <- get_exons(gtf = test_gtf, ucsc_chr = TRUE, ignore.strand = TRUE)

test_ers_delta <- get_ers_delta(ers = test_ers, opt_exons = exons_no_overlap)

test_opt_ers <- get_opt_ers(ers = test_ers, ers_delta = test_ers_delta)

print(test_opt_ers)


opt_ers <- ODER(
     bw_paths = test_bw_path, auc_raw = auc_example,
     auc_target = 40e6 * 100, chrs = c("chr21", "chr22"),
    genome = "hg38", mccs = c(5, 10), mrgs = c(10, 20),
     gtf = test_gtf, ucsc_chr = TRUE, ignore.strand = TRUE,
     exons_no_overlap = NULL, bw_chr = "nochr"
 )

