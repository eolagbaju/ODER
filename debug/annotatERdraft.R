library(dasper)

load(file="/data/recount/GTEx_SRP012682/gtex_split_read_table_annotated_rda/lung_split_read_table_annotated.rda")

lung_junc <- gtex_split_read_table_annotated_only_junc_coverage

# junctions <-
#   junction_load(
#     junction_paths = lung_junc,
#     controls = "fibroblasts",
#     chrs = c("chr21", "chr22")
#   )

#keep <- c("junID","chr","start","stop","strand","countsSamples")
keep <- c("chr","start","stop","strand")
lung_trim <- lung_junc[keep]
colnames(lung_trim) <- c("chr","start","end","strand")

# RSE needs colnames to be numeric, converts True False to 1s and 0s
testlung <- lung_junc %>% 
  dplyr::select(chr,start,stop,strand,precBoundDonor) %>%
  dplyr::mutate(precBoundDonor = ifelse(precBoundDonor == "TRUE",1,0))


test_junction <- SummarizedExperiment::makeSummarizedExperimentFromDataFrame(testlung)



gtf_url <- "http://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz"
gtf_path <- ODER:::.file_cache(gtf_url)
junction_annot(junctions = test_junction,ref = gtf_path)


lungrange <-  GenomicRanges::makeGRangesFromDataFrame(head(lung_junc,500))
annotated_junctions <- junction_annot(junctions = lungrange,ref = gtf_path)