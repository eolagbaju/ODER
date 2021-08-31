gtex_url <- "https://storage.googleapis.com/gtex_analysis_v6p/rna_seq_data/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct.gz"
gtex_path <- ODER:::.file_cache(gtex_url)
library(data.table)
gtex_data = fread(gtex_path)
outpath <- "/home/eolagbaju/projects/data"

data <- gtex_data %>% dplyr::mutate(Name = stringr::str_split_fixed(Name, "\\.", n=2)[,1])

save(data, file = stringr::str_c(outpath, "/data.Rdata"))

names <- colnames(data)
names <- names %>%
  stringr::str_replace_all("\\.+", "_") %>%
  stringr::str_replace("_$", "") %>%
  tolower()

for(i in 3:ncol(data)){
  
  file_name <- names[i]
  print(stringr::str_c(Sys.time(), "-", i, "-", file_name))
  df <- as.data.frame(data)[,c(1,i)] %>% dplyr::filter(.[[2]] > 0.1)
  save(df, file = stringr::str_c(outpath, "/", file_name, ".Rdata"))
  #write.table(df, str_c(outpath, "/", file_name, ".txt"), quote=FALSE, row.names=FALSE, sep="\t")
  
  rm(df)
}

############## expressed gene stuff ###################

gtex_metadata <- recount::all_metadata("gtex")
gtex_metadata <- gtex_metadata %>%
  as.data.frame() %>%
  dplyr::filter(project == "SRP012682")

url <- recount::download_study(
  project = "SRP012682",
  type = "samples",
  download = FALSE
)

bw_path <- ODER:::.file_cache(url[84])
gtf_url <- "http://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz"
gtf_path <- ODER:::.file_cache(gtf_url)

opt_ers <- ODER(
  bw_paths = bw_path, auc_raw = gtex_metadata[["auc"]][84],
  auc_target = 40e6 * 100, chrs = c("chr21","chr22"),
  genome = "hg38", mccs = c(2, 4, 6, 8, 10), mrgs = c(10, 20, 30),
  gtf = gtf_path, ucsc_chr = TRUE, ignore.strand = TRUE,
  exons_no_overlap = NULL, bw_chr = "chr"
)

test_juncs <- SummarizedExperiment::rowRanges(dasper::junctions_example)

aers <- annotatERs(opt_ers[["opt_ers"]], junc_data = test_juncs,
                   gtf_path = gtf_path, ensembl = TRUE)

fibro_df <- get_tissue(tissue = "cells - transformed fibroblasts")
expressed_genes <- get_expressed_genes(gtf_path = gtf_path, tissue_df = fibro_df)
aers_ng <- get_nearest_expressed_genes(annot_ers = aers, exp_genes = expressed_genes, gtf_path = gtf_path)
