gtex_metadata <- recount::all_metadata("gtex")
gtex_metadata <- gtex_metadata %>%
  as.data.frame() %>%
  dplyr::filter(project == "SRP012682")


url <- recount::download_study(
  project = "SRP012682",
  type = "samples",
  download = FALSE
)

run <- gtex_metadata[["run"]][[1]]

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
lung_junc_21_22 <- dplyr::filter(lung_junc,chr == 21 | chr == 22)
save(lung_junc_21_22, file = "lung_junc_example.RData")

aers <- annotatERs(opt_ers[["opt_ers"]], junc_data = lung_junc,
                   gtf_path = gtf_path, ensembl = TRUE)

aers2 <- aers
genesource <- character(length(aers2))

mcols(aers2)$gene_source <- genesource
mcols(aers2[lengths(mcols(aers2)[["genes"]])>0])$gene_source <- "junction(s)"

#aers[lengths(mcols(aers)[["genes"]])==0] #getting annotated junctions with no associated gene
ngaers <- aers[lengths(mcols(aers)[["genes"]])==0]
ngaers2 <- ngaers
gtf_gr <- rtracklayer::import(gtf_path)
genes_gr <- gtf_gr[gtf_gr$type == "gene"]
GenomeInfoDb::seqlevelsStyle(genes_gr) <- "UCSC"
#?GenomicRanges::distanceToNearest() function to use? or ?GenomicRanges::nearest()
nearest_genes <- GenomicRanges::nearest(ngaers,genes_gr)
ngdist <- GenomicRanges::distanceToNearest(aers,genes_gr)
overmaxdist <- ngdist[mcols(ngdist)[["distance"]]>10000]
mcols(aers2[queryHits(overmaxdist)])[["genes"]] <- ""
mcols(aers2[queryHits(overmaxdist)])[["gene_source"]] <- "Too far"
#GenomicRanges::distanceToNearest(ngaers,genes_gr)
missing_genes <- mcols(genes_gr[nearest_genes])[["gene_id"]]
#mcols(ngaers2)[["genes"]] <- missing_genes

mcols(aers2[mcols(ngaers2)[["og_index"]]])[["gene"]] <- missing_genes
mcols(aers2[mcols(ngaers2)[["og_index"]]])[["gene_source"]] <- "nearest gtf genes"


ref_results <- refine_ERs(aers)

refined_ers <- ref_results[[1]]

changes <- ref_results[[2]]

###### generate gene count matrix ######

annot_bed <- rtracklayer::export.bed(aers,"annotation.bed")

mean_coverage <- megadepth::get_coverage(bigwig_file = bw_path,op = "mean", annotation = annot_bed)

#mclist <- as.list(round(mcols(mean_coverage)[["score"]])) %>% S4Vectors::SimpleList()

mcmat <- as.matrix(mcols(mean_coverage)[["score"]])

#se1 <- SummarizedExperiment::SummarizedExperiment(assays = mcmat, rowRanges = aers, metadata = mcols(aers))

#test_juncs <- SummarizedExperiment::rowRanges(dasper::junctions_example)
test_grs <- GenomicRanges::GRanges(
  seqnames = S4Vectors::Rle(c("chr21", "chr22"), c(10, 1)),
  ranges = IRanges::IRanges(
    start = c(5306423, 24738, 2218, 5033895, 17554, 5080044, 3080053, 16570, 508007, 15005, 20312),
    end = c(5323718, 24891, 3407, 5033980, 17728, 5080091, 3080081, 16723, 508009, 15038, 20582),
    names = head(letters, 11)
  ),
  strand = S4Vectors::Rle((c("+", "-")), c(6, 5)),
  score = 1:11,
  GC = seq(1, 0, length = 11)
)
test_grs2 <- GenomicRanges::GRanges( # explain that ths shouldn't overlap
  seqnames = S4Vectors::Rle(c("chr1", "chr2"), c(10, 1)),
  ranges = IRanges::IRanges(
    start = c(5506469, 36738, 218, 8033895, 18554, 5180044, 2080053, 14570, 508007, 19005, 23312),
    end = c(5523703, 37891, 407, 8033980, 18728, 5180091, 2080081, 14723, 508009, 19038, 23582),
    names = head(letters, 11)
  ),
  strand = S4Vectors::Rle((c("+", "-")), c(6, 5))
)

test_annot_ers <- annotatERs(
  opt_ers = test_grs, junc_data = lung_junc,
  gtf_path = gtf_path, chrs_to_keep = c("chr21", "chr22")
)

test_annot_ers2 <- annotatERs(
  opt_ers = test_grs2, junc_data = lung_junc,
  gtf_path = gtf_path, chrs_to_keep = c("chr1", "chr2")
)

annot_bed1 <- rtracklayer::export.bed(test_annot_ers,"annotation1.bed")
annot_bed2 <- rtracklayer::export.bed(test_annot_ers2,"annotation2.bed")

chr_info <- get_chr_info(chrs = c("chr1","chr2","chr21","chr22"), genome = "hg38")
chrseqlen1 <- chr_info[["size"]][3:4]
names(chrseqlen1) <- c("chr21","chr22")
chrseqlen2 <- chr_info[["size"]][1:2]
names(chrseqlen2) <- c("chr1","chr2")


GenomeInfoDb::seqlengths(test_grs) <- chrseqlen1
GenomeInfoDb::seqlengths(test_grs2) <- chrseqlen2

big_file1 <- tempfile()

#test_bigwig <- rtracklayer::export(object = test_grs, format = "bigwig", con = "test1.bw") #, con = big_file1)
rtracklayer::export(object = test_grs, format = "bigwig", con = "test1.bw") #, con = big_file1)
#test_bigwig2 <- rtracklayer::export(object = test_grs2, format = "bigwig", con =  "test2.bw")
rtracklayer::export(object = test_grs2, format = "bigwig", con =  "test2.bw")

mean_coverage1 <- megadepth::get_coverage(bigwig_file = "test1.bw", op = "mean", annotation = annot_bed1)
mean_coverage2 <- megadepth::get_coverage(bigwig_file = "test2.bw", op = "mean", annotation = annot_bed2)

mcmat1 <- as.matrix(mcols(mean_coverage1)[["score"]])
mcmat2 <- as.matrix(mcols(mean_coverage2)[["score"]])

# mcmat <- as.matrix(mcols(mean_coverage)[["score"]])
# #se1 <- SummarizedExperiment::SummarizedExperiment(assays = mcmat, rowRanges = aers, metadata = mcols(aers))


# gtf_url <- "http://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz"
# gtf_path <- ODER:::.file_cache(gtf_url)
#https://storage.googleapis.com/gtex_analysis_v6p/rna_seq_data/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct.gz
gtex_url <- "https://storage.googleapis.com/gtex_analysis_v6p/rna_seq_data/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct.gz"
gtex_path <- ODER:::.file_cache(gtex_url)
library(data.table)
dt = fread(gtex_path)