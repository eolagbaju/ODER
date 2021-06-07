#' Get junctions that overlap the optimally defined ERs
#'
#' @param opt_ers optimally defined ERs (the product of the ODER function)
#' @param junc_data junction data that should match the ERs passed into opt_ers
#' @param gtf_path a gtf file with exon data
#'
#' @return optimally defined ers annotated with junction and gene information
#' @export
#'
#' @examples
#'
#' #'
#' \dontshow{
#' gtf_url <- paste0(
#'     "http://ftp.ensembl.org/pub/release-103/gtf/",
#'     "homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz"
#' )
#' # .file_cache is an internal function to download a bigwig file from a link
#' # if the file has been downloaded recently, it will be retrieved from a cache
#' gtf_path <- ODER:::.file_cache(gtf_url)
#' }
#'
#' example_ers <- GenomicRanges::GRanges(
#'     seqnames = S4Vectors::Rle(c("chr21", "chr22"), c(10, 1)),
#'     ranges = IRanges::IRanges(
#'         start = c(
#'             5026423, 24738, 5032218, 5033895, 17554, 50800446,
#'             50800539, 16570, 50800790, 15005, 20312
#'         ),
#'         end = c(
#'             5323718, 24891, 5033407, 5033980, 17728, 50800910,
#'             50800817, 16723, 50800910, 15038, 20582
#'         ),
#'         names = head(letters, 11)
#'     ),
#'     strand = S4Vectors::Rle((c("+", "-")), c(6, 5)),
#'     score = 1:11,
#'     GC = seq(1, 0, length = 11)
#' )
#'
#' example_junctions <- SummarizedExperiment::rowRanges(dasper::junctions_example)
#'
#' example_er_juncs <- get_junctions(
#'     opt_ers = example_ers,
#'     junc_data = example_junctions,
#'     gtf_path = gtf_path
#' )
#'
#' print(example_er_juncs)
get_junctions <- function(opt_ers, junc_data, gtf_path) {
    if (methods::is(junc_data, "data.frame")) {
        junc_data <- GenomicRanges::makeGRangesFromDataFrame(junc_data)
    }
    annotated_junctions <- dasper::junction_annot(junctions = junc_data, ref = gtf_path)
    GenomeInfoDb::seqlevelsStyle(annotated_junctions) <- "UCSC" # to match opt_ers

    hits <- GenomicRanges::findOverlaps(opt_ers, annotated_junctions)
    ann_junc_hits <- annotated_junctions[S4Vectors::subjectHits(hits)]
    ann_junc_hits$er_index <- S4Vectors::queryHits(hits)

    gene_id_list <- GenomicRanges::mcols(ann_junc_hits)[["gene_id_junction"]] %>%
        S4Vectors::split(f = ann_junc_hits$er_index)

    gene_id_list <- lapply(gene_id_list, function(x) {
        return(unique(unlist(x)))
    })

    ann_junc_hits <- ann_junc_hits %>% S4Vectors::split(f = ann_junc_hits$er_index)

    er_indices <- unique(S4Vectors::queryHits(hits))

    miss_ers <- numeric(0)
    for (i in 1:length(opt_ers)) {
        if (!(i %in% er_indices)) {
            miss_ers <- c(miss_ers, i)
        }
    }

    empty_grl <- GenomicRanges::GRangesList(
        lapply(miss_ers, function(x) {
            return(GenomicRanges::GRanges())
        })
    )
    names(empty_grl) <- miss_ers
    combi_ajh <- suppressWarnings(c(ann_junc_hits, empty_grl))
    sorted_combi_ajh <- combi_ajh[order(as.integer(names(combi_ajh)))]

    empty_gil <- vector("list", length = length(miss_ers))
    names(empty_gil) <- miss_ers

    combi_gil <- c(gene_id_list, empty_gil)
    sorted_combi_gil <- combi_gil[order(as.integer(names(combi_gil)))]
    sorted_combi_gil <- sorted_combi_gil %>% IRanges::CharacterList()

    GenomicRanges::mcols(opt_ers)$grl <- sorted_combi_ajh
    GenomicRanges::mcols(opt_ers)$genes <- sorted_combi_gil

    return(opt_ers)
}











# find_intersect <- function(erstart,erend,junctions){
#
#   for (i in 1:length(ranges(junctions))){
#     jlist <- list()
#     if (erstart<=start(ranges(junctions)[i]) & erend>=start(ranges(junctions)[i])){
#       append(jlist,stringr::str_c(as.character(start(ranges(junctions)[i])),
#                                   " - ",
#                                   as.character(end(ranges(junctions)[i])))
#       )
#     }else if (erstart<=end(ranges(junctions)[i]) & erend>=end(ranges(junctions)[i])){
#       append(jlist,stringr::str_c(as.character(start(ranges(junctions)[i])),
#                                   " - ",
#                                   as.character(end(ranges(junctions)[i])))
#       )
#     }
#   }
#
#   return(jlist)
# }
#
# .gr_convert <- function(bplist){
#   #str_extract(er_junc_hits[["junction_coords"]][[1]][[1]],":(.*?):") - gets ranges
#   #str_extract(er_junc_hits[["junction_coords"]][[1]][[1]],"(.*?):") - gets chr/seqname
#   #str_sub(er_junc_hits[["junction_coords"]][[1]][[1]],-1,-1) - get the strand
#   grlist <- GenomicRanges::GRangesList()
#
#   for (i in seq_along(bplist)){
#     seqname <- gsub(":","",str_extract(bplist[i],"(.*?):"))
#     ranges <- gsub(":","",str_extract(bplist[i],":(.*?):"))
#     strand <- str_sub(bplist[i],-1,-1)
#     gr <- GenomicRanges::GRanges(
#       seqnames = seqname,
#       ranges = ranges,
#       strand <- strand
#     )
#
#     grlist[[i]] <- gr
#   }
#
#   return(grlist)
#
# }
#
# extract_gene_ids <- function(grange){
#   genes <- unlist(mcols(grange)[["gene_id_junction"]])
#   return(genes)
# }
#
# gene_id_extract <- function(grange){
#   return(unique(unlist(mcols(grange)[["gene_id_junction"]])))
# }
