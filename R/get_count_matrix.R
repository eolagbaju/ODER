#' Generate the count matrix
#' 
#' Scores the mean coverage of the expressed regions as a count matrix
#'
#' @param bw_paths Vector containing the bigwig file paths to read in
#' @param annot_ers GRangesList containing the annotated ERs (product of annotatERs())
#' @param cols A dataframe containing the information to be used as colData for the
#' output. If NULL then the bw_paths will be used for the colData
#'
#' @return A Ranged Summarized Experiment containing the gene counts as an assay
#' @export
#'
#' @examples
#' \dontshow{
#' library(magrittr)
#' gtex_metadata <- recount::all_metadata("gtex")
#' gtex_metadata <- gtex_metadata %>%
#'     as.data.frame() %>%
#'     dplyr::filter(project == "SRP012682")
#'
#' gtf_url <- "http://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz"
#' gtf_path <- ODER:::.file_cache(gtf_url)
#'
#' url <- recount::download_study(
#'     project = "SRP012682",
#'     type = "samples",
#'     download = FALSE
#' ) # .file_cache is an internal function to download a bigwig file from a link
#' # if the file has been downloaded recently, it will be retrieved from a cache
#'
#' bw_path <- ODER:::.file_cache(url[84])
#' bw_path2 <- ODER:::.file_cache(url[91])
#'
#' opt_ers <- suppressWarnings(ODER(
#'     bw_paths = bw_path, auc_raw = gtex_metadata[["auc"]][84],
#'     auc_target = 40e6 * 100, chrs = c("chr21", "chr22"),
#'     genome = "hg38", mccs = c(2, 4, 6, 8, 10), mrgs = c(10, 20, 30),
#'     gtf = gtf_path, ucsc_chr = TRUE, ignore.strand = TRUE,
#'     exons_no_overlap = NULL, bw_chr = "chr"
#' ))
#'
#' opt_ers2 <- suppressWarnings(ODER(
#'     bw_paths = bw_path, auc_raw = gtex_metadata[["auc"]][91],
#'     auc_target = 40e6 * 100, chrs = c("chr21", "chr22"),
#'     genome = "hg38", mccs = c(2, 4, 6, 8, 10), mrgs = c(10, 20, 30),
#'     gtf = gtf_path, ucsc_chr = TRUE, ignore.strand = TRUE,
#'     exons_no_overlap = NULL, bw_chr = "chr"
#' ))
#'
#' test_juncs <- SummarizedExperiment::rowRanges(dasper::junctions_example)
#' aers <- suppressWarnings(annotatERs(opt_ers[["opt_ers"]],
#'     junc_data = test_juncs,
#'     gtf_path = gtf_path, ensembl = TRUE
#' ))
#' aers2 <- suppressWarnings(annotatERs(opt_ers2[["opt_ers"]],
#'     junc_data = test_juncs,
#'     gtf_path = gtf_path, ensembl = TRUE
#' ))
#' }
#' example_bw_paths <- c(bw_path, bw_path2)
#' annot_ersl <- GenomicRanges::GRangesList(aers, aers2)
#'
#' example_cm <- get_count_matrix(bw_paths = example_bw_paths, annot_ers = annot_ersl)
#' example_cm
get_count_matrix <- function(bw_paths, annot_ers, cols = NULL) {
    megadepth::install_megadepth()

    if (is.null(cols)) {
        cols <- as.data.frame(bw_paths)
    }

    if (methods::is(annot_ers, "GRanges")) {
        rranges <- GenomicRanges::granges(annot_ers)
        fil <- tempfile("annotation.bed")
        annot_bed <- rtracklayer::export.bed(annot_ers, fil)
        mean_coverage <- megadepth::get_coverage(bigwig_file = bw_paths, op = "mean", annotation = annot_bed)
        gene_counts <- as.matrix(S4Vectors::mcols(mean_coverage)[["score"]])
        se <- SummarizedExperiment::SummarizedExperiment(
            assays = gene_counts,
            rowRanges = rranges,
            colData = cols
        )

        return(se)
    }

    rranges <- GenomicRanges::granges(annot_ers[[1]])
    gene_counts <- matrix(nrow = length(annot_ers[[1]]), ncol = length(bw_paths))

    for (i in seq_along(bw_paths)) {
        fil <- tempfile(stringr::str_c("annotation", as.character(i), ".bed"))
        annot_bed <- rtracklayer::export.bed(annot_ers[[i]], fil)
        mean_coverage <- megadepth::get_coverage(bigwig_file = bw_paths[i], op = "mean", annotation = annot_bed)
        gene_counts[, i] <- S4Vectors::mcols(mean_coverage)[["score"]]
    }


    se <- SummarizedExperiment::SummarizedExperiment(
        assays = gene_counts,
        rowRanges = rranges,
        colData = cols
    )

    return(se)
}
