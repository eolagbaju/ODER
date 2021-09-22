#' Generate the count matrix
#'
#' Scores the mean coverage of the expressed regions as a count matrix
#'
#' @param bw_paths Vector containing the bigwig file paths to read in
#' @param annot_ers GRangesList containing the annotated ERs
#' (product of \code{annotatERs})
#' @param cols A dataframe containing the information to be used as colData for
#'  the output. If NULL then the bw_paths will be used for the colData
#'
#' @return A Ranged Summarized Experiment containing the gene counts as an assay
#' @export
#'
#' @examples
#' megadepth::install_megadepth()
#'
#' if (!exists("rec_url")) {
#'     rec_url <- recount::download_study(
#'         project = "SRP012682",
#'         type = "samples",
#'         download = FALSE
#'     )
#' }
#' # file_cache is an internal function to download a bigwig file from a link
#' # if the file has been downloaded recently, it will be retrieved from a cache
#' bw_path <- file_cache(rec_url[1])
#'
#' ex_opt_ers <- GenomicRanges::GRanges(
#'     seqnames = S4Vectors::Rle(c("chr1", "chr2"), c(4, 1)),
#'     ranges = IRanges::IRanges(
#'         start = c(1:5),
#'         end = seq(100, 500, 100)
#'     )
#' )
#'
#' example_cm <- get_count_matrix(
#'     bw_paths = c(bw_path, bw_path),
#'     annot_ers = ex_opt_ers
#' )
#' example_cm
get_count_matrix <- function(bw_paths, annot_ers, cols = NULL) {
    if (is.null(cols)) {
        cols <- data.frame(BigWig_paths = bw_paths)
    }

    if (!methods::is(annot_ers, "GRanges")) {
        stop("annot_ers must be Granges")
    }
    # empty matrix to store the mean coverage for each bigwig
    gene_counts <- matrix(nrow = length(annot_ers), ncol = length(bw_paths))
    # converting the annotated ERs into a bed file
    fil <- tempfile("er.bed")
    er_bed <- rtracklayer::export.bed(object = annot_ers, con = fil)

    for (i in seq_along(bw_paths)) {
        mean_coverage <- megadepth::get_coverage(
            bigwig_file = bw_paths[i],
            op = "mean",
            annotation = er_bed
        )
        gene_counts[, i] <- S4Vectors::mcols(mean_coverage)[["score"]]
    }


    se <- SummarizedExperiment::SummarizedExperiment(
        assays = gene_counts,
        rowRanges = annot_ers,
        colData = cols
    )

    return(se)
}
