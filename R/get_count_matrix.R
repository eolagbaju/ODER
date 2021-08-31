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
#' if (!exists("gtex_metadata")) {
#'     gtex_metadata <- recount::all_metadata("gtex")
#'     gtex_metadata <- gtex_metadata %>%
#'         as.data.frame() %>%
#'         dplyr::filter(project == "SRP012682")
#' }
#' if (!exists("gtf_path")) {
#'     gtf_url <- "http://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz"
#'     gtf_path <- .file_cache(gtf_url)
#' }
#' if (!exists("rec_url")) {
#'     rec_url <- recount::download_study(
#'         project = "SRP012682",
#'         type = "samples",
#'         download = FALSE
#'     ) # .file_cache is an internal function to download a bigwig file from a link
#'     # if the file has been downloaded recently, it will be retrieved from a cache
#' }
#' bw_path <- .file_cache(rec_url[1])
#' bw_path2 <- .file_cache(rec_url[6])
#'
#' if (!exists("opt_ers1")) {
#'     opt_ers1 <- ODER(
#'         bw_paths = bw_path, auc_raw = auc_example,
#'         auc_target = 40e6 * 100, chrs = c("chr21"),
#'         genome = "hg38", mccs = c(5, 10), mrgs = c(10, 20),
#'         gtf = gtf_path, ucsc_chr = TRUE, ignore.strand = TRUE,
#'         exons_no_overlap = NULL, bw_chr = "chr"
#'     )
#' }
#' if (!exists("opt_ers6")) {
#'     opt_ers6 <- ODER(
#'         bw_paths = bw_path2, auc_raw = gtex_metadata[["auc"]][6],
#'         auc_target = 40e6 * 100, chrs = c("chr21"),
#'         genome = "hg38", mccs = c(5, 10), mrgs = c(10, 20),
#'         gtf = gtf_path, ucsc_chr = TRUE, ignore.strand = TRUE,
#'         exons_no_overlap = NULL, bw_chr = "chr"
#'     )
#' }
#' test_juncs <- lung_junc_21_22
#'
#' if (!exists("genom_state")) {
#'     genom_state <- generate_genomic_state(
#'         gtf = gtf_path,
#'         chrs_to_keep = c("21"), ensembl = TRUE
#'     )
#' }
#' if (!exists("aers1fcm")) {
#'     aers1fcm <- suppressWarnings(annotatERs(head(opt_ers1[["opt_ers"]], 500),
#'         junc_data = test_juncs,
#'         gtf_path = gtf_path, ensembl = TRUE, genom_state = genom_state
#'     ))
#' }
#' if (!exists("aers6fcm")) {
#'     aers6fcm <- suppressWarnings(annotatERs(head(opt_ers6[["opt_ers"]], 500),
#'         junc_data = test_juncs, genom_state = genom_state,
#'         gtf_path = gtf_path, ensembl = TRUE
#'     ))
#' }
#' }
#' example_bw_paths <- c(bw_path, bw_path2)
#' annot_ersl <- GenomicRanges::GRangesList(aers1fcm, aers6fcm)
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
