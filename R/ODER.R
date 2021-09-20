#' Generates the optimum expressed regions
#'
#' Returns the optimum definition of the expressed regions by finding the ideal
#' MCC (Mean Coverage Cutoff) and MRG (Max Region Gap). The combination of MCC
#' and MRG that returns the expressed region with the smallest exon delta is the
#' most ideal.
#'
#' @param exons_no_overlap Optimum set of exons to help calculate deltas
#' @param file_type Describes if the BigWigs are stranded or not. Either
#' "stranded" or non-stranded
#' @inheritParams get_coverage
#' @inheritParams get_ers
#' @inheritParams get_exons
#' @inheritParams get_strand_ers
#'
#' @return list containing optimised ERs, optimal pair of MCC/MRGs and
#' \code{delta_df}
#' @export
#'
#' @examples
#' \dontshow{
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
#' if (!exists("gtf_path")) {
#'     gtf_url <- paste0(
#'         "http://ftp.ensembl.org/pub/release-103/gtf/",
#'         "homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz"
#'     )
#'     gtf_path <- file_cache(gtf_url)
#' }
#' # As of rtracklayer 1.25.16, BigWig is not supported on Windows.
#' }
#' if (!xfun::is_windows()) {
#'     opt_ers <- ODER(
#'         bw_paths = bw_path, auc_raw = gtex_lung_auc_1,
#'         auc_target = 40e6 * 100, chrs = c("chr21", "chr22"),
#'         genome = "hg38", mccs = c(5, 10), mrgs = c(10, 20),
#'         gtf = gtf_path, ucsc_chr = TRUE, ignore.strand = TRUE,
#'         exons_no_overlap = NULL, bw_chr = "chr"
#'     )
#'
#'     opt_ers
#' }
ODER <- function(bw_paths, auc_raw, auc_target, chrs = "", genome = "hg38",
    mccs, mrgs, gtf = NULL, ucsc_chr, ignore.strand, exons_no_overlap = NULL,
    biotype = "Non-overlapping", bw_chr = "chr", file_type = "non-stranded",
    bw_pos = NULL, bw_neg = NULL, auc_raw_pos = NULL, auc_raw_neg = NULL) {
    if (is.null(gtf)) stop("gtf must be provided.")

    if (file_type == "non-stranded") {
        coverage <- get_coverage(
            bw_paths = bw_paths, auc_raw = auc_raw, auc_target = auc_target,
            chrs = chrs, genome = "hg38", bw_chr = bw_chr
        )
        ers <- get_ers(coverage = coverage, mccs = mccs, mrgs = mrgs)
        if (!is.null(gtf)) {
            exons_no_overlap <- get_exons(
                gtf = gtf, ucsc_chr = ucsc_chr,
                ignore.strand = ignore.strand, biotype = biotype
            )
        }
        ers_delta <- get_ers_delta(ers = ers, opt_exons = exons_no_overlap)

        opt_ers <- get_opt_ers(ers = ers, ers_delta = ers_delta)

        return(opt_ers)
    } else if (file_type == "stranded") {
        stranded_ers <- get_strand_ers(
            bw_pos = bw_pos, bw_neg = bw_neg, auc_raw_pos = auc_raw_pos, auc_raw_neg = auc_raw_neg,
            auc_target = auc_target, chrs = chrs, mccs = mccs, mrgs = mrgs
        )
        if (!is.null(gtf)) {
            exons_no_overlap <- get_exons(
                gtf = gtf,
                ucsc_chr = ucsc_chr,
                ignore.strand = ignore.strand
            )
        }
        ers_delta <- get_ers_delta(ers = stranded_ers, opt_exons = exons_no_overlap)

        opt_ers <- get_opt_ers(ers = stranded_ers, ers_delta = ers_delta)

        return(opt_ers)
    }
}
