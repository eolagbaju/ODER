#' Generates the optimum expressed regions
#'
#' @param exons_no_overlap Optimum set of exons to help calculate deltas
#' @inheritParams get_coverage
#' @inheritParams get_ers
#' @inheritParams get_exons
#' @inheritParams get_coverage
#'
#' @return list containing optimised ERs, optimal pair of MCC/MRGs and
#' \code{delta_df}
#' @export
#'
#' @examples
#' \dontshow{
#' url <- recount::download_study(
#'     project = "SRP012682",
#'     type = "samples",
#'     download = FALSE
#' ) # .file_cache is an internal function to download a bigwig file from a link
#' # if the file has been downloaded recently, it will be retrieved from a cache
#'
#' bw_path <- ODER:::.file_cache(url[1])
#' gtf_url <- paste0(
#'     "http://ftp.ensembl.org/pub/release-103/gtf/",
#'     "homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz"
#' )
#' gtf_path <- ODER:::.file_cache(gtf_url)
#' }
#'
#' opt_ers <- ODER(
#'     bw_paths = bw_path, auc_raw = auc_example,
#'     auc_target = 40e6 * 100, chrs = c("chr21", "chr22"),
#'     genome = "hg38", mccs = c(5, 10), mrgs = c(10, 20),
#'     gtf = gtf_path, ucsc_chr = TRUE, ignore.strand = TRUE,
#'     exons_no_overlap = NULL, bw_chr = "chr"
#' )
#'
#' opt_ers
ODER <- function(bw_paths, auc_raw, auc_target, chrs = "", genome = "hg38",
    mccs, mrgs,
    gtf = NULL, ucsc_chr, ignore.strand,
    exons_no_overlap = NULL,
    bw_chr = "chr") {
    if (is.null(gtf) && is.null(opt_ers)) stop("One of gtf OR opt_gr must be provided.")

    coverage <- get_coverage(
        bw_paths = bw_paths, auc_raw = auc_raw, auc_target = auc_target,
        chrs = chrs, genome = "hg38", bw_chr = bw_chr
    )

    ers <- get_ers(coverage = coverage, mccs = mccs, mrgs = mrgs)

    if (!is.null(gtf)) {
        exons_no_overlap <- get_exons(
            gtf = gtf,
            ucsc_chr = ucsc_chr,
            ignore.strand = ignore.strand
        )
    }

    ers_delta <- get_ers_delta(ers = ers, opt_exons = exons_no_overlap)

    opt_ers <- get_opt_ers(ers = ers, ers_delta = ers_delta)

    return(opt_ers)
}
