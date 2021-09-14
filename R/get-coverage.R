#' Obtain the mean coverage across multiple BigWig files
#'
#' \code{get_coverage} returns the mean coverage of the BigWig files passed in.
#' Internally, this operates through `derfinder::loadCoverage`.
#'
#' @param bw_paths path(s) to bigwig file(s) with the RNA-seq data that you want
#'   the #'   coverage of.
#' @param auc_raw vector containing AUCs(Area Under Coverage) matching the order
#'   of bigwig path(s).
#' @param auc_target total AUC to normalise all samples to e.g. 40e6 * 100 would
#'   be the estimated total auc for sample sequenced to 40 million reads of
#'   100bp in length.
#' @param chrs chromosomes to obtain mean coverage for, default is "" giving
#'   every chromosome. Can take UCSC format(chrs = "chr1") or just the
#'   chromosome i.e. chrs = c(1,X)
#' @param genome the UCSC genome you want to use, the default is hg38.
#' @param bw_chr specifies whether the bigwig files has the chromosomes labelled
#'   with a "chr" preceding the chromosome i.e. "chr1" vs "1". Can be either
#'   "chr" or "nochr" with "chr" being the default.
#'
#' @return a list of Rles detailing the mean coverage per chromosome passed in.
#' @export
#'
#' @examples
#' \dontshow{
#' if (!exists("rec_url")) {
#'     rec_url <- recount::download_study(
#'         project = "SRP012682",
#'         type = "samples",
#'         download = FALSE
#'     ) # .file_cache is an internal function to download a bigwig file from a link
#'     # if the file has been downloaded recently, it will be retrieved from a cache
#'
#'     bw_path <- .file_cache(rec_url[1])
#' }
#' }
#' # As of rtracklayer 1.25.16, BigWig is not supported on Windows.
#' if (!xfun::is_windows()) {
#'     eg_coverage <- get_coverage(
#'         bw_paths = bw_path,
#'         auc_raw = 11872688252,
#'         auc_target = 40e6 * 100,
#'         chrs = c("chr21", "chr22")
#'     )
#'     print(eg_coverage)
#' }
get_coverage <- function(bw_paths, auc_raw, auc_target, chrs = "", genome = "hg38", bw_chr = "chr") {
    if (!is.numeric(auc_raw) | !is.numeric(auc_target)) {
        stop("Please enter a valid number for the auc values")
    } else if (missing(bw_paths)) {
        stop("bw_paths is empty")
    } else if (bw_check(bw_paths) == FALSE) {
        stop("Please check your bigwig file paths are correct")
    }

    chrs <- informatting(chrs)

    chr_info <- get_chr_info(chrs = chrs, genome = genome)

    all_chrs_mean_cov <- list()

    print(stringr::str_c(Sys.time(), " - Obtaining mean coverage across ", length(bw_paths), " samples"))

    for (i in 1:seq_along(chr_info)) {
        print(stringr::str_c(Sys.time(), " - ", chr_info[["chrom"]][i]))
        # loading coverage information for designated chromosomes and merging them into a dataframe
        chr_mean_cov <-
            derfinder::loadCoverage(
                files = bw_paths,
                totalMapped = auc_raw, # normalise by auc here as for bws, more accurate since Rail-RNA clips reads
                targetSize = auc_target, # these come from derfinder::filterData
                chr = chr_formatting(chr_info$chrom[i], bw_chr),
                chrlen = chr_info$size[i],
                inputType = "BigWig",
                returnMean = TRUE,
                returnCoverage = FALSE,
                verbose = FALSE,
                cutoff = NULL
            ) # setting cutoff as NULL here and instead to be applied in findRegions()
        # storing the mean coverage in a list
        all_chrs_mean_cov[[chr_info[["chrom"]][i]]] <- chr_mean_cov["meanCoverage"]
    }
    return(all_chrs_mean_cov)
}

#' Get information from UCSC about the chromosomes passed in
#'
#' Download information about each of the chromosomes passed in, most importantly
#' the size.
#'
#' @param chrs chromosomes to look up (must match UCSC format)
#' @param genome the UCSC genome to look at see \url{https://genome.ucsc.edu/}.
#'
#' @return a dataframe with data on the passed in chromosomes
#' @export
#'
#' @examples
#' eg_info <- get_chr_info(chrs = c("chr21", "chr22"), genome = "hg38")
#'
#' print(eg_info)
get_chr_info <- function(chrs, genome) {
    chrs <- informatting(chrs)

    all_UCSC_chr <- GenomeInfoDb::getChromInfoFromUCSC(genome)[["chrom"]]
    if (all_UCSC_chr %contain% chrs) {
        chr_info <- GenomeInfoDb::getChromInfoFromUCSC(genome) %>%
            dplyr::filter(chrom %in% chrs)
    } else if (chrs == "") {
        chr_info <- GenomeInfoDb::getChromInfoFromUCSC(genome)
    } else {
        stop("Non-UCSC chromosome format was entered")
    }

    return(chr_info)
}
