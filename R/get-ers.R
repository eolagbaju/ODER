#' Define sets of ERs
#'
#' \code{get_ers} defines expressed regions across an inputted range of mean
#' coverage cut-offs (MCCs) and max region gaps (MRGs) from the coverage.
#'
#' @param coverage the coverage of the bigwig files passed into
#'   \code{\link{get_coverage}}.
#' @param mccs numeric vector containing the mean coverage cut-offs to apply.
#' MCCs are the mininmum number of reads that a base has to have to be
#' considered as expressed.
#' @param mrgs numeric vector containing the max region gaps to apply. MRGs are
#' the maximum number of bases between ERs that can be below the MCC and the
#' adjacent ERs will be merged.
#'
#' @return list containing sets of ERs, each generated using a particular
#'   combination of MCC and MRG.
#' @export
#'
#' @examples
#' # gtex_SRP012682_SRX222703_lung_coverage_1 is from the package data folder
#' eg_ers <- get_ers(
#'     coverage = gtex_SRP012682_SRX222703_lung_coverage_1,
#'     mccs = c(5, 10),
#'     mrgs = c(10, 20)
#' )
#'
#' eg_ers
get_ers <- function(coverage, mccs, mrgs) {
    er_entry_check(coverage, mccs, mrgs)
    ers <- list()
    for (j in seq_along(mccs)) { # getting the various MCC-MRG pairings ready
        mcc_label <- stringr::str_c("mcc_", mccs[j])
        for (k in seq_along(mrgs)) {
            mrg_label <- stringr::str_c("mrg_", mrgs[k])
            ers[[mcc_label]][[mrg_label]] <- list()
        }
    }
    for (i in seq_along(coverage)) { ##### Generate ERs #####
        message(stringr::str_c(
            Sys.time(), " - Generating ERs for ", names(coverage)[i]
        ))
        for (j in seq_along(mccs)) {
            mcc_label <- stringr::str_c("mcc_", mccs[j])
            er_mcc <- derfinder::findRegions( # generate ERs at particular mcc
                position = S4Vectors::Rle(
                    TRUE, length(coverage[[i]][["meanCoverage"]])
                ),
                fstats = coverage[[i]][["meanCoverage"]],
                chr = names(coverage)[i], cutoff = mccs[j], maxRegionGap = 0L,
                # setting maxClusterGap to chr length (ignore clusters)
                # to reduce run time. No impact on ERs
                maxClusterGap = length(coverage[[i]][["meanCoverage"]]),
                verbose = FALSE
            ) # collapse ERs with less than a MRG apart
            for (k in seq_along(mrgs)) {
                mrg_label <- stringr::str_c("mrg_", mrgs[k])
                ers[[mcc_label]][[mrg_label]][[
                names(coverage)[i]]] <- er_mcc %>%
                    GenomicRanges::reduce(min.gapwidth = mrgs[k])
            }
        }
    } ##### Merge ERs across chromosomes #####
    for (j in seq_along(mccs)) {
        mcc_label <- stringr::str_c("mcc_", mccs[j])
        for (k in seq_along(mrgs)) {
            mrg_label <- stringr::str_c("mrg_", mrgs[k])
            ers[[mcc_label]][[mrg_label]] <-
                ers[[mcc_label]][[mrg_label]] %>%
                GenomicRanges::GRangesList() %>%
                unlist() %>%
                sort()
            names(ers[[mcc_label]][[mrg_label]]) <- NULL
        } # enables conversion into dataframe otherwise all rows will have chr
    }
    return(ers)
}

#' Define sets of ERs for stranded Bigwigs
#'
#' \code{get_strand_ers} defines ERs across an inputted range of mean coverage
#' cut-offs (MCCs) and max region gaps (MRGs) from the coverage.
#'
#' @param bw_pos positive strand bigwig file
#' @param bw_neg negative strand bigwig file
#' @param auc_raw_pos vector containing AUCs(Area Under Coverage) matching the
#' order of the positive bigwig paths.
#' @param auc_raw_neg vector containing AUCs(Area Under Coverage) matching the
#' order of the negative bigwig paths.
#' @param auc_target total AUC to normalise all samples to. E.g. 40e6 * 100
#'   would be the estimated total auc for sample sequenced to 40 million reads
#'   of 100bp in length.
#' @param chrs chromosomes to obtain mean coverage for, default is "" giving
#'   every chromosome. Can take UCSC format(chrs = "chr1")
#'   or just the chromosome i.e. chrs = c(1,X)
#' @param mccs mean coverage cut-offs to apply.
#' @param mrgs max region gaps to apply.
#' @param bw_chr specifies whether the bigwig files has the chromosomes labelled
#'   with a "chr" preceding the chromosome i.e. "chr1" vs "1". Can be either
#'   "chr" or "nochr" with "chr" being the default.
#'
#' @return list containing sets of stranded ERs, each generated using a
#' particular combination of MCC and MRG.
#' @export
#'
#' @examples
#' library("magrittr")
#' if (!exists("gtex_metadata")) {
#'     gtex_metadata <- recount::all_metadata("gtex")
#'     gtex_metadata <- gtex_metadata %>%
#'         as.data.frame() %>%
#'         dplyr::filter(project == "SRP012682")
#' }
#' if (!exists("rec_url")) {
#'     rec_url <- recount::download_study(
#'         project = "SRP012682",
#'         type = "samples",
#'         download = FALSE
#'     )
#' }
#' # file_cache is an internal function to download a bigwig file from a link
#' # if the file has been downloaded recently, it will be retrieved from a cache
#' bw_plus <- file_cache(rec_url[58])
#' bw_minus <- file_cache(rec_url[84])
#'
#' # As of rtracklayer 1.25.16, BigWig is not supported on Windows.
#' if (!xfun::is_windows()) {
#'     stranded_ers <- get_strand_ers(
#'         bw_pos = bw_plus, bw_neg = bw_minus,
#'         auc_raw_pos = gtex_metadata[["auc"]][58],
#'         auc_raw_neg = gtex_metadata[["auc"]][84], auc_target = 40e6 * 100,
#'         chrs = "chr21", mccs = c(5, 10), mrgs = c(10, 20)
#'     )
#'     stranded_ers
#' }
get_strand_ers <- function(bw_pos, bw_neg, auc_raw_pos, auc_raw_neg, auc_target,
    chrs, mccs, mrgs, bw_chr = "chr") {
    plus_coverage <- get_coverage(
        bw_paths = bw_pos, auc_raw = auc_raw_pos,
        auc_target = auc_target, chrs = chrs,
        bw_chr = bw_chr
    )
    minus_coverage <- get_coverage(
        bw_paths = bw_neg, auc_raw = auc_raw_neg,
        auc_target = auc_target, chrs = chrs,
        bw_chr = bw_chr
    )

    ers_plus <- get_ers(coverage = plus_coverage, mccs = mccs, mrgs = mrgs)
    ers_minus <- get_ers(coverage = minus_coverage, mccs = mccs, mrgs = mrgs)

    sublist <- vector("list", length(ers_plus[[1]]))
    ers_combi <- vector("list", length(ers_plus))

    for (i in seq_along(ers_combi)) {
        ers_combi[[i]] <- sublist
    }
    # combining the positive and negative strands into one combined ER
    for (i in seq_along(ers_plus)) {
        names(ers_combi) <- names(ers_plus)
        for (j in seq_along(ers_plus[[i]])) {
            names(ers_combi[[i]]) <- names(ers_plus[[j]])
            BiocGenerics::strand(ers_plus[[i]][[j]]) <- "+"
            BiocGenerics::strand(ers_minus[[i]][[j]]) <- "-"
            ers_combi[[i]][[j]] <- c(ers_plus[[i]][[j]], ers_minus[[i]][[j]])
            GenomeInfoDb::sortSeqlevels(ers_combi[[i]][[j]])
            BiocGenerics::sort(ers_combi[[i]][[j]])
        }
    }

    return(ers_combi)
}

#' Checks for valid input toget_ers
#'
#' @param coverage bigwig coverage
#' @param mccs Mean Cutoff Coverages
#' @param mrgs Max Region Gaps
#'
#' @keywords internal
#' @noRd
er_entry_check <- function(coverage, mccs, mrgs) {
    if (missing(coverage)) {
        stop("Coverage is missing")
    } else if (missing(mccs)) {
        stop("Mean Coverage Cutoff is empty")
    } else if (missing(mrgs)) {
        stop("Max Region Gap is empty")
    }
}
