#' Generating the mean coverage of the expressed regions
#'
#' \code{get_coverage} returns the coverage of the BigWig data passed in
#'
#' @param bw_paths paths to bigwig files with the RNA-seq data that you want the coverage of.
#' @param auc_raw vector containing aucs matching the order of bigwig paths.
#' @param auc_target total auc to normalise all samples to. E.g. 40e6 * 100
#'   would be the estimated total auc for sample sequenced to 40 million reads
#'   of 100bp in length.
#' @param chrs chromosomes to obtain mean coverage for, default is "" giving every chromosome
#' @param genome the UCSC genome you want to use, the default is hg38
#'
#' @return a list of Rles per chromosome passed in
#' @export
#'
#' @examples
#' eg_coverage <- get_coverage(
#'     # 	bw_paths = bw_path,
#'     auc_raw = gtex_metadata[["auc"]][1],
#'     auc_target = 40e6 * 100,
#'     chrs = test_chrs
#' )
get_coverage <- function(bw_paths, auc_raw, auc_target, chrs = "", genome = "hg38") {
    if (!is.numeric(auc_raw) | !is.numeric(auc_target)) {
        stop("Please enter a valid number for the auc values")
    } else if (missing(bw_paths)) {
        stop("bw_paths is empty")
    } else if (bw_check(bw_paths) == FALSE) {
        stop("Please check your bigwig file paths are correct")
    }

    chr_info <- get_chr_info(chrs = chrs, genome = genome)

    all_chrs_mean_cov <- list()

    print(stringr::str_c(Sys.time(), " - Obtaining mean coverage across ", length(bw_paths), " samples"))

    for (i in 1:nrow(chr_info)) {
        print(stringr::str_c(Sys.time(), " - ", chr_info[["chrom"]][i]))
        # loading coverage information for designated chromosomes and merging them into a dataframe
        chr_mean_cov <-
            derfinder::loadCoverage(
                files = bw_paths,
                totalMapped = auc_raw, # normalise by auc here as for bws, more accurate since Rail-RNA clips reads
                targetSize = auc_target, # these come from derfinder::filterData
                chr = chr_info$chrom[i],
                chrlen = chr_info$size[i],
                inputType = "BigWig",
                returnMean = T,
                returnCoverage = F,
                verbose = F,
                cutoff = NULL
            ) # setting cutoff as NULL here and instead to be applied in findRegions()
        # storing the mean coverage in a list
        all_chrs_mean_cov[[chr_info[["chrom"]][i]]] <- chr_mean_cov["meanCoverage"]
    }
    return(all_chrs_mean_cov)
}

#' Get information from UCSC about the chromosomes passed in
#'
#' @param chrs chromosomes to look up (must match UCSC format)
#' @param genome the UCSC genome to look at
#'
#' @return a dataframe with data on the passed in chromosomes
#' @export
#'
#' @examples
#' eg_info <- get_chr_info(chrs = chrs, genome = "hg38")
get_chr_info <- function(chrs, genome) {
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

#' Define sets of ERs
#'
#' \code{gen_ERs} defines ERs across an inputted range of mean coverage cut-offs
#' (MCCs) and max region gaps (MRGs).
#'
#' @param coverage the coverage of the bigwig files passed into \code{\link{get_coverage}}.
#' @param mccs mean coverage cut-offs to apply.
#' @param mrgs max region gaps to apply.
#'
#' @return list containing sets of ERs, each generated using a particular
#'   combination of MCC and MRG.
#' @export
#'
#' @examples
#' eg_ers <- get_ers(coverage = coverage, mccs = c(5, 10), mrgs = c(10, 20))
get_ers <- function(coverage, mccs, mrgs) {
    if (missing(coverage)) {
        stop("Coverage is missing")
    } else if (missing(mccs)) {
        stop("Mean Coverage Cutoff is empty")
    } else if (missing(mrgs)) {
        stop("Max Region Gap is empty")
    }

    ers <- list()

    for (j in 1:length(mccs)) {
        mcc_label <- stringr::str_c("mcc_", mccs[j])

        for (k in 1:length(mrgs)) {
            mrg_label <- stringr::str_c("mrg_", mrgs[k])

            ers[[mcc_label]][[mrg_label]] <- list()
        }
    }

    ##### Generate ERs #####

    for (i in 1:length(coverage)) {
        print(stringr::str_c(Sys.time(), " - Generating ERs for ", names(coverage)[i]))

        for (j in 1:length(mccs)) {
            mcc_label <- stringr::str_c("mcc_", mccs[j])

            # generate ERs at particular mcc
            suppressMessages(
                er_mcc <- derfinder::findRegions(
                    position = S4Vectors::Rle(TRUE, length(coverage[[i]][["meanCoverage"]])),
                    fstats = coverage[[i]][["meanCoverage"]],
                    chr = names(coverage)[i],
                    cutoff = mccs[j],
                    maxRegionGap = 0L,
                    maxClusterGap = length(coverage[[i]][["meanCoverage"]]), # setting this to chr length (ignore clusters) to reduce run time. No impact on ERs
                    verbose = F
                )
            )

            for (k in 1:length(mrgs)) {
                mrg_label <- stringr::str_c("mrg_", mrgs[k])

                # collapse ERs with less than a MRG apart
                ers[[mcc_label]][[mrg_label]][[names(coverage)[i]]] <- er_mcc %>%
                    GenomicRanges::reduce(min.gapwidth = mrgs[k])
            }
        }
    }

    ##### Merge ERs across chromosomes #####

    for (j in 1:length(mccs)) {
        mcc_label <- stringr::str_c("mcc_", mccs[j])

        for (k in 1:length(mrgs)) {
            mrg_label <- stringr::str_c("mrg_", mrgs[k])

            ers[[mcc_label]][[mrg_label]] <-
                ers[[mcc_label]][[mrg_label]] %>%
                GenomicRanges::GRangesList() %>%
                unlist() %>%
                sort()
        }
    }

    return(ers)
}
