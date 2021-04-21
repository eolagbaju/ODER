
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
#' # coverage_example <- load(file = "data/coverage_example.rda")
#'
#' eg_ers <- get_ers(coverage = coverage_example, mccs = c(5, 10), mrgs = c(10, 20))
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
