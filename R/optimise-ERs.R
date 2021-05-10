#' Obtain set of non-overlapping exons
#'
#' @param gtf Either a string containg the path to a .gtf file or a pre-imported
#'   gtf using \code{\link[rtracklayer]{import}}.
#' @param ucsc_chr logical scalar, determining whether to add "chr" prefix to
#'   the seqnames of non-overlapping exons and change "chrMT" -> "chrM". Note,
#'   if set to TRUE and seqnames already have "chr", it will not add another.
#' @param ignore.strand logical value for input into
#'   \code{\link[GenomicRanges]{findOverlaps}}, default is True.
#'
#' @return GRanges object containing non-overlapping exons.
#' @export
#'
#' @examples
#' gtf_url <- paste0(
#'     "http://ftp.ensembl.org/pub/release-103/gtf/",
#'     "homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz"
#' )
#' gtf_path <- ODER:::.file_cache(gtf_url)
#' eg_exons_no_overlap <- get_exons(
#'     gtf = gtf_path,
#'     ucsc_chr = TRUE,
#'     ignore.strand = TRUE
#' )
#'
#' print(eg_exons_no_overlap)
get_exons <- function(gtf, ucsc_chr, ignore.strand = TRUE) {
    if (is.character(gtf)) {
        if (!xor(
            stringr::str_sub(gtf, -4, -1) == ".gtf",
            stringr::str_sub(gtf, -7, -1) == ".gtf.gz"
        )) {
            stop("Please check your gtf file path")
        }
        print(stringr::str_c(Sys.time(), " - Loading in GTF..."))

        gtf_gr <- rtracklayer::import(gtf)
    } else {
        gtf_gr <- gtf
    }

    print(stringr::str_c(Sys.time(), " - Obtaining non-overlapping exons"))

    exons_gr <- gtf_gr[gtf_gr$type == "exon"]
    exons_gr <- exons_gr[!duplicated(exons_gr$exon_id)]

    exons_hits <- GenomicRanges::findOverlaps(exons_gr,
        drop.self = TRUE,
        ignore.strand = ignore.strand
    )

    exons_no_overlap_gr <- exons_gr[-c(S4Vectors::queryHits(exons_hits) %>% unique())]

    # check - no overlaps

    if (ucsc_chr) {
        GenomeInfoDb::seqlevels(exons_no_overlap_gr) <-
            GenomeInfoDb::seqlevels(exons_no_overlap_gr) %>%
            stringr::str_replace("chr", "") %>%
            stringr::str_c("chr", .) %>%
            stringr::str_replace("chrMT", "chrM")
    }

    return(exons_no_overlap_gr)
}


#' Calculates delta for sets of ERs
#'
#' @param ers Sets of ERs across various MCCs/MRGs - output of
#' \code{\link{get_ers}}.
#' @param opt_exons GRanges object that contains the regions that ideally, you
#' want to the ER definitions to match - output of \code{\link{get_exons}}.
#' @param delta_fun Function that calculates the delta between ERs and
#'   \code{opt_exons}. Takes as input a set of ERs from \code{ers} and
#'   \code{opt_exons}. Then outputs a tibble/dataframe containing the summarised
#'   delta scores for that set of one set of ERs.
#'
#' @return tibble/dataframe containing summarised delta values. One row per set
#'   of ERs.
#' @export
#'
#' @examples
#' gtf_url <- paste0(
#'     "http://ftp.ensembl.org/pub/release-103/gtf/",
#'     "homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz"
#' )
#' gtf_path <- ODER:::.file_cache(gtf_url)
#'
#' eg_opt_exons <- get_exons(
#'     gtf = gtf_path,
#'     ucsc_chr = TRUE,
#'     ignore.strand = TRUE
#' )
#'
#' eg_ers_delta <- get_ers_delta(
#'     ers = ers_example, # ers_example is from the package data folder
#'     opt_exons = eg_opt_exons,
#'     delta_fun = ODER:::.delta
#' ) # .delta is ODER's default, you can pass in your own if you have one
#'
#' print(eg_ers_delta)
get_ers_delta <- function(ers, opt_exons, delta_fun = ODER:::.delta) {
    if (missing(ers)) {
        stop("No ERs were entered")
    } else if (missing(opt_exons)) {
        stop("No opt_exons were entered")
    }

    print(stringr::str_c(Sys.time(), " - Calculating delta for ERs..."))

    mcc_labels <- names(ers)

    delta_df <- dplyr::tibble()

    for (i in 1:length(mcc_labels)) {
        mrg_labels <- names(ers[[mcc_labels[i]]])

        for (j in 1:length(mrg_labels)) {
            delta_summarised <-
                delta_fun(
                    query = ers[[mcc_labels[i]]][[mrg_labels[j]]],
                    subject = opt_exons
                )

            delta_df <- delta_df %>%
                dplyr::bind_rows(delta_summarised %>%
                    dplyr::mutate(
                        mcc = stringr::str_remove(
                            mcc_labels[i],
                            stringr::fixed("mcc_")
                        ),
                        mrg = stringr::str_remove(
                            mrg_labels[j],
                            stringr::fixed("mrg_")
                        )
                    ))
        }
    }

    delta_df <- delta_df %>%
        dplyr::select(mcc, mrg, dplyr::everything())

    delta_df[["mcc"]] <- as.numeric(as.character(delta_df[["mcc"]]))
    delta_df[["mrg"]] <- as.numeric(as.character(delta_df[["mrg"]]))

    return(delta_df)
}


#' Obtains optimised set of ERs
#'
#' @param ers Sets of ERs across various MCCs/MRGs - output of
#' \code{\link{get_ers}}.
#' @param ers_delta tibble/dataframe containing summarised delta values. One row per set
#'   of ERs.
#'
#' @return list containing optimised ERs, optimal pair of MCC/MRGs and
#' \code{delta_df}
#' @export
#'
#' @examples
#'
#' gtf_url <- paste0(
#'     "http://ftp.ensembl.org/pub/release-103/gtf/",
#'     "homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz"
#' )
#' gtf_path <- ODER:::.file_cache(gtf_url)
#' exons_no_overlap <- get_exons(
#'     gtf = gtf_path,
#'     ucsc_chr = TRUE,
#'     ignore.strand = TRUE
#' )
#'
#' ers_delta <- get_ers_delta(
#'     ers = ers_example, # ers_example is from the package data folder
#'     opt_exons = exons_no_overlap
#' )
#'
#' opt_ers <- get_opt_ers(
#'     ers = ers_example,
#'     ers_delta = ers_delta
#' )
get_opt_ers <- function(ers, ers_delta) {
    if (missing(ers)) {
        stop("No ERs were entered")
    } else if (missing(ers_delta)) {
        stop("No ers_delta were entered")
    }

    print(stringr::str_c(Sys.time(), " - Obtaining optimal set of ERs..."))

    delta_opt <-
        ers_delta %>%
        dplyr::filter(median == min(median)) %>% # with the lowest median ER delta
        dplyr::filter(n_eq_0 == max(n_eq_0)) # and highest num of delta equal to 0

    mcc_label <- stringr::str_c("mcc_", as.character(delta_opt[["mcc"]]))
    mrg_label <- stringr::str_c("mrg_", as.character(delta_opt[["mrg"]]))

    opt_ers <-
        list(
            opt_ers = ers[[mcc_label]][[mrg_label]],
            opt_mcc_mrg = c(mcc_label, mrg_label),
            deltas = ers_delta
        )

    return(opt_ers)
}


#' Default delta functions
#'
#' @param query Set of ERs
#' @param subject Optimum exons
#'
#' @return summarised delta scores
#'
#' @keywords internal
#' @noRd
.delta <- function(query, subject) {
    # finding ovelaps of exons and expressed regions
    hits <- GenomicRanges::findOverlaps(query = query, subject = subject)

    # obtain situation where 1 ER overlaps multiple exons...
    n_dis_exons_ab_1 <-
        hits %>%
        as.data.frame() %>%
        dplyr::group_by(queryHits) %>%
        dplyr::summarise(n_dis_exons = dplyr::n_distinct(subjectHits)) %>%
        dplyr::filter(n_dis_exons > 1)

    # ...and remove them
    hits <- hits[!(S4Vectors::queryHits(hits) %in% n_dis_exons_ab_1$queryHits)]

    delta_raw <-
        dplyr::bind_cols(
            query[S4Vectors::queryHits(hits)] %>%
                as.data.frame(row.names = NULL) %>%
                dplyr::select(seqnames, start, end),
            subject[S4Vectors::subjectHits(hits)] %>%
                as.data.frame(row.names = NULL) %>%
                dplyr::select(seqnames1 = seqnames, start1 = start, end1 = end)
        ) %>%
        dplyr::mutate(
            start_diff = start - start1,
            end_diff = end - end1,
            delta = abs(start_diff) + abs(end_diff)
        )

    delta_summarised <-
        dplyr::tibble(
            sum = sum(delta_raw$delta),
            mean = mean(delta_raw$delta),
            median = median(delta_raw$delta),
            n_eq_0 = sum(delta_raw$delta == 0),
            propor_eq_0 = mean(delta_raw$delta == 0)
        )

    return(delta_summarised)
}
