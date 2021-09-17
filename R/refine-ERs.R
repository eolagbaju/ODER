#' Refines the ERs start and end points
#'
#' Uses the junctions added by \code{\link{annotatERs}} to modify the starts and
#' ends of the expressed regions. When a junction intersects an expressed region
#' depending on whether it is the start or end or both, the regions corresponding
#' starts and ends will be modified.
#'
#' As junctions mark intron boundaries, the expressed region will be changed to
#' either being one less or one more than the junction end.
#'
#' @param annot_ers ERs that have been annotated (result of annotatER)
#'
#' @return Genomic ranges with refined base pair starts and ends
#' @export
#'
#' @examples
#' # example set of ers to save runtime
#' ex_annot_ers <- GenomicRanges::GRanges( # this is created to not overlap
#'     seqnames = S4Vectors::Rle(c("chr21"), c(3)),
#'     ranges = IRanges::IRanges(
#'         start = c(5093576, 5097663, 5162182),
#'         end = c(5093833, 5097762, 5162257)
#'     ),
#'     grl = GenomicRanges::GRangesList(
#'         GenomicRanges::GRangesList(
#'             GenomicRanges::GRanges(
#'                 seqnames = S4Vectors::Rle(c("chr21"), c(1)),
#'                 ranges = IRanges::IRanges(
#'                     start = c(5093712),
#'                     end = c(5093744)
#'                 )
#'             ),
#'             GenomicRanges::GRanges(
#'                 seqnames = S4Vectors::Rle(c("chr21"), c(1)),
#'                 ranges = IRanges::IRanges(
#'                     start = c(5097642),
#'                     end = c(5097669)
#'                 )
#'             ),
#'             GenomicRanges::GRanges(
#'                 seqnames = S4Vectors::Rle(c("chr21"), c(1)),
#'                 ranges = IRanges::IRanges(
#'                     start = c(5162249),
#'                     end = c(5162287)
#'                 )
#'             )
#'         )
#'     ),
#'     annotation = c("intron", "intron", "intron")
#' )
#' refined_ers <- refine_ERs(ex_annot_ers)
#'
#' refined_ers
refine_ERs <- function(annot_ers) {

    # isolating the ers that are labelled as introns or intergenic
    annot_ers <- annot_ers[S4Vectors::mcols(annot_ers)[["annotation"]] %in% c("intron", "intergenic")]
    # getting the ers that only have one or two overlapping junctions
    annot_ers <- annot_ers[lengths(S4Vectors::mcols(annot_ers)[["grl"]]) == 2 |
        lengths(S4Vectors::mcols(annot_ers)[["grl"]]) == 1]
    # getting rid of ers with 2 junctions that overlap each other - seems to take forever tho
    # annot_ers <- annot_ers[unlist(lapply(S4Vectors::mcols(annot_ers)[["grl"]],inv_colgrs))]

    changes <- logical(0) # logical vector to record changes

    message(stringr::str_c(Sys.time(), " - Refining the Expressed regions..."))

    refined_results <- modify_ers(annot_ers)

    return(refined_results[[1]][refined_results[[2]]])
}


#' Changes the ers starts and ends based on the junctions
#'
#' @param annot_ers annotated ers with overlapping junctions
#'
#' @return ers with refined starts and ends and a logical vector listing the
#' changes
#' @keywords internal
#' @noRd
modify_ers <- function(annot_ers) {
    changes <- logical(0) # logical vector to record changes

    for (a in seq_along(annot_ers)) {
        er_start <- BiocGenerics::start(IRanges::ranges(annot_ers)[a])
        er_end <- BiocGenerics::end(IRanges::ranges(annot_ers)[a])
        change <- FALSE
        # seeing if single junction overlaps
        if (length(IRanges::ranges(S4Vectors::mcols(annot_ers)[["grl"]])[[a]]) == 1) {
            j_start <- as.integer(BiocGenerics::start(IRanges::ranges(S4Vectors::mcols(annot_ers)[["grl"]])[a]))
            j_end <- as.integer(BiocGenerics::end(IRanges::ranges(S4Vectors::mcols(annot_ers)[["grl"]])[a]))
            if (inbetween(value = j_start, rstart = er_start, rend = er_end) &
                inbetween(value = j_end, rstart = er_start, rend = er_end)) {
                er_start <- j_start + 1
                er_end <- j_end - 1
                change <- TRUE
            } else if (inbetween(value = j_start, rstart = er_start, rend = er_end)) {
                er_end <- j_start - 1
                change <- TRUE
            } else if (inbetween(value = j_end, rstart = er_start, rend = er_end)) {
                er_start <- j_end + 1
                change <- TRUE
            }

            # seeing if one or two junctions overlap
        } else if (length(IRanges::ranges(S4Vectors::mcols(annot_ers)[["grl"]])[[a]]) == 2) {
            if (colgrs(IRanges::ranges(S4Vectors::mcols(annot_ers)[["grl"]])[[a]])) {
                changes <- c(changes, change)
                next
            }
            if (IRanges::ranges(S4Vectors::mcols(annot_ers)[["grl"]])[[a]][1] >
                IRanges::ranges(S4Vectors::mcols(annot_ers)[["grl"]])[[a]][2]) {
                er_start <- as.integer(BiocGenerics::end(IRanges::ranges(S4Vectors::mcols(annot_ers)[["grl"]])[[a]][2])) + 1
                er_end <- as.integer(BiocGenerics::start(IRanges::ranges(S4Vectors::mcols(annot_ers)[["grl"]])[[a]][1])) - 1
                change <- TRUE
            } else {
                er_start <- as.integer(BiocGenerics::end(IRanges::ranges(S4Vectors::mcols(annot_ers)[["grl"]])[[a]][1])) + 1
                er_end <- as.integer(BiocGenerics::start(IRanges::ranges(S4Vectors::mcols(annot_ers)[["grl"]])[[a]][2])) - 1
                change <- TRUE
            }
        }

        BiocGenerics::start(annot_ers[a]) <- er_start
        BiocGenerics::end(annot_ers[a]) <- er_end
        changes <- c(changes, change)
    }
    return(list(annot_ers, changes))
}
