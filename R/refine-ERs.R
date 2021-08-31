#' Refines the ERs start and end points
#'
#' Uses the junctions added by \code{\link{annotatERs}} to modify the starts and
#' ends of the expressed regions. When a junction intersects an expressed region
#' depending on whether it is the start or end or both the regions corresponding
#' starts and ends will be modified.
#'
#' As junctions mark intron boundaries, the expressed region will be changed to
#' either being one less or more than the junction end.
#'
#' @param annot_ers ERs that have been annotated (result of annotatER)
#'
#' @return Genomic ranges with refined base pair starts and ends and a logical
#' vector listing changes
#' @export
#'
#' @examples
#' \dontshow{
#'
#' if (!exists("rec_url")) {
#'     rec_url <- recount::download_study(
#'         project = "SRP012682",
#'         type = "samples",
#'         download = FALSE
#'     ) # .file_cache is an internal function to download a bigwig file from a link
#'     # if the file has been downloaded recently, it will be retrieved from a cache
#' }
#' bw_path <- .file_cache(rec_url[1])
#' if (!exists("gtf_path")) {
#'     gtf_url <- paste0(
#'         "http://ftp.ensembl.org/pub/release-103/gtf/",
#'         "homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz"
#'     )
#'     gtf_path <- .file_cache(gtf_url)
#' }
#' }
#' if (!exists("opt_ers1")) {
#'     opt_ers1 <- ODER(
#'         bw_paths = bw_path, auc_raw = auc_example,
#'         auc_target = 40e6 * 100, chrs = c("chr21"),
#'         genome = "hg38", mccs = c(5, 10), mrgs = c(10, 20),
#'         gtf = gtf_path, ucsc_chr = TRUE, ignore.strand = TRUE,
#'         exons_no_overlap = NULL, bw_chr = "chr"
#'     )
#' }
#' if (!exists("genom_state")) {
#'     genom_state <- generate_genomic_state(
#'         gtf = gtf_path,
#'         chrs_to_keep = c("21"), ensembl = TRUE
#'     )
#' }
#'
#' junctions <- lung_junc_21_22
#' if (!exists("annot_ers1")) {
#'     annot_ers1 <- annotatERs(
#'         opt_ers = head(opt_ers1[["opt_ers"]], 100), junc_data = junctions,
#'         gtf_path = gtf_path, chrs_to_keep = c("21"), ensembl = TRUE,
#'         genom_state = genom_state
#'     )
#' }
#'
#' annot_ers <- add_expressed_genes(tissue = "lung", gtf_path = gtf_path, annot_ers = annot_ers1)
#' refined_ers <- refine_ERs(annot_ers)
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

    print(stringr::str_c(Sys.time(), " - Refining the Expressed regions..."))

    refined_results <- modify_ers(annot_ers)

    return(refined_results[[1]][refined_results[[2]]])
}


#' Changes the ers starts and ends based on the junctions
#'
#' @param annot_ers annotated ers with overlapping junctions
#'
#' @return ers with refined starts and ends
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
