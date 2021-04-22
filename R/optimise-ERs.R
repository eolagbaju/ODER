#' Obtain set of non-overlapping exons
#'
#' @param gtf Either a string containg the path to a .gtf file or a pre-imported gtf using
#'  \code{\link[rtracklayer]{import}}.
#' @param ucsc_chr logical scalar, determining whether to add "chr" prefix to the
#'  seqnames of non-overlapping exons and change "chrMT" -> "chrM". Note, if
#'  set to TRUE and seqnames already have "chr", it will not add another.
#' @param ignore.strand logical value for input into
#'  \code{\link[GenomicRanges]{findOverlaps}}, default is True.
#'
#' @return GRanges object containing non-overlapping exons.
#' @export
#'
#' @examples
#' gtf_url <- "http://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz"
#' gtf_tempfile <- file.path(tempdir(), "gtf_file.gtf")
#' download.file(url = gtf_url, destfile = gtf_tempfile)
#'
#' eg_exons_no_overlap <- get_exons(gtf = gtf_tempfile, ucsc_chr = T, ignore.strand = T)
get_exons <- function(gtf, ucsc_chr, ignore.strand = T) {
    if (is.character(gtf) & stringr::str_sub(gtf, -4, -1) != ".gtf") {
        stop("Please check your gtf file path")
    }

    if (is.character(gtf)) {
        print(stringr::str_c(Sys.time(), " - Loading in GTF..."))

        gtf_gr <- rtracklayer::import(gtf)
    } else {
        gtf_gr <- gtf
    }

    print(stringr::str_c(Sys.time(), " - Obtaining non-overlapping exons"))

    exons_gr <- gtf_gr[gtf_gr$type == "exon"]
    exons_gr <- exons_gr[!duplicated(exons_gr$exon_id)]

    exons_hits <- GenomicRanges::findOverlaps(exons_gr,
        drop.self = T,
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
