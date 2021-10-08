#' Obtain set of non-overlapping exons
#'
#' Downloads a well-defined set of exons to be used in obtaining the optimum set
#' of Expressed regions. These exons are used in calculating the exon deltas.
#'
#' @param gtf Either a string containg the path to a .gtf file or a pre-imported
#'   gtf using `rtracklayer::import` .
#' @param ucsc_chr logical scalar, determining whether to add "chr" prefix to
#'   the seqnames of non-overlapping exons and change "chrMT" -> "chrM". Note,
#'   if set to TRUE and seqnames already have "chr", it will not add another.
#' @param ignore.strand logical value for input into
#'   \code{\link[GenomicRanges]{findOverlaps}}, default is True.
#' @param biotype Filters the GTF file passed in to what would be considered the
#' "Gold Standard" exons. The Default is "Non-overlapping" but the options are:
#' "Non-overlapping" (exons that don't intersect each other),
#' "Three Prime" (3' UTR), "Five Prime" (5' UTR), "Internal" (Internal coding),
#' "lncRNA" (Long Non-Coding RNA), "ncRNA" (Non-Coding RNA) and "Pseudogene"
#'
#'
#' @return GRanges object containing non-overlapping exons.
#' @export
#'
#' @describeIn get_opt_ers Filter for the exons to calculate the deltas
#' against
#'
#' @examples
#'
#' gtf_url <- paste0(
#'     "http://ftp.ensembl.org/pub/release-103/gtf/",
#'     "homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz"
#' )
#' gtf_path <- file_cache(gtf_url)
#'
#' gtf_gr <- rtracklayer::import(gtf_path)
#'
#' eg_opt_exons <- get_exons(
#'     gtf = gtf_gr,
#'     ucsc_chr = TRUE,
#'     ignore.strand = TRUE
#' )
#'
#'
#' eg_opt_exons
get_exons <- function(gtf, ucsc_chr, ignore.strand = TRUE,
    biotype = "Non-overlapping") {
    if (is.character(gtf)) {
        if (!xor(
            stringr::str_sub(gtf, -4, -1) == ".gtf",
            stringr::str_sub(gtf, -7, -1) == ".gtf.gz"
        )) {
            stop("Please check your gtf file path")
        }
        message(stringr::str_c(Sys.time(), " - Loading in GTF..."))
        gtf_gr <- rtracklayer::import(gtf)
    } else {
        gtf_gr <- gtf
    }
    if (biotype == "Non-overlapping") {
        return(no_exons(gtf_gr, ignore.strand, ucsc_chr))
    }
    gtf.df.pc.tsl1 <- gtf_proc(gtf_gr) # GTF Processing
    gtf.df <- as.data.frame(gtf_gr, stringsAsFactor = FALSE)

    all_data_gr <- gold_exons(gtf.df)

    if (biotype == "Three Prime" | biotype == "3 Prime" | biotype == "3'") {
        return(tp_exons(gtf.df.pc.tsl1, all_data_gr))
    }
    if (biotype == "Five Prime" | biotype == "5 Prime" | biotype == "5'") {
        return(fp_exons(gtf.df.pc.tsl1, all_data_gr))
    }
    if (biotype == "Internal") {
        return(int_exons(gtf.df.pc.tsl1, all_data_gr))
    }
    if (biotype == "lncRNA" | biotype == "LNCRNA" | biotype == "lncrna") {
        return(lnc_exons(all_data_gr, gtf.df))
    }
    if (biotype == "ncRNA" | biotype == "NCRNA" | biotype == "ncrna") {
        return(nc_exons(all_data_gr, gtf.df))
    }
    if (biotype == "pseudo" | biotype == "Pseudo" |
        biotype == "pseudogene" | biotype == "Pseudogene") {
        return(pseudo_exons(all_data_gr, gtf.df))
    }
}

#' Calculates delta for sets of ERs
#'
#' Calculates the median exon delta and the number of ERs with an exon delta of
#' 0 by comparing each combination of MCC and MRG with the optimum exons from
#' the ensembl database.
#'
#' @param ers Sets of ERs across various MCCs/MRGs - output of
#' \code{\link{get_ers}}.
#' @param opt_exons GRanges object that contains the regions that ideally, you
#' want the ER definitions to match - output of \code{\link{get_exons}}.
#' @param delta_fun Function that calculates the delta between ERs and
#'   \code{opt_exons}. Takes as input a set of ERs from \code{ers} and
#'   \code{opt_exons}. Then outputs a tibble/dataframe containing the summarised
#'   delta scores for that set of one set of ERs.
#'
#' @return tibble/dataframe containing summarised delta values. One row per set
#'   of ERs.
#' @export
#'
#' @describeIn get_opt_ers Method to get ers delta to help determine the
#'  optimum ers
#'
#' @examples
#' data(gtex_SRP012682_SRX222703_lung_ers_1, package = "ODER")
#'
#' if (!exists("eg_ers_delta")) {
#'     eg_ers_delta <- get_ers_delta(
#'         ers = gtex_SRP012682_SRX222703_lung_ers_1,
#'         # gtex_SRP012682_SRX222703_lung_ers_1 is from the package data folder
#'         opt_exons = eg_opt_exons
#'     ) # .delta is ODER's default and is used if delta_fun is left NULL
#'     # you can pass in your own if you have one
#' }
#' eg_ers_delta
get_ers_delta <- function(ers, opt_exons, delta_fun = NULL) {
    if (missing(ers)) {
        stop("No ERs were entered")
    } else if (missing(opt_exons)) {
        stop("No opt_exons were entered")
    }

    if (is.null(delta_fun)) delta_fun <- .delta

    message(stringr::str_c(Sys.time(), " - Calculating delta for ERs..."))

    mcc_labels <- names(ers)

    delta_df <- dplyr::tibble()

    for (i in seq_along(mcc_labels)) {
        mrg_labels <- names(ers[[mcc_labels[i]]])

        for (j in seq_along(mrg_labels)) {
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
#' Uses a delta calculating function and a well defined set of exons to find
#' which combination of MCC and MRG gives the best definition of the Expressed
#' regions.
#'
#' @param ers Sets of ERs across various MCCs/MRGs - output of
#' \code{\link{get_ers}}.
#' @param ers_delta tibble/dataframe containing summarised delta values. One row
#' per set of ERs.
#'
#' @return list containing optimised ERs, optimal pair of MCC/MRGs and
#' \code{delta_df}
#' @export
#' @examples
#' data(gtex_SRP012682_SRX222703_lung_ers_1, package = "ODER")
#' opt_ers <- get_opt_ers(
#'     ers = gtex_SRP012682_SRX222703_lung_ers_1,
#'     ers_delta = eg_ers_delta
#' )
#' opt_ers
get_opt_ers <- function(ers, ers_delta) {
    if (missing(ers)) {
        stop("No ERs were entered")
    } else if (missing(ers_delta)) {
        stop("No ers_delta were entered")
    }

    message(stringr::str_c(Sys.time(), " - Obtaining optimal set of ERs..."))

    delta_opt <-
        ers_delta %>%
        dplyr::filter(median == min(median)) %>%
        dplyr::filter(n_eq_0 == max(n_eq_0)) # with the lowest median ER delta
    # and highest num of delta equal to 0

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


#' Default delta function
#'
#' Calculates the difference between the starts and ends of the ERs and the set
#' of exons. Default function used when using \code{\link{get_ers_delta}}.
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

#' Filters GRanges to non-overlapping exons
#'
#' @param gtf_gr GRanges
#' @param ignore.strand logical argument
#' @param ucsc_chr logical scalar
#'
#' @return exons_no_overlap_gr genomic ranges
#' @keywords internal
#' @noRd
no_exons <- function(gtf_gr, ignore.strand, ucsc_chr) {
    message(stringr::str_c(
        Sys.time(),
        " - Obtaining non-overlapping exons"
    ))
    exons_gr <- gtf_gr[gtf_gr$type == "exon"]
    exons_gr <- exons_gr[!duplicated(exons_gr$exon_id)]

    exons_hits <- GenomicRanges::findOverlaps(exons_gr,
        drop.self = TRUE,
        ignore.strand = ignore.strand
    )
    exons_no_overlap_gr <- exons_gr[-c(
        S4Vectors::queryHits(exons_hits) %>% unique()
    )]
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

#' Pre-processes gtf for use in exon location filtering
#'
#' @param gtf_gr GRanges
#'
#' @return exons_no_overlap_gr genomic ranges
#' @keywords internal
#' @noRd
gtf_proc <- function(gtf_gr) {
    # GTF Processing
    gtf_gr <- GenomeInfoDb::keepStandardChromosomes(gtf_gr,
        species = "Homo_sapiens",
        pruning.mode = "coarse"
    )
    gtf.df <- as.data.frame(gtf_gr, stringsAsFactor = FALSE)
    # Select protein-coding genes and transcripts
    gtf.df.pc <- gtf.df %>% dplyr::filter(
        gene_biotype %in% "protein_coding",
        transcript_biotype %in% "protein_coding"
    )
    # Select TSL level 1
    gtf.df.pc$transcript_support_level <- gsub(
        "\\s*\\([^\\)]+\\)",
        "",
        as.numeric(gtf.df.pc$transcript_support_level)
    )
    gtf.df.pc.tsl1 <- gtf.df.pc %>% dplyr::filter(
        transcript_support_level %in% 1
    )
    return(gtf.df.pc.tsl1)
}

#' Getting the gold standard exons
#'
#' Select protein-coding genes and transcripts
#'
#' @param gtf_gr GRanges
#'
#' @return all_data_gr gold standard genomic ranges
#' @keywords internal
#' @noRd
gold_exons <- function(gtf.df) {
    ##############################################################
    ################# FINDING GOLD STANDARD #######################
    ##############################################################
    # Prepare the GTF for gold standard analysis
    # we need all exon types and utrs for comparison
    # need a label whether terminal exon or internal

    gtf.df$exon_number <- as.numeric(gtf.df$exon_number)

    all_exons <- gtf.df %>%
        dplyr::filter(type %in% "exon") %>%
        dplyr::group_by(transcript_id) %>%
        dplyr::mutate(
            internal_cds = ifelse(
                exon_number == max(exon_number) | exon_number == min(
                    exon_number
                ), "NO", "YES"
            )
        ) %>%
        as.data.frame()
    all_utrs <- gtf.df %>% dplyr::filter(type %in% c(
        "five_prime_utr",
        "three_prime_utr"
    ))
    all_utrs$internal_cds <- "NO"
    # combine all GTF data together
    all_data <- rbind(all_exons, all_utrs)
    all_data_gr <- GenomicRanges::makeGRangesFromDataFrame(
        all_data,
        keep.extra.columns = TRUE
    )
    return(all_data_gr)
}


#' Getting the three prime coding exons
#'
#' Concept:
#' 1. If 3'UTR overlaps a 3'UTR of a transcript from the SAME gene,
#' then this should be retained.
#' 2. If 3'UTR overlaps a coding exon of a transcript from the SAME gene,
#' then this should be removed.
#' 3. If 3'UTR overlaps any part of a transcript from a
#' DIFFERENT gene (any strand), then this should be removed.
#'
#' @param gtf.df.pc.tsl1 GRanges
#' @param all_data_gr gold standar granges
#'
#' @return threeprime_exons 3' UTR exons
#' @keywords internal
#' @noRd
tp_exons <- function(gtf.df.pc.tsl1, all_data_gr) {
    message(stringr::str_c(Sys.time(), " - Obtaining Three Prime exons"))
    utr_all <- gtf.df.pc.tsl1 %>% dplyr::filter(type %in% "three_prime_utr")
    # Collapsing the 3'UTRs among the transcripts for each gene
    utr_all_collapse <- collapse_gtf(utr_all, "three_prime_utr_id")
    utr.gr <- GenomicRanges::makeGRangesFromDataFrame(utr_all_collapse,
        keep.extra.columns = TRUE
    )
    # Compute the overlap
    x <- IRanges::findOverlapPairs(utr.gr,
        all_data_gr,
        ignore.strand = TRUE
    ) %>%
        as.data.frame() %>%
        dplyr::mutate(gold = ifelse(
            first.X.group_name %in% second.X.gene_id,
            ifelse(second.X.type %in% "three_prime_utr", "YES",
                ifelse(second.X.type %in% "exon" &
                    second.internal_cds %in% "NO", "YES", "NO")
            ), "NO"
        ))
    # Extract the 3'UTRs which failed our conditions
    utr.failed <- x %>%
        dplyr::filter(gold == "NO") %>%
        dplyr::select(first.X.three_prime_utr_id) %>%
        dplyr::distinct(first.X.three_prime_utr_id) %>%
        plyr::rename(c("first.X.three_prime_utr_id" = "three_prime_utr_id"))
    # remove the failed ones
    utr.gs <- dplyr::anti_join(
        utr_all_collapse, utr.failed,
        by = "three_prime_utr_id"
    )
    threeprime_exons <- GenomicRanges::makeGRangesFromDataFrame(
        utr.gs,
        keep.extra.columns = TRUE
    )
    GenomeInfoDb::seqlevelsStyle(threeprime_exons) <- "UCSC"
    return(threeprime_exons)
}

#' Getting the five prime coding exons
#'
#' @param gtf.df.pc.tsl1 GRanges
#' @param all_data_gr gold standar granges
#'
#' @return fiveprime_exons 5' UTR exons
#' @keywords internal
#' @noRd
fp_exons <- function(gtf.df.pc.tsl1, all_data_gr) {
    message(stringr::str_c(Sys.time(), " - Obtaining Five Prime exons"))
    five_prime <- gtf.df.pc.tsl1 %>% dplyr::filter(
        type %in% "five_prime_utr"
    )
    # Collapsing the 5'UTRs among the transcripts for each gene
    five_prime_collapse <- collapse_gtf(five_prime, "five_prime_utr_id")
    five.gr <- GenomicRanges::makeGRangesFromDataFrame(
        five_prime_collapse,
        keep.extra.columns = TRUE
    )
    # Compute the overlap
    y <- IRanges::findOverlapPairs(five.gr,
        all_data_gr,
        ignore.strand = TRUE
    ) %>%
        as.data.frame() %>%
        dplyr::mutate(gold = ifelse(
            first.X.group_name %in% second.X.gene_id,
            ifelse(second.X.type %in% "five_prime_utr", "YES",
                ifelse(second.X.type %in% "exon" &
                    second.internal_cds %in% "NO", "YES", "NO")
            ), "NO"
        ))
    five.failed <- y %>%
        dplyr::filter(gold == "NO") %>%
        dplyr::select(first.X.five_prime_utr_id) %>%
        dplyr::distinct(first.X.five_prime_utr_id) %>%
        plyr::rename(c("first.X.five_prime_utr_id" = "five_prime_utr_id"))
    five.gs <- dplyr::anti_join(
        five_prime_collapse, five.failed,
        by = "five_prime_utr_id"
    )
    fiveprime_exons <- GenomicRanges::makeGRangesFromDataFrame(
        five.gs,
        keep.extra.columns = TRUE
    )
    GenomeInfoDb::seqlevelsStyle(fiveprime_exons) <- "UCSC"
    return(fiveprime_exons)
}

#' Getting the Internal exons
#'
#' Concept
#' 1. If ICE overlaps an ICE of a transcript from the SAME gene,
#' then this should be retained.
#' 2. If ICE overlaps a terminal exon or a UTR of a transcript from
#' the SAME gene, then this should be removed.
#' 3. If ICE overlaps any part of a transcript from a DIFFERENT gene,
#' then this should be removed.
#'
#' @param gtf.df.pc.tsl1 GRanges
#' @param all_data_gr gold standard granges
#'
#' @return int_exons internal exons
#' @keywords internal
#' @noRd
int_exons <- function(gtf.df.pc.tsl1, all_data_gr) {
    internal.cds <- .int_prep(gtf.df.pc.tsl1)
    # Collapsing the ICEs amongst the transcripts for each gene
    internal_cds_collapse <- collapse_gtf(internal.cds, "cds_id")
    internal.cds.gr <- GenomicRanges::makeGRangesFromDataFrame(
        internal_cds_collapse,
        keep.extra.columns = TRUE
    )
    y <- IRanges::findOverlapPairs( # Compute the overlap
        internal.cds.gr, all_data_gr,
        ignore.strand = TRUE
    ) %>%
        as.data.frame() %>%
        dplyr::mutate(
            gold = ifelse(
                first.X.group_name %in% second.X.gene_id,
                ifelse(second.internal_cds %in% "NO", "NO", "YES"), "NO"
            )
        )
    # Extract the ICEs which failed our conditions
    cds.failed <- y %>%
        dplyr::filter(gold == "NO") %>%
        dplyr::select(first.X.cds_id) %>%
        dplyr::distinct(first.X.cds_id) %>%
        plyr::rename(c("first.X.cds_id" = "cds_id"))
    internal.cds.gs <- dplyr::anti_join( # remove the failed ones
        internal_cds_collapse, cds.failed,
        by = "cds_id"
    )
    internal_exons <- GenomicRanges::makeGRangesFromDataFrame(
        internal.cds.gs,
        keep.extra.columns = TRUE
    )
    GenomeInfoDb::seqlevelsStyle(internal_exons) <- "UCSC"
    return(internal_exons)
}

#' Preprocessing for the internal exons
#'
#' @param gtf.df.pc.tsl1 GRanges
#'
#' @return internal.cds data frame
#' @keywords internal
#' @noRd
.int_prep <- function(gtf.df.pc.tsl1) {
    message(stringr::str_c(Sys.time(), " - Obtaining Internal coding exons"))
    # Extract all the coding exons
    cds_all <- gtf.df.pc.tsl1 %>%
        dplyr::filter(type %in% "CDS") %>%
        dplyr::mutate(CDS_id = paste(gene_id,
            seqnames,
            start, end, strand,
            sep = ":"
        ))
    # number of coding exons per transcript
    cds_count <- table(cds_all$transcript_id) %>% as.data.frame()
    # Only the transcripts with >2 coding exons will contain internal exons,
    # so I remove the transcripts with < 3 coding exons
    trans_to_remove <- cds_count %>%
        dplyr::filter(Freq < 3) %>%
        plyr::rename(c("Var1" = "transcript_id"))
    cds.filt <- dplyr::anti_join(cds_all, trans_to_remove, by = "transcript_id")
    cds.filt$exon_number <- as.numeric(cds.filt$exon_number)
    # extract INTERNAL CODING EXONS (ICEs):
    # remove the first and the last coding exon
    internal.cds <- cds.filt %>%
        dplyr::group_by(transcript_id) %>%
        dplyr::filter(
            exon_number < max(exon_number),
            exon_number > min(exon_number)
        ) %>%
        as.data.frame()
    return(internal.cds)
}

#' Getting the Long Non-coding RNAs
#'
# 1. If lncrna overlaps another transcript (ANY PART) from
# the SAME gene, then this should be retained.
# 2. If lncrna overlaps any part of a transcript from a DIFFERENT gene,
# then this should be removed.
#'
#' @param all_data_gr gold standard granges
#' @param gtf.df GRanges
#'
#' @return lncrna_exons Long Non-coding RNAs
#' @keywords internal
#' @noRd
lnc_exons <- function(all_data_gr, gtf.df) {
    lncRNA <- c(
        "non_coding", "3prime_overlapping_ncRNA", "antisense", "lincRNA",
        "sense_intronic", "sense_overlapping", "macro_lncRNA", "lncRNA"
    )
    gtf.df.pc <- gtf.df %>% dplyr::filter(
        gene_biotype %in% lncRNA,
        transcript_biotype %in% lncRNA
    )
    gtf.df.pc$transcript_support_level <- gsub( # Select TSL level 1
        "\\s*\\([^\\)]+\\)", "", as.numeric(gtf.df.pc$transcript_support_level)
    )
    gtf.df.pc.tsl1 <- gtf.df.pc %>% dplyr::filter(
        transcript_support_level %in% 1
    )
    lncrna.gtf <- gtf.df.pc.tsl1 %>% dplyr::filter(type %in% "exon")
    # Collapsing the transcripts for each gene
    lncrna_all_collapse <- collapse_gtf(lncrna.gtf, "lncrna_id")
    lncrna.gr <- GenomicRanges::makeGRangesFromDataFrame(
        lncrna_all_collapse,
        keep.extra.columns = TRUE
    )
    # Compute the overlap
    y <- IRanges::findOverlapPairs(lncrna.gr,
        all_data_gr,
        ignore.strand = TRUE
    ) %>%
        as.data.frame() %>%
        dplyr::mutate(gold = ifelse(
            first.X.group_name %in% second.X.gene_id, "YES", "NO"
        ))
    lncrna.failed <- y %>%
        dplyr::filter(gold == "NO") %>%
        dplyr::select(first.X.lncrna_id) %>%
        dplyr::distinct(first.X.lncrna_id) %>%
        plyr::rename(c("first.X.lncrna_id" = "lncrna_id"))
    lncrna.gs <- dplyr::anti_join(
        lncrna_all_collapse, lncrna.failed,
        by = "lncrna_id"
    )
    lncrna_exons <- GenomicRanges::makeGRangesFromDataFrame(
        lncrna.gs,
        keep.extra.columns = TRUE
    )
    GenomeInfoDb::seqlevelsStyle(lncrna_exons) <- "UCSC"
    return(lncrna_exons)
}


#' Getting the  Non-coding RNAs
#'
#' 1. If ncrna overlaps another transcript (ANY PART) from the SAME gene,
#' then this should be retained.
#' 2. If ncrna overlaps any part of a transcript from a DIFFERENT gene,
#' then this should be removed.
#'
#' @param all_data_gr gold standard granges
#' @param gtf.df GRanges
#'
#' @return ncrna_exons Non-coding exons
#' @keywords internal
#' @noRd
nc_exons <- function(all_data_gr, gtf.df) {
    message(stringr::str_c(Sys.time(), " - Obtaining Non-Coding RNA"))
    ncRNA <- c(
        "miRNA", "misc_RNA", "rRNA", "snRNA", "snoRNA", "vaultRNA"
    )
    ncrna.gtf <- gtf.df %>% dplyr::filter(
        gene_biotype %in% ncRNA, type %in% "exon"
    )
    # Collapsing among the transcripts for each gene
    ncrna_all_collapse <- collapse_gtf(ncrna.gtf, "ncrna_id")
    ncrna.gr <- GenomicRanges::makeGRangesFromDataFrame(
        ncrna_all_collapse,
        keep.extra.columns = TRUE
    )
    # Compute the overlap
    y <- IRanges::findOverlapPairs(ncrna.gr,
        all_data_gr,
        ignore.strand = TRUE
    ) %>%
        as.data.frame() %>%
        dplyr::mutate(
            gold = ifelse(
                first.X.group_name %in% second.X.gene_id, "YES", "NO"
            )
        )
    ncrna.failed <- y %>%
        dplyr::filter(gold == "NO") %>%
        dplyr::select(first.X.ncrna_id) %>%
        dplyr::distinct(first.X.ncrna_id) %>%
        plyr::rename(c("first.X.ncrna_id" = "ncrna_id"))
    ncrna.gs <- dplyr::anti_join(
        ncrna_all_collapse, ncrna.failed,
        by = "ncrna_id"
    )
    ncrna_exons <- GenomicRanges::makeGRangesFromDataFrame(
        ncrna.gs,
        keep.extra.columns = TRUE
    )
    GenomeInfoDb::seqlevelsStyle(ncrna_exons) <- "UCSC"
    return(ncrna_exons)
}


#' Getting the pseudogene exons
#'
#' 1. If pseudogene overlaps another transcript (ANY PART) from the
#' SAME gene, then this should be retained.
#' 2. If pseudogene overlaps any part of a transcript from a
#' DIFFERENT gene, then this should be removed.
#'
#' @param all_data_gr gold standard granges
#' @param gtf.df GRanges
#'
#' @return lncrna_exons Long Non-coding RNAs
#' @keywords internal
#' @noRd
pseudo_exons <- function(all_data_gr, gtf.df) {
    message(stringr::str_c(Sys.time(), " - Obtaining Pseudogene"))
    data(pseudogene, package = "ODER")
    gtf.df.pc <- gtf.df %>% dplyr::filter(
        gene_biotype %in% pseudogene,
        transcript_biotype %in% pseudogene
    ) # Select TSL level 1
    gtf.df.pc$transcript_support_level <- gsub(
        "\\s*\\([^\\)]+\\)", "", as.numeric(gtf.df.pc$transcript_support_level)
    )
    tsl1 <- gtf.df.pc %>% dplyr::filter(transcript_support_level %in% 1)
    pseudo.gtf <- tsl1 %>% dplyr::filter( # Select pseudoGenes
        gene_biotype %in% pseudogene, type %in% "exon"
    )
    pseudo_all_collapse <- collapse_gtf(pseudo.gtf, "pseudoGene_id")
    pseudogene.gr <- GenomicRanges::makeGRangesFromDataFrame(
        pseudo_all_collapse,
        keep.extra.columns = TRUE
    ) # Compute the overlap
    y <- IRanges::findOverlapPairs(pseudogene.gr, all_data_gr,
        ignore.strand = TRUE
    ) %>%
        as.data.frame() %>%
        dplyr::mutate(gold = ifelse(
            first.X.group_name %in% second.X.gene_id, "YES", "NO"
        ))
    pseudogene.failed <- y %>%
        dplyr::filter(gold == "NO") %>%
        dplyr::select(first.X.pseudoGene_id) %>%
        dplyr::distinct(first.X.pseudoGene_id) %>%
        plyr::rename(c("first.X.pseudoGene_id" = "pseudoGene_id"))
    pseudogene.gs <- dplyr::anti_join(
        pseudo_all_collapse, pseudogene.failed,
        by = "pseudoGene_id"
    )
    pseudogene_exons <- GenomicRanges::makeGRangesFromDataFrame(
        pseudogene.gs,
        keep.extra.columns = TRUE
    )
    GenomeInfoDb::seqlevelsStyle(pseudogene_exons) <- "UCSC"
    return(pseudogene_exons)
}

#' Collapsing among the transcripts for each gene
#'
#' @param a gtf in GRanges
#' @param exon_type string of exon location/type to collapse on
#'
#' @return all_collapse
#' @keywords internal
#' @noRd
collapse_gtf <- function(gtf, exon_type) {
    # Collapsing among the transcripts for each gene
    all_grList <- GenomicRanges::makeGRangesListFromDataFrame(
        gtf,
        split.field = "gene_id", names.field = "transcript_id"
    )
    all_collapse <- IRanges::reduce(
        all_grList,
        with.revmap = TRUE
    ) %>%
        as.data.frame() %>%
        dplyr::mutate(
            elements_collapsed = lengths(revmap),
            "{exon_type}" := paste(
                group_name, seqnames, start, end,
                strand, elements_collapsed,
                sep = ":"
            )
        )
    return(all_collapse)
}
