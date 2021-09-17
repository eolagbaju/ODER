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
#' @examples
#' \dontshow{
#' if (!exists("gtf_path")) {
#'     gtf_url <- paste0(
#'         "http://ftp.ensembl.org/pub/release-103/gtf/",
#'         "homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz"
#'     )
#'     gtf_path <- file_cache(gtf_url)
#' }
#'
#' if (!exists("gtf_gr")) {
#'     gtf_gr <- rtracklayer::import(gtf_path)
#' }
#' }
#' if (!exists("eg_opt_exons")) {
#'     eg_opt_exons <- get_exons(
#'         gtf = gtf_gr,
#'         ucsc_chr = TRUE,
#'         ignore.strand = TRUE
#'     )
#' }
#'
#' message(eg_opt_exons)
get_exons <- function(gtf, ucsc_chr, ignore.strand = TRUE, biotype = "Non-overlapping") {
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
        message(stringr::str_c(Sys.time(), " - Obtaining non-overlapping exons"))

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

    # GTF Processing

    gtf_gr <- GenomeInfoDb::keepStandardChromosomes(gtf_gr, species = "Homo_sapiens", pruning.mode = "coarse")
    gtf.df <- as.data.frame(gtf_gr, stringsAsFactor = FALSE)

    # Select protein-coding genes and transcripts
    gtf.df.pc <- gtf.df %>% dplyr::filter(gene_biotype %in% "protein_coding", transcript_biotype %in% "protein_coding")


    # Select TSL level 1
    gtf.df.pc$transcript_support_level <- gsub("\\s*\\([^\\)]+\\)", "", as.numeric(gtf.df.pc$transcript_support_level))
    gtf.df.pc.tsl1 <- gtf.df.pc %>% dplyr::filter(transcript_support_level %in% 1)

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
        dplyr::mutate(internal_cds = ifelse(exon_number == max(exon_number) | exon_number == min(exon_number), "NO", "YES")) %>%
        as.data.frame()

    all_utrs <- gtf.df %>% dplyr::filter(type %in% c("five_prime_utr", "three_prime_utr"))
    all_utrs$internal_cds <- "NO"

    # combine all GTF data together
    all_data <- rbind(all_exons, all_utrs)
    all_data_gr <- GenomicRanges::makeGRangesFromDataFrame(all_data, keep.extra.columns = TRUE)

    if (biotype == "Three Prime" | biotype == "3 Prime" | biotype == "3'") {
        message(stringr::str_c(Sys.time(), " - Obtaining Three Prime exons"))
        ##############################################################
        ################# Extract 3 prime UTRs #######################
        ##############################################################

        utr_all <- gtf.df.pc.tsl1 %>% dplyr::filter(type %in% "three_prime_utr")

        # Collapsing the 3'UTRs among the transcripts for each gene
        utr_all_grList <- GenomicRanges::makeGRangesListFromDataFrame(utr_all, split.field = "gene_id", names.field = "transcript_id")

        utr_all_collapse <- IRanges::reduce(utr_all_grList, with.revmap = TRUE) %>%
            as.data.frame() %>%
            dplyr::mutate(elements_collapsed = lengths(revmap), three_prime_utr_id = paste(group_name, seqnames, start, end, strand, elements_collapsed, sep = ":"))


        ##############################################################
        ################# 3 prime UTRs #######################
        ##############################################################

        ## concept
        # 1. If 3'UTR overlaps a 3'UTR of a transcript from the SAME gene, then this should be retained.
        # 2. If 3'UTR overlaps a coding exon of a transcript from the SAME gene, then this should be removed.
        # 3. If 3'UTR overlaps any part of a transcript from a DIFFERENT gene (any strand), then this should be removed.


        utr.gr <- GenomicRanges::makeGRangesFromDataFrame(utr_all_collapse, keep.extra.columns = TRUE)

        # Compute the overlap
        x <- IRanges::findOverlapPairs(utr.gr, all_data_gr, ignore.strand = TRUE) %>%
            as.data.frame() %>%
            dplyr::mutate(gold = ifelse(first.X.group_name %in% second.X.gene_id,
                ifelse(second.X.type %in% "three_prime_utr",
                    "YES",
                    ifelse(second.X.type %in% "exon" & second.internal_cds %in% "NO",
                        "YES",
                        "NO"
                    )
                ),
                "NO"
            ))

        # Extract the 3'UTRs which failed our conditions
        utr.failed <- x %>%
            dplyr::filter(gold == "NO") %>%
            dplyr::select(first.X.three_prime_utr_id) %>%
            dplyr::distinct(first.X.three_prime_utr_id) %>%
            plyr::rename(c("first.X.three_prime_utr_id" = "three_prime_utr_id"))


        # remove the failed ones
        utr.gs <- dplyr::anti_join(utr_all_collapse, utr.failed, by = "three_prime_utr_id")


        threeprime_exons <- GenomicRanges::makeGRangesFromDataFrame(utr.gs, keep.extra.columns = TRUE)
        GenomeInfoDb::seqlevelsStyle(threeprime_exons) <- "UCSC"

        return(threeprime_exons)
    }

    if (biotype == "Five Prime" | biotype == "5 Prime" | biotype == "5'") {
        message(stringr::str_c(Sys.time(), " - Obtaining Five Prime exons"))
        #######################################################################
        ################# Extract five prime  #################################
        #######################################################################

        # extract 5' UTRs
        five_prime <- gtf.df.pc.tsl1 %>% dplyr::filter(type %in% "five_prime_utr")

        # Collapsing the 5'UTRs among the transcripts for each gene
        five_prime_grList <- GenomicRanges::makeGRangesListFromDataFrame(five_prime, split.field = "gene_id", names.field = "transcript_id")

        five_prime_collapse <- IRanges::reduce(five_prime_grList, with.revmap = TRUE) %>%
            as.data.frame() %>%
            dplyr::mutate(elements_collapsed = lengths(revmap), five_prime_utr_id = paste(group_name, seqnames, start, end, strand, elements_collapsed, sep = ":"))


        # five = five_prime_collapse %>% dplyr::filter(width >=40)

        ##############################################################
        ################# five prime #######################
        ##############################################################

        ## concept
        # same as three prime

        five.gr <- GenomicRanges::makeGRangesFromDataFrame(five_prime_collapse, keep.extra.columns = TRUE)

        # Compute the overlap
        y <- IRanges::findOverlapPairs(five.gr, all_data_gr, ignore.strand = TRUE) %>%
            as.data.frame() %>%
            dplyr::mutate(gold = ifelse(first.X.group_name %in% second.X.gene_id,
                ifelse(second.X.type %in% "five_prime_utr",
                    "YES",
                    ifelse(second.X.type %in% "exon" & second.internal_cds %in% "NO",
                        "YES",
                        "NO"
                    )
                ),
                "NO"
            ))



        five.failed <- y %>%
            dplyr::filter(gold == "NO") %>%
            dplyr::select(first.X.five_prime_utr_id) %>%
            dplyr::distinct(first.X.five_prime_utr_id) %>%
            plyr::rename(c("first.X.five_prime_utr_id" = "five_prime_utr_id"))


        five.gs <- dplyr::anti_join(five_prime_collapse, five.failed, by = "five_prime_utr_id")

        fiveprime_exons <- GenomicRanges::makeGRangesFromDataFrame(five.gs, keep.extra.columns = TRUE)
        GenomeInfoDb::seqlevelsStyle(fiveprime_exons) <- "UCSC"
        return(fiveprime_exons)
    }

    if (biotype == "Internal") {
        message(stringr::str_c(Sys.time(), " - Obtaining Internal coding exons"))

        #######################################################################
        ################# Extract Internal coding exons (ICE) #######################
        #######################################################################

        # Extract all the coding exons
        cds_all <- gtf.df.pc.tsl1 %>%
            dplyr::filter(type %in% "CDS") %>%
            dplyr::mutate(CDS_id = paste(gene_id, seqnames, start, end, strand, sep = ":"))

        # number of coding exons per transcript
        cds_count <- table(cds_all$transcript_id) %>% as.data.frame()


        # Only the transcripts with >2 coding exons will contain internal exons, so I remove the transcripts with < 3 coding exons
        trans_to_remove <- cds_count %>%
            dplyr::filter(Freq < 3) %>%
            plyr::rename(c("Var1" = "transcript_id"))
        cds.filt <- dplyr::anti_join(cds_all, trans_to_remove, by = "transcript_id")
        cds.filt$exon_number <- as.numeric(cds.filt$exon_number)

        # extract INTERNAL CODING EXONS (ICEs): remove the first and the last coding exon
        internal.cds <- cds.filt %>%
            dplyr::group_by(transcript_id) %>%
            dplyr::filter(exon_number < max(exon_number), exon_number > min(exon_number)) %>%
            as.data.frame()

        # Collapsing the ICEs amongst the transcripts for each gene
        internal_cds_grList <- GenomicRanges::makeGRangesListFromDataFrame(internal.cds, split.field = "gene_id", names.field = "transcript_id")

        internal_cds_collapse <- IRanges::reduce(internal_cds_grList, with.revmap = TRUE) %>%
            as.data.frame() %>%
            dplyr::mutate(elements_collapsed = lengths(revmap), cds_id = paste(group_name, seqnames, start, end, strand, elements_collapsed, sep = ":"))

        ##############################################################
        ################# ICE #######################
        ##############################################################

        ## concept
        # 1. If ICE overlaps an ICE of a transcript from the SAME gene, then this should be retained.
        # 2. If ICE overlaps a terminal exon or a UTR of a transcript from the SAME gene, then this should be removed.
        # 3. If ICE overlaps any part of a transcript from a DIFFERENT gene, then this should be removed.

        internal.cds.gr <- GenomicRanges::makeGRangesFromDataFrame(internal_cds_collapse, keep.extra.columns = TRUE)

        # Compute the overlap
        y <- IRanges::findOverlapPairs(internal.cds.gr, all_data_gr, ignore.strand = TRUE) %>%
            as.data.frame() %>%
            dplyr::mutate(gold = ifelse(first.X.group_name %in% second.X.gene_id, ifelse(second.internal_cds %in% "NO", "NO", "YES"), "NO"))

        # Extract the ICEs which failed our conditions
        cds.failed <- y %>%
            dplyr::filter(gold == "NO") %>%
            dplyr::select(first.X.cds_id) %>%
            dplyr::distinct(first.X.cds_id) %>%
            plyr::rename(c("first.X.cds_id" = "cds_id"))

        # remove the failed ones
        internal.cds.gs <- dplyr::anti_join(internal_cds_collapse, cds.failed, by = "cds_id")

        internal_exons <- GenomicRanges::makeGRangesFromDataFrame(internal.cds.gs, keep.extra.columns = TRUE)
        GenomeInfoDb::seqlevelsStyle(internal_exons) <- "UCSC"
        return(internal_exons)
    }

    if (biotype == "lncRNA" | biotype == "LNCRNA" | biotype == "lncrna") {
        message(stringr::str_c(Sys.time(), " - Obtaining Long Non-Coding RNA"))
        #######################################################################
        #################  Extract lncRNAs  ###################################
        #######################################################################
        lncRNA <- c(
            "non_coding",
            "3prime_overlapping_ncRNA",
            "antisense",
            "lincRNA",
            "sense_intronic",
            "sense_overlapping",
            "macro_lncRNA",
            "lncRNA"
        )

        gtf.df.pc <- gtf.df %>% dplyr::filter(gene_biotype %in% lncRNA, transcript_biotype %in% lncRNA)

        # Select TSL level 1
        gtf.df.pc$transcript_support_level <- gsub("\\s*\\([^\\)]+\\)", "", as.numeric(gtf.df.pc$transcript_support_level))
        gtf.df.pc.tsl1 <- gtf.df.pc %>% dplyr::filter(transcript_support_level %in% 1)

        lncrna.gtf <- gtf.df.pc.tsl1 %>% dplyr::filter(type %in% "exon")

        # Collapsing the transcripts for each gene
        lncrna_all_grList <- GenomicRanges::makeGRangesListFromDataFrame(lncrna.gtf, split.field = "gene_id", names.field = "transcript_id")

        lncrna_all_collapse <- IRanges::reduce(lncrna_all_grList, with.revmap = TRUE) %>%
            as.data.frame() %>%
            dplyr::mutate(elements_collapsed = lengths(revmap), lncrna_id = paste(group_name, seqnames, start, end, strand, elements_collapsed, sep = ":"))

        ##############################################################
        ################# lncrna #######################
        ##############################################################

        ## concept
        # 1. If lncrna overlaps another transcript (ANY PART) from the SAME gene, then this should be retained.
        # 2. If lncrna overlaps any part of a transcript from a DIFFERENT gene, then this should be removed.

        lncrna.gr <- GenomicRanges::makeGRangesFromDataFrame(lncrna_all_collapse, keep.extra.columns = TRUE)

        # Compute the overlap
        y <- IRanges::findOverlapPairs(lncrna.gr, all_data_gr, ignore.strand = TRUE) %>%
            as.data.frame() %>%
            dplyr::mutate(gold = ifelse(first.X.group_name %in% second.X.gene_id, "YES", "NO"))


        lncrna.failed <- y %>%
            dplyr::filter(gold == "NO") %>%
            dplyr::select(first.X.lncrna_id) %>%
            dplyr::distinct(first.X.lncrna_id) %>%
            plyr::rename(c("first.X.lncrna_id" = "lncrna_id"))


        lncrna.gs <- dplyr::anti_join(lncrna_all_collapse, lncrna.failed, by = "lncrna_id")

        lncrna_exons <- GenomicRanges::makeGRangesFromDataFrame(lncrna.gs, keep.extra.columns = TRUE)
        GenomeInfoDb::seqlevelsStyle(lncrna_exons) <- "UCSC"
        return(lncrna_exons)
    }

    if (biotype == "ncRNA" | biotype == "NCRNA" | biotype == "ncrna") {
        message(stringr::str_c(Sys.time(), " - Obtaining Non-Coding RNA"))

        ####################################################################
        #################  Extract ncRNAs  #################################
        ####################################################################

        ncRNA <- c(
            "miRNA",
            "misc_RNA",
            "rRNA",
            "snRNA",
            "snoRNA",
            "vaultRNA"
        )

        # Select ncRNA
        # Comments: Transcript support level on ncRNA gene biotype has not been performed by Ensembl i.e. TSL category == NA. Therefore, all transcripts were taken into account.
        ncrna.gtf <- gtf.df %>% dplyr::filter(gene_biotype %in% ncRNA, type %in% "exon")


        # Collapsing among the transcripts for each gene
        ncrna_all_grList <- GenomicRanges::makeGRangesListFromDataFrame(ncrna.gtf, split.field = "gene_id", names.field = "transcript_id")

        ncrna_all_collapse <- IRanges::reduce(ncrna_all_grList, with.revmap = TRUE) %>%
            as.data.frame() %>%
            dplyr::mutate(elements_collapsed = lengths(revmap), ncrna_id = paste(group_name, seqnames, start, end, strand, elements_collapsed, sep = ":"))

        # ncrna = ncrna_all_collapse %>% dplyr::filter(width >=40)

        ##############################################################
        ################# ncrna #######################
        ##############################################################

        ## concept
        # 1. If ncrna overlaps another transcript (ANY PART) from the SAME gene, then this should be retained.
        # 2. If ncrna overlaps any part of a transcript from a DIFFERENT gene, then this should be removed.

        ncrna.gr <- GenomicRanges::makeGRangesFromDataFrame(ncrna_all_collapse, keep.extra.columns = TRUE)

        # Compute the overlap
        y <- IRanges::findOverlapPairs(ncrna.gr, all_data_gr, ignore.strand = TRUE) %>%
            as.data.frame() %>%
            dplyr::mutate(gold = ifelse(first.X.group_name %in% second.X.gene_id, "YES", "NO"))


        ncrna.failed <- y %>%
            dplyr::filter(gold == "NO") %>%
            dplyr::select(first.X.ncrna_id) %>%
            dplyr::distinct(first.X.ncrna_id) %>%
            plyr::rename(c("first.X.ncrna_id" = "ncrna_id"))


        ncrna.gs <- dplyr::anti_join(ncrna_all_collapse, ncrna.failed, by = "ncrna_id")
        ncrna_exons <- GenomicRanges::makeGRangesFromDataFrame(ncrna.gs, keep.extra.columns = TRUE)
        GenomeInfoDb::seqlevelsStyle(ncrna_exons) <- "UCSC"
        return(ncrna_exons)
    }

    if (biotype == "pseudo" | biotype == "Pseudo" | biotype == "pseudogene" | biotype == "Pseudogene") {
        message(stringr::str_c(Sys.time(), " - Obtaining Pseudogene"))

        ####################################################################
        #################  Extract pseudogenes  #################################
        ####################################################################

        pseudogene <- c(
            "pseudogene",
            "processed_pseudogene",
            "unprocessed_pseudogene",
            "transcribed_processed_pseudogene",
            "transcribed_unitary_pseudogene",
            "transcribed_unprocessed_pseudogene",
            "translated_processed_pseudogene",
            "unitary_pseudogene",
            "unprocessed_pseudogene",
            "TR_V_pseudogene",
            "TR_J_pseudogene",
            "rRNA_pseudogene",
            "polymorphic_pseudogene",
            "IG_V_pseudogene",
            "IG_pseudogene",
            "IG_J_pseudogene",
            "IG_C_pseudogene"
        )

        gtf.df.pc <- gtf.df %>% dplyr::filter(gene_biotype %in% pseudogene, transcript_biotype %in% pseudogene)


        # Select TSL level 1
        gtf.df.pc$transcript_support_level <- gsub("\\s*\\([^\\)]+\\)", "", as.numeric(gtf.df.pc$transcript_support_level))
        gtf.df.pc.tsl1 <- gtf.df.pc %>% dplyr::filter(transcript_support_level %in% 1)


        # Select pseudoGenes
        pseudo.gtf <- gtf.df.pc.tsl1 %>% dplyr::filter(gene_biotype %in% pseudogene, type %in% "exon")


        # Collapsing among the transcripts for each gene
        pseudo_all_grList <- GenomicRanges::makeGRangesListFromDataFrame(pseudo.gtf, split.field = "gene_id", names.field = "transcript_id")

        pseudo_all_collapse <- IRanges::reduce(pseudo_all_grList, with.revmap = TRUE) %>%
            as.data.frame() %>%
            dplyr::mutate(elements_collapsed = lengths(revmap), pseudoGene_id = paste(group_name, seqnames, start, end, strand, elements_collapsed, sep = ":"))


        # pseudogene = pseudo_all_collapse %>% dplyr::filter(width >=40)

        ##############################################################
        ################# pseudogene #######################
        ##############################################################

        ## concept
        # 1. If pseudogene overlaps another transcript (ANY PART) from the SAME gene, then this should be retained.
        # 2. If pseudogene overlaps any part of a transcript from a DIFFERENT gene, then this should be removed.

        pseudogene.gr <- GenomicRanges::makeGRangesFromDataFrame(pseudo_all_collapse, keep.extra.columns = TRUE)

        # Compute the overlap
        y <- IRanges::findOverlapPairs(pseudogene.gr, all_data_gr, ignore.strand = TRUE) %>%
            as.data.frame() %>%
            dplyr::mutate(gold = ifelse(first.X.group_name %in% second.X.gene_id, "YES", "NO"))


        pseudogene.failed <- y %>%
            dplyr::filter(gold == "NO") %>%
            dplyr::select(first.X.pseudoGene_id) %>%
            dplyr::distinct(first.X.pseudoGene_id) %>%
            plyr::rename(c("first.X.pseudoGene_id" = "pseudoGene_id"))


        pseudogene.gs <- dplyr::anti_join(pseudo_all_collapse, pseudogene.failed, by = "pseudoGene_id")
        pseudogene_exons <- GenomicRanges::makeGRangesFromDataFrame(pseudogene.gs, keep.extra.columns = TRUE)
        GenomeInfoDb::seqlevelsStyle(pseudogene_exons) <- "UCSC"
        return(pseudogene_exons)
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
#' @examples
#' \dontshow{
#' if (!exists("gtf_path")) {
#'     gtf_url <- paste0(
#'         "http://ftp.ensembl.org/pub/release-103/gtf/",
#'         "homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz"
#'     )
#'     # file_cache is an internal function to download a bigwig file from a
#'     # link
#'     # if the file has been downloaded recently, it will be retrieved from a
#'     # cache gtf_path <- file_cache(gtf_url)
#' }
#' if (!exists("eg_opt_exons")) {
#'     eg_opt_exons <- get_exons(
#'         gtf = gtf_path,
#'         ucsc_chr = TRUE,
#'         ignore.strand = TRUE
#'     )
#' }
#' }
#' if (!exists("eg_ers_delta")) {
#'     eg_ers_delta <- get_ers_delta(
#'         ers = ers_example, # ers_example is from the package data folder
#'         opt_exons = eg_opt_exons,
#'         delta_fun = ODER:::.delta
#'     ) # .delta is ODER's default, you can pass in your own if you have one
#' }
#' message(eg_ers_delta)
get_ers_delta <- function(ers, opt_exons, delta_fun = ODER:::.delta) {
    if (missing(ers)) {
        stop("No ERs were entered")
    } else if (missing(opt_exons)) {
        stop("No opt_exons were entered")
    }

    message(stringr::str_c(Sys.time(), " - Calculating delta for ERs..."))

    mcc_labels <- names(ers)

    delta_df <- dplyr::tibble()

    for (i in 1:seq_along(mcc_labels)) {
        mrg_labels <- names(ers[[mcc_labels[i]]])

        for (j in 1:seq_along(mrg_labels)) {
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
#' \dontshow{
#' if (!exists("gtf_path")) {
#'     gtf_url <- paste0(
#'         "http://ftp.ensembl.org/pub/release-103/gtf/",
#'         "homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz"
#'     )
#'     gtf_path <- file_cache(gtf_url)
#' }
#' if (!exists("eg_opt_exons")) {
#'     eg_opt_exons <- get_exons(
#'         gtf = gtf_path,
#'         ucsc_chr = TRUE,
#'         ignore.strand = TRUE
#'     )
#' }
#' if (!exists("eg_ers_delta")) {
#'     eg_ers_delta <- get_ers_delta(
#'         ers = ers_example, # ers_example is from the package data folder
#'         opt_exons = eg_opt_exons
#'     )
#' }
#' }
#' opt_ers <- get_opt_ers(
#'     ers = ers_example,
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
