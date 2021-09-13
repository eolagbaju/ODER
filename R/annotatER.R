#' Connects ERs to genes using junction data, then classifies ERs into "exonic",
#' "intronic", "intergenic", or a combination of these categories
#'
#' Finds the overlap between junctions and ERs, then adds gene info and junction
#' info as metadata columns. Then, uses a `gtf` file or a `Txdb` passed in to
#' generate a genomic state used to label each ER as to whether they are exonic,
#' intronic, intergenic or none.
#'
#' @inheritParams get_junctions
#' @param genom_state a genomic state object
#' @param gtf gtf in a GRanges object, pre-imported using
#' `rtracklayer::import` . This is used to provide the gene information
#' for annotation.
#' @param txdb [TxDb-class][GenomicFeatures::TxDb-class] (txdb object) to create
#' genomic state. This is used to annotate the expressed regions as exonic, intronic
#' or intergenic.
#' @return annotated ERs
#' @export
#'
#' @examples
#' \dontshow{
#' if (!exists("gtf_path")) {
#'     gtf_url <- paste0(
#'         "http://ftp.ensembl.org/pub/release-103/gtf/",
#'         "homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz"
#'     )
#'     gtf_path <- .file_cache(gtf_url)
#' }
#'
#' if (!exists("gtf_gr")) {
#'     gtf_gr <- rtracklayer::import(gtf_path)
#' }
#' }
#' ex_opt_ers <- GenomicRanges::GRanges(
#'     seqnames = S4Vectors::Rle(c("chr21"), c(5)),
#'     ranges = IRanges::IRanges(
#'         start = c(5032176, 5033408, 5034717, 5035188, 5036577),
#'         end = c(5032217, 5033425, 5034756, 5035189, 5036581)
#'     )
#' )
#'
#' junctions <- SummarizedExperiment::rowRanges(dasper::junctions_example)
#'
#' chrs_to_keep <- c("21", "22")
#' #### preparing the txdb and genomstate object(s)
#' hg38_chrominfo <- GenomeInfoDb::getChromInfoFromUCSC("hg38")
#' new_info <- hg38_chrominfo$size[match(
#'     chrs_to_keep,
#'     GenomeInfoDb::mapSeqlevels(hg38_chrominfo$chrom, "Ensembl")
#' )]
#' names(new_info) <- chrs_to_keep
#' gtf_gr_tx <- GenomeInfoDb::keepSeqlevels(gtf_gr,
#'     chrs_to_keep,
#'     pruning.mode = "tidy"
#' )
#' GenomeInfoDb::seqlengths(gtf_gr_tx) <- new_info
#' GenomeInfoDb::seqlevelsStyle(gtf_gr_tx) <- "UCSC"
#' rtracklayer::genome(gtf_gr_tx) <- "hg38"
#'
#' ucsc_txdb <- GenomicFeatures::makeTxDbFromGRanges(gtf_gr_tx)
#' genom_state <- derfinder::makeGenomicState(txdb = ucsc_txdb)
#' ens_txdb <- ucsc_txdb
#' GenomeInfoDb::seqlevelsStyle(ens_txdb) <- "Ensembl"
#' ################### end of genomstate creation
#'
#' if (!exists("annot_ers1")) {
#'     annot_ers1 <- annotatERs(
#'         opt_ers = ex_opt_ers, junc_data = junctions,
#'         gtf = gtf_gr, txdb = ens_txdb, genom_state = genom_state
#'     )
#' }
#'
#' annot_ers1
annotatERs <- function(opt_ers,
    junc_data,
    genom_state,
    gtf,
    txdb) {
    ann_opt_ers <- get_junctions(
        opt_ers = opt_ers, junc_data = junc_data,
        txdb = txdb
    )
    print(stringr::str_c(Sys.time(), " - Annotating the Expressed regions..."))

    annot_ers <- derfinder::annotateRegions(
        regions = ann_opt_ers,
        genomicState = genom_state[["fullGenome"]],
        maxgap = -1L, minoverlap = 1L
    )
    annot_table <- annotate_table(annot_ers[["countTable"]])
    GenomicRanges::mcols(ann_opt_ers)$annotation <- annot_table[["region_annot"]]
    GenomicRanges::mcols(ann_opt_ers)$og_index <- annot_table[["ER_index"]]

    ## add missing junction genes
    genesource <- character(length(ann_opt_ers))
    GenomicRanges::mcols(ann_opt_ers)$gene_source <- genesource
    ## getting the gtf in to a Granges format
    gtf_gr <- gtf_load(gtf)

    # if all of none of the ers are matched up with a gene, this is run
    if (all(identical(unique(unlist(GenomicRanges::mcols(ann_opt_ers)$genes)), character(0)))) {
        genes_gr <- gtf_gr[gtf_gr$type == "gene"]
        GenomeInfoDb::seqlevelsStyle(genes_gr) <- "UCSC"
        nearest_genes <- GenomicRanges::nearest(ann_opt_ers, genes_gr)
        missing_genes <- GenomicRanges::mcols(genes_gr[nearest_genes])[["gene_id"]]
        GenomicRanges::mcols(ann_opt_ers)[["genes"]] <- missing_genes
        GenomicRanges::mcols(ann_opt_ers[GenomicRanges::mcols(ann_opt_ers)[["og_index"]]])[["genes"]] <- missing_genes
        GenomicRanges::mcols(ann_opt_ers[GenomicRanges::mcols(ann_opt_ers)[["og_index"]]])[["gene_source"]] <- "nearest gtf genes"

        ngdist <- GenomicRanges::distanceToNearest(ann_opt_ers, genes_gr)
        overmaxdist <- ngdist[GenomicRanges::mcols(ngdist)[["distance"]] > 10000]
        if (length(overmaxdist) > 0) {
            GenomicRanges::mcols(ann_opt_ers[queryHits(overmaxdist)])[["genes"]] <- ""
            GenomicRanges::mcols(ann_opt_ers[queryHits(overmaxdist)])[["gene_source"]] <- "Too far"
        }
        print(stringr::str_c(Sys.time(), " - done!"))

        return(ann_opt_ers)
    }
    GenomicRanges::mcols(ann_opt_ers[lengths(GenomicRanges::mcols(ann_opt_ers)[["genes"]]) > 0])$gene_source <- "junction(s)"


    ng_ann_opt_ers <- ann_opt_ers[lengths(GenomicRanges::mcols(ann_opt_ers)[["genes"]]) == 0]
    genes_gr <- gtf_gr[gtf_gr$type == "gene"]
    GenomeInfoDb::seqlevelsStyle(genes_gr) <- "UCSC"
    nearest_genes <- GenomicRanges::nearest(ng_ann_opt_ers, genes_gr)
    missing_genes <- GenomicRanges::mcols(genes_gr[nearest_genes])[["gene_id"]]
    GenomicRanges::mcols(ng_ann_opt_ers)[["genes"]] <- missing_genes
    GenomicRanges::mcols(ann_opt_ers[GenomicRanges::mcols(ng_ann_opt_ers)[["og_index"]]])[["genes"]] <- missing_genes
    GenomicRanges::mcols(ann_opt_ers[GenomicRanges::mcols(ng_ann_opt_ers)[["og_index"]]])[["gene_source"]] <- "nearest gtf genes"

    ngdist <- GenomicRanges::distanceToNearest(ann_opt_ers, genes_gr)
    overmaxdist <- ngdist[GenomicRanges::mcols(ngdist)[["distance"]] > 10000]
    # GenomicRanges::mcols(ann_opt_ers[queryHits(overmaxdist)])[["genes"]] <- ""
    # GenomicRanges::mcols(ann_opt_ers[queryHits(overmaxdist)])[["gene_source"]] <- "Too far"
    if (length(overmaxdist) > 0) {
        GenomicRanges::mcols(ann_opt_ers[queryHits(overmaxdist)])[["genes"]] <- ""
        GenomicRanges::mcols(ann_opt_ers[queryHits(overmaxdist)])[["gene_source"]] <- "Too far"
    }

    print(stringr::str_c(Sys.time(), " - done!"))

    return(ann_opt_ers)
}




#' Get junctions that overlap the optimally defined ERs
#'
#' The junctions will be added to the ERs passed in as a metadata column, a
#' `GRangeslist` with each GRanges element of the list corresponding to it's
#' associated ER. If there are no overlaps, the `GRangeslist` will be empty. Any
#' ERs that do not have an overlapping junction will have an empty GRanges. The
#' respective genes that each ER could be associated with will also be passed in
#' as a metadata column, a character list.
#'
#' @seealso dasper::junction_annot
#'
#' @param opt_ers optimally defined ERs (the product of the ODER function)
#' @param junc_data junction data that should match the ERs passed into opt_ers
#' @param txdb a gtf file or txdb ([TxDb-class][GenomicFeatures::TxDb-class])
#' with exon data. This will ultimately be used to annotate the junction data
#' passed in.
#'
#' @return optimally defined ers annotated with junction and gene information
#' @export
#'
#' @examples
#' \dontshow{
#' if (!exists("gtf_path")) {
#'     gtf_url <- paste0(
#'         "http://ftp.ensembl.org/pub/release-103/gtf/",
#'         "homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz"
#'     )
#'     # .file_cache is an internal function to download a bigwig file from a link
#'     # if the file has been downloaded recently, it will be retrieved from a cache
#'     gtf_path <- .file_cache(gtf_url)
#' }
#'
#' if (!exists("gtf_gr")) {
#'     gtf_gr <- rtracklayer::import(gtf_path)
#' }
#'
#' chrs_to_keep <- c("21", "22")
#' if (!exists("ens_txdb")) {
#'     #### preparing the txdb object
#'     hg38_chrominfo <- GenomeInfoDb::getChromInfoFromUCSC("hg38")
#'     new_info <- hg38_chrominfo$size[match(
#'         chrs_to_keep,
#'         GenomeInfoDb::mapSeqlevels(hg38_chrominfo$chrom, "Ensembl")
#'     )]
#'     names(new_info) <- chrs_to_keep
#'     gtf_gr_tx <- GenomeInfoDb::keepSeqlevels(gtf_gr,
#'         chrs_to_keep,
#'         pruning.mode = "tidy"
#'     )
#'     GenomeInfoDb::seqlengths(gtf_gr_tx) <- new_info
#'     GenomeInfoDb::seqlevelsStyle(gtf_gr_tx) <- "UCSC"
#'     rtracklayer::genome(gtf_gr_tx) <- "hg38"
#'
#'     ens_txdb <- GenomicFeatures::makeTxDbFromGRanges(gtf_gr_tx)
#'     GenomeInfoDb::seqlevelsStyle(ens_txdb) <- "Ensembl"
#'     ################### end of txdb creation
#' }
#' }
#' example_ers <- GenomicRanges::GRanges(
#'     seqnames = S4Vectors::Rle(c("chr21", "chr22"), c(10, 1)),
#'     ranges = IRanges::IRanges(
#'         start = c(
#'             5026423, 24738, 5032218, 5033895, 17554, 50800446,
#'             50800539, 16570, 50800790, 15005, 20312
#'         ),
#'         end = c(
#'             5323718, 24891, 5033407, 5033980, 17728, 50800910,
#'             50800817, 16723, 50800910, 15038, 20582
#'         ),
#'         names = head(letters, 11)
#'     ),
#'     strand = S4Vectors::Rle((c("+", "-")), c(6, 5)),
#'     score = 1:11,
#'     GC = seq(1, 0, length = 11)
#' )
#'
#' example_junctions <- SummarizedExperiment::rowRanges(dasper::junctions_example)
#'
#' example_er_juncs <- get_junctions(
#'     opt_ers = example_ers,
#'     junc_data = example_junctions,
#'     txdb = ens_txdb # can either be a gtf file or txdb in the Ensembl format
#' )
#'
#' print(example_er_juncs)
get_junctions <- function(opt_ers, junc_data, txdb) {
    if (methods::is(junc_data, "data.frame")) {
        junc_data <- GenomicRanges::makeGRangesFromDataFrame(junc_data)
    }
    GenomeInfoDb::seqlevelsStyle(junc_data) <- "NCBI"
    annotated_junctions <- dasper::junction_annot(junctions = junc_data, ref = txdb)
    GenomeInfoDb::seqlevelsStyle(annotated_junctions) <- "UCSC" # to match opt_ers

    print(stringr::str_c(Sys.time(), " - Finding junctions overlapping ers..."))

    hits <- GenomicRanges::findOverlaps(opt_ers, annotated_junctions)
    ann_junc_hits <- annotated_junctions[S4Vectors::subjectHits(hits)]
    ann_junc_hits$er_index <- S4Vectors::queryHits(hits)

    gene_id_list <- GenomicRanges::mcols(ann_junc_hits)[["gene_id_junction"]] %>%
        S4Vectors::split(f = ann_junc_hits$er_index)

    gene_id_list <- lapply(gene_id_list, function(x) {
        return(unique(unlist(x)))
    })

    ann_junc_hits <- ann_junc_hits %>% S4Vectors::split(f = ann_junc_hits$er_index)

    er_indices <- unique(S4Vectors::queryHits(hits))

    miss_ers <- numeric(0)
    for (i in 1:length(opt_ers)) {
        if (!(i %in% er_indices)) {
            miss_ers <- c(miss_ers, i)
        }
    }

    empty_grl <- GenomicRanges::GRangesList(
        lapply(miss_ers, function(x) {
            return(GenomicRanges::GRanges())
        })
    )
    names(empty_grl) <- miss_ers
    combi_ajh <- suppressWarnings(c(ann_junc_hits, empty_grl))
    sorted_combi_ajh <- combi_ajh[order(as.integer(names(combi_ajh)))]

    empty_gil <- vector("list", length = length(miss_ers))
    names(empty_gil) <- miss_ers

    combi_gil <- c(gene_id_list, empty_gil)
    sorted_combi_gil <- combi_gil[order(as.integer(names(combi_gil)))]
    sorted_combi_gil <- sorted_combi_gil %>% IRanges::CharacterList()

    GenomicRanges::mcols(opt_ers)$grl <- sorted_combi_ajh
    GenomicRanges::mcols(opt_ers)$genes <- sorted_combi_gil

    return(opt_ers)
}


#' Converting count table output of derfinder to region annotation
#'
#' \code{convert_annot_count_table_to_region_annot} takes as input the derfinder
#' output and converts to useful region annotation - "intron", "exon",
#' "intergenic"... etc
#'
#' @param count_table count table output from
#'   \code{\link{derfinder::annotateRegions}}
#'
#' @return df with the notation "intron", "exon", "intergenic"
#'
#' @keywords internal
#' @noRd
annotate_table <- function(count_table) {
    count_table_tmp <-
        count_table %>%
        dplyr::mutate(
            exon = ifelse(exon > 0, 1, 0),
            intergenic = ifelse(intergenic > 0, 1, 0),
            intron = ifelse(intron > 0, 1, 0),
            annot_code_tmp = stringr::str_c(exon, intergenic, intron)
        ) %>%
        tibble::as_tibble()

    count_table[["region_annot"]] <-
        count_table_tmp[["annot_code_tmp"]] %>%
        purrr::map(count_conversion) %>%
        unlist()

    count_table_w_region_annot <-
        count_table %>%
        dplyr::mutate(ER_index = dplyr::row_number())

    return(count_table_w_region_annot)
}

#' Converting the count table to annotations
#'
#' @param annot_code
#'
#' @return table with annotation added on
#' @keywords internal
#' @noRd
count_conversion <- function(annot_code) {
    region_annot <-
        switch(annot_code,
            `100` = "exon",
            `110` = "exon, intergenic",
            `101` = "exon, intron",
            `111` = "exon, intergenic, intron",
            `010` = "intergenic",
            `011` = "intron, intergenic",
            `001` = "intron",
            `000` = "none"
        )

    return(region_annot)
}
