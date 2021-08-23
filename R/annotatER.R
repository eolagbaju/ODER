#' Finds overlapping junctions and annotates ERs as exonic, intronic, intergenic,
#' some combination or none of those
#'
#' Looks at junctions passed in to find any overlaps and adds them in along with
#' other information as metadata columns. Then uses a gtf file or a Txdb passed
#' in to generate a genomic state and then labels each ER as to whether they are
#' exonic, intronic, intergenic on none
#'
#'
#' @inheritParams get_junctions
#' @inheritParams generate_genomic_state
#'
#' @return annotated ERs
#' @export
#'
#' @examples
#' \dontshow{
#' url <- recount::download_study(
#'     project = "SRP012682",
#'     type = "samples",
#'     download = FALSE
#' ) # .file_cache is an internal function to download a bigwig file from a link
#' # if the file has been downloaded recently, it will be retrieved from a cache
#'
#' bw_path <- ODER:::.file_cache(url[1])
#' gtf_url <- paste0(
#'     "http://ftp.ensembl.org/pub/release-103/gtf/",
#'     "homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz"
#' )
#' gtf_path <- ODER:::.file_cache(gtf_url)
#' }
#'
#' opt_ers <- ODER(
#'     bw_paths = bw_path, auc_raw = auc_example,
#'     auc_target = 40e6 * 100, chrs = c("chr21", "chr22"),
#'     genome = "hg38", mccs = c(5, 10), mrgs = c(10, 20),
#'     gtf = gtf_path, ucsc_chr = TRUE, ignore.strand = TRUE,
#'     exons_no_overlap = NULL, bw_chr = "chr"
#' )
#'
#' junctions <- SummarizedExperiment::rowRanges(dasper::junctions_example)
#' annot_ers <- annotatERs(
#'     opt_ers = opt_ers[["opt_ers"]], junc_data = junctions,
#'     gtf_path = gtf_path, chrs_to_keep = c("21", "22"), ensembl = TRUE
#' )
#'
#' annot_ers
annotatERs <- function(opt_ers, junc_data, gtf_path, txdb = NULL,
    chrs_to_keep = c(1:22, "X", "Y", "MT"),
    ensembl = TRUE) {
    ann_opt_ers <- get_junctions(
        opt_ers = opt_ers, junc_data = junc_data,
        gtf_path = gtf_path
    )


    genom_state <- generate_genomic_state(
        gtf = gtf_path, txdb = txdb,
        chrs_to_keep = informatting2(chrs_to_keep),
        ensembl = ensembl
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
    GenomicRanges::mcols(ann_opt_ers[lengths(GenomicRanges::mcols(ann_opt_ers)[["genes"]]) > 0])$gene_source <- "junction(s)"


    ng_ann_opt_ers <- ann_opt_ers[lengths(GenomicRanges::mcols(ann_opt_ers)[["genes"]]) == 0]
    gtf_gr <- rtracklayer::import(gtf_path)
    genes_gr <- gtf_gr[gtf_gr$type == "gene"]
    GenomeInfoDb::seqlevelsStyle(genes_gr) <- "UCSC"
    nearest_genes <- GenomicRanges::nearest(ng_ann_opt_ers, genes_gr)
    missing_genes <- GenomicRanges::mcols(genes_gr[nearest_genes])[["gene_id"]]
    GenomicRanges::mcols(ng_ann_opt_ers)[["genes"]] <- missing_genes
    GenomicRanges::mcols(ann_opt_ers[GenomicRanges::mcols(ng_ann_opt_ers)[["og_index"]]])[["genes"]] <- missing_genes
    GenomicRanges::mcols(ann_opt_ers[GenomicRanges::mcols(ng_ann_opt_ers)[["og_index"]]])[["gene_source"]] <- "nearest gtf genes"

    ngdist <- GenomicRanges::distanceToNearest(ann_opt_ers, genes_gr)
    overmaxdist <- ngdist[GenomicRanges::mcols(ngdist)[["distance"]] > 10000]
    GenomicRanges::mcols(ann_opt_ers[queryHits(overmaxdist)])[["genes"]] <- ""
    GenomicRanges::mcols(ann_opt_ers[queryHits(overmaxdist)])[["gene_source"]] <- "Too far"


    print(stringr::str_c(Sys.time(), " - done!"))

    return(ann_opt_ers)
}




#' Get junctions that overlap the optimally defined ERs
#'
#' The junctions will be added to the ERs passed in as a metadata column, a
#' GRangeslist with each GRanges element of the list corresponding to it's
#' associated ER. If there are no overlaps, the GRangeslist will be empty.
#' Any ERs that do not have an overlapping junction will have an empty GRanges. 
#' The respective genes that each ER could be associated with will also be passed 
#' in as a metadata column, a character list
#'
#' Each Granges of the GRangeslist will have the metadata columns of "in_ref",
#' "gene_id_start", "tx_name_start", "exon_name_start", "strand_start",
#' "exon_width_start","gene_id_end", "tx_name_end", "exon_name_end", "strand_end",
#' "exon_width_end", "gene_id_junction", "strand_junction", "type" and "er_index
#'  added by dasper's junction_annot. 
#'
#' @seealso dasper::junction_annot
#'
#' @param opt_ers optimally defined ERs (the product of the ODER function)
#' @param junc_data junction data that should match the ERs passed into opt_ers
#' @param gtf_path a gtf file with exon data
#'
#' @return optimally defined ers annotated with junction and gene information
#' @export
#'
#' @examples
#' \dontshow{
#' gtf_url <- paste0(
#'     "http://ftp.ensembl.org/pub/release-103/gtf/",
#'     "homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz"
#' )
#' # .file_cache is an internal function to download a bigwig file from a link
#' # if the file has been downloaded recently, it will be retrieved from a cache
#' gtf_path <- ODER:::.file_cache(gtf_url)
#' }
#'
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
#'     gtf_path = gtf_path
#' )
#'
#' print(example_er_juncs)
get_junctions <- function(opt_ers, junc_data, gtf_path) {
    if (methods::is(junc_data, "data.frame")) {
        junc_data <- GenomicRanges::makeGRangesFromDataFrame(junc_data)
    }
    GenomeInfoDb::seqlevelsStyle(junc_data) <- "NCBI"
    annotated_junctions <- dasper::junction_annot(junctions = junc_data, ref = gtf_path)
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


#' Generating a genomic state object from Txdb or gtf
#'
#' \code{generate_genomic_state} takes txdb object or a gtf file and makes a
#' genomic state to be used in the annotation of the expressed regions 
#'
#' @param txdb txdb object, if one is not entered a gtf file needs to be
#' @param chrs_to_keep chromosomes to keep in genomic state (in NCBI format i.e.
#' 1,2,3...22,X,Y,MT)
#' @param gtf gtf file, will be converted into a txdb and then genomic state, if
#' no gtf is entered a txdb must be
#' @param ensembl logical variable to say whether the gtf file entered is from
#' ensembl, the default is true
#'
#' @return a genomic state
#'
#' @export
#'
#' @examples
#' \dontshow{
#' gtf_url <- paste0(
#'     "http://ftp.ensembl.org/pub/release-103/gtf/",
#'     "homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz"
#' )
#' # .file_cache is an internal function to download a bigwig file from a link
#' # if the file has been downloaded recently, it will be retrieved from a cache
#' gtf_path <- ODER:::.file_cache(gtf_url)
#' }
#' genom_state <- generate_genomic_state(
#'     gtf = gtf_path,
#'     chrs_to_keep = c("1", "2", "X"), ensembl = TRUE
#' )
generate_genomic_state <- function(gtf = NULL, txdb = NULL,
    chrs_to_keep = c(1:22, "X", "Y", "MT"),
    ensembl = TRUE) {
    if (is.null(gtf) && is.null(txdb)) {
        stop("Neither a gtf or a txdb was entered")
    } else if (!(is.null(txdb))) {
        return(derfinder::makeGenomicState(txdb = txdb))
    } else if (ensembl) {
        print(stringr::str_c(Sys.time(), " - Generating a genomic state..."))

        gtf_txdb <- GenomicFeatures::makeTxDbFromGFF(file = gtf, format = "gtf")
        gtf_txdb <- GenomeInfoDb::keepSeqlevels(gtf_txdb, chrs_to_keep)

        hg38_chrominfo <- GenomeInfoDb::getChromInfoFromUCSC("hg38")
        new_info <- hg38_chrominfo$size[match(
            GenomeInfoDb::seqlevels(gtf_txdb),
            GenomeInfoDb::mapSeqlevels(hg38_chrominfo$chrom, "Ensembl")
        )]
        names(new_info) <- chrs_to_keep

        gtf_gr <- rtracklayer::import(gtf)
        gtf_gr <- GenomeInfoDb::keepSeqlevels(gtf_gr,
            chrs_to_keep,
            pruning.mode = "tidy"
        )
        GenomeInfoDb::seqlengths(gtf_gr) <- new_info
        GenomeInfoDb::seqlevelsStyle(gtf_gr) <- "UCSC"
        rtracklayer::genome(gtf_gr) <- "hg38"

        pro_txdb <- GenomicFeatures::makeTxDbFromGRanges(gtf_gr)

        return(derfinder::makeGenomicState(txdb = pro_txdb))
    }
}

#' Converting count table output of derfinder to region annotation
#'
#' \code{convert_annot_count_table_to_region_annot} takes as input the derfinder
#'  output and converts to useful region annotation - "intron", "exon", "intergenic"... etc
#'
#' @param count_table count table output from \code{\link{derfinder::annotateRegions}}
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
