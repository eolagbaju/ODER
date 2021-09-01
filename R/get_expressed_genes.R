#' Adding the nearest expressed genes
#'
#' Updating expressed regions with the expressed gene that is closest to it.
#' After entering the tissue that has been sequenced, the nearest gene and
#' nearest expressed gene will be added to the metadata columns of the annotated
#' ERs.
#'
#' @inheritParams get_tissue
#' @inheritParams get_expressed_genes
#' @inheritParams get_nearest_expressed_genes
#' @return Granges with annotated ERs and details of their nearest expressed
#' genes
#' @export
#' @examples
#' \dontshow{
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
#'
#' if (!exists("opt_ers1")) {
#'     opt_ers1 <- ODER(
#'         bw_paths = bw_path, auc_raw = auc_example,
#'         auc_target = 40e6 * 100, chrs = c("chr21", "chr22"),
#'         genome = "hg38", mccs = c(5, 10), mrgs = c(10, 20),
#'         gtf = gtf_path, ucsc_chr = TRUE, ignore.strand = TRUE,
#'         exons_no_overlap = NULL, bw_chr = "chr"
#'     )
#' }
#' junctions <- SummarizedExperiment::rowRanges(dasper::junctions_example)
#' if (!exists("genom_state")) {
#'     genom_state <- generate_genomic_state(
#'         gtf = gtf_path,
#'         chrs_to_keep = c("1", "2", "X"), ensembl = TRUE
#'     )
#' }
#' if (!exists("annot_ers1")) {
#'     annot_ers1 <- annotatERs(
#'         opt_ers = head(opt_ers1[["opt_ers"]], 100), junc_data = junctions,
#'         gtf_path = gtf_path, chrs_to_keep = c("21", "22"), ensembl = TRUE,
#'         genom_state = genom_state
#'     )
#' }
#' annot_ers <- add_expressed_genes(
#'     tissue = "lung", gtf_path = gtf_path,
#'     annot_ers = annot_ers1
#' )
add_expressed_genes <- function(input_file = NULL, tissue, gtf_path,
    species = "Homo_sapiens", annot_ers) {
    tissue_df <- get_tissue(input_file = input_file, tissue = tissue)

    expressed_genes <- get_expressed_genes(gtf_path = gtf_path, species = species, tissue_df = tissue_df)

    full_annot_ers <- get_nearest_expressed_genes(annot_ers = annot_ers, exp_genes = expressed_genes, gtf_path = gtf_path)

    return(full_annot_ers)
}



#' Get gene data for a tissue
#'
#' Generate a dataframe for the tissue entered in containing a list of expressed
#' genes that appear above a threshold of RPKM>0.1.
#' (RPKM - Reads Per Kilobase of transcript, per Million mapped reads (RPKM) is
#' a normalized unit of transcript expression.)
#'
#' Can take a GTEX median expression file as an input or use the default (v6)
#'
#' @param input_file GTEX median expression file, if left as NULL
#' \code{\link{get_tissue}} will use the default file
#' @param tissue Tissue to filter for. See tissue_options for options
#'
#' @return Dataframe containing expressed genes
#' @keywords internal
#' @noRd
get_tissue <- function(input_file = NULL, tissue) {
    if (is.null(input_file)) {
        gtex_url <- "https://storage.googleapis.com/gtex_analysis_v6p/rna_seq_data/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct.gz"
        gtex_path <- .file_cache(gtex_url)
        gtex_data <- data.table::fread(gtex_path)
    } else {
        gtex_data <- data.table::fread(input_file)
    }

    data <- gtex_data %>% dplyr::mutate(Name = stringr::str_split_fixed(Name, "\\.", n = 2)[, 1])

    names <- colnames(data) %>%
        stringr::str_replace_all("\\.+", "_") %>%
        stringr::str_replace("_$", "") %>%
        tolower()

    tissue_index <- match(tissue, tissue_options) + 2

    df <- data[, c(1, tissue_index)] %>% dplyr::filter(.[[2]] > 0.1)

    return(df)
}

#' Get the expressed genes
#'
#' Takes in a gtf file and a dataframe with expressed genes from a specific
#' tissue. Filters the gtf file for genes and keeps those that are present in the
#' tissue dataframe passed in.
#'
#' @param gtf_path gtf file path
#' @param species character string containing the species to filter for,
#' Homo sapiens is the default
#' @param tissue_df dataframe containing the expressed genes for a particular
#' tissue
#'
#' @return GRanges with the expressed genes for a specific tissue
#' @keywords internal
#' @noRd
get_expressed_genes <- function(gtf_path, species = "Homo_sapiens", tissue_df) {
    gtf <- rtracklayer::import(gtf_path)
    gtf <- GenomeInfoDb::keepStandardChromosomes(gtf, species = species, pruning.mode = "coarse")
    GenomeInfoDb::seqlevelsStyle(gtf) <- "UCSC" # add chr to seqnames
    genesgtf <- gtf[S4Vectors::mcols(gtf)[["type"]] == "gene"]
    gtf.gene <- as.data.frame(genesgtf)

    gtf.exp.gr <- dplyr::semi_join(gtf.gene, tissue_df, by = c("gene_id" = "Name")) %>%
        GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = TRUE)

    return(gtf.exp.gr)
}

#' Get the expressed genes
#'
#' Adds the overall nearest gene and the nearest expressed gene to the annotated
#' expressed regions as metadata columns. Takes in annotated ers produced by
#' \code{\link{annotatERs}}, expressed genes from \code{\link{get_expressed_genes}}
#' and a gtf file.
#'
#' @param gtf_path gtf file path
#' @param annot_ers annotated ERs, should have an mcols column called "annotation"
#' @param exp_genes GRanges containing the expressed genes of a particular tissue
#'
#' @return GRanges with the expressed genes for a specific tissue
#' @keywords internal
#' @noRd
get_nearest_expressed_genes <- function(annot_ers, exp_genes, gtf_path) {
    gtf <- rtracklayer::import(gtf_path)
    gtf <- GenomeInfoDb::keepStandardChromosomes(gtf, species = "Homo_sapiens", pruning.mode = "coarse")
    GenomeInfoDb::seqlevelsStyle(gtf) <- "UCSC" # add chr to seqnames
    genesgtf <- gtf[S4Vectors::mcols(gtf)[["type"]] == "gene"]

    # annot_ers <- annot_ers[S4Vectors::mcols(annot_ers)[["annotation"]] %in% c("intron", "intergenic")]

    nearest_hit <- GenomicRanges::nearest(annot_ers, genesgtf, select = c("arbitrary"), ignore.strand = FALSE)
    S4Vectors::mcols(annot_ers)[["nearest_gene_v94_name"]] <- S4Vectors::mcols(genesgtf[nearest_hit])[["gene_id"]]

    exp_hit <- GenomicRanges::nearest(annot_ers, exp_genes, select = c("arbitrary"), ignore.strand = FALSE)
    S4Vectors::mcols(annot_ers)[["nearest_expressed_gene_v94_name"]] <- S4Vectors::mcols(exp_genes[exp_hit])[["gene_id"]]

    return(annot_ers)
}
