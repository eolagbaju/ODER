#' Adding the nearest expressed genes
#'
#' Updating expressed regions with the expressed gene that is closest to it.
#' After entering the tissue that has been sequenced, the nearest gene and
#' nearest expressed gene will be added to the metadata columns of the annotated
#' ERs.
#'
#' @param input_file GTEX median expression file, if left as NULL the default
#' file will be used.
#' @param tissue Tissue to filter for. See \code{\link{tissue_options}} for
#' options
#' @param gtf Either a string containg the path to a .gtf file or a pre-imported
#'   gtf using `rtracklayer::import` . Provides gene data to help determine
#'   the nearest gene and nearest expressed gene.
#' @param species character string containing the species to filter for,
#' Homo sapiens is the default
#' @param annot_ers annotated ERs i.e. the product of \code{\link{annotatERs}},
#' should have an mcols column called "annotation"
#' @return Granges with annotated ERs and details of their nearest expressed
#'    genes
#'
#' @export
#' @examples
#'
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
#'
#' ex_opt_ers <- GenomicRanges::GRanges(
#'     seqnames = S4Vectors::Rle(c("chr21", "chr22"), c(2, 2)),
#'     ranges = IRanges::IRanges(
#'         start = c(5116369, 5118691, 5125879, 5128214),
#'         end = c(5117231, 5118847, 5125988, 5128403)
#'     )
#' )
#'
#' ex_opt_ers_w_exp_genes <- add_expressed_genes(
#'     tissue = "lung", gtf = gtf_gr,
#'     annot_ers = ex_opt_ers
#' )
#'
#' ex_opt_ers_w_exp_genes
add_expressed_genes <- function(input_file = NULL, tissue, gtf,
    species = "Homo_sapiens", annot_ers) {
    tissue_df <- get_tissue(input_file = input_file, tissue = tissue)

    expressed_genes <- get_expressed_genes(gtf = gtf, species = species, tissue_df = tissue_df)

    full_annot_ers <- get_nearest_expressed_genes(annot_ers = annot_ers, exp_genes = expressed_genes, gtf = gtf)

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
        gtex_path <- file_cache(gtex_url)
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
#' tissue. Filters the gtf file for genes and keeps those that are present in
#' the tissue dataframe passed in.
#'
#' @param gtf gtf file path or gtf GRanges
#' @param species character string containing the species to filter for,
#' Homo sapiens is the default
#' @param tissue_df dataframe containing the expressed genes for a particular
#' tissue
#'
#' @return GRanges with the expressed genes for a specific tissue
#' @keywords internal
#' @noRd
get_expressed_genes <- function(gtf, species = "Homo_sapiens", tissue_df) {
    # gtf <- rtracklayer::import(gtf_path)
    gtf <- gtf_load(gtf)
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
#' \code{\link{annotatERs}}, expressed genes and a gtf file.
#'
#' @param gtf gtf file path or gtf GRanges
#' @param annot_ers annotated ERs, should have an mcols column called
#' "annotation"
#' @param exp_genes GRanges containing the expressed genes of a particular
#' tissue
#'
#' @return GRanges with the expressed genes for a specific tissue
#' @keywords internal
#' @noRd
get_nearest_expressed_genes <- function(annot_ers, exp_genes, gtf) {
    gtf <- gtf_load(gtf)
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
