#' An example AUC value
#'
#' An Area Under Coverage (AUC) value for a user to try out the package and to
#' pass in for tests. From the GTEX data set and project SRP012682, the actual
#' value is 11872688252.
#'
#' @format A numeric value
#' @usage data(gtex_SRP012682_SRX222703_lung_auc_1)
#' @source See example.R in data-raw
#'
"gtex_SRP012682_SRX222703_lung_auc_1"
#' An example object containing coverage
#'
#' Coverage generated for a user to try out the package and
#' to pass in for tests. Coverage of chromosomes 21 and 22 from the project
#' SRP012682.
#'
#' @format A list of length 2 containing 2 Rles for chromosomes 21 and
#' 22 respectively
#' @usage data(gtex_SRP012682_SRX222703_lung_coverage_1)
#' @source See example.R in data-raw
#'
"gtex_SRP012682_SRX222703_lung_coverage_1"
#' An example set of Expressed Regions
#'
#' An example set of Expressed Regions generated for a user to try out the
#' package and to pass in for tests. Generated using
#' gtex_SRP012682_SRX222703_lung_coverage_1 and MCCs of 5 and 10 and MRGs of
#' 10 and 20.
#'
#' @format A list containing two lists (for each mcc) each with a set of genomic
#' ranges for the different combinations of mcc and mrg
#' @usage data(gtex_SRP012682_SRX222703_lung_ers_1)
#' @source See example.R in data-raw
#'
"gtex_SRP012682_SRX222703_lung_ers_1"
#' An example set of ER deltas
#'
#' This set of deltas was calculated using gtex_lung_ers_1 and exons from
#' ensembl.
#'
#' @format A tibble/dataframe with the sums, means, medians, n_eq_0 and
#' propor_eq_0 for each combination of mccs (5 & 10) and mrgs (10 & 20)
#' @usage data(gtex_SRP012682_SRX222703_lung_erdelta_1)
#' @source See example.R in data-raw
#'
"gtex_SRP012682_SRX222703_lung_erdelta_1"
#' Junction data of chromosomes 21 and 22 from a lung tissue sample
#'
#' These junctions were sampled from a local junction file.
#'
#' @format A dataframe with the junction ID, chromosome, start and ends, strand,
#' number of samples, acceptor and donor
#' @usage data(lung_junc_21_22)
#' @source GTEx
#'
"lung_junc_21_22"
#' The different tissues that can be filtered on for gene expression
#'
#' These options were derived from the contents of the GTEx analysis gene median
#' RPKM file.
#'
#' @format A character vector with all of the tissue options available to filter
#' on. These are to be used in conjunction with the
#' \code{\link{add_expressed_genes}} function.
#' @usage data(tissue_options)
#' @source local data
#'
"tissue_options"
#' Different transcript biotypes that count as pseudogene
#'
#' These are the various transcript biotypes typically found in the transcript
#' biotype column of a gtf file.
#'
#' @format A character vector with all of the different pseudogene categories
#' \code{\link{get_exons}} function.
#' @usage data(pseudogene)
#' @source See exon_biotypes.R in data-raw
#'
"pseudogene"
