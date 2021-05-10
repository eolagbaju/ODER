#' An example AUC value
#'
#' An AUC(Area Under Coverage) value for a user to try out the package and to 
#' pass in for tests. From the GTEX data set and project SRP012682, the actual
#' value is 11872688252.
#'
#' @format A numeric value
#'
#' @source See example.R in data-raw
#'
"auc_example"
#' An example coverage
#'
#' An example coverage generated for a user to try out the package and
#' to pass in for tests. Coverage of chromosomes 21 and 22 from the project
#' SRP012682.
#'
#' @format A list of length 2 containing 2 Rles for chromosomes 21 and 22 respectively
#'
#' @source See example.R in data-raw
#'
"coverage_example"
#' An example set of Expressed Regions
#'
#' An example set of Expressed Regions generated for a user to try out the 
#' package and to pass in for tests. Generated using coverage_example and MCCs 
#' of 5 and 10 and MRGs of 10 and 20. 
#'
#' @format A list containing two lists (for each mcc) each with a set of genomic ranges
#' for the different combinations of mcc and mrg
#'
#' @source See example.R in data-raw
#'
"ers_example"
#' An example set of ER deltas
#'
#' This set of deltas was calculated using ers_example and exons from ensembl.
#'
#' @format A tibble/dataframe with the sums, means, medians, n_eq_0 and
#' propor_eq_0 for each combination of mccs (5 & 10) and mrgs (10 & 20)
#'
#' @source See example.R in data-raw
#'
"ers_delta_example"
