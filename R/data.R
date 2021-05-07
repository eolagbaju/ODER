#' An example AUC value
#'
#' An AUC(Area Under Coverage) value for a user to try out the package and to pass in for tests
#'
#' @format A numeric value
#'
#' @source See example.R in data-raw
#'
"auc_example"
#' An example coverage
#'
#' An example coverage generated for a user to try out the package and
#' to pass in for tests
#'
#' @format A list of length 2 containing 2 Rles for chromosomes 21 and 22 respectively
#'
#' @source See example.R in data-raw
#'
"coverage_example"
#' An example set of Expressed Regions
#'
#' An example set of Expressed Regions generated for a user to try out the package and
#' to pass in for tests
#'
#' @format A list containing two lists (for each mcc) each with a set of genomic ranges
#' for the different combinations of mcc and mrg
#'
#' @source See example.R in data-raw
#'
"ers_example"
#' An example set of ER deltas
#'
#' This set of deltas was calculated using Expressed Regions from chromosomes 21
#' and 22 and exons from the ensemble hg38 set
#'
#' @format A tibble/dataframe with the sums, means, medians, n_eq_0 and
#' propor_eq_0 for each combination of mccs (5 & 10) and mrgs (10 & 20)
#'
#' @source See example.R in data-raw
#'
"ers_delta_example"
