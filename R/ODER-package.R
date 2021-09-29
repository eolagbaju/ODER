#' ODER: Optimising the Definition of Expressed Regions
#'
#' The aim of ODER is to identify previously unannotated expressed regions (ERs)
#' using RNA-sequencing data. For this purpose, ODER defines and optimises the
#' definition of ERs, then connected these ERs to genes using junction data. In
#' this way, ODER improves gene annotation. Gene annotation is a staple input of
#' many bioinformatic pipelines and a more complete gene annotation can enable
#' more accurate interpretation of disease associated variants.
#'
#' @docType package
#' @name ODER
NULL

#' @importFrom data.table :=
#' @importFrom IRanges ranges
#' @importFrom BiocGenerics start end
#' @importFrom S4Vectors mcols
NULL
