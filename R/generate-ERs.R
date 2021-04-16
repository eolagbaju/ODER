
# usethis::use_pipe(export = TRUE)

#' Generating the mean coverage of the expressed regions
#'
#' @param bw_paths
#' @param auc_raw
#' @param auc_target
#' @param chrs
#' @param genome
#'
#' @return
#' @export
#'
#' @examples
get_coverage <- function(bw_paths,auc_raw,auc_target,chrs,genome="hg38"){
  chr_info <- GenomeInfoDb::getChromInfoFromUCSC(genome) %>%
    dplyr::filter(chrom %in% chrs)

  all_chrs_mean_cov <- list()

  print(stringr::str_c(Sys.time(), " - Obtaining mean coverage across ", length(bw_paths), " samples")) # sub for a rutils function maybe?

  for(i in 1:nrow(chr_info)){

    print(stringr::str_c(Sys.time(), " - ", chr_info[["chrom"]][i]))# sub for a rutils function maybe?
    #loading coverage information for designated chromosomes and merging them into a dataframe
    chr_mean_cov <-
      derfinder::loadCoverage(files = bw_paths,
                              totalMapped = auc_raw, # normalise by auc here as for bws, more accurate since Rail-RNA clips reads
                              targetSize = auc_target,
                              chr = chr_info$chrom[i],
                              chrlen = chr_info$size[i], 
                              inputType = "BigWig",
                              returnMean = T,
                              returnCoverage = F,
                              verbose = F,
                              cutoff = NULL) # setting cutoff as NULL here and instead to be applied in findRegions()
    #storing the mean coverage in a list
    all_chrs_mean_cov[[chr_info[["chrom"]][i]]] <- chr_mean_cov["meanCoverage"]

  }

  return(all_chrs_mean_cov)

}
