#' Cache a file if it is not found locally
#'
#' `.file_cache` will use [BiocFileCache][BiocFileCache::BiocFileCache-class]
#' will cache the file for faster repeated retrival, if it not found locally
#' (i.e. a URL).
#'
#' @param file_path a path to file of interest.
#'
#' @return file_path of cached file or unchanged file_path if found locally.
#'
#' @keywords internal
#' @noRd
.file_cache <- function(file_path) {
  if (!file.exists(file_path)) {

    # suppress warning for tidyverse deprecated funs (select_() used over select()) in BiocFileCache
    # exact = TRUE means exact match required, if F then uses regex search
    suppressWarnings(
      file_path <- BiocFileCache::bfcrpath(BiocFileCache::BiocFileCache(ask = FALSE),
                                           file_path,
                                           exact = TRUE
      )
    )
  }

  return(file_path)
}

#' Check if larger vector contains a smaller vector
#'
#' @param values bigger vector to search
#' @param x smaller vector to search for 
#'
#' @return
#'
#' @keywords internal
#' @noRd
"%contain%" <- function(values,x) {
  tx <- table(x)
  tv <- table(values)
  z <- tv[names(tx)] - tx
  all(z >= 0 & !is.na(z))
} #https://stackoverflow.com/questions/34445106/test-if-vector-is-contained-in-another-vector-including-repetitions
