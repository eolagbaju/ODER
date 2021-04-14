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
