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
#' @seealso \url{https://stackoverflow.com/questions/34445106/test-if-vector-is-contained-in-another-vector-including-repetitions}
#'
#' @keywords internal
#' @noRd
"%contain%" <- function(values, x) {
    tx <- table(x)
    tv <- table(values)
    z <- tv[names(tx)] - tx
    all(z >= 0 & !is.na(z))
}

#' Check if valid bigwig path(s) is passed in
#'
#' @param bws bigwig path or paths passed in
#'
#' @return TRUE unless an invalid path(s) is passed in
#'
#' @keywords internal
#' @noRd
bw_check <- function(bws) {
    if (length(bws) == 1 & stringr::str_sub(bws, -3, -1) != ".bw") {
        return(FALSE)
    } else if (length(bws) >= 1) {
        for (fp in bws) {
            if (stringr::str_sub(fp, -3, -1) != ".bw") {
                return(FALSE)
            }
        }
    }
    return(TRUE)
}

#' Modify chromosome format for get-coverage
#'
#' @param chr chromosome of interest
#' @param chr_format chromosome format, either chr or nochr
#'
#' @return string with unmodified chromosome format or modified
#' @keywords internal
#' @noRd
chr_formatting <- function(chr, chr_format) {
    if (chr_format == "chr") {
        return(chr)
    } else if (chr_format == "nochr") {
        mod_chr <- stringr::str_replace(chr, "chr", "")
        return(mod_chr)
    }
}
