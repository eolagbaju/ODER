#' Cache a file if it is not found locally
#'
#' `file_cache` will use [BiocFileCache][BiocFileCache::BiocFileCache-class]
#' will cache the file for faster repeated retrival, if it is not found locally
#' (i.e. a URL).
#'
#' @param file_path a path to file of interest.
#'
#' @return file_path of cached file or unchanged file_path if found locally.
#'
#' @export
#'
#' @examples
#' if (!exists("rec_url")) {
#'     rec_url <- recount::download_study(
#'         project = "SRP012682",
#'         type = "samples",
#'         download = FALSE
#'     )
#' }
#' eg_bwfile <- file_cache(rec_url[1])
#' eg_bwfile
file_cache <- function(file_path) {
    if (!file.exists(file_path)) {

        # suppress warning for tidyverse deprecated funs (select_() used over
        # select()) in BiocFileCache
        # exact = TRUE means exact match required, if F then uses regex search
        # suppressWarnings(
        file_path <- BiocFileCache::bfcrpath(BiocFileCache::BiocFileCache(
            ask = FALSE
        ),
        file_path,
        exact = TRUE
        )
        # )
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
#' @seealso \url{https://stackoverflow.com/questions/34445106/test-if-vector-is-
#' contained-in-another-vector-including-repetitions}
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
    if (length(bws) == 1) {
        if (stringr::str_sub(bws, -3, -1) != ".bw") {
            return(FALSE)
        }
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
    } else if (chr_format == "nochr") { # MT for mitochondrial
        if (grepl("M", chr)) {
            return("MT")
        } else {
            mod_chr <- stringr::str_replace(chr, "chr", "")
            return(mod_chr)
        }
    }
}

#' Modify chromosome format for get_chr_info
#'
#' @param chrs vector of desired chromosomes
#'
#' @return appropriately formatted chromosomes for
#' @keywords internal
#' @noRd
informatting <- function(chrs) {
    default <- c(
        "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
        "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
        "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"
    )

    if (length(chrs) == 1) {
        if (chrs == "") {
            return(default)
        }
        if (stringr::str_detect(chrs, "chr")) {
            return(chrs)
        }
        if (chrs == "M" | chrs == "MT") {
            return("chrM")
        } else {
            return(stringr::str_c("chr", as.character(chrs)))
        }
    }
    if (length(chrs) > 1) {
        if (stringr::str_detect(chrs[1], "chr")) {
            return(chrs)
        }
    }

    mod_chrs <- vector(mode = "character", length = length(chrs))
    counter <- 1
    for (chr in chrs) {
        chr <- as.character(chr)
        if (chr == "M" | chr == "MT") {
            mod_chrs[counter] <- "chrM"
        } else if (stringr::str_sub(chr, 1, 3) != "chr") {
            mod_chrs[counter] <- stringr::str_c("chr", chr)
        } else {
            mod_chrs[counter] <- chr
        }
        counter <- counter + 1
    }
    return(mod_chrs)
}

#' Modify chromosome format for annotatERs
#'
#' @param chrs vector of desired chromosomes
#'
#' @return appropriately formatted chromosomes for
#' @keywords internal
#' @noRd
informatting2 <- function(chrs) {
    default <- c(
        "1", "2", "3", "4", "5", "6", "7", "8", "9",
        "10", "11", "12", "13", "14", "15", "16", "17",
        "18", "19", "20", "21", "22", "X", "Y", "MT"
    )


    if (length(chrs) == 1) {
        if (chrs == "") {
            return(default)
        }
        if (!(stringr::str_detect(chrs, "chr"))) {
            return(chrs)
        }
        if (chrs == "M" | chrs == "MT") {
            return("MT")
        }
    }
    if (length(chrs) > 1) {
        if (!(stringr::str_detect(chrs[1], "chr"))) {
            return(chrs)
        }
    }

    mod_chrs <- vector(mode = "character", length = length(chrs))
    counter <- 1
    for (chr in chrs) {
        chr <- as.character(chr)
        if (chr == "chrM" | chr == "chrMT") {
            mod_chrs[counter] <- "MT"
        } else if (stringr::str_sub(chr, 1, 3) == "chr") {
            mod_chrs[counter] <- stringr::str_remove(chr, "chr")
        } else {
            mod_chrs[counter] <- chr
        }
        counter <- counter + 1
    }
    return(mod_chrs)
}

#' Check if value is between two other values or is within a range
#'
#' @param value value to check
#' @param rstart start of range
#' @param rend end of range
#'
#' @return TRUE or FALSE
#' @keywords internal
#' @noRd
inbetween <- function(value, rstart, rend) {
    if ((value > rstart) & (value < rend)) {
        return(TRUE)
    } else if (value == rstart) {
        return(TRUE)
    } else if (value == rend) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

#' Finds overlapping genomic ranges (specifically junctions)
#'
#' @param x metadata containing one or two genomic ranges
#'
#' @return TRUE or FALSE
#' @keywords internal
#' @noRd
colgrs <- function(x) {
    if (length(x) == 1) {
        return(FALSE)
    } else if (GenomicRanges::countOverlaps(x[1], x[2]) > 0) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

#' Finds non-overlapping genomic ranges (specifically junctions)
#'
#' @param x metadata containing one or two genomic ranges
#'
#' @return TRUE or FALSE
#' @keywords internal
#' @noRd
inv_colgrs <- function(x) {
    if (length(x) == 1) {
        return(TRUE)
    } else if (GenomicRanges::countOverlaps(x[1], x[2]) > 0) {
        return(FALSE)
    } else {
        return(TRUE)
    }
}

#' Checks for and returns a GTF file in the form of a Genomic Ranges object
#'
#' @param a gtf in GRanges form or a gtf file path
#'
#' @return gtf genomic ranges
#' @keywords internal
#' @noRd
gtf_load <- function(gtf) {
    if (methods::is(gtf, "GenomicRanges")) {
        return(gtf)
    } else if (is.character(gtf)) {
        return(rtracklayer::import(gtf))
    } else {
        stop("Invalid gtf argument passed in")
    }
}
