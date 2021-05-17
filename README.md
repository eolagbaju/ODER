
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ODER

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![BioC
status](http://www.bioconductor.org/shields/build/release/bioc/ODER.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/ODER)
[![R-CMD-check-bioc](https://github.com/eolagbaju/ODER/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/eolagbaju/ODER/actions)
[![Codecov test
coverage](https://codecov.io/gh/eolagbaju/ODER/branch/master/graph/badge.svg)](https://codecov.io/gh/eolagbaju/ODER?branch=master)
<!-- badges: end -->

The goal of `ODER` is to well define expressed regions from
RNA-squencing experiments so that they can be confidently identified as
either currently annotated exons or previously unannotated exons. If the
exons are unannotated, `ODER` also enables you to annotate these
unannotated expressed regions.

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `ODER` using from
[Bioconductor](http://bioconductor.org/) the following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("ODER")
```

And the development version from
[GitHub](https://github.com/eolagbaju/ODER) with:

``` r
BiocManager::install("eolagbaju/ODER")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library("ODER")
url <- recount::download_study(
    project = "SRP012682",
    type = "samples",
    download = FALSE
)
#> Setting options('download.file.method.GEOquery'='auto')
#> Setting options('GEOquery.inmemory.gpl'=FALSE)
# .file_cache is an internal function to download a bigwig file from a link if the file has been downloaded recently, it will be retrieved from a cache
bw_path <- ODER:::.file_cache(url[1])
gtf_url <- "http://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz"
gtf_path <- ODER:::.file_cache(gtf_url)

opt_ers <- ODER(
    bw_paths = bw_path, auc_raw = auc_example,
    auc_target = 40e6 * 100, chrs = c("chr20", "chr21", "chr22"),
    genome = "hg38", mccs = c(2, 4, 6, 8, 10), mrgs = c(10, 20, 30),
    gtf = gtf_path, ucsc_chr = TRUE, ignore.strand = TRUE,
    exons_no_overlap = NULL, bw_chr = "chr"
)
#> Loading required package: BiocGenerics
#> Loading required package: parallel
#> 
#> Attaching package: 'BiocGenerics'
#> The following objects are masked from 'package:parallel':
#> 
#>     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
#>     clusterExport, clusterMap, parApply, parCapply, parLapply,
#>     parLapplyLB, parRapply, parSapply, parSapplyLB
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from 'package:base':
#> 
#>     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
#>     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
#>     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
#>     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
#>     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
#>     union, unique, unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> Loading required package: stats4
#> 
#> Attaching package: 'S4Vectors'
#> The following object is masked from 'package:base':
#> 
#>     expand.grid
#> [1] "2021-05-17 16:25:43 - Obtaining mean coverage across 1 samples"
#> [1] "2021-05-17 16:25:43 - chr20"
#> [1] "2021-05-17 16:25:44 - chr21"
#> [1] "2021-05-17 16:25:45 - chr22"
#> [1] "2021-05-17 16:25:46 - Generating ERs for chr20"
#> [1] "2021-05-17 16:25:52 - Generating ERs for chr21"
#> [1] "2021-05-17 16:25:55 - Generating ERs for chr22"
#> [1] "2021-05-17 16:26:00 - Loading in GTF..."
#> [1] "2021-05-17 16:26:47 - Obtaining non-overlapping exons"
#> [1] "2021-05-17 16:26:49 - Calculating delta for ERs..."
#> [1] "2021-05-17 16:26:54 - Obtaining optimal set of ERs..."

opt_ers
#> $opt_ers
#> GRanges object with 25393 ranges and 0 metadata columns:
#>           seqnames            ranges strand
#>              <Rle>         <IRanges>  <Rle>
#>       [1]    chr20     167593-167621      *
#>       [2]    chr20     191726-191774      *
#>       [3]    chr20     269685-269730      *
#>       [4]    chr20     271105-271130      *
#>       [5]    chr20     271210-271267      *
#>       ...      ...               ...    ...
#>   [25389]    chr22 50798779-50798835      *
#>   [25390]    chr22 50798884-50799149      *
#>   [25391]    chr22 50799209-50799284      *
#>   [25392]    chr22 50799669-50799744      *
#>   [25393]    chr22 50800460-50800587      *
#>   -------
#>   seqinfo: 3 sequences from an unspecified genome; no seqlengths
#> 
#> $opt_mcc_mrg
#> [1] "mcc_10" "mrg_30"
#> 
#> $deltas
#> # A tibble: 15 x 7
#>      mcc   mrg     sum  mean median n_eq_0 propor_eq_0
#>    <dbl> <dbl>   <int> <dbl>  <dbl>  <int>       <dbl>
#>  1     2    10 4718631  870.  149      795       0.147
#>  2     2    20 4250197  830.  145      805       0.157
#>  3     2    30 3839111  794.  138      814       0.168
#>  4     4    10 4395515  901.  132      814       0.167
#>  5     4    20 3868780  847.  123      835       0.183
#>  6     4    30 3426418  793.  116      847       0.196
#>  7     6    10 4271226  953.  130      835       0.186
#>  8     6    20 3583476  868.  113      852       0.206
#>  9     6    30 3139405  806.  107      870       0.223
#> 10     8    10 4069944  981.  126      839       0.202
#> 11     8    20 3406090  899.  111      858       0.227
#> 12     8    30 2917974  822.  100      884       0.249
#> 13    10    10 3891345 1011.  124      798       0.207
#> 14    10    20 3181382  909.  106      827       0.236
#> 15    10    30 2790475  854.   94.5    853       0.261
```

<img src="man/figures/README-plot_example-1.png" width="100%" />

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date.

You can also embed plots, for example:

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub\!

## Citation

Below is the citation output from using `citation('ODER')` in R. Please
run this yourself to check for any updates on how to cite **ODER**.

``` r
# print(citation('ODER'), bibtex = TRUE)
```

Please note that the `ODER` was only made possible thanks to many other
R and bioinformatics software authors, which are cited either in the
vignettes and/or the paper(s) describing this package.

## Code of Conduct

Please note that the `ODER` project is released with a [Contributor Code
of Conduct](http://bioconductor.org/about/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.

## Development tools

  - Continuous code testing is possible thanks to [GitHub
    actions](https://www.tidyverse.org/blog/2020/04/usethis-1-6-0/)
    through *[usethis](https://CRAN.R-project.org/package=usethis)*,
    *[remotes](https://CRAN.R-project.org/package=remotes)*, and
    *[rcmdcheck](https://CRAN.R-project.org/package=rcmdcheck)*
    customized to use [Bioconductor’s docker
    containers](https://www.bioconductor.org/help/docker/) and
    *[BiocCheck](https://bioconductor.org/packages/3.12/BiocCheck)*.
  - Code coverage assessment is possible thanks to
    [codecov](https://codecov.io/gh) and
    *[covr](https://CRAN.R-project.org/package=covr)*.
  - The [documentation website](http://eolagbaju.github.io/ODER) is
    automatically updated thanks to
    *[pkgdown](https://CRAN.R-project.org/package=pkgdown)*.
  - The code is styled automatically thanks to
    *[styler](https://CRAN.R-project.org/package=styler)*.
  - The documentation is formatted thanks to
    *[devtools](https://CRAN.R-project.org/package=devtools)* and
    *[roxygen2](https://CRAN.R-project.org/package=roxygen2)*.

For more details, check the `dev` directory.

This package was developed using
*[biocthis](https://bioconductor.org/packages/3.12/biocthis)*.
