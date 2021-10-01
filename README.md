
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

The goal of `ODER` is to **O**ptimise the **D**efinition of
**E**xpressed **R**egions. `ODER` is a packaged form of the method
developed in the Zhang et al. 2020 publication: [Incomplete annotation
has a disproportionate impact on our understanding of Mendelian and
complex neurogenetic
disorders](https://www.science.org/doi/10.1126/sciadv.aay8299). For a
more detailed explanation of using `ODER`, please see the
[vignette](https://eolagbaju.github.io/ODER/articles/ODERflow.html). For
more explanation of the methodology behind `ODER`, see the mehtods
section of the [original
publication](https://www.science.org/doi/10.1126/sciadv.aay8299).

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

## Citation

Below is the citation output from using `citation('ODER')` in R. Please
run this yourself to check for any updates on how to cite **ODER**.

``` r
message(citation("ODER"), bibtex = TRUE)
#> list(title = "Optimising the Definition of Expressed Regions", author = list(list(given = NULL, family = "eolagbaju", role = NULL, email = NULL, comment = NULL)), year = "2021", url = "http://www.bioconductor.org/packages/ODER", note = "https://github.com/eolagbaju/ODER/ODER - R package version 0.99.23", doi = "10.18129/B9.bioc.ODER")list(title = "Optimising the Definition of Expressed Regions", author = list(list(given = NULL, family = "eolagbaju", role = NULL, email = NULL, comment = NULL)), year = "2021", journal = "bioRxiv", doi = "10.1101/TODO", url = "https://www.biorxiv.org/content/10.1101/TODO")TRUE
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
