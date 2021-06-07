library(dasper)

megadepth::install_megadepth()
ref <- GenomicState::GenomicStateHub(version = "31", genome = "hg38", filetype = "TxDb")[[1]]

# Obtain the urls to the remotely hosted GTEx BigWig files
url <- recount::download_study(
  project = "SRP012682",
  type = "samples",
  download = FALSE
)

# cache the file using BiocFileCache for faster retrieval
bw_path <- dasper:::.file_cache(url[1])

junctions_example_1_path <-
  system.file(
    "extdata",
    "junctions_example_1.txt",
    package = "dasper",
    mustWork = TRUE
  )
junctions_example_2_path <-
  system.file(
    "extdata",
    "junctions_example_2.txt",
    package = "dasper",
    mustWork = TRUE
  )

# only keep chromosomes 21 + 22 for speed
junctions <-
  junction_load(
    junction_paths = c(
      junctions_example_1_path,
      junctions_example_2_path
    ),
    controls = "fibroblasts",
    chrs = c("chr21", "chr22")
  )

junctions <- junctions[, 1:5]

junctions <- junction_annot(junctions, ref)

head(SummarizedExperiment::rowData(junctions))
#SummarizedExperiment::rowRanges(junctions)
#SummarizedExperiment::assay(junctions)
#SummarizedExperiment::colData(junctions)