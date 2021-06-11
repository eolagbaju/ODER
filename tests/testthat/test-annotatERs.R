# SummarizedExperiment::rowRanges(dasper::junctions_example)
# mcols(mcols(test_er_aj[1])[["grl"]][[1]])[["type"]]
gtf_url <- paste0(
    "http://ftp.ensembl.org/pub/release-103/",
    "gtf/homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz"
)

gtf_path <- ODER:::.file_cache(gtf_url)

test_grs <- GenomicRanges::GRanges(
    seqnames = S4Vectors::Rle(c("chr21", "chr22"), c(10, 1)),
    ranges = IRanges::IRanges(
        start = c(5026423, 24738, 5032218, 5033895, 17554, 50800446, 50800539, 16570, 50800790, 15005, 20312),
        end = c(5323718, 24891, 5033407, 5033980, 17728, 50800910, 50800817, 16723, 50800910, 15038, 20582),
        names = head(letters, 11)
    ),
    strand = S4Vectors::Rle((c("+", "-")), c(6, 5)),
    score = 1:11,
    GC = seq(1, 0, length = 11)
)
test_grs2 <- GenomicRanges::GRanges(
    seqnames = S4Vectors::Rle(c("chr1", "chr2"), c(10, 1)),
    ranges = IRanges::IRanges(
        start = c(5026423, 24738, 5032218, 5033895, 17554, 50800446, 50800539, 16570, 50800790, 15005, 20312),
        end = c(5323718, 24891, 5033407, 5033980, 17728, 50800910, 50800817, 16723, 50800910, 15038, 20582),
        names = head(letters, 11)
    ),
    strand = S4Vectors::Rle((c("+", "-")), c(6, 5)),
    score = 1:11,
    GC = seq(1, 0, length = 11)
)

test_juncs <- SummarizedExperiment::rowRanges(dasper::junctions_example)
# GenomeInfoDb::seqlevelsStyle(test_juncs) <- "NCBI"

test_er_aj <- suppressWarnings(get_junctions(
    opt_ers = test_grs,
    junc_data = test_juncs,
    gtf_path = gtf_path
))

test_er_aj2 <- suppressWarnings(get_junctions(
    opt_ers = test_grs2,
    junc_data = test_juncs,
    gtf_path = gtf_path
))

test_gene <- unique(unlist(GenomicRanges::mcols(GenomicRanges::mcols(test_er_aj[1])[["grl"]][[1]])[["gene_id_junction"]]))
# test if the specific ranges are caught and geneid
# dasper example junctions don't appear to have genes labelled?
# ranges(mcols(test_er_aj[3])[["grl"]][[1]])[1]

test_that("get_junctions works", {
    expect_true(methods::is(test_er_aj, "GenomicRanges"))
    expect_equal(length(IRanges::ranges(GenomicRanges::mcols(test_er_aj[1])[["grl"]][[1]])), 53)
    expect_equal(length(IRanges::ranges(GenomicRanges::mcols(test_er_aj[2])[["grl"]][[1]])), 0)
    expect_equal(length(IRanges::ranges(GenomicRanges::mcols(test_er_aj2[1:10])[["grl"]][[1]])), 0)
    expect_equal(IRanges::ranges(GenomicRanges::mcols(test_er_aj[3])[["grl"]][[1]])[1], IRanges::IRanges(start = 5026423, end = 5323718))
    expect_equal(test_gene, "ENSG00000277117")
})
