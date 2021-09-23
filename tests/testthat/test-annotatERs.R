if (!exists("gtf_path")) {
    gtf_url <- paste0(
        "http://ftp.ensembl.org/pub/release-103/",
        "gtf/homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz"
    )
    gtf_path <- ODER::file_cache(gtf_url)
}
if (!exists("gtf_grs")) {
    gtf_grs <- rtracklayer::import(gtf_path)
}
gtf_grs_ann <- gtf_grs
GenomeInfoDb::seqlevelsStyle(gtf_grs) <- "UCSC"
gtf_grs <- GenomeInfoDb::keepSeqlevels(gtf_grs, c("chr21", "chr22"), pruning.mode = "coarse")
exons_gr <- gtf_grs[gtf_grs$type == "exon"]
genes_gr <- gtf_grs[gtf_grs$type == "gene"]
if (!exists("rec_url")) {
    rec_url <- recount::download_study(
        project = "SRP012682",
        type = "samples",
        download = FALSE
    )
}

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
test_grs2 <- GenomicRanges::GRanges( # this is created to not overlap
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



chrs_to_keep <- c("21", "22")
#### preparing the txdb object
hg38_chrominfo <- GenomeInfoDb::getChromInfoFromUCSC("hg38")
new_info <- hg38_chrominfo$size[match(
    chrs_to_keep,
    GenomeInfoDb::mapSeqlevels(hg38_chrominfo$chrom, "Ensembl")
)]
names(new_info) <- chrs_to_keep
gtf_gr_tx <- GenomeInfoDb::keepSeqlevels(gtf_grs_ann,
    chrs_to_keep,
    pruning.mode = "tidy"
)
GenomeInfoDb::seqlengths(gtf_gr_tx) <- new_info
GenomeInfoDb::seqlevelsStyle(gtf_gr_tx) <- "UCSC"
GenomeInfoDb::genome(gtf_gr_tx) <- "hg38"

ucsc_txdb <- suppressWarnings(GenomicFeatures::makeTxDbFromGRanges(gtf_gr_tx))
test_gstate <- derfinder::makeGenomicState(txdb = ucsc_txdb)
ens_txdb <- ucsc_txdb
GenomeInfoDb::seqlevelsStyle(ens_txdb) <- "Ensembl"
################### end of txdb creation

test_juncs <- SummarizedExperiment::rowRanges(dasper::junctions_example)

test_er_aj <- suppressWarnings(get_junctions(
    opt_ers = test_grs,
    junc_data = test_juncs,
    txdb = ens_txdb
))

test_er_aj2 <- suppressWarnings(get_junctions(
    opt_ers = test_grs2,
    junc_data = test_juncs,
    txdb = ens_txdb
))

test_gene <- unique(unlist(GenomicRanges::mcols(
    GenomicRanges::mcols(test_er_aj[1])[["grl"]][[1]]
)[["gene_id_junction"]]))


test_annot_ers <- suppressWarnings(annotatERs(
    opt_ers = test_grs, junc_data = test_juncs,
    genom_state = test_gstate,
    gtf = gtf_grs_ann, txdb = ens_txdb
))

test_optgrs2 <- GenomicRanges::GRanges( # this is created to not overlap
    seqnames = S4Vectors::Rle(c("chr21"), c(6)),
    ranges = IRanges::IRanges(
        start = c(5040739, 46665406, 5038740, 46663060, 5032176, 46661574),
        end = c(5040744, 46665562, 5038775, 46663135, 5032217, 46661583),
    )
) # first 2 ranges intergenic, 2 introns, 2exons
test_annopt_ers <- suppressWarnings(annotatERs(
    opt_ers = test_optgrs2, junc_data = test_juncs,
    genom_state = test_gstate,
    gtf = gtf_grs_ann, txdb = ens_txdb
))


test_that("get_junctions works", {
    expect_true(methods::is(test_er_aj, "GenomicRanges"))
    expect_equal(length(IRanges::ranges(GenomicRanges::mcols(test_er_aj[1])[["grl"]][[1]])), 53)
    expect_equal(length(IRanges::ranges(GenomicRanges::mcols(test_er_aj[2])[["grl"]][[1]])), 0)
    expect_equal(length(IRanges::ranges(GenomicRanges::mcols(test_er_aj2[1:11])[["grl"]][[1]])), 0)
    expect_equal(
        IRanges::ranges(GenomicRanges::mcols(test_er_aj[3])[["grl"]][[1]])[1],
        IRanges::IRanges(start = 5026423, end = 5323718)
    )
    expect_equal(test_gene, "ENSG00000277117")
    # add test to make sure the junctions actually overlap
    expect_true(GenomicRanges::countOverlaps(test_er_aj[3], mcols(test_er_aj)[["grl"]][[3]]) > 0)
    # test when there is no overlaps the findoverlaps agrees
    expect_equal(unname(GenomicRanges::countOverlaps(test_er_aj[7], test_juncs)), 0)
    # test index that appears later to make sure indexing is consistent
})


test_that("annotatERs works", {
    expect_true("exon, intergenic, intron" %in% S4Vectors::mcols(test_annot_ers)[["annotation"]])
    expect_true("intron" %in% S4Vectors::mcols(test_annot_ers)[["annotation"]])
    expect_true("intergenic" %in% S4Vectors::mcols(test_annot_ers)[["annotation"]])
    expect_true("none" %in% S4Vectors::mcols(test_annot_ers)[["annotation"]])
    expect_true("genes" %in% names(S4Vectors::mcols(test_annot_ers)))
    expect_true("og_index" %in% names(S4Vectors::mcols(test_annot_ers)))
    expect_equal(unname(GenomicRanges::countOverlaps(test_annot_ers[1], S4Vectors::mcols(test_annot_ers[1])[["grl"]][[1]])), 53)
    # test that they have the right annotation
    # check if intergenic overlap with gtf (they shouldnt), see if exons have gene id associated?

    expect_true(all(GenomicRanges::countOverlaps(test_annopt_ers[S4Vectors::mcols(test_annopt_ers)[["annotation"]] == "intergenic"], exons_gr) == 0))
    expect_true(all(GenomicRanges::countOverlaps(test_annopt_ers[S4Vectors::mcols(test_annopt_ers)[["annotation"]] == "intergenic"], gtf_grs) == 0))
    expect_true(all(GenomicRanges::countOverlaps(test_annopt_ers[S4Vectors::mcols(test_annopt_ers)[["annotation"]] == "intron"], exons_gr) == 0))
    expect_true(all(GenomicRanges::countOverlaps(test_annopt_ers[S4Vectors::mcols(test_annopt_ers)[["annotation"]] == "intron"], gtf_grs) > 0))
    expect_true(all(GenomicRanges::countOverlaps(test_annopt_ers[S4Vectors::mcols(test_annopt_ers)[["annotation"]] == "exon"], exons_gr) > 0))

    expect_equal(IRanges::ranges(test_optgrs2[3]), IRanges::ranges(test_annopt_ers[3]))

    expect_equal(GenomicRanges::mcols(test_annot_ers)[["genes"]][[4]], "ENSG00000277117")
    expect_equal(GenomicRanges::mcols(test_annot_ers)[["gene_source"]][[4]], "nearest gtf genes")
    expect_true(GenomicRanges::mcols(GenomicRanges::distanceToNearest(test_grs, genes_gr))[4, 1] < 10000)
    expect_equal(GenomicRanges::mcols(test_annot_ers)[["genes"]][[5]], "")
    expect_equal(GenomicRanges::mcols(test_annot_ers)[["gene_source"]][[5]], "Too far")
    expect_true(GenomicRanges::mcols(GenomicRanges::distanceToNearest(test_grs, genes_gr))[5, 1] > 10000)
    expect_false("" %in% GenomicRanges::mcols(test_annot_ers)[["gene_source"]])
})
