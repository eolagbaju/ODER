if (!exists("gtex_metadata")) {
    gtex_metadata <- recount::all_metadata("gtex")
    gtex_metadata <- gtex_metadata %>%
        as.data.frame() %>%
        dplyr::filter(project == "SRP012682")
}
if (!exists("gtf_path")) {
    gtf_url <- paste0(
        "http://ftp.ensembl.org/pub/release-103/",
        "gtf/homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz"
    )
    gtf_path <- ODER:::.file_cache(gtf_url)
}
gtf_grs <- rtracklayer::import(gtf_path)
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
test_bw_path <- ODER:::.file_cache(rec_url[1])
if (!exists("test_opt_ers1")) {
    test_opt_ers1 <- suppressWarnings(ODER(
        bw_paths = test_bw_path, auc_raw = gtex_metadata[["auc"]][1],
        auc_target = 40e6 * 100, chrs = c("chr21", "chr22"),
        genome = "hg38", mccs = c(5, 10), mrgs = c(10, 20),
        gtf = gtf_path, ucsc_chr = TRUE, ignore.strand = TRUE,
        exons_no_overlap = NULL, bw_chr = "chr"
    ))
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

test_juncs <- SummarizedExperiment::rowRanges(dasper::junctions_example)

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

test_gene <- unique(unlist(GenomicRanges::mcols(
    GenomicRanges::mcols(test_er_aj[1])[["grl"]][[1]]
)[["gene_id_junction"]]))

test_gstate <- suppressWarnings(generate_genomic_state(
    gtf = gtf_path,
    chrs_to_keep = c("21", "22"),
    ensembl = TRUE
))

test_annot_ers <- suppressWarnings(annotatERs(
    opt_ers = test_grs, junc_data = test_juncs,
    gtf_path = gtf_path, chrs_to_keep = c("chr21", "chr22")
))
if (!exists("test_annot_opters1")) {
    test_annot_opters1 <- suppressWarnings(annotatERs(
        opt_ers = test_opt_ers1[["opt_ers"]], junc_data = test_juncs,
        gtf_path = gtf_path, chrs_to_keep = c("chr21", "chr22")
    ))
}


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

test_that("generate_genomic_state works", {
    expect_true(methods::is(test_gstate[[1]], "GenomicRanges"))
    expect_true(methods::is(test_gstate[[2]], "GenomicRanges"))
    expect_equal(length(test_gstate), 2)
    expect_equal(unique(seqnames(test_gstate[[1]])), factor(c("chr21", "chr22")))
    expect_equal(names(test_gstate), c("fullGenome", "codingGenome"))
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
    expect_true(all(GenomicRanges::countOverlaps(test_annot_opters1[S4Vectors::mcols(test_annot_opters1)[["annotation"]] == "intergenic"], exons_gr) == 0))
    expect_true(all(GenomicRanges::countOverlaps(test_annot_opters1[S4Vectors::mcols(test_annot_opters1)[["annotation"]] == "intergenic"], gtf_grs) == 0))
    expect_true(all(GenomicRanges::countOverlaps(test_annot_opters1[S4Vectors::mcols(test_annot_opters1)[["annotation"]] == "intron"], exons_gr) == 0))
    expect_true(all(GenomicRanges::countOverlaps(test_annot_opters1[S4Vectors::mcols(test_annot_opters1)[["annotation"]] == "intron"], gtf_grs) > 0))
    expect_true(all(GenomicRanges::countOverlaps(test_annot_opters1[S4Vectors::mcols(test_annot_opters1)[["annotation"]] == "exon"], exons_gr) > 0))
    expect_equal(IRanges::ranges(test_opt_ers1[["opt_ers"]][5479]), IRanges::ranges(test_annot_opters1[5479]))
    expect_equal(GenomicRanges::mcols(test_annot_ers)[["genes"]][[4]], "ENSG00000277117")
    expect_equal(GenomicRanges::mcols(test_annot_ers)[["gene_source"]][[4]], "nearest gtf genes")
    expect_true(GenomicRanges::mcols(GenomicRanges::distanceToNearest(test_grs, genes_gr))[4, 1] < 10000)
    expect_equal(GenomicRanges::mcols(test_annot_ers)[["genes"]][[5]], "")
    expect_equal(GenomicRanges::mcols(test_annot_ers)[["gene_source"]][[5]], "Too far")
    expect_true(GenomicRanges::mcols(GenomicRanges::distanceToNearest(test_grs, genes_gr))[5, 1] > 10000)
    expect_false("" %in% GenomicRanges::mcols(test_annot_ers)[["gene_source"]])
})
