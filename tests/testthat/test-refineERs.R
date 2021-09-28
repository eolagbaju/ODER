gr1 <- GenomicRanges::GRanges("chr21", IRanges::IRanges(5026423, 5323718), strand = "+")
gr2 <- GenomicRanges::GRanges("chr21", IRanges::IRanges(c(13261708, 13482787), c(13482707, 13851281)), strand = "-")
gr3 <- GenomicRanges::GRanges("chr21", IRanges::IRanges(c(18365503, 18365938), c(18365600, 18857780)), strand = "-")
gr4 <- GenomicRanges::GRanges("chr21", IRanges::IRanges(c(18858010, 19351750), c(19350050, 19351960)), strand = "+")
gr5 <- GenomicRanges::GRanges("chr21", IRanges::IRanges(19408735, 19772936), strand = "+")
gr6 <- GenomicRanges::GRanges("chr21", IRanges::IRanges(20091736, 20179922), strand = "-")
gr7 <- GenomicRanges::GRanges("chr21", IRanges::IRanges(21091736, 21091922), strand = "-")
gr8 <- GenomicRanges::GRanges("chr21", IRanges::IRanges(c(22136878, 22150426), c(22150523, 22161511)), strand = "+")
gr9 <- GenomicRanges::GRanges("chr21", IRanges::IRanges(c(25085652, 25086122), c(25086438, 25087197)), strand = "+")
grl <- GenomicRanges::GRangesList(gr1, gr2, gr3, gr4, gr5, gr6, gr7, gr8, gr9)
test_annot_ers <- GenomicRanges::GRanges(
    seqnames = S4Vectors::Rle(c("chr21"), c(9)),
    ranges = IRanges::IRanges(
        start = c(5026400, 13482700, 18365500, 19350000, 19772900, 20179900, 21091800, 22150500, 25086100),
        end = c(5026500, 13482800, 18366000, 19352000, 19773000, 20180000, 211875305, 22150540, 25086150),
        names = c(1:9)
    ),
    strand = S4Vectors::Rle((c("+", "-", "-", "+", "+", "-", "-", "+", "+"))),
    grl = grl,
    annotation = c("intron", "intron", "intron", "intergenic", "intergenic", "intergenic", "exon", "intergenic", "intron")
)

test_annot_ers2 <- test_annot_ers[S4Vectors::mcols(test_annot_ers)[["annotation"]] %in% c("intron", "intergenic")]
# getting the ers that only have one or two overlapping junctions
test_annot_ers2 <- test_annot_ers2[lengths(S4Vectors::mcols(test_annot_ers2)[["grl"]]) == 2 |
    lengths(S4Vectors::mcols(test_annot_ers2)[["grl"]]) == 1]

test_refine_results <- modify_ers(test_annot_ers2)

test_parref_ers <- test_refine_results[[1]]

test_changes <- test_refine_results[[2]]

test_refined_ers <- refine_ERs(test_annot_ers)

test_that("refine_ERs works", {
    expect_true(methods::is(test_refined_ers, "GenomicRanges"))
    expect_true(methods::is(test_changes, "logical"))
    expect_equal(length(test_refine_results), 2)
    expect_equal(length(test_refined_ers), 6)
    expect_true(all(lengths(S4Vectors::mcols(test_parref_ers)[["grl"]]) <= 2))
    expect_true("intron" %in% S4Vectors::mcols(test_refined_ers)[["annotation"]])
    expect_true("intergenic" %in% S4Vectors::mcols(test_refined_ers)[["annotation"]])
    expect_false("exon" %in% S4Vectors::mcols(test_refined_ers)[["annotation"]])
    expect_equal(BiocGenerics::end(test_refined_ers[1]), BiocGenerics::start(S4Vectors::mcols(test_annot_ers)[["grl"]][[1]]) - 1)
    expect_equal(BiocGenerics::start(test_refined_ers[2]), BiocGenerics::end(S4Vectors::mcols(test_annot_ers)[["grl"]][[2]])[1] + 1)
    expect_equal(BiocGenerics::end(test_refined_ers[2]), BiocGenerics::start(S4Vectors::mcols(test_annot_ers)[["grl"]][[2]])[2] - 1)
    expect_equal(BiocGenerics::end(test_annot_ers["9"]), BiocGenerics::end(test_parref_ers["9"]))
})
