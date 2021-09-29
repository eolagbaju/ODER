if (!exists("gtf_path")) {
    gtf_url <- paste0(
        "http://ftp.ensembl.org/pub/release-103/gtf/",
        "homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz"
    )
    gtf_path <- file_cache(gtf_url)
}

if (!exists("gtf_grs")) {
    gtf_grs <- rtracklayer::import(gtf_path)
}

ex_opt_ers <- GenomicRanges::GRanges(
    seqnames = S4Vectors::Rle(c("chr21", "chr22"), c(2, 2)),
    ranges = IRanges::IRanges(
        start = c(5116369, 5118691, 5125879, 5128214),
        end = c(5117231, 5118847, 5125988, 5128403)
    )
)

liver_tissue <- get_tissue(tissue = "liver")
stomach_tissue <- get_tissue(tissue = "stomach")

livexpr_genes <- get_expressed_genes(gtf = gtf_grs, tissue_df = liver_tissue)

full_annot_liver_ers <- get_nearest_expressed_genes(annot_ers = ex_opt_ers, exp_genes = livexpr_genes, gtf = gtf_grs)

test_that("get_tissue works", {
    expect_equal(colnames(liver_tissue)[2], "liver")
    expect_equal(colnames(stomach_tissue)[2], "stomach")
})

test_that("get_expressed_genes works", {
    expect_true(methods::is(livexpr_genes, "GenomicRanges"))
    expect_true(all(grepl(pattern = "chr", x = GenomicRanges::seqnames(livexpr_genes))))
    expect_true("gene_id" %in% colnames(S4Vectors::mcols(livexpr_genes)))
})

test_that("get_nearest_expressed_genes works", {
    expect_equal(
        S4Vectors::mcols(full_annot_liver_ers)[["nearest_gene_v94_name"]],
        c("ENSG00000276612", "ENSG00000276612", "ENSG00000277248", "ENSG00000277248")
    )
    expect_equal(
        S4Vectors::mcols(full_annot_liver_ers)[["nearest_expressed_gene_v94_name"]],
        c("ENSG00000264462", "ENSG00000264462", "ENSG00000225480", "ENSG00000225480")
    )
})
