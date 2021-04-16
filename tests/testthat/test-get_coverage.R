test_that("get_chr_info works", {
  all_UCSC_chr <- GenomeInfoDb::getChromInfoFromUCSC("hg38")[["chrom"]]
  chrs <- c("chr1","chr2","chr3")
  test_info <- get_chr_info(chrs=chrs,genome="hg38")
  expect_true(all_UCSC_chr %contain% "chr1")
  expect_true(all_UCSC_chr %contain% chrs)
  expect_error(get_chr_info(chrs="chromononsense",genome = "hg38"),"Non-UCSC chromosome format was entered")
  expect_equal(nrow(test_info),length(chrs))
  expect_type(test_info[[1]],"character")
  expect_type(test_info[[2]],"integer")
  expect_type(test_info[[3]],"logical")
  expect_type(test_info[[4]],"logical")

})
