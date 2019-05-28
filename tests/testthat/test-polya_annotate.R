context("annotation of polya predictions table")

test_that("annotation with annotables works",{

  expect_error(annotate_with_annotables())
  expect_error(annotate_with_annotables(empty_polya_data_table))

  expect_error(annotate_with_annotables(example_valid_polya_table))
  expect_error(annotate_with_annotables(empty_polya_data_table,"grcm38"))

})


test_that("annotation with biomart works",{



  number_of_reads_per_sample=20000
  number_of_transcripts_per_sample=50
  #mouse_genes_ensembl_mart <- biomaRt::useMart("ensembl",dataset="mmusculus_gene_ensembl")



  expect_error(annotate_with_biomart())
  expect_error(annotate_with_biomart(empty_polya_data_table))

  expect_error(annotate_with_biomart(example_valid_polya_table))
  expect_error(annotate_with_biomart(empty_polya_data_table))


})
