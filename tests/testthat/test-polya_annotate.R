context("annotation of polya predictions table")

test_that("annotation with annotables works",{

  expect_error(annotate_with_annotables())
  expect_error(annotate_with_annotables(empty_polya_data_table))

  expect_error(annotate_with_annotables(example_valid_polya_table))
  expect_error(annotate_with_annotables(empty_polya_data_table,"grcm38"))

})
