context("test functions for data import and filtering")

library(assertthat)
library(testthat)
library(assertive)

empty_sample_table = data.frame()

sample_table_with_correct_columns_wrong_data = data.frame(sample_name = c("iorem","ipsum"),polya_path=c("iorem","ipsum"))
sample_table_with_wrong_columns = data.frame(sample_name_wrong = c("iorem","ipsum"),polya_path_wrong=c("iorem","ipsum"))

test_that("read_polya_single correctly parses provided parameters",{
  expect_error(read_polya_single(""))
  expect_error(read_polya_single("/foo/bar"))
  expect_error(read_polya_single(" ",gencode="TRUE"))
})




test_that("Sample table is correctly provided",{
  expect_error(read_polya_multiple())
  expect_error(read_polya_multiple(samples_table = empty_sample_table))
  expect_error(read_polya_multiple(samples_table = sample_table_with_wrong_columns))
  expect_error(read_polya_multiple(samples_table = sample_table_with_correct_columns_wrong_data))
})

test_that("Filtering works",{
  expect_error(remove_failed_reads())
  expect_error(remove_failed_reads(empty_sample_table))
})
