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
  expect_error(read_polya_single(" ",gencode="FALSE"))
  expect_error(read_polya_single(" ",gencode=FALSE,sample_name="sample1"))
  expect_error(read_polya_single(sample_tempfile1,gencode="FALSE",sample_name="sample1"))
})

test_that("polya data are correctly read from files",{
  empty_tempfile = tempfile()
  expect_error(read_polya_single(empty_tempfile))
  expect_message(read_polya_single(sample_tempfile1),"Loading data from",fixed=TRUE)
  expect_message(read_polya_single(sample_tempfile1,gencode=TRUE),"Loading data from",fixed=TRUE)
  expect_warning(read_polya_single(sample_tempfile1,gencode=TRUE,sample_name="sample1"),"sample_name was provided in the input file",fixed=TRUE,all=FALSE)
  expect_message(read_polya_multiple(example_sample_table),"Loading data from",fixed=TRUE)
  expect_length(rownames(read_polya_multiple(example_sample_table)),number_of_reads_per_sample*2)
})



test_that("Sample table is correctly provided",{
  expect_error(read_polya_multiple())
  expect_error(read_polya_multiple(samples_table = empty_sample_table))
  expect_error(read_polya_multiple(samples_table = sample_table_with_wrong_columns))
  expect_error(read_polya_multiple(samples_table = sample_table_with_correct_columns_wrong_data))
})

test_that("Filtering works",{
  example_polya_table_read_from_disk <- read_polya_multiple(samples_table = example_sample_table)
  expect_error(remove_failed_reads())
  expect_error(remove_failed_reads(empty_sample_table))
  expect_silent(remove_failed_reads(example_valid_polya_table))
})
