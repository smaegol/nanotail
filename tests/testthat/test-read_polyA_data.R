context("test functions for data import")

library(assertthat)
library(testthat)
library(assertive)

empty_sample_table = data.frame()

test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("path is correctly provided",{
  expect_error(read_polya_single(""))
  expect_error(read_polya_single("/foo/bar"))
  expect_error(read_polya_single(" ",ensembl="TRUE"))
})



test_that("Sample table is correctly provided",{
  expect_error(read_polya_multiple())
  expect_error(read_polya_multiple(empty_sample_table))
})

test_that("Filtering works",{
  expect_error(remove_failed_reads())
  expect_error(remove_failed_reads(empty_sample_table))
})
