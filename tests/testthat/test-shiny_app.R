context("Test shiny app")

library(assertthat)
library(testthat)
library(assertive)

test_that("valid parameters are provided for nanoTailApp()",{


  expect_error(nanoTailApp())
  expect_error(nanoTailApp(empty_polya_data_table))
  #expect_silent(nanoTailApp(example_valid_polya_table))
  expect_error(nanoTailApp(example_valid_polya_table %>% dplyr::select(-polya_length)))
  expect_error(nanoTailApp(example_valid_polya_table %>% dplyr::select(-transcript)))
  expect_error(nanoTailApp(example_valid_polya_table %>% dplyr::select(-sample_name)))
})
