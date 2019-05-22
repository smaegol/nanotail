context("Test plotting functions")

library(assertthat)
library(testthat)
library(assertive)


test_that("valid parameters are provided for plot_polya_distribution()",{

  expect_error(plot_polya_distribution(example_valid_polya_table))
  expect_error(plot_polya_distribution(example_valid_polya_table,groupingFactor = "non-existent_group"))
  expect_silent(plot_polya_distribution(example_valid_polya_table,groupingFactor = "group"))
  expect_error(plot_polya_distribution(example_valid_polya_table,groupingFactor = "group",scale_x_limit_low = 0))
  expect_silent(plot_polya_distribution(example_valid_polya_table,groupingFactor = "group",scale_x_limit_low = 0,scale_x_limit_high = 200))
  expect_warning(plot_polya_distribution(example_valid_polya_table,groupingFactor = "group",scale_x_limit_low = 0,scale_x_limit_high = 200,color_palette = "non_existent_palette"),"Unknown palette",all = FALSE,fixed=TRUE)
})
