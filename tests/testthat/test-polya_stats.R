context("test functions for statistics calculation")


library(assertthat)
library(testthat)
library(assertive)



test_that("Valid input parameters provided for calculate_polya_stats()",{
  expect_error(calculate_polya_stats())
  expect_error(calculate_polya_stats(empty_polya_data_table))

  expect_silent(calculate_polya_stats(example_valid_polya_table,grouping_factor="group"))
  expect_silent(calculate_polya_stats(example_valid_polya_table,grouping_factor="group"))
  expect_silent(calculate_polya_stats(example_valid_polya_table,grouping_factor="group"))
  expect_silent(calculate_polya_stats(example_valid_polya_table,grouping_factor="group"))
  expect_error(calculate_polya_stats(example_valid_polya_table,grouping_factor="group",stat_test = "ttest"),"Please provide one of available statistical tests")

  expect_silent(calculate_polya_stats(example_valid_polya_table,grouping_factor="group",stat_test="glm",custom_glm_formula = "polya_length ~ group"))
  expect_silent(calculate_polya_stats(example_valid_polya_table,grouping_factor="group",stat_test="glm",custom_glm_formula = polya_length ~ group))
  expect_error(calculate_polya_stats(example_valid_polya_table,grouping_factor="group",stat_test="glm",custom_glm_formula = polya_length ~ group + group2))
  expect_error(calculate_polya_stats(example_valid_polya_table,grouping_factor="group",stat_test="glm",custom_glm_formula = polya_length ~ NULL))

  expect_silent(calculate_polya_stats(example_valid_polya_table_3levels,grouping_factor="group",condition1 = "group1",condition2 = "group2"))
  expect_error(calculate_polya_stats(example_valid_polya_table_3levels,grouping_factor="group",condition1 = "group1",condition2 = "nonexistent_group"))
  expect_error(calculate_polya_stats(example_valid_polya_table_3levels,grouping_factor="group",condition2 = "group2",condition1 = "nonexistent_group"))
  expect_error(calculate_polya_stats(example_valid_polya_table_3levels,grouping_factor="group"))
  expect_error(calculate_polya_stats(example_polya_table_sample1,grouping_factor="group"))

  expect_error(calculate_polya_stats(example_valid_polya_table,grouping_factor="non-existent-group"))
  expect_error(calculate_polya_stats(example_valid_polya_table,grouping_factor="group",use_dwell_time = "TRUE"))
  expect_error(calculate_polya_stats(example_valid_polya_table,grouping_factor="group",min_reads = "2"))
  expect_error(calculate_polya_stats(example_valid_polya_table,grouping_factor="group",min_reads = NULL))
})


test_that("Valid input parameters provided for calculate_pca()",{

  summarized_example_valid_polya_table <- summarize_polya(example_valid_polya_table,summary_factors = "sample_name")

  expect_error(calculate_pca())
  expect_error(calculate_pca(empty_polya_data_table))
  expect_error(calculate_pca(example_valid_polya_table))
  expect_error(calculate_pca(example_valid_polya_table,parameter = "polya_imagined_summary"))
  expect_named(calculate_pca(summarized_example_valid_polya_table),c("pca","sample_names"),ignore.order = TRUE)
  expect_error(calculate_pca(summarized_example_valid_polya_table %>% dplyr::select(-sample_name)))
  expect_error(calculate_pca(summarized_example_valid_polya_table %>% dplyr::ungroup() %>% dplyr::select(-transcript)))
})


test_that("summarize polyA is working",{
  expect_error(summarize_polya())
  expect_error(summarize_polya(empty_polya_data_table))
  expect_error(summarize_polya(example_valid_polya_table,summary_factors = 1))
  expect_error(summarize_polya(example_valid_polya_table_3levels,summary_factors = "non-existent-factor"))
  expect_silent(summarize_polya(example_valid_polya_table_3levels))
})


test_that("Valid input parameters provided for calculate_diff_exp_binom()",{
  expect_error(calculate_diff_exp_binom())
  expect_error(calculate_diff_exp_binom(empty_polya_data_table))

  expect_error(calculate_diff_exp_binom(example_valid_polya_table))
  expect_error(calculate_diff_exp_binom(example_valid_polya_table,grouping_factor = "group"))

  expect_error(calculate_diff_exp_binom(example_valid_polya_table,grouping_factor = "group",condition1 = "group1", condition2 = "non-existent-factor-level"))
  expect_error(calculate_diff_exp_binom(example_valid_polya_table,grouping_factor = "group",condition2 = "group1", condition1 = "non-existent-factor-level"))

  expect_length(rownames(calculate_diff_exp_binom(example_valid_polya_table,grouping_factor = "group",condition1 = "group1",condition2="group2")),50)
})

