context("Test plotting functions")

library(assertthat)
library(testthat)
library(assertive)


test_that("valid parameters are provided for plot_polya_distribution()",{

  expect_error(plot_polya_distribution())
  expect_error(plot_polya_distribution(example_valid_polya_table))
  expect_error(plot_polya_distribution(example_valid_polya_table,groupingFactor = "non-existent_group"))
  expect_silent(plot_polya_distribution(example_valid_polya_table,groupingFactor = "group"))
  expect_silent(plot_polya_distribution(example_valid_polya_table,groupingFactor = "group",condition1="group1",condition2="group2"))
  expect_error(plot_polya_distribution(example_valid_polya_table,groupingFactor = "group",condition1="group1",condition2="non-existent_group"))
  expect_error(plot_polya_distribution(example_valid_polya_table,groupingFactor = "group",condition2="group1",condition1="non-existent_group"))
  expect_error(plot_polya_distribution(example_valid_polya_table,groupingFactor = "group",condition1="group1",condition2="group1"))
  expect_error(plot_polya_distribution(example_valid_polya_table,groupingFactor = "group",scale_x_limit_low = 0))
  expect_silent(plot_polya_distribution(example_valid_polya_table,groupingFactor = "group",scale_x_limit_low = 0,scale_x_limit_high = 200))
  expect_error(plot_polya_distribution(example_valid_polya_table,groupingFactor = "group",scale_x_limit_low = "0",scale_x_limit_high = 200))
  expect_error(plot_polya_distribution(example_valid_polya_table,groupingFactor = "group",scale_x_limit_low = 0,scale_x_limit_high = "200"))
  expect_warning(plot_polya_distribution(example_valid_polya_table,groupingFactor = "group",scale_x_limit_low = 0,scale_x_limit_high = 200,color_palette = "non_existent_palette"),"Please provide valid color palette",all = FALSE,fixed=TRUE)
})

test_that("valid parameters are provided for plot_polya_boxplot()",{

  expect_error(plot_polya_boxplot())
  expect_error(plot_polya_boxplot(example_valid_polya_table))
  expect_error(plot_polya_boxplot(example_valid_polya_table,groupingFactor = "non-existent_group"))
  expect_silent(plot_polya_boxplot(example_valid_polya_table,groupingFactor = "group"))
  expect_silent(plot_polya_boxplot(example_valid_polya_table,groupingFactor = "group",condition1="group1",condition2="group2"))
  expect_error(plot_polya_boxplot(example_valid_polya_table,groupingFactor = "group",condition1="group1",condition2="non-existent_group"))
  expect_error(plot_polya_boxplot(example_valid_polya_table,groupingFactor = "group",condition2="group1",condition1="non-existent_group"))
  expect_error(plot_polya_boxplot(example_valid_polya_table,groupingFactor = "group",condition1="group1",condition2="group1"))
  expect_error(plot_polya_boxplot(example_valid_polya_table,groupingFactor = "group",scale_y_limit_low = 0))
  expect_silent(plot_polya_boxplot(example_valid_polya_table,groupingFactor = "group",scale_y_limit_low = 0,scale_y_limit_high = 200))
  expect_error(plot_polya_boxplot(example_valid_polya_table,groupingFactor = "group",scale_y_limit_low = "0",scale_y_limit_high = 200))
  expect_error(plot_polya_boxplot(example_valid_polya_table,groupingFactor = "group",scale_y_limit_low = 0,scale_y_limit_high = "200"))
  expect_warning(plot_polya_boxplot(example_valid_polya_table,groupingFactor = "group",scale_y_limit_low = 0,scale_y_limit_high = 200,color_palette = "non_existent_palette"),"Please provide valid color palette",all = FALSE,fixed=TRUE)
})


test_that("valid parameters are provided for plot_counts_scatter()",{

  summarized_example_valid_polya_table <- summarize_polya(example_valid_polya_table,summary_factors = "group")

  expect_error(plot_counts_scatter())
  expect_error(plot_counts_scatter(example_valid_polya_table))
  expect_error(plot_counts_scatter(summarized_example_valid_polya_table,groupingFactor = "group"))
  expect_warning(plot_counts_scatter(summarized_example_valid_polya_table,groupingFactor = "group",condition1="group1",condition2="group2"),"Ignoring unknown aesthetics: text",fixed=TRUE)
  expect_error(plot_counts_scatter(summarized_example_valid_polya_table,groupingFactor = "group",condition1="group1",condition2="non-existent_group"))
  expect_error(plot_counts_scatter(summarized_example_valid_polya_table,groupingFactor = "group",condition2="group1",condition1="non-existent_group"))
  expect_error(plot_counts_scatter(summarized_example_valid_polya_table,groupingFactor = "group",condition1="group1",condition2="group1"))
  expect_error(plot_counts_scatter(summarized_example_valid_polya_table,groupingFactor = "group",max_counts = 0))
  expect_warning(plot_counts_scatter(summarized_example_valid_polya_table,groupingFactor = "group",condition1="group1",condition2="group2",max_counts=200),"Ignoring unknown aesthetics: text",fixed=TRUE)
  expect_error(plot_counts_scatter(summarized_example_valid_polya_table,groupingFactor = "group",condition1="group1",condition2="group2",max_counts="200"))
  expect_warning(plot_counts_scatter(summarized_example_valid_polya_table,groupingFactor = "group",condition1="group1",condition2="group2",max_counts=200,color_palette = "non_existent_palette"),"Please provide valid color palette",all = FALSE,fixed=TRUE)

})


test_that("valid parameters are provided for plot_nanopolish_qc()",{

  nanopolish_qc_of_example_valid_polya_table <- get_nanopolish_processing_info(example_valid_polya_table)
  grouped_nanopolish_qc_of_example_valid_polya_table <- get_nanopolish_processing_info(example_valid_polya_table,grouping_factor = "sample_name")

  expect_error(plot_nanopolish_qc())
  expect_error(plot_nanopolish_qc(example_valid_polya_table))
  expect_error(plot_nanopolish_qc(empty_polya_data_table))
  expect_silent(plot_nanopolish_qc(nanopolish_qc_of_example_valid_polya_table))
    expect_silent(plot_nanopolish_qc(nanopolish_qc_of_example_valid_polya_table,frequency = FALSE))
  expect_error(plot_nanopolish_qc(nanopolish_qc_of_example_valid_polya_table,frequency = "FALSE"))
  expect_warning(plot_nanopolish_qc(grouped_nanopolish_qc_of_example_valid_polya_table,color_palette = "non_existent_palette"),"Please provide valid color palette",all = FALSE,fixed=TRUE)
})



test_that("valid parameters are provided for plot_volcano()",{

  binom_test_output_of_example_valid_polya_table <- calculate_diff_exp_binom(example_valid_polya_table,grouping_factor = "group",condition1 = "group1",condition2 = "group2")
  expect_error(plot_volcano())
  expect_error(plot_volcano(example_valid_polya_table))
  expect_error(plot_volcano(empty_polya_data_table))
  expect_warning(plot_volcano(binom_test_output_of_example_valid_polya_table),"Ignoring unknown aesthetics: text",fixed=TRUE)
  expect_error(plot_volcano(binom_test_output_of_example_valid_polya_table %>% dplyr::select(-transcript)))
  expect_error(plot_volcano(binom_test_output_of_example_valid_polya_table %>% dplyr::select(-fold_change)))
  expect_error(plot_volcano(binom_test_output_of_example_valid_polya_table %>% dplyr::select(-padj)))
  expect_error(plot_volcano(binom_test_output_of_example_valid_polya_table %>% dplyr::select(-significance)))
  expect_warning(plot_volcano(binom_test_output_of_example_valid_polya_table,color_palette = "non_existent_palette"),"Please provide valid color palette",all = FALSE,fixed=TRUE)
})

test_that("valid parameters are provided for plot_MA()",{

  binom_test_output_of_example_valid_polya_table <- calculate_diff_exp_binom(example_valid_polya_table,grouping_factor = "group",condition1 = "group1",condition2 = "group2")
  expect_error(plot_MA())
  expect_error(plot_MA(example_valid_polya_table))
  expect_error(plot_MA(empty_polya_data_table))
  expect_warning(plot_MA(binom_test_output_of_example_valid_polya_table),"Ignoring unknown aesthetics: text",fixed=TRUE)
  expect_error(plot_MA(binom_test_output_of_example_valid_polya_table %>% dplyr::select(-transcript)))
  expect_error(plot_MA(binom_test_output_of_example_valid_polya_table %>% dplyr::select(-fold_change)))
  expect_error(plot_MA(binom_test_output_of_example_valid_polya_table %>% dplyr::select(-mean_expr)))
  expect_error(plot_MA(binom_test_output_of_example_valid_polya_table %>% dplyr::select(-significance)))
  expect_warning(plot_MA(binom_test_output_of_example_valid_polya_table,color_palette = "non_existent_palette"),"Please provide valid color palette",all = FALSE,fixed=TRUE)
})
