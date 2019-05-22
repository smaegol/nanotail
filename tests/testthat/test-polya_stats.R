context("test functions for statistics calculation")


library(assertthat)
library(testthat)
library(assertive)


empty_polya_data_table = data.frame()

one_line_polya_data_table = data.frame(sample_name="wt1",group="wt",read_id="1234567890",polya_length=100,position=1,contig="contig",leader_start=1,adapter_start=1,polya_start=1,qc_tag="PASS",transcript="transcript")
two_line_polya_data_table = data.frame(sample_name=c("wt1","mut1"),group=c("wt","mut"),read_id=c("1234567890","1234567891"),polya_length=c(100,20),position=rep(1,2),contig=rep("contig1",2),leader_start=rep(1,2),adapter_start=rep(100,2),polya_start=rep(1,2),qc_tag=rep("PASS",2),transcript=rep("transcript",2))
correct_polya_data_table = data.frame(sample_name=c(rep("wt1",100),rep("wt2",100),rep("mut1",100),rep("mut2",100)),group=c(rep("wt",200),rep("mut",200)),read_id=paste0(seq(1,400),""),transcript=rep(paste0("transcript",seq(1,100)),4),polya_length=sample(300,400,replace=TRUE),qc_tag=rep("PASS",400),dwell_time=sample(10000,400,replace=TRUE))

temp_output_file <- tempfile()

test_that("Valid input parameters provided for calculate_polya_stats()",{
  expect_error(calculate_polya_stats())
  expect_error(calculate_polya_stats(empty_polya_data_table))
  expect_error(calculate_polya_stats(empty_polya_data_table))
  expect_error(calculate_polya_stats(one_line_polya_data_table,min_reads = "1"))
  expect_error(calculate_polya_stats(two_line_polya_data_table,min_reads = 0),NA)
  expect_known_output(calculate_polya_stats(correct_polya_data_table,grouping_factor="group"),temp_output_file,print=TRUE)
  expect_known_output(calculate_polya_stats(correct_polya_data_table,grouping_factor="group",stat_test = "Wilcoxon"),temp_output_file,print=TRUE)
  expect_known_output(calculate_polya_stats(correct_polya_data_table,grouping_factor="group",stat_test = "Wilcoxon",alpha = 0.05),temp_output_file,print=TRUE)
  expect_known_output(calculate_polya_stats(correct_polya_data_table,grouping_factor="group",stat_test = "Wilcoxon",alpha = 0.05),temp_output_file,print=TRUE)
  expect_known_output(calculate_polya_stats(correct_polya_data_table,grouping_factor="group",stat_test = "Wilcoxon",alpha = 0.05,condition1 = "wt",condition2="mut"),temp_output_file,print=TRUE)

  expect_named(calculate_polya_stats(correct_polya_data_table,grouping_factor="group"),c("summary","summary_short","stats"),ignore.order = TRUE)
  expect_named(calculate_polya_stats(correct_polya_data_table,grouping_factor="group",stat_test="Wilcoxon"),c("summary","summary_short","stats"),ignore.order = TRUE)
  expect_named(calculate_polya_stats(correct_polya_data_table,grouping_factor="group",stat_test="KS"),c("summary","summary_short","stats"),ignore.order = TRUE)
  expect_named(calculate_polya_stats(correct_polya_data_table,grouping_factor="group",stat_test="glm"),c("summary","summary_short","stats"),ignore.order = TRUE)
  expect_error(calculate_polya_stats(correct_polya_data_table,grouping_factor="group",stat_test = "ttest"),"Please provide one of available statistical tests")

})



test_that("summarize polyA",{
  expect_error(summarize_polya())
  expect_error(summarize_polya(empty_polya_data_table))
  expect_error(summarize_polya(one_line_polya_data_table),NA)
})
