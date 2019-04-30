context("test functions for statistics calculation")


library(assertthat)
library(testthat)
library(assertive)


empty_polya_data_table = data.frame()

one_line_polya_data_table = data.frame(sample_name="wt1",group="wt",read_id="1234567890",polya_length=100,position=1,contig="contig",leader_start=1,adapter_start=1,polya_start=1,qc_tag="PASS",transcript="transcript")
two_line_polya_data_table = data.frame(sample_name=c("wt1","mut1"),group=c("wt","mut"),read_id=c("1234567890","1234567891"),polya_length=c(100,20),position=rep(1,2),contig=rep("contig1",2),leader_start=rep(1,2),adapter_start=rep(100,2),polya_start=rep(1,2),qc_tag=rep("PASS",2),transcript=rep("transcript",2))



test_that("Valid input provided",{
  expect_error(calculate_polya_stats())
  expect_error(calculate_polya_stats(empty_polya_data_table))
  expect_error(calculate_polya_stats(one_line_polya_data_table,min_reads = "1"))
  expect_error(calculate_polya_stats(two_line_polya_data_table,min_reads = 0),NA)
})

test_that("summarize polyA",{
  expect_error(summarize_polya())
  expect_error(summarize_polya(empty_polya_data_table))
  expect_error(summarize_polya(one_line_polya_data_table),NA)
})
