library(testthat)
library(nanotail)


#define testing environment for all tests:
empty_polya_data_table = data.frame()
qc_values <- c(rep("PASS",10),"ADAPTER","READ_FAILED_LOAD","SUFFCLIP","NOREGION")
one_line_polya_data_table = data.frame(sample_name="wt1",group="wt",read_id="1234567890",polya_length=100,position=1,contig="contig",leader_start=1,adapter_start=1,polya_start=1,qc_tag="PASS",transcript="transcript")
two_line_polya_data_table = data.frame(sample_name=c("wt1","mut1"),group=c("wt","mut"),read_id=c("1234567890","1234567891"),polya_length=c(100,20),position=rep(1,2),contig=rep("contig1",2),leader_start=rep(1,2),adapter_start=rep(100,2),polya_start=rep(1,2),qc_tag=rep("PASS",2),transcript=rep("transcript",2))

# generation of valid multi-sample polyA table (resembling output of read_polya_multiple())
number_of_reads_per_sample=20000
number_of_transcripts_per_sample=50
example_polya_table_sample1 = data.frame(sample_name="sample1",group="group1",read_id=paste0("sample1_",seq(1,number_of_reads_per_sample)),transcript=rep(paste0("transcript",seq(1,number_of_transcripts_per_sample)),number_of_reads_per_sample/number_of_transcripts_per_sample),polya_length=c(rlnorm(number_of_reads_per_sample,meanlog=4,sdlog=1)),qc_tag=qc_values[sample(14,20000,replace = TRUE)],position=sample(1000,number_of_reads_per_sample,replace=TRUE),leader_start=sample(seq(1,100),number_of_reads_per_sample,replace=TRUE),adapter_start=sample(seq(100,200),number_of_reads_per_sample,replace=TRUE),polya_start=sample(seq(200,500),number_of_reads_per_sample,replace=TRUE),transcript_start=sample(seq(1000,10000),number_of_reads_per_sample,replace=TRUE),read_rate=rnorm(n = number_of_reads_per_sample,mean = 70,sd=20))
example_polya_table_sample2 = data.frame(sample_name="sample2",group="group2",read_id=paste0("sample2_",seq(1,number_of_reads_per_sample)),transcript=rep(paste0("transcript",seq(1,number_of_transcripts_per_sample)),number_of_reads_per_sample/number_of_transcripts_per_sample),polya_length=c(rlnorm(number_of_reads_per_sample,meanlog=3,sdlog=1)),qc_tag=qc_values[sample(14,20000,replace = TRUE)],position=sample(1000,number_of_reads_per_sample,replace=TRUE),leader_start=sample(seq(1,100),number_of_reads_per_sample,replace=TRUE),adapter_start=sample(seq(100,200),number_of_reads_per_sample,replace=TRUE),polya_start=sample(seq(200,500),number_of_reads_per_sample,replace=TRUE),transcript_start=sample(seq(1000,10000),number_of_reads_per_sample,replace=TRUE),read_rate=rnorm(n = number_of_reads_per_sample,mean = 70,sd=20))
example_polya_table_sample3 = data.frame(sample_name="sample3",group="group3",read_id=paste0("sample2_",seq(1,number_of_reads_per_sample)),transcript=rep(paste0("transcript",seq(1,number_of_transcripts_per_sample)),number_of_reads_per_sample/number_of_transcripts_per_sample),polya_length=c(rlnorm(number_of_reads_per_sample,meanlog=3,sdlog=1)),qc_tag=qc_values[sample(14,20000,replace = TRUE)],position=sample(1000,number_of_reads_per_sample,replace=TRUE),leader_start=sample(seq(1,100),number_of_reads_per_sample,replace=TRUE),adapter_start=sample(seq(100,200),number_of_reads_per_sample,replace=TRUE),polya_start=sample(seq(200,500),number_of_reads_per_sample,replace=TRUE),transcript_start=sample(seq(1000,10000),number_of_reads_per_sample,replace=TRUE),read_rate=rnorm(n = number_of_reads_per_sample,mean = 70,sd=20))
example_valid_polya_table = rbind(example_polya_table_sample1,example_polya_table_sample2)
example_valid_polya_table_3levels = rbind(example_polya_table_sample1,example_polya_table_sample2,example_polya_table_sample3)

example_polya_table_mouse = data.frame(sample_name="sample1",group="group1",read_id=paste0("sample1_",seq(1,number_of_reads_per_sample)),transcript=rep(paste0("ENSMUST000001049",seq(1,number_of_transcripts_per_sample)),number_of_reads_per_sample/number_of_transcripts_per_sample),ensembl_transcript_id_short=rep(paste0("ENSMUST000001049",seq(1,number_of_transcripts_per_sample)),number_of_reads_per_sample/number_of_transcripts_per_sample),polya_length=c(rlnorm(number_of_reads_per_sample,meanlog=4,sdlog=1)),qc_tag=qc_values[sample(14,20000,replace = TRUE)],position=sample(1000,number_of_reads_per_sample,replace=TRUE),leader_start=sample(seq(1,100),number_of_reads_per_sample,replace=TRUE),adapter_start=sample(seq(100,200),number_of_reads_per_sample,replace=TRUE),polya_start=sample(seq(200,500),number_of_reads_per_sample,replace=TRUE),transcript_start=sample(seq(1000,10000),number_of_reads_per_sample,replace=TRUE),read_rate=rnorm(n = number_of_reads_per_sample,mean = 70,sd=20))

# generation of valid polyA tables to be saved to temp files (for testing or data import functions)
example_polya_table_sample1 = data.frame(read_id=paste0("sample1_",seq(1,number_of_reads_per_sample)),contig=rep(paste0("transcript",seq(1,number_of_transcripts_per_sample)),number_of_reads_per_sample/number_of_transcripts_per_sample),polya_length=c(rlnorm(number_of_reads_per_sample,meanlog=4,sdlog=1)),qc_tag=qc_values[sample(14,20000,replace = TRUE)],position=sample(1000,number_of_reads_per_sample,replace=TRUE),leader_start=sample(seq(1,100),number_of_reads_per_sample,replace=TRUE),adapter_start=sample(seq(100,200),number_of_reads_per_sample,replace=TRUE),polya_start=sample(seq(200,500),number_of_reads_per_sample,replace=TRUE),transcript_start=sample(seq(1000,10000),number_of_reads_per_sample,replace=TRUE),read_rate=rnorm(n = number_of_reads_per_sample,mean = 70,sd=20))
example_polya_table_sample2 = data.frame(read_id=paste0("sample2_",seq(1,number_of_reads_per_sample)),contig=rep(paste0("transcript",seq(1,number_of_transcripts_per_sample)),number_of_reads_per_sample/number_of_transcripts_per_sample),polya_length=c(rlnorm(number_of_reads_per_sample,meanlog=3,sdlog=1)),qc_tag=qc_values[sample(14,20000,replace = TRUE)],position=sample(1000,number_of_reads_per_sample,replace=TRUE),leader_start=sample(seq(1,100),number_of_reads_per_sample,replace=TRUE),adapter_start=sample(seq(100,200),number_of_reads_per_sample,replace=TRUE),polya_start=sample(seq(200,500),number_of_reads_per_sample,replace=TRUE),transcript_start=sample(seq(1000,10000),number_of_reads_per_sample,replace=TRUE),read_rate=rnorm(n = number_of_reads_per_sample,mean = 70,sd=20))

sample_tempfile1 <- tempfile()
sample_tempfile2 <- tempfile()

write.table(example_polya_table_sample1,sample_tempfile1,quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)
write.table(example_polya_table_sample2,sample_tempfile2,quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)

example_sample_table <- data.frame(polya_path=c(sample_tempfile1,sample_tempfile2),sample_name=c("sample1","sample2"),group=c("group1","group2"))


test_check("nanotail")
