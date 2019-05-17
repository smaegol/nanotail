#' Read Single Nanopolish polyA preditions from file
#'
#' This is the basic function used to import output from \code{nanopolish polya} to R
#'
#' @param polya_path path to nanopolish output file
#' @param sample_name sample name (optional), provided as a string.
#' If specified will be included as an additional column sample_name.
#' @param gencode are contig names GENCODE-compliant.
#' Can get transcript names and ensembl_transcript IDs if reads were mapped for example to Gencode reference transcriptome
#'
#' @family read_data
#'
#' @export
#'
#' @return a [tibble][tibble::tibble-package] with polya predictions
#'
read_polya_single <- function(polya_path, gencode = TRUE, sample_name = NA) {
    # required asserts

    #check if parameters are provided
    if (missing(polya_path)) {
      stop("The path to polyA predictions (argument polya_path) is missing",
           call. = FALSE)
    }
    assert_that(is_a_non_missing_nor_empty_string(polya_path),msg = "Empty string provided as an input. Please provide a polya_path as a string")
    assert_that(is_existing_file(polya_path),msg=paste("File ",polya_path," not exists",sep=""))
    assert_that(is_non_empty_file(polya_path),msg=paste("File ",polya_path," is empty",sep=""))
    assert_that(is_a_bool(gencode),msg="Please provide TRUE/FALSE values for gencode parameter")

    message(paste0("Loading data from ",polya_path))

    #integer64 set to "numeric" to avoid inconsistences when called from read_polya_multiple
    polya_data <- data.table::fread(polya_path, integer64 = "numeric", data.table = F,header=TRUE,stringsAsFactors = FALSE,check.names = TRUE) %>% dplyr::as_tibble()
    polya_data %<>% dplyr::mutate(polya_length = round(polya_length))
    # change first column name
    colnames(polya_data)[1] <- "read_id"
    # ENSMUST00000103410.2|ENSMUSG00000076609.2|OTTMUSG00000053470.1|OTTMUST00000133385.1|RP23-435I4.10-001|Igkc|532|IG_C_gene|
    # ENSMUST00000052902.8|Gm9797|k|5db8d4b0-ba94-40b5-a3f7-54afd85f51e1_0_ENSMUSG00000045455.8 if ENSEMBL based contig names - read ensembl IDS,
    # transcript names
    if (gencode == TRUE) {
        transcript_names <- gsub(".*?\\|.*?\\|.*?\\|.*?\\|.*?\\|(.*?)\\|.*", "\\1", polya_data$contig)
        polya_data$transcript <- transcript_names
        ensembl_transcript_ids <- gsub("^(.*?)\\|.*\\|.*", "\\1", polya_data$contig)
        ensembl_transcript_ids_short <- gsub("(.*)\\..", "\\1", ensembl_transcript_ids)
        polya_data$ensembl_transcript_id_full <- ensembl_transcript_ids
        polya_data$ensembl_transcript_id_short <- ensembl_transcript_ids_short
    }

    if(!is.na(sample_name)) {
      message("sample_name is not NA")
      polya_data$sample_name = sample_name
    }

    return(polya_data)
}


# TODO Benchmark na lapply/for loop

#' Reads multiple nanopolish polyA predictions at once
#'
#' This function can be used to load any number of files with polyA predictions with single invocation,
#' allowing for metadata specification.
#'
#'
#' @param samples_table data.frame or tibble containing samples metadata and paths to files.
#' Should have at least two columns: \itemize{
#' \item polya_path - containing path to the polya predictions file
#' \item sample_name - unique name of the sample
#' \item additional columns can provide metadata which will be included in the final table
#' }
#' @param ... - additional parameters to pass to read_polya_single(), like gencode=[TRUE/FALSE]
#'
#' @return a [tibble][tibble::tibble-package] containing polyA predictions for all specified samples, with metadata provided in samples_table
#' stored as separate columns
#'
#' @family read_data
#'
#' @export
#'
read_polya_multiple <- function(samples_table,...) {

  if (missing(samples_table)) {
    stop("Samples table argument is missing",
         call. = FALSE)
  }

  assert_that(has_rows(samples_table),msg = "Empty data frame provided as an input (samples_table). Please provide samples_table describing data to load")

  samples_data <- samples_table %>% as.tbl() %>% dplyr::mutate(polya_path = as.character(polya_path)) %>% dplyr::group_by(sample_name) %>% dplyr::mutate(polya_contents=purrr::map(polya_path, function(x) read_polya_single(x))) %>% dplyr::ungroup() %>% dplyr::select(-polya_path)
  polya_data <- tidyr::unnest(samples_data)
  return(polya_data)
}


#' Removes reads which failed during Nanopolish polya processing
#'
#' Convenient function to quickly remove all reads failing during nanopolish polya processing
#'
#' @param polya_data output table from \code{read_polya_single} or \code{read_polya_multiple}
#'
#' @return a [tibble][tibble::tibble-package] with only reads having qc_tag=='PASS'
#'
#' @export
#'
remove_failed_reads <- function(polya_data) {

  if (missing(polya_data)) {
    stop("Please provide data.frame with polyA predictions as an input.",
         call. = FALSE)
  }

  assert_that(has_rows(polya_data),msg = "Empty data frame provided as an input (polya_data). Please provide valid input")

  filtered_polya_data <- polya_data %>% dplyr::filter(qc_tag=='PASS')
  return(filtered_polya_data)
}


