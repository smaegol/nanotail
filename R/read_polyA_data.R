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
#' @seealso \link{read_polya_multiple}
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
    assertthat::assert_that(assertive::is_a_non_missing_nor_empty_string(polya_path),msg = "Empty string provided as an input. Please provide a polya_path as a string")
    assertthat::assert_that(assertive::is_existing_file(polya_path),msg=paste("File ",polya_path," not exists",sep=""))
    assertthat::assert_that(assertive::is_non_empty_file(polya_path),msg=paste("File ",polya_path," is empty",sep=""))
    assertthat::assert_that(assertive::is_a_bool(gencode),msg="Please provide TRUE/FALSE values for gencode parameter")

    message(paste0("Loading data from ",polya_path))

    #integer64 set to "numeric" to avoid inconsistences when called from read_polya_multiple
    polya_data <- data.table::fread(polya_path, integer64 = "numeric", data.table = F,header=TRUE,stringsAsFactors = FALSE,check.names = TRUE,showProgress = FALSE) %>% dplyr::as_tibble()
    polya_data <- polya_data %>% dplyr::mutate(polya_length = round(polya_length),dwell_time=transcript_start-polya_start)
    # change first column name
    colnames(polya_data)[1] <- "read_id"
    # transcript names, if mapping to gencode transcriptome
    if (gencode == TRUE) {
        transcript_names <- gsub(".*?\\|.*?\\|.*?\\|.*?\\|.*?\\|(.*?)\\|.*", "\\1", polya_data$contig)
        polya_data$transcript <- transcript_names
        ensembl_transcript_ids <- gsub("^(.*?)\\|.*\\|.*", "\\1", polya_data$contig)
        ensembl_transcript_ids_short <- gsub("(.*)\\..", "\\1", ensembl_transcript_ids) # without version number
        polya_data$ensembl_transcript_id_full <- ensembl_transcript_ids
        polya_data$ensembl_transcript_id_short <- ensembl_transcript_ids_short
    }
    else {
      polya_data <- polya_data %>% dplyr::rename(transcript = contig)
    }

    if(!is.na(sample_name)) {
      # set sample_name (if was set)
      if (! "sample_name" %in% colnames(polya_data)) {
        warning("sample_name was provided in the input file. Overwriting with the provided one")
      }
      polya_data$sample_name = sample_name
      polya_data$sample_name <- as.factor(polya_data$sample_name)
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
#' }
#' Additional columns can provide metadata which will be included in the final table
#' @param ... - additional parameters to pass to read_polya_single(), like gencode=(TRUE/FALSE)
#'
#' @return a [tibble][tibble::tibble-package] containing polyA predictions for all specified samples, with metadata provided in samples_table
#' stored as separate columns
#'
#' @seealso \link{read_polya_single}
#'
#' @export
#'
read_polya_multiple <- function(samples_table,...) {

  if (missing(samples_table)) {
    stop("Samples table argument is missing",
         call. = FALSE)
  }

  assertthat::assert_that(assertive::has_rows(samples_table),msg = "Empty data frame provided as an input (samples_table). Please provide samples_table describing data to load")
  assertthat::assert_that("polya_path" %in% colnames(samples_table),msg = "Samples table should contain at least polya_path and sample_name columns")
  assertthat::assert_that("sample_name" %in% colnames(samples_table),msg = "Samples table should contain at least polya_path and sample_name columns")

  samples_data <- samples_table %>% dplyr::as.tbl() %>% dplyr::mutate_if(is.character,as.factor) %>% dplyr::mutate(polya_path = as.character(polya_path)) %>% dplyr::group_by(sample_name) %>% dplyr::mutate(polya_contents=purrr::map(polya_path, function(x) read_polya_single(x))) %>% dplyr::ungroup() %>% dplyr::select(-polya_path)
  polya_data <- tidyr::unnest(samples_data)

  return(polya_data)
}


#' Removes reads which failed during Nanopolish polya processing
#'
#' Convenient function to quickly remove all reads failing during nanopolish polya processing
#'
#' @param polya_data output table from \link{read_polya_single} or \link{read_polya_multiple}
#'
#' @return a [tibble][tibble::tibble-package] with only reads having qc_tag=='PASS'
#'
#' @export
#'
#' @seealso \link{read_polya_single}, \link{read_polya_multiple}
#'
remove_failed_reads <- function(polya_data) {

  if (missing(polya_data)) {
    stop("Please provide data.frame with polyA predictions as an input.",
         call. = FALSE)
  }

  assertthat::assert_that(assertive::has_rows(polya_data),msg = "Empty data frame provided as an input (polya_data). Please provide valid input")

  filtered_polya_data <- polya_data %>% dplyr::filter(qc_tag=='PASS')
  return(filtered_polya_data)
}


