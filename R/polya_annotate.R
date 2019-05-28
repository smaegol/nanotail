#' Annotate polyA predictions using annotables
#'
#' @param polya_data polya data table to annotate
#' @param genome valid genome from annotables to use for annotation
#'
#' @return a \link[tibble]{tibble}
#' @export
#'
annotate_with_annotables <- function(polya_data,genome) {

  if ( !requireNamespace('annotables',quietly = TRUE) ) {
    stop("NanoTail requires 'annotables'. Please install it using
         install.packages('devtools')
        devtools::install_github('stephenturner/annotables')")
  }
  require(annotables)


  if (missing(polya_data)) {
    stop("Please provide data.frame with polyA predictions as an input.",
         call. = FALSE)
  }

  if (missing(genome)) {
    stop("Please provide valid genome from annotables package to use for annotation.",
         call. = FALSE)
  }

  assertthat::assert_that(assertive::has_rows(polya_data),msg = "Empty data frame provided as an input (polya_data). Please provide valid input")

  tx_to_gene_table = paste0(genome,"_tx2gene")
  polya_data_annotated <-  polya_data %>% dplyr::left_join( eval(as.symbol(tx_to_gene_table)),by=c("ensembl_transcript_id_short"  = "enstxp")) %>% dplyr::inner_join(eval(as.name(genome)))

  return(polya_data_annotated)
  }


#' Title
#'
#' @param polya_data polya data table to annotate
#' @param attributes_to_get what annotations should be retrieved. Default = c('external_gene_name','description','transcript_biotype')
#' @param transcript_id which column should be matched in the target mart
#' @param mart_to_use mart object created with \link[biomaRt]{useMart} or \link[biomaRt]{useEnsembl}
#'
#' @return a \link[tibble]{tibble}
#' @export
annotate_with_biomart <- function(polya_data,attributes_to_get=c('ensembl_transcript_id','external_gene_name','description','transcript_biotype'),filters='ensembl_transcript_id',mart_to_use=NA) {

  if (missing(polya_data)) {
    stop("Please provide data.frame with polyA predictions as an input.",
         call. = FALSE)
  }

  if (missing(mart_to_use)) {
    stop("Please provide valid mart object",
         call. = FALSE)
  }

  assertthat::assert_that(assertive::has_rows(polya_data),msg = "Empty data frame provided as an input (polya_data). Please provide valid input")
  assertthat::assert_that(class(mart_to_use)=="Mart",msg="Please provide valid mart object")
  assertthat::assert_that(length(attributes)>0,msg="please provide attributes")

  ensembl_ids = unique(polya_data$ensembl_transcript_id_short)
  ensembl_ids <- ensembl_ids[!is.na(ensembl_ids)]

  polya_data <- polya_data %>% dplyr::rename(ensembl_transcript_id = ensembl_transcript_id_short)

  # if using biomaRt version older than from Bioconductor 3.9, it cannot process more than 500 values at once
  if (packageVersion("biomaRt")<"2.40.0") {
    number_of_items <- length(ensembl_ids)
    annotation_data=data.frame()
    for (z in seq(1,number_of_items,by = 500)) {

      annotation_data_temp<-getBM(attributes=attributes_to_get, filters =filters, values = ensembl_ids[z:(z+499)], mart = mart_to_use)
      #print(annotation_data_temp)
      annotation_data<-rbind(annotation_data,annotation_data_temp)
    }
  }
  # since biomaRt 2.40 batch submission is possible
  else {
    annotation_data<-biomaRt::getBM(attributes=attributes_to_get, filters = filters, values = ensembl_ids, mart = mart_to_use)
  }
  polya_data_annotated <-  polya_data %>% dplyr::left_join(annotation_data)

  return(polya_data_annotated)
}
