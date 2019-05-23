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

