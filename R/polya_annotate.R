#' Title
#'
#' @param polya_data
#' @param genome
#'
#' @return
#' @export
#'
#' @examples
annotate_with_annotables <- function(polya_data,genome) {

  if ( !requireNamespace('annotables',quietly = TRUE) ) {
    stop("NanoTail requires 'annotables'. Please install it using
         install.packages('devtools')
        devtools::install_github('stephenturner/annotables')")
  }


  if (missing(polya_data)) {
    stop("Please provide data.frame with polyA predictions as an input.",
         call. = FALSE)
  }

  if (missing(genome)) {
    stop("Please provide valid genome from annotables package to use for annotation.",
         call. = FALSE)
  }

  assert_that(has_rows(polya_data),msg = "Empty data frame provided as an input (polza?data). Please provide valid input")

  tx_to_gene_table = paste0(genome,"_tx2gene")
  print(tx_to_gene_table)
  polya_data_annotated <-  polya_data %>% dplyr::inner_join(eval(as.symbol(tx_to_gene_table)),by=c("ensembl_transcript_id_short"  = "enstxp")) %>% dplyr::inner_join(eval(as.name(genome)))

  return(polya_data_annotated)
  }


genome_test = "grch38"
genome_test2<-paste0("annotables::",genome_test)
