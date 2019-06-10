#' Calculates geometric mean
#'
#' @param x input vector
#' @param na.rm should NA values be removed?
#'
#' @return geometric mean of values provided as an input
#' @export
#'
#' @examples
#' a <- rnorm(100,33,5)
#' gm_mean(a)
gm_mean = function(x, na.rm=TRUE){

  assertthat::assert_that(is.vector(x),msg = "Please provide numeric vector as an input for gm_mean")
  assertthat::assert_that(length(x)>0,msg = "Empty vector provided as input")
  assertthat::assert_that(assertive::is_a_bool(na.rm),msg = "Please provide boolean value for na.rm option")
  assertthat::assert_that(assertive::is_numeric(x),msg = "Please provide numeric vector as input")
  if (length(x)==1) {
    gm_mean=x[1]
  }
  else {
    gm_mean = exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  return(gm_mean)
}



#' Subsample a date frame
#'
#' Uses base subsetting and \link{sample} or dplyr \link[dplyr]{sample_n} or \link[dplyr{sample_frac}] to get the subset of the bigger data.frame or tibble
#'
#' @param input_table input table for subsampling
#' @param groupingFactor grouping factor(s)
#' @param reads_to_subsample specify absolute number of rows to subsample from the data frame (group-wise)
#' @param fraction_to_subsample specify fraction of rows to subsample from the data frame (group-wise)
#'
#' @return \link{tibble}
#' @export
#'
subsample_table <- function(input_table,groupingFactor=NA,reads_to_subsample=NA,fraction_to_subsample=NA)
{
  if (missing(input_table)) {
    stop("PolyA predictions are missing. Please provide a valid polya_data argument",
         call. = FALSE)
  }

  #assertthat::assert_that(is.numeric(reads_to_subsample),"Non-numeric argument for reads_to_subsample")

  assertthat::assert_that((!is.na(reads_to_subsample)) || (!is.na(fraction_to_subsample)),msg = "Please provide either reads_to_subsample or fraction_to_subsample")
  assertthat::assert_that((!is.na(reads_to_subsample) && is.na(fraction_to_subsample)) || (is.na(reads_to_subsample) && (!is.na(fraction_to_subsample))),msg = "Please provide either reads_to_subsample or fraction_to_subsample (not both of them)")

  #if set to 0 - do not subsample - return input table))
  if ((!is.na(reads_to_subsample) && reads_to_subsample==0) || (!is.na(fraction_to_subsample) && fraction_to_subsample==0.0)) {
    return(input_table)
  }
  else {
    if(!is.na(groupingFactor)) {
      # group, if required
      assertthat::assert_that(groupingFactor %in% colnames(input_table),msg=paste0(groupingFactor," is not a column of input dataset"))
      input_table <- input_table %>% group_by(.dots = groupingFactor)
      if (!is.na(reads_to_subsample)) {
        input_table <- dplyr::sample_n(input_table,reads_to_subsample)
      }
      else if (!is.na(fraction_to_subsample)) {
        input_table <- dplyr::sample_frac(input_table,fraction_to_subsample)
      }
    }
    else {
      if (any(class(polya_test_lymph2)=="grouped_df")) {
        grouping_var = dplyr::group_vars(input_table)
        #input_table %>% ungroup(input_table)
      }
      if (!is.na(reads_to_subsample)) {
        input_table <- input_table[sample(nrow(input_table),reads_to_subsample),]
      }
      else if (!is.na(fraction_to_subsample)) {
        input_table <- dplyr::sample_frac(input_table,fraction_to_subsample)
      }
    }

    return(input_table)
  }
}




#' Default theme for ggplot2-based plots in the NanoTail package
#'
axis_elements_size=15
axis_titles_size=18
nanotail_ggplot2_theme <- ggplot2::theme(axis.title = ggplot2::element_text(size=axis_titles_size),axis.text = ggplot2::element_text(size=axis_elements_size),legend.text = ggplot2::element_text(size=axis_elements_size),legend.title = ggplot2::element_text(size=axis_titles_size))

