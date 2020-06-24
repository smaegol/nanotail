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
#' Uses base subsetting and \link{sample} or dplyr \link[dplyr]{sample_n} or \link[dplyr]{sample_frac} to get the subset of the bigger data.frame or tibble
#'
#' @param input_table input table for subsampling
#' @param groupingFactor grouping factor(s)
#' @param subsample specify absolute number of rows or fraction to subsample from the data frame (group-wise)
#'
#' @return \link{tibble}
#' @export
#'
subsample_table <- function(input_table,groupingFactor=NA,subsample=NA)
{
  if (missing(input_table)) {
    stop("PolyA predictions are missing. Please provide a valid polya_data argument",
         call. = FALSE)
  }

  #assertthat::assert_that(is.numeric(reads_to_subsample),"Non-numeric argument for reads_to_subsample")

  assertthat::assert_that(!is.na(subsample),msg = "Please provide subsample option as an integer or fraction")

  if (isTRUE(round(subsample)==subsample)) {
    subsample_number = TRUE
  }
  else {
    subsample_number = FALSE
  }

  #if set to 0 - do not subsample - return input table))
  if (subsample==0) {
    return(input_table)
  }
  else {
    if(!is.na(groupingFactor)) {
      # group, if required
      assertthat::assert_that(groupingFactor %in% colnames(input_table),msg=paste0(groupingFactor," is not a column of input dataset"))
      input_table <- input_table %>% group_by(.dots = groupingFactor)
      if (subsample_number) {
        input_table <- dplyr::sample_n(input_table,subsample)
      }
      else {
        input_table <- dplyr::sample_frac(input_table,subsample)
      }
    }
    else {
      if (any(class(polya_test_lymph2)=="grouped_df")) {
        grouping_var = dplyr::group_vars(input_table)
        #input_table %>% ungroup(input_table)
      }
      if (subsample_number) {
        input_table <- input_table[sample(nrow(input_table),subsample),]
      }
      else {
        input_table <- dplyr::sample_frac(input_table,subsample)
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




# based on https://community.rstudio.com/t/spread-with-multiple-value-columns/5378/2
#' Spread multiple columns
#'
#' @param df data frame to apply spread on
#' @param key as in \link{spread}
#' @param value vector of columns to be taken as value for \link{spread}
#'
#' @return \link{tibble}
#' @export
#'
spread_multiple <- function(df, key, value) {
  # quote key
  keyq <- rlang::enquo(key)
  # break value vector into quotes
  valueq <- rlang::enquo(value)
  s <- rlang::quos(!!valueq)
  df %>% tidyr::gather(variable, value, !!!s) %>%
    tidyr::unite(temp, !!keyq, variable) %>%
    tidyr::spread(temp, value)
}


#' Calculates scaling vector for virtual gel plotting
#'
#' @param input_data input polyA table for calculation of scaling factor (count of reads)
#' @param groupingFactor for which factor calculate counts
#'
#' @return named vector
#' @export
#'
calculate_scaling_vector_for_virutal_gel <- function(input_data,groupingFactor) {
  scaling_vector <- summarize_polya(input_data,transcript_id_column = groupingFactor) %>% dplyr::select(!!rlang::sym(groupingFactor),counts) %>% tibble::deframe()
  return(scaling_vector)
  }



StatMedianLine <- ggplot2::ggproto("StatMedianLine", ggplot2::Stat,
                          compute_group = function(data, scales) {
                            transform(data, yintercept=median(y))
                          },
                          required_aes = c("x", "y")
)

#' Helper function for calculating median stat for violin/boxpolot ggplot plots
#'
#' @param mapping 
#' @param data 
#' @param geom 
#' @param position 
#' @param na.rm 
#' @param show.legend 
#' @param inherit.aes 
#' @param ... 
#'
#' @return
#' @export
#'
stat_median_line <- function(mapping = NULL, data = NULL, geom = "hline",
                             position = "identity", na.rm = FALSE, show.legend = NA, 
                             inherit.aes = TRUE, ...) {
  ggplot2::layer(
    stat = StatMedianLine, data = data, mapping = mapping, geom = geom, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}