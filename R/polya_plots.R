#' PCA biplot
#'
#' Plots PCA biplot using \link[ggbiplot]{ggbiplot}
#'
#' @param pca_object pca object
#' @param samples_names names of samples to be shown on the plot
#'
#' @return \link[ggplot2]{ggplot} object

#' @export
#'
plot_polyA_PCA <- function(pca_object,samples_names) {

  ggbiplot::ggbiplot(pca_object,var.axes = F,labels = samples_names)

}


#' Plots polya distribution as ndensity plot
#'
#' @param polya_data input data with polyA predictions
#' @param groupingFactor how to group
#' @param scale_x_limit_low lower limit of x axis
#' @param scale_x_limit_high upper limit of x axis
#' @param color_palette RcolorBrewer palette
#' @param reverse_palette should color palette be reversed? (TRUE or FALSE)
#' @param plot_title (optional) plot title (parameter for ggtitle)
#' @param parameter_to_plot what to plot on x scale (defaults to polya_length)
#' @param condition1 First condition to include on the plot
#' @param condition2 Second condition to include on the plot
#'
#' @return \link[ggplot2]{ggplot} object
#' @export
#'

plot_polya_distribution <- function(polya_data, groupingFactor=NA, parameter_to_plot = "polya_length", scale_x_limit_low=NA, scale_x_limit_high=NA, color_palette = "Set1", reverse_palette = FALSE, plot_title = NA, condition1=NA,condition2=NA) {


  if (missing(polya_data)) {
    stop("PolyA predictions are missing. Please provide a valid polya_data argument",
         call. = FALSE)
  }

  assertthat::assert_that(groupingFactor %in% colnames(polya_data),msg=paste0(groupingFactor," is not a column of input dataset"))
  assertthat::assert_that(assertive::is_a_bool(reverse_palette),msg="Please provide boolen value for reverse_palette option")


  if (!is.na(condition1)) {
    if(!is.na(condition2)) {
      assertthat::assert_that(condition1 %in% levels(polya_data[[groupingFactor]]),msg=paste0(condition1," is not a level of ",grouping_factor," (groupingFactor)"))
      assertthat::assert_that(condition2 %in% levels(polya_data[[groupingFactor]]),msg=paste0(condition2," is not a level of ",grouping_factor," (groupingFactor)"))
      assertthat::assert_that(condition2 != condition1,msg="condition2 should be different than condition1")
      polya_data <- polya_data %>% dplyr::filter(!!rlang::sym(groupingFactor) %in% c(condition1,condition2))
    }
  }


  if (!is.na(groupingFactor)) {
    distribution_plot <- ggplot2::ggplot(polya_data,ggplot2::aes_string(x=parameter_to_plot,color=groupingFactor)) + ggplot2::geom_density(size=1,ggplot2::aes(y=..ndensity..)) + ggplot2::theme_bw() + ggplot2::ylab("normalized frequency")
    if (reverse_palette) {
      distribution_plot <- distribution_plot + ggplot2::scale_colour_brewer(palette = color_palette,direction=-1)
    }
    else {
      distribution_plot <- distribution_plot + ggplot2::scale_colour_brewer(palette = color_palette)
    }
  }
  else {
    distribution_plot <- ggplot2::ggplot(polya_data,ggplot2::aes_string(x=parameter_to_plot)) + ggplot2::geom_density(size=1,ggplot2::aes(y=..ndensity..)) + ggplot2::theme_bw() + ggplot2::xlab("normalized frequency")
  }
  if (!is.na(scale_x_limit_low)) {
    if (is.na(scale_x_limit_high)) {
      stop("Please provide upper limit for x scale")
    }
    assertthat::assert_that(assertive::is_numeric(scale_x_limit_low),msg="Please provide numeric value for scale_x_limit_low")
    assertthat::assert_that(assertive::is_numeric(scale_x_limit_high),msg="Please provide numeric value for scale_x_limit_high")
    distribution_plot <- distribution_plot + ggplot2::scale_x_continuous(limits=c(scale_x_limit_low,scale_x_limit_high))
  }


  if(!is.na(plot_title)){
    distribution_plot <- distribution_plot + ggplot2::ggtitle(plot_title)
  }

  return(distribution_plot)

}





#' Plots boxplot of estimated polya lengths
#'
#' @param polya_data input table with polyA predictions
#' @param groupingFactor which factor to use for grouping
#' @param scale_y_limit_low lower limit of Y axis
#' @param scale_y_limit_high upper limit of Y axis
#' @param color_palette RColorBrewer palette to use
#' @param reverse_palette Should the palette be reversed?
#' @param plot_title Plot title (optional)
#' @param condition1 First condition to include on the plot
#' @param condition2 Second condition to include on the plot
#'
#' @return \link[ggplot2]{ggplot} object
#' @export
#'
plot_polya_boxplot <- function(polya_data, groupingFactor, scale_y_limit_low=NA, scale_y_limit_high=NA, color_palette = "Set1", reverse_palette = FALSE, plot_title = NA,condition1=NA,condition2=NA) {


  if (missing(polya_data)) {
    stop("PolyA predictions are missing. Please provide a valid polya_data argument",
         call. = FALSE)
  }

  assertthat::assert_that(groupingFactor %in% colnames(polya_data),msg=paste0(groupingFactor," is not a column of input dataset"))
  assertthat::assert_that(assertive::is_a_bool(reverse_palette),msg="Please provide boolen value for reverse_palette option")

  if (!is.na(condition1)) {
    if(!is.na(condition2)) {
      assertthat::assert_that(condition1 %in% levels(polya_data[[groupingFactor]]),msg=paste0(condition1," is not a level of ",grouping_factor," (groupingFactor)"))
      assertthat::assert_that(condition2 %in% levels(polya_data[[groupingFactor]]),msg=paste0(condition2," is not a level of ",grouping_factor," (groupingFactor)"))
      assertthat::assert_that(condition2 != condition1,msg="condition2 should be different than condition1")
      polya_data <- polya_data %>% dplyr::filter(!!rlang::sym(groupingFactor) %in% c(condition1,condition2))
    }
  }

  transcripts_boxplot <- ggplot2::ggplot(polya_data,ggplot2::aes_string(x=groupingFactor,y="polya_length")) + ggplot2::geom_boxplot()
  if (!is.na(scale_y_limit_low)) {
    if (is.na(scale_y_limit_high)) {
      stop("Please provide both limits for y scale")
    }
    assertthat::assert_that(assertive::is_numeric(scale_y_limit_low),msg="Please provide numeric value for scale_y_limit_low")
    assertthat::assert_that(assertive::is_numeric(scale_y_limit_high),msg="Please provide numeric value for scale_y_limit_high")
    transcripts_boxplot <- transcripts_boxplot + ggplot2::scale_y_continuous(limits=c(scale_y_limit_low,scale_y_limit_high))
  }
  if (reverse_palette) {
    transcripts_boxplot <- transcripts_boxplot + ggplot2::scale_colour_brewer(palette = color_palette,direction=-1)
  }
  else {
    transcripts_boxplot <- transcripts_boxplot + ggplot2::scale_colour_brewer(palette = color_palette)
  }

  if(!is.na(plot_title)){
    transcripts_boxplot <- transcripts_boxplot + ggplot2::ggtitle(plot_title)
  }

  return(transcripts_boxplot)
}


#' Title
#'
#' @param polya_data_summarized polyA predictions table, summarized using \link{summarize_polya}
#' @param groupingFactor name of column used for grouping
#' @param color_palette RColorBrewer color palette to be used
#' @param reverse_palette should the palette be reversed
#' @param condition1 first condition to use for plotting
#' @param condition2 second condition to use for plotting
#' @param min_counts minimum number of counts to be shown
#' @param max_counts maximum number of counts to be shown
#'
#' @return \link[ggplot2]{ggplot} object
#' @export
#'
plot_counts_scatter <- function(polya_data_summarized, groupingFactor = NA, color_palette = "Set1", reverse_palette = FALSE,condition1 = NA, condition2 = NA,min_counts=0,max_counts=NA) {


  if (missing(polya_data_summarized)) {
    stop("Summarized PolyA predictions are missing. Please provide a valid polya_data argument",
         call. = FALSE)
  }


  assertthat::assert_that(!is.na(condition1),msg = "Please specify conditions for comparison")
  assertthat::assert_that(!is.na(condition2),msg = "Please specify conditions for comparison")
  assertthat::assert_that(!is.na(groupingFactor),msg = "Please specify groupingFactor")
  assertthat::assert_that(condition1 %in% levels(polya_data_summarized[[groupingFactor]]),msg=paste0(condition1," is not a level of ",grouping_factor," (groupingFactor)"))
  assertthat::assert_that(condition2 %in% levels(polya_data_summarized[[groupingFactor]]),msg=paste0(condition2," is not a level of ",grouping_factor," (groupingFactor)"))
  assertthat::assert_that(condition2 != condition1,msg="condition2 should be different than condition1")
  assertthat::assert_that("counts" %in% colnames(polya_data_summarized),msg = "Please provide summarized polya table (using summarize_polya()) as an input")
  assertthat::assert_that(assertive::is_a_bool(reverse_palette),msg="Please provide boolen value for reverse_palette option")

  polya_data_summarized_counts_xy<-polya_data_summarized %>% dplyr::group_by(transcript,!!rlang::sym(groupingFactor)) %>% dplyr::summarize(counts_sum=sum(counts)) %>% tidyr::spread_(groupingFactor,"counts_sum")

  polya_data_summarized_counts_xy[is.na(polya_data_summarized_counts_xy)] <- 0

  polya_data_summarized_counts_xy <- polya_data_summarized_counts_xy %>% dplyr::filter(!!rlang::sym(condition1)>=min_counts,!!rlang::sym(condition2)>=min_counts)

  if (!is.na(max_counts)) {
    assertthat::assert_that(assertive::is_numeric(max_counts),msg="Please provide numeric value for max_counts")
    polya_data_summarized_counts_xy <- polya_data_summarized_counts_xy %>% dplyr::filter(!!rlang::sym(condition1)<=max_counts,!!rlang::sym(condition2)<=max_counts)
  }

  counts_scatter_plot<-ggplot2::ggplot(polya_data_summarized_counts_xy,ggplot2::aes(x=!!rlang::sym(condition1),y=!!rlang::sym(condition2))) + ggplot2::geom_point(ggplot2::aes(text=transcript),alpha=0.7)

  if (reverse_palette) {
    counts_scatter_plot <- counts_scatter_plot + ggplot2::scale_colour_brewer(palette = color_palette,direction=-1)
  }
  else {
    counts_scatter_plot <- counts_scatter_plot + ggplot2::scale_colour_brewer(palette = color_palette)
  }

  return(counts_scatter_plot)
}



#' Plot Nanopolish polya QC
#'
#' @param nanopolish_processing_info output of \link{get_nanopolish_processing_info}
#' @param color_palette RcolorBrewer palette to use
#' @param reverse_palette Should the palette be reversed
#' @param frequency show frequency plot instead of counts plot
#'
#' @return \link[ggplot2]{ggplot} object
#' @export
#'
plot_nanopolish_qc <- function(nanopolish_processing_info,color_palette = "Set1", reverse_palette = FALSE,frequency=TRUE) {


  if (missing(nanopolish_processing_info)) {
    stop("nanopolish processing info is missing. Please provide a valid nanopolish_processing_info argument",
         call. = FALSE)
  }

  assertthat::assert_that(assertive::has_rows(nanopolish_processing_info),msg = "Empty data.frame provided as an input")
  basic_colnames = c("qc_tag","n")
  assertthat::assert_that(basic_colnames[1] %in% colnames(nanopolish_processing_info),msg="qc_tag column is missing in the input. Is that valid output of get_nanopolish_processing_info()?")
  assertthat::assert_that(basic_colnames[2] %in% colnames(nanopolish_processing_info),msg="n column is missing in the input. Is that valid output of get_nanopolish_processing_info()?")

  assertthat::assert_that(assertive::is_a_bool(frequency),msg="Non-boolean value provided for option frequency")
  assertthat::assert_that(assertive::is_a_bool(reverse_palette),msg="Please provide boolen value for reverse_palette option")

  # if there were multiple samples compared
  if (ncol(nanopolish_processing_info)>2) {
    grouping_colname = setdiff(colnames(nanopolish_processing_info),basic_colnames)
    nanopolish_qc_plot <- ggplot2::ggplot(nanopolish_processing_info,ggplot2::aes(x=!!rlang::sym(grouping_colname),fill=qc_tag,y=n))
    if(frequency) {
      nanopolish_qc_plot <- nanopolish_qc_plot + ggplot2::geom_bar(stat="identity",position="fill") + ggplot2::ylab("frequency")
    }
    else{
      nanopolish_qc_plot <- nanopolish_qc_plot + ggplot2::geom_bar(stat="identity",position="stack") + ggplot2::ylab("count")
    }
    nanopolish_qc_plot <- nanopolish_qc_plot + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,vjust = 0.7))
  }
  else {
    nanopolish_qc_plot <- ggplot2::ggplot(nanopolish_processing_info,ggplot2::aes(x=qc_tag,y=n)) + ggplot2::geom_bar(stat="identity")
  }

  if (reverse_palette) {
    nanopolish_qc_plot <- nanopolish_qc_plot + ggplot2::scale_fill_brewer(palette = color_palette,direction=-1)
  }
  else {
    nanopolish_qc_plot <- nanopolish_qc_plot + ggplot2::scale_fill_brewer(palette = color_palette)
  }

  return(nanopolish_qc_plot)
}



#' Plots volcano plot of differential expression analysis
#'
#' @param input_data a table with output from \link{calculate_diff_exp_binom} or \link{calculate_polya_stats}
#' @param color_palette RColorBrewer palette to use
#' @param reverse_palette should the palette be reserved
#'
#' @return \link[ggplot2]{ggplot} object
#' @export
#'
plot_volcano <- function(input_data,color_palette = "Set1", reverse_palette = FALSE) {


  if (missing(input_data)) {
    stop("nanopolish processing info is missing. Please provide a valid nanopolish_processing_info argument",
         call. = FALSE)
  }

  assertthat::assert_that(assertive::has_rows(input_data),msg = "Empty data.frame provided as an input")
  assertthat::assert_that("fold_change" %in% colnames(input_data),msg = "Input table is not a valid input for plot_volcano(). fold_change column is missing.")
  assertthat::assert_that("padj" %in% colnames(input_data),msg = "Input table is not a valid input for plot_volcano(). padj column is missing.")
  assertthat::assert_that("significance" %in% colnames(input_data),msg = "Input table is not a valid input for plot_volcano(). significance column is missing.")
  assertthat::assert_that("transcript" %in% colnames(input_data),msg = "Input table is not a valid input for plot_volcano(). transcript column is missing.")
  assertthat::assert_that(assertive::is_a_bool(reverse_palette),msg="Please provide boolen value for reverse_palette option")

  volcano_plot <- ggplot2::ggplot(input_data,ggplot2::aes(x=log2(fold_change),y=-log10(padj),col=significance)) + ggplot2::geom_point(ggplot2::aes(text=transcript))

  if (reverse_palette) {
    volcano_plot <- volcano_plot + ggplot2::scale_color_brewer(palette = color_palette,direction=-1)
  }
  else {
    volcano_plot <- volcano_plot + ggplot2::scale_color_brewer(palette = color_palette)
  }
  return(volcano_plot)

}


#' Plots MA plot of differential expression analysis
#'
#' Crates simple MA plot, with log10(mean expression) on the X-axis and log2(fold_change) on the Y-axis
#'
#'
#' @param input_data a table with output from \link{calculate_diff_exp_binom}
#' @param color_palette RColorBrewer palette to use
#' @param reverse_palette should the palette be reserved
#'
#' @return \link[ggplot2]{ggplot} object
#' @export
#'

plot_MA <- function(input_data,color_palette = "Set1", reverse_palette = FALSE) {


  if (missing(input_data)) {
    stop("nanopolish processing info is missing. Please provide a valid nanopolish_processing_info argument",
         call. = FALSE)
  }

  assertthat::assert_that(assertive::has_rows(input_data),msg = "Empty data.frame provided as an input")
  assertthat::assert_that("fold_change" %in% colnames(input_data),msg = "Input table is not a valid input for plot_volcano(). fold_change column is missing.")
  assertthat::assert_that("mean_expr" %in% colnames(input_data),msg = "Input table is not a valid input for plot_volcano(). padj column is missing.")
  assertthat::assert_that("significance" %in% colnames(input_data),msg = "Input table is not a valid input for plot_volcano(). significance column is missing.")
  assertthat::assert_that("transcript" %in% colnames(input_data),msg = "Input table is not a valid input for plot_volcano(). transcript column is missing.")
  assertthat::assert_that(assertive::is_a_bool(reverse_palette),msg="Please provide boolen value for reverse_palette option")

  MA_plot <- ggplot2::ggplot(input_data,ggplot2::aes(x=log10(mean_expr),y=log2(fold_change),col=significance)) + ggplot2::geom_point(ggplot2::aes(text=transcript))

  if (reverse_palette) {
    MA_plot <- MA_plot + ggplot2::scale_color_brewer(palette = color_palette,direction=-1)
  }
  else {
    MA_plot <- MA_plot + ggplot2::scale_color_brewer(palette = color_palette)
  }
  return(MA_plot)

}

