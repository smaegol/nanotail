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
#' @param parameter_to_plot what to plot on x scale (defaults to polya_length)
#' @param condition1 First condition to include on the plot
#' @param condition2 Second condition to include on the plot
#' @param show_center_values Show center values as vertical line. Possible values: "none","median","mean"
#' @param subsample Subsample input table, provide either absolute number or fraction
#' @param ndensity Should ndensity (scaled density) be plotted instead of normal denisty (default = TRUE)
#' @param ... parameters passed to .basic_aesthetics function (scale_x_limit_low = NA, scale_x_limit_high = NA, scale_y_limit_low = NA, scale_y_limit_high = NA, color_palette = "Set1",plot_title=NA)
#'
#' @return \link[ggplot2]{ggplot} object
#' @export
#'
plot_polya_distribution <- function(polya_data, groupingFactor=NA, parameter_to_plot = "polya_length", condition1=NA,condition2=NA,show_center_values="none",subsample=NA,ndensity=TRUE,...) {


  if (missing(polya_data)) {
    stop("PolyA predictions are missing. Please provide a valid polya_data argument",
         call. = FALSE)
  }


  assertthat::assert_that(show_center_values %in% c("none","median","mean","gm_mean"))







  if (!is.na(groupingFactor)) {
    assertthat::assert_that(groupingFactor %in% colnames(polya_data),msg=paste0(groupingFactor," is not a column of input dataset"))
    if (!is.na(condition1)) {
      if(!is.na(condition2)) {
        assertthat::assert_that(condition1 %in% levels(polya_data[[groupingFactor]]),msg=paste0(condition1," is not a level of ",grouping_factor," (groupingFactor)"))
        assertthat::assert_that(condition2 %in% levels(polya_data[[groupingFactor]]),msg=paste0(condition2," is not a level of ",grouping_factor," (groupingFactor)"))
        assertthat::assert_that(condition2 != condition1,msg="condition2 should be different than condition1")
        polya_data <- polya_data %>% dplyr::filter(!!rlang::sym(groupingFactor) %in% c(condition1,condition2))
      }
    }
    if(!is.na(subsample)) {
      polya_data <- subsample_table(polya_data,groupingFactor = groupingFactor,subsample=subsample)
    }
    distribution_plot <- ggplot2::ggplot(polya_data,ggplot2::aes_string(x=parameter_to_plot,color=groupingFactor))
    if (ndensity) {
      distribution_plot <- distribution_plot + ggplot2::geom_line(stat="density",size=1,ggplot2::aes(y=..ndensity..)) + ggplot2::xlab("normalized frequency")
    }
    else {
      distribution_plot <- distribution_plot + ggplot2::geom_line(stat="density",size=1,ggplot2::aes(y=..density..)) + ggplot2::xlab("densityy")
    }
    distribution_plot <- distribution_plot + ggplot2::theme_bw()

    if(show_center_values!="none") {
      center_values = polya_data %>% dplyr::group_by(!!rlang::sym(groupingFactor)) %>% dplyr::summarize(median_value = median(polya_length,na.rm = TRUE),mean_value=mean(polya_length,na.rm=TRUE),gm_mean_value=gm_mean(polya_length,na.rm=TRUE))
      if(show_center_values=='median') {
            distribution_plot <- distribution_plot + ggplot2::geom_vline(data=center_values,ggplot2::aes_string(xintercept="median_value",color=groupingFactor),linetype="longdash")
      }
      else if(show_center_values=='mean') {
        distribution_plot <- distribution_plot + ggplot2::geom_vline(data=center_values,ggplot2::aes_string(xintercept="mean_value",color=groupingFactor),linetype="longdash")
      }
      else if(show_center_values=='gm_mean') {
        distribution_plot <- distribution_plot + ggplot2::geom_vline(data=center_values,ggplot2::aes_string(xintercept="gm_mean_value",color=groupingFactor),linetype="longdash")
      }
    }
  }
  else {
    if(!is.na(subsample)) {
      polya_data <- subsample_table(polya_data,subsample=subsample)
    }
    distribution_plot <- ggplot2::ggplot(polya_data,ggplot2::aes_string(x=parameter_to_plot))
    if (ndensity) {
      distribution_plot <- distribution_plot + ggplot2::geom_line(stat="density",size=1,ggplot2::aes(y=..ndensity..)) + ggplot2::xlab("normalized frequency")
    }
    else {
      distribution_plot <- distribution_plot + ggplot2::geom_line(stat="density",size=1,ggplot2::aes(y=..density..)) + ggplot2::xlab("density")
    }

    distribution_plot <- distribution_plot  + ggplot2::theme_bw()
    if(show_center_values=="median") {
      distribution_plot <- distribution_plot + ggplot2::geom_vline(aes(xintercept=median(polya_length)),linetype="longdash")
    }
    else if (show_center_values=="mean") {
      distribution_plot <- distribution_plot + ggplot2::geom_vline(aes(xintercept=mean(polya_length)),linetype="longdash")
    }
  }

  distribution_plot <- .basic_aesthetics(distribution_plot,...)




  distribution_plot <- distribution_plot + nanotail_ggplot2_theme

  return(distribution_plot)

}




#' Plots boxplot of estimated polya lengths
#'
#' @param polya_data input table with polyA predictions
#' @param groupingFactor which factor to use for grouping
#' @param condition1 First condition to include on the plot
#' @param condition2 Second condition to include on the plot
#' @param violin Should violin plot be plotted instead of boxplot?
#' @param ... parameters passed to .basic_aesthetics function (scale_x_limit_low = NA, scale_x_limit_high = NA, scale_y_limit_low = NA, scale_y_limit_high = NA, color_palette = "Set1",plot_title=NA)
#'
#' @return \link[ggplot2]{ggplot} object
#' @export
#'
plot_polya_boxplot <- function(polya_data, groupingFactor,condition1=NA,condition2=NA,violin=FALSE,...) {


  if (missing(polya_data)) {
    stop("PolyA predictions are missing. Please provide a valid polya_data argument",
         call. = FALSE)
  }

  assertthat::assert_that(groupingFactor %in% colnames(polya_data),msg=paste0(groupingFactor," is not a column of input dataset"))


  if (!is.na(condition1)) {
    if(!is.na(condition2)) {
      assertthat::assert_that(condition1 %in% levels(polya_data[[groupingFactor]]),msg=paste0(condition1," is not a level of ",grouping_factor," (groupingFactor)"))
      assertthat::assert_that(condition2 %in% levels(polya_data[[groupingFactor]]),msg=paste0(condition2," is not a level of ",grouping_factor," (groupingFactor)"))
      assertthat::assert_that(condition2 != condition1,msg="condition2 should be different than condition1")
      polya_data <- polya_data %>% dplyr::filter(!!rlang::sym(groupingFactor) %in% c(condition1,condition2))
    }
  }


  transcripts_boxplot <- ggplot2::ggplot(polya_data,ggplot2::aes_string(x=groupingFactor,y="polya_length"))
  if(violin) {
    transcripts_boxplot <- transcripts_boxplot + ggplot2::geom_violin()
  }
  else{
    transcripts_boxplot <- transcripts_boxplot + ggplot2::geom_boxplot()
  }

  transcripts_boxplot <- .basic_aesthetics(transcripts_boxplot,...)

  return(transcripts_boxplot)
}


#' Title
#'
#' @param polya_data_summarized polyA predictions table, summarized using \link{summarize_polya}
#' @param groupingFactor name of column used for grouping
#' @param condition1 first condition to use for plotting
#' @param condition2 second condition to use for plotting
#' @param ... parameters passed to .basic_aesthetics function (scale_x_limit_low = NA, scale_x_limit_high = NA, scale_y_limit_low = NA, scale_y_limit_high = NA, color_palette = "Set1",plot_title=NA)
#' @param min_counts minimum number of counts to be shown
#' @param max_counts maximum number of counts to be shown
#' #'
#' @return \link[ggplot2]{ggplot} object
#' @export
#'
plot_counts_scatter <- function(polya_data_summarized, groupingFactor = NA, condition1 = NA, condition2 = NA,min_counts = 0, max_counts = 0,points_coloring_factor =NA, repel_elements=NA,repel_group=NA,...) {


  if (missing(polya_data_summarized)) {
    stop("Summarized PolyA predictions are missing. Please provide a valid polya_data argument",
         call. = FALSE)
  }


  assertthat::assert_that(!is.na(condition1),msg = "Please specify conditions for comparison")
  assertthat::assert_that(!is.na(condition2),msg = "Please specify conditions for comparison")
  assertthat::assert_that(!is.na(groupingFactor),msg = "Please specify groupingFactor")
  assertthat::assert_that(condition1 %in% levels(polya_data_summarized[[groupingFactor]]),msg=paste0(condition1," is not a level of ",groupingFactor," (groupingFactor)"))
  assertthat::assert_that(condition2 %in% levels(polya_data_summarized[[groupingFactor]]),msg=paste0(condition2," is not a level of ",groupingFactor," (groupingFactor)"))
  assertthat::assert_that(condition2 != condition1,msg="condition2 should be different than condition1")
  assertthat::assert_that("counts" %in% colnames(polya_data_summarized),msg = "Please provide summarized polya table (using summarize_polya()) as an input")



  if(!is.na(points_coloring_factor)){
    assertthat::assert_that(points_coloring_factor %in% colnames(polya_data_summarized),msg = "Please provide valid points_coloring_factor as an input")
    polya_data_summarized_counts_xy<-polya_data_summarized %>% dplyr::group_by(transcript,!!rlang::sym(groupingFactor),!!rlang::sym(points_coloring_factor)) %>% dplyr::summarize(counts_sum=sum(counts)) %>% tidyr::spread_(groupingFactor,"counts_sum")
  }
  else {
    polya_data_summarized_counts_xy<-polya_data_summarized %>% dplyr::group_by(transcript,!!rlang::sym(groupingFactor)) %>% dplyr::summarize(counts_sum=sum(counts)) %>% tidyr::spread_(groupingFactor,"counts_sum")
  }
  polya_data_summarized_counts_xy[is.na(polya_data_summarized_counts_xy)] <- 0

  polya_data_summarized_counts_xy <- polya_data_summarized_counts_xy %>% dplyr::filter(!!rlang::sym(condition1)>=min_counts,!!rlang::sym(condition2)>=min_counts,!!rlang::sym(condition1)<=max_counts,!!rlang::sym(condition2)<=max_counts)

  if (!is.na(max_counts)) {
    assertthat::assert_that(assertive::is_numeric(max_counts),msg="Please provide numeric value for max_counts")
    polya_data_summarized_counts_xy <- polya_data_summarized_counts_xy %>% dplyr::filter(!!rlang::sym(condition1)<=max_counts,!!rlang::sym(condition2)<=max_counts)
  }

  if(!is.na(points_coloring_factor)){
    counts_scatter_plot<-ggplot2::ggplot(polya_data_summarized_counts_xy,ggplot2::aes(x=!!rlang::sym(condition1),y=!!rlang::sym(condition2),colour=!!rlang::sym(points_coloring_factor))) + ggplot2::geom_point(ggplot2::aes(text=transcript),alpha=0.7)
  }
  else
    {
    counts_scatter_plot<-ggplot2::ggplot(polya_data_summarized_counts_xy,ggplot2::aes(x=!!rlang::sym(condition1),y=!!rlang::sym(condition2))) + ggplot2::geom_point(ggplot2::aes(text=transcript),alpha=0.7)
  }

  if (!is.na(repel_elements)) {
    assertthat::assert_that(assertive::is_numeric(repel_elements),msg="Please provide numeric paraemter for repel_elements")
    if(!is.na(repel_group)) {
    assertthat::assert_that(assertive::is_numeric(repel_elements),msg="Please provide numeric paraemter for repel_elements")
    counts_scatter_plot <- counts_scatter_plot + ggrepel::geom_text_repel(data=polya_data_summarized_counts_xy %>% dplyr::ungroup() %>% dplyr::filter(!!rlang::sym(points_coloring_factor) == repel_group) %>% dplyr::arrange(dplyr::desc(!!rlang::sym(condition1)))[1:repel_elements,], ggplot2::aes(label=polya_data_summarized_counts_xy %>% dplyr::ungroup() %>% dplyr::filter(!!rlang::sym(points_coloring_factor) == repel_group) %>% dplyr::arrange(dplyr::desc(!!rlang::sym(condition1))) %>% dplyr::select(transcript) %>% as.vector()[1:repel_elements]))
    }
    else {
      counts_scatter_plot <- counts_scatter_plot + ggrepel::geom_text_repel(data=polya_data_summarized_counts_xy %>% dplyr::ungroup() %>% dplyr::arrange(dplyr::desc(!!rlang::sym(condition1))) %>% dplyr::slice(1:repel_elements), ggplot2::aes(label=polya_data_summarized_counts_xy %>% dplyr::ungroup() %>% dplyr::arrange(dplyr::desc(!!rlang::sym(condition1))) %>% dplyr::slice(1:repel_elements) %>% dplyr::select(transcript) %>% as.vector()))
    }
  }
  counts_scatter_plot <- .basic_aesthetics(counts_scatter_plot,...)


  return(counts_scatter_plot)
}



#' Plot Nanopolish polya QC
#'
#' @param nanopolish_processing_info output of \link{get_nanopolish_processing_info}
#' @param frequency show frequency plot instead of counts plot
#' @param ... parameters passed to .basic_aesthetics function (scale_x_limit_low = NA, scale_x_limit_high = NA, scale_y_limit_low = NA, scale_y_limit_high = NA, color_palette = "Set1",plot_title=NA)
#'
#' @return \link[ggplot2]{ggplot} object
#' @export
#'
plot_nanopolish_qc <- function(nanopolish_processing_info, frequency=TRUE,...) {


  if (missing(nanopolish_processing_info)) {
    stop("nanopolish processing info is missing. Please provide a valid nanopolish_processing_info argument",
         call. = FALSE)
  }

  assertthat::assert_that(assertive::has_rows(nanopolish_processing_info),msg = "Empty data.frame provided as an input")
  basic_colnames = c("qc_tag","n")
  assertthat::assert_that(basic_colnames[1] %in% colnames(nanopolish_processing_info),msg="qc_tag column is missing in the input. Is that valid output of get_nanopolish_processing_info()?")
  assertthat::assert_that(basic_colnames[2] %in% colnames(nanopolish_processing_info),msg="n column is missing in the input. Is that valid output of get_nanopolish_processing_info()?")

  assertthat::assert_that(assertive::is_a_bool(frequency),msg="Non-boolean value provided for option frequency")


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

  nanopolish_qc_plot <- .basic_aesthetics(nanopolish_qc_plot,color_mode = "fill",...)
  return(nanopolish_qc_plot)
}



#' Plots volcano plot of differential expression analysis
#'

#' @param input_data a table with output from \link{calculate_diff_exp_binom} or \link{calculate_polya_stats}
#' @param ... parameters passed to .basic_aesthetics function (scale_x_limit_low = NA, scale_x_limit_high = NA, scale_y_limit_low = NA, scale_y_limit_high = NA, color_palette = "Set1",plot_title=NA)
#' #'
#' @return \link[ggplot2]{ggplot} object
#' @export
#'
plot_volcano <- function(input_data,...) {


  if (missing(input_data)) {
    stop("nanopolish processing info is missing. Please provide a valid nanopolish_processing_info argument",
         call. = FALSE)
  }

  assertthat::assert_that(assertive::has_rows(input_data),msg = "Empty data.frame provided as an input")
  assertthat::assert_that("fold_change" %in% colnames(input_data),msg = "Input table is not a valid input for plot_volcano(). fold_change column is missing.")
  assertthat::assert_that("padj" %in% colnames(input_data),msg = "Input table is not a valid input for plot_volcano(). padj column is missing.")
  assertthat::assert_that("significance" %in% colnames(input_data),msg = "Input table is not a valid input for plot_volcano(). significance column is missing.")
  assertthat::assert_that("transcript" %in% colnames(input_data),msg = "Input table is not a valid input for plot_volcano(). transcript column is missing.")


  volcano_plot <- ggplot2::ggplot(input_data,ggplot2::aes(x=log2(fold_change),y=-log10(padj),col=significance)) + ggplot2::geom_point(ggplot2::aes(text=transcript))


  volcano_plot <- .basic_aesthetics(volcano_plot,...)

  return(volcano_plot)

}


#' Plots MA plot of differential expression analysis
#'
#' Crates simple MA plot, with log10(mean expression) on the X-axis and log2(fold_change) on the Y-axis
#'
#'
#' @param input_data a table with output from \link{calculate_diff_exp_binom}
#' @param ... parameters passed to .basic_aesthetics function (scale_x_limit_low = NA, scale_x_limit_high = NA, scale_y_limit_low = NA, scale_y_limit_high = NA, color_palette = "Set1",plot_title=NA)
#'
#' @return \link[ggplot2]{ggplot} object
#' @export
#'

plot_MA <- function(input_data,...) {


  if (missing(input_data)) {
    stop("nanopolish processing info is missing. Please provide a valid nanopolish_processing_info argument",
         call. = FALSE)
  }

  assertthat::assert_that(assertive::has_rows(input_data),msg = "Empty data.frame provided as an input")
  assertthat::assert_that("fold_change" %in% colnames(input_data),msg = "Input table is not a valid input for plot_MA(). fold_change column is missing.")
  assertthat::assert_that("mean_expr" %in% colnames(input_data),msg = "Input table is not a valid input for plot_MA(). mean_expr column is missing.")
  assertthat::assert_that("significance" %in% colnames(input_data),msg = "Input table is not a valid input for plot_MA(). significance column is missing.")
  assertthat::assert_that("transcript" %in% colnames(input_data),msg = "Input table is not a valid input for plot_MA(). transcript column is missing.")


  MA_plot <- ggplot2::ggplot(input_data,ggplot2::aes(x=log10(mean_expr),y=log2(fold_change),col=significance)) + ggplot2::geom_point(ggplot2::aes(text=transcript))

  MA_plot <- .basic_aesthetics(MA_plot,...)

  return(MA_plot)

}

#' Title
#'
#' @param annotated_polya_data data frame(or tibble) with polyA predictions and associated annotations
#' @param grouping_factor column in polya_data_table specifing factor grouping samples
#' @param annotation_factor column specifying factor grouping transcripts by annotation
#' @param condition1 if only 2 conditions to show, choose which one is first
#' @param condition2 if only 2 conditions to show, choose which one is second
#' @param annotation_levels vector specifying selected annotation levels from annotation_factor
#' @param ... parameters passed to .basic_aesthetics function (scale_x_limit_low = NA, scale_x_limit_high = NA, scale_y_limit_low = NA, scale_y_limit_high = NA, color_palette = "Set1",plot_title=NA)
#'
#' @return
#' @export
#'
plot_annotations_comparison_boxplot <- function(annotated_polya_data,annotation_factor = NA,grouping_factor = NA,condition1=NA,condition2=NA,annotation_levels=c(),violin=FALSE,...) {

  if (missing(annotated_polya_data)) {
    stop("Annotated polya data table is missing. Please provide a valid input",
         call. = FALSE)
  }

  assertthat::assert_that(assertive::has_rows(annotated_polya_data),msg = "Empty data.frame provided as an input")


  if(!is.na(annotation_factor)) {
    assertthat::assert_that(annotation_factor %in% colnames(annotated_polya_data),msg=paste0(annotation_factor," is not a column of input data.frame"))
  }

  if (length(annotation_levels)>0) {
    #assertthat::assert_that(all(annotation_levels %in% levels(annotated_polya_data[[annotation_factor]])),msg="non-existing factor levels specified for annotation_levels parameter")
    annotated_polya_data <- annotated_polya_data %>% dplyr::filter(!!rlang::sym(annotation_factor) %in% annotation_levels)
    print(annotation_levels)
  }


  if (!is.na(condition1)) {
    if(!is.na(condition2)) {
      assertthat::assert_that(condition1 %in% levels(polya_data[[groupingFactor]]),msg=paste0(condition1," is not a level of ",grouping_factor," (groupingFactor)"))
      assertthat::assert_that(condition2 %in% levels(polya_data[[groupingFactor]]),msg=paste0(condition2," is not a level of ",grouping_factor," (groupingFactor)"))
      assertthat::assert_that(condition2 != condition1,msg="condition2 should be different than condition1")
      annotated_polya_data <- annotated_polya_data %>% dplyr::filter(!!rlang::sym(groupingFactor) %in% c(condition1,condition2))
    }
  }

  if (!is.na(grouping_factor)) {
    assertthat::assert_that(grouping_factor %in% colnames(annotated_polya_data))
    transcripts_boxplot <- ggplot2::ggplot(annotated_polya_data,ggplot2::aes_string(x=annotation_factor,y="polya_length",color=grouping_factor))
  }
  else {
    transcripts_boxplot <- ggplot2::ggplot(annotated_polya_data,ggplot2::aes_string(x=annotation_factor,y="polya_length"))
  }
  if(violin) {
    transcripts_boxplot <- transcripts_boxplot + ggplot2::geom_violin(position = ggplot2::position_dodge())
  }
  else{
    transcripts_boxplot <- transcripts_boxplot + ggplot2::geom_boxplot(position=ggplot2::position_dodge())
  }

  transcripts_boxplot <- .basic_aesthetics(transcripts_boxplot,...)
  transcripts_boxplot <- transcripts_boxplot + ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, vjust = 1, hjust = 1))
  return(transcripts_boxplot)

}


#' Title
#'
#' @param ggplot_object ggplot2 object to manipulate asesthetics
#' @param scale_x_limit_low lower limit of x continuous scale
#' @param scale_x_limit_high upper limit of x continuous scale
#' @param scale_y_limit_low lower limit of y continuous scale
#' @param scale_y_limit_high upper limit of y continuous scale
#' @param color_palette color palette (one from RColorBrewer of ggsci packages)
#' @param plot_title Title of the plot
#'
#' @return \link[ggplot2]{ggplot} object
#'
.basic_aesthetics <- function(ggplot_object,scale_x_limit_low = NA, scale_x_limit_high = NA, scale_y_limit_low = NA, scale_y_limit_high = NA, color_palette = "Set1",plot_title=NA,color_mode="color",axis_titles_size=16)
{


  if(missing(ggplot_object)) {
    stop("ggplot object required as the input")
  }

  if (!is.na(scale_x_limit_low)) {
    if (is.na(scale_x_limit_high)) {
      stop("Please provide both limits for x scale")
    }
    assertthat::assert_that(assertive::is_numeric(scale_x_limit_low),msg="Please provide numeric value for scale_x_limit_low")
    assertthat::assert_that(assertive::is_numeric(scale_x_limit_high),msg="Please provide numeric value for scale_x_limit_high")
    ggplot_object <- ggplot_object + ggplot2::scale_x_continuous(limits=c(scale_x_limit_low,scale_x_limit_high))
  }


  if (!is.na(scale_y_limit_low)) {
    if (is.na(scale_y_limit_high)) {
      stop("Please provide both limits for y scale")
    }
    assertthat::assert_that(assertive::is_numeric(scale_y_limit_low),msg="Please provide numeric value for scale_y_limit_low")
    assertthat::assert_that(assertive::is_numeric(scale_y_limit_high),msg="Please provide numeric value for scale_y_limit_high")
    ggplot_object <- ggplot_object + ggplot2::scale_y_continuous(limits=c(scale_y_limit_low,scale_y_limit_high))
  }


  valid_color_palettes_ggsci = c("npg", "aaas", "nejm", "lancet", "jama",
                          "jco", "ucscgb", "d3", "locuszoom",
                          "igv", "uchicago", "startrek", "tron",
                          "futurama", "rickandmorty", "simpsons")
  valid_color_palettes_RColorBrewer <- rownames(RColorBrewer::brewer.pal.info)

  if(!is.na(color_palette)) {
    if(!color_palette %in% c(valid_color_palettes_RColorBrewer,valid_color_palettes_ggsci)) {
      warning("Please provide valid color palette from RColorBrewer or ggsci packages")
    }
  }

  if(color_mode %in% c("color","both")) {
    if(color_palette %in% c(valid_color_palettes_ggsci)) {
      ggplot_object <- ggplot_object + eval(parse(text = paste0("ggsci::scale_color_",color_palette,"()")))
    }
    else if (color_palette %in% valid_color_palettes_RColorBrewer) {
      ggplot_object <- ggplot_object + ggplot2::scale_colour_brewer(palette = color_palette)
    }
  }
  else if (color_mode %in% c("fill","both")) {
    if(color_palette %in% c(valid_color_palettes_ggsci)) {
      ggplot_object <- ggplot_object + eval(parse(text = paste0("ggsci::scale_fill_",color_palette,"()")))
    }
    else if (color_palette %in% valid_color_palettes_RColorBrewer) {
      ggplot_object <- ggplot_object + ggplot2::scale_fill_brewer(palette = color_palette)
    }
  }

  if(!is.na(plot_title)){
    ggplot_object <- ggplot_object + ggplot2::ggtitle(plot_title)
  }

  return(ggplot_object)

}
